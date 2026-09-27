//! Benchmark harness: `bench [h] [d] [leaf]` (defaults: 64 4 0).

use rand::{thread_rng, RngCore};
use std::time::Instant;
use tlots::hashing::N;
use tlots::hypertree::{self, Dealer, HyperParams};
use tlots::lots_cff::{self, LotsParams};
use tlots::wotsplus;

fn ms(d: std::time::Duration) -> f64 {
    d.as_secs_f64() * 1000.0
}

/// Reconstruct the CFF-Lamport signature the way a per-coordinate Lagrange
/// interpolation does, recomputing the quorum's coefficients for each of the
/// `kappa * R` coordinates. `lots_cff::combine` computes them once instead;
/// this is here only to measure what that saves.
fn reconstruct_per_coordinate(
    partials: &[(u64, Vec<[u64; tlots::shamir::R]>)],
    t: usize,
    kappa: usize,
) -> Vec<[u8; 32]> {
    let quorum = &partials[..t];
    (0..kappa)
        .map(|pos| {
            let mut blocks = [0u64; tlots::shamir::R];
            for (r, out) in blocks.iter_mut().enumerate() {
                let points: Vec<(u64, u64)> =
                    quorum.iter().map(|(x, shares)| (*x, shares[pos][r])).collect();
                *out = tlots::shamir::reconstruct(&points);
            }
            tlots::shamir::blocks_to_secret(&blocks)
        })
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let h: u32 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(64);
    let d: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);
    let leaf: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);
    let (e, kappa) = (261, 123);

    let mut rng = thread_rng();
    let (mut id, mut master) = ([0u8; N], [0u8; N]);
    rng.fill_bytes(&mut id);
    rng.fill_bytes(&mut master);
    let params = HyperParams::new(h, d, LotsParams::new(e, kappa), id, master).expect("invalid (h, d)");
    let lh = params.layer_height();
    println!("h={h}, d={d} (subtree height {lh}), e={e}, kappa={kappa}, leaf={leaf}");

    let t0 = Instant::now();
    let mut dealer = Dealer::key_gen(params);
    println!("KeyGen (top subtree): {:.1} s", t0.elapsed().as_secs_f64());

    let t0 = Instant::now();
    let (kp, cert, kts, tts) = dealer.resolve_leaf_timed(leaf);
    println!("resolve_leaf: {:.1} s", t0.elapsed().as_secs_f64());
    for (k, (kt, tt)) in kts.iter().zip(&tts).enumerate() {
        let what = if k == 0 { "LOTS ".to_string() } else { format!("WOTS+ r={k}") };
        println!("  {what:<10} keygen {:>10.1} ms, tree {:>8.1} ms", ms(*kt), ms(*tt));
    }

    println!("\nCache demonstration (issue_next_leaf, rebuilt? bottom first):");
    for _ in 0..3 {
        let t0 = Instant::now();
        let (_, _, rebuilt) = dealer.issue_next_leaf_verbose().unwrap();
        println!("  {:>10.3} ms  {:?}", ms(t0.elapsed()), rebuilt);
    }

    // Sizes, grouped as in the paper's size table.
    let wots_sig = wotsplus::LEN * N;
    let path = lh as usize * N;
    let index_bytes = 8; // the global LOTS index, a u64: h/8 B at h = 64
    let rows = [
        ("pk^B = (root, id_MT)", 2 * N),
        ("sigma^L, pk_cmpl", e * N),
        ("id_msg, salt, MT_0.Path, index", (lh as usize + 2) * N + index_bytes),
        ("Sigma: d-1 pairs (sig^W_r, Path_r)", (d as usize - 1) * (wots_sig + path)),
    ];
    println!("\nSizes:");
    for (name, b) in rows {
        println!("  {name:<36} {b:>7} B");
    }
    let total: usize = rows[1..].iter().map(|r| r.1).sum();
    println!("  {:<36} {:>7} B", "total |sigma|", total);
    println!("  {:<36} {:>7} B", "share of one leaf (sh_j)", e * tlots::shamir::R * tlots::shamir::ELEM_BYTES);
    println!("  {:<36} {:>7} B", "sig proof pi (e - kappa hashes)", (e - kappa) * N);

    println!(
        "\n{:>6} {:>6} {:>10} {:>10} {:>10} {:>12} {:>10} {:>10} {:>14}",
        "N", "T", "Dist", "(of vk)", "PartSign", "Reconstruct", "Verify", "Total", "(per-coord.)"
    );
    let msg = b"benchmark message";
    for &(n, t) in &[(10usize, 5usize), (50, 25), (100, 25), (100, 50), (400, 100), (1000, 200)] {
        // Dist = Algorithm 4's IssueShares: Shamir sharing plus the
        // verification-key commitments that make the abort identifiable.
        let t0 = Instant::now();
        let (shares, vk) = lots_cff::dist(&dealer.pp, leaf, &kp.sk, t, n, &mut rng);
        let dist = t0.elapsed();
        let t0 = Instant::now();
        let vk2 = lots_cff::verification_keys(&dealer.pp, leaf, &shares);
        let vk_time = t0.elapsed(); // the vk share of Dist, measured on its own
        assert_eq!(vk, vk2);

        let t0 = Instant::now();
        let partials: Vec<_> = shares
            .iter()
            .take(t)
            .map(|sh| hypertree::part_sign(&dealer.params.lots, &dealer.pp, leaf, msg, sh))
            .collect();
        let part = t0.elapsed(); // all T signers together

        let t0 = Instant::now();
        let sig = hypertree::combine(&dealer.params.lots, &dealer.pp, msg, t, leaf, cert.clone(), &partials).unwrap();
        let comb = t0.elapsed();

        let t0 = Instant::now();
        let ok = hypertree::verify(&dealer.params, &dealer.id, &dealer.pk, msg, &sig);
        let ver = t0.elapsed();
        assert!(ok, "signature failed to verify");
        assert_eq!(sig.size_bytes(), total);

        // The same interpolation with the Lagrange coefficients recomputed for
        // every coordinate, which is what sharing them across the kappa * R
        // coordinates of one quorum saves.
        let t0 = Instant::now();
        let per_coord = reconstruct_per_coordinate(&partials, t, kappa);
        let comb_naive = t0.elapsed();
        assert_eq!(per_coord, sig.sigma_cff);

        println!(
            "{:>6} {:>6} {:>10.3} {:>10.3} {:>10.3} {:>12.3} {:>10.3} {:>10.3} {:>14.3}",
            n,
            t,
            ms(dist),
            ms(vk_time),
            ms(part),
            ms(comb),
            ms(ver),
            ms(dist + part + comb + ver),
            ms(comb_naive)
        );
    }

    // Identifiable abort (Section 3): one participant submits a malformed
    // partial signature, so combining fails; every participant opens its
    // verification key, which names the culprit, and the honest remainder
    // still produces a valid signature.
    let (n, t) = (50usize, 25usize);
    let (shares, vk) = lots_cff::dist(&dealer.pp, leaf, &kp.sk, t, n, &mut rng);
    let mut partials: Vec<_> = shares
        .iter()
        .map(|sh| hypertree::part_sign(&dealer.params.lots, &dealer.pp, leaf, msg, sh))
        .collect();
    let culprit = 7usize;
    partials[culprit].1[0][0] ^= 1;

    let bad = hypertree::combine(
        &dealer.params.lots,
        &dealer.pp,
        msg,
        t,
        leaf,
        cert.clone(),
        &partials[..t],
    )
    .unwrap();
    let bad_ok = hypertree::verify(&dealer.params, &dealer.id, &dealer.pk, msg, &bad);

    let t0 = Instant::now();
    let proofs: Vec<_> = shares
        .iter()
        .map(|sh| lots_cff::sig_proof(&dealer.params.lots, &dealer.pp, leaf, msg, sh))
        .collect();
    let proof_time = t0.elapsed();

    let t0 = Instant::now();
    let blamed: Vec<u32> = shares
        .iter()
        .zip(&partials)
        .zip(&proofs)
        .filter(|((sh, (_, partial)), proof)| {
            !lots_cff::verify_partial(
                &dealer.params.lots,
                &dealer.pp,
                leaf,
                sh.id,
                &vk[sh.id as usize],
                msg,
                partial,
                proof,
            )
        })
        .map(|((sh, _), _)| sh.id)
        .collect();
    let check_time = t0.elapsed();

    let honest: Vec<_> = partials
        .iter()
        .enumerate()
        .filter(|(j, _)| !blamed.contains(&(*j as u32)))
        .take(t)
        .map(|(_, p)| p.clone())
        .collect();
    let good = hypertree::combine(&dealer.params.lots, &dealer.pp, msg, t, leaf, cert, &honest).unwrap();
    let good_ok = hypertree::verify(&dealer.params, &dealer.id, &dealer.pk, msg, &good);

    println!("\nIdentifiable abort (N={n}, T={t}, participant {culprit} corrupted):");
    println!("  signature from the corrupted quorum verifies: {bad_ok}");
    println!("  blamed by their verification keys:            {blamed:?}");
    println!("  signature after dropping them verifies:       {good_ok}");
    println!(
        "  SigProof (all {n} signers): {:.3} ms, checking all {n} proofs: {:.3} ms",
        ms(proof_time),
        ms(check_time)
    );
    assert!(!bad_ok && good_ok && blamed == vec![culprit as u32]);
}
