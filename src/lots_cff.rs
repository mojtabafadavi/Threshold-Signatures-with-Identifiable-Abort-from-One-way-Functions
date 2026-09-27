//! Tweaked CFF-based Lamport one-time signature (Section 2.1) and its
//! (T,N)-threshold version with identifiable abort (Sections 3 and 5).
//!
//! For a key with tweak `t` (its global index) and public parameters
//! `P_LOTS`, `P_msg` (`PubParams::lots`, `PubParams::lots_msg`, both keyed
//! on the long-term public key as Algorithm 4 specifies):
//!
//! * `LKeyGen(t)`: `x_i <- {0,1}^lambda`, `y_i = H(P_LOTS, t||i, x_i)`;
//! * `LSign(sk, msg, t)`: `B_msg = g(H(P_msg, t, msg))`, output `(x_i)_{i in B_msg}`;
//! * `LDerivePK(msg, sigma, pk_cmpl, t)`: `y'_i = H(P_LOTS, t||i, x_i)` for
//!   `i in B_msg`, and `y'_i = y_i` from `pk_cmpl` otherwise;
//! * `LVerify(pk, msg, sigma, t)`: accept iff `y'_i = y_i` for all `i in B_msg`.
//!
//! The threshold half (Section 3) Shamir-shares each `x_i` (`dist`), lets
//! each participant reveal only its shares of `B_msg` (`part_sign`), and
//! reconstructs `(x_i)_{i in B_msg}` from any `T` of them (`combine`).
//! Algorithm 4's verification keys (`party_vk`) and Section 5's
//! `TLOTS.SigProof` / `TLOTS.SignFailed` check (`sig_proof`,
//! `verify_partial`) give the identifiable abort.
//!
//! Components are indexed `0..e` here (the paper writes `[e]`).

use crate::cff::Cff;
use crate::hashing::{thash, tweak, tweak::LotsTweak, Digest32, PubParams};
use crate::shamir;
use rand::RngCore;
use std::sync::Arc;

#[derive(Clone)]
pub struct LotsParams {
    pub e: usize,
    pub kappa: usize,
    pub cff: Arc<Cff>,
}

impl LotsParams {
    pub fn new(e: usize, kappa: usize) -> Self {
        LotsParams {
            e,
            kappa,
            cff: Arc::new(Cff::new(e, kappa)),
        }
    }

    /// `B_msg = g_lambda(H(P_msg, t, msg))`, sorted ascending.
    pub fn block_for(&self, pp: &PubParams, t: LotsTweak, msg: &[u8]) -> Vec<usize> {
        let digest = thash(pp.lots_msg(), &tweak::lots_msg(t), msg);
        self.cff.digest_to_block(&digest)
    }
}

/// `y_i = H(P_LOTS, t||i, x_i)`.
pub fn pk_component(pp: &PubParams, t: LotsTweak, i: usize, x: &[u8; 32]) -> Digest32 {
    thash(pp.lots(), &tweak::lots_key(t, i as u32), x)
}

pub struct LotsKeyPair {
    /// The key's tweak `t` (its global leaf index).
    pub tweak: LotsTweak,
    pub sk: Vec<[u8; 32]>, // (x_0, ..., x_{e-1})
    pub pk: Vec<Digest32>, // (y_0, ..., y_{e-1})
}

/// `LKeyGen(1^lambda, t)`.
pub fn key_gen(params: &LotsParams, pp: &PubParams, t: LotsTweak, rng: &mut impl RngCore) -> LotsKeyPair {
    let mut sk = Vec::with_capacity(params.e);
    let mut pk = Vec::with_capacity(params.e);
    for i in 0..params.e {
        let mut x = [0u8; 32];
        rng.fill_bytes(&mut x);
        pk.push(pk_component(pp, t, i, &x));
        sk.push(x);
    }
    LotsKeyPair { tweak: t, sk, pk }
}

/// `LSign(sk, msg, t)`.
pub fn sign(params: &LotsParams, pp: &PubParams, t: LotsTweak, sk: &[[u8; 32]], msg: &[u8]) -> Vec<[u8; 32]> {
    params.block_for(pp, t, msg).iter().map(|&i| sk[i]).collect()
}

/// `pk_cmpl = (y_i)_{i not in B_msg}`, in ascending order of `i`.
pub fn complementary_pk(params: &LotsParams, pp: &PubParams, t: LotsTweak, pk: &[Digest32], msg: &[u8]) -> Vec<Digest32> {
    let block = params.block_for(pp, t, msg);
    (0..params.e).filter(|i| block.binary_search(i).is_err()).map(|i| pk[i]).collect()
}

/// `LDerivePK(msg, sigma, pk_cmpl, t)`: the full candidate key, or `None`
/// if the signature or `pk_cmpl` has the wrong length.
pub fn derive_pk(
    params: &LotsParams,
    pp: &PubParams,
    t: LotsTweak,
    msg: &[u8],
    sig: &[[u8; 32]],
    pk_cmpl: &[Digest32],
) -> Option<Vec<Digest32>> {
    let block = params.block_for(pp, t, msg);
    if sig.len() != block.len() || pk_cmpl.len() != params.e - block.len() {
        return None;
    }
    let (mut s, mut c) = (sig.iter(), pk_cmpl.iter());
    let mut out = Vec::with_capacity(params.e);
    for i in 0..params.e {
        if block.binary_search(&i).is_ok() {
            out.push(pk_component(pp, t, i, s.next()?));
        } else {
            out.push(*c.next()?);
        }
    }
    Some(out)
}

/// `LVerify(pk, msg, sigma, t)`.
pub fn verify(params: &LotsParams, pp: &PubParams, t: LotsTweak, pk: &[Digest32], msg: &[u8], sig: &[[u8; 32]]) -> bool {
    let block = params.block_for(pp, t, msg);
    pk.len() == params.e
        && block.len() == sig.len()
        && block.iter().zip(sig).all(|(&i, x)| pk_component(pp, t, i, x) == pk[i])
}

// ---------------------------------------------------------------------
// Threshold version: TCDist / TCPartSign / TCCombine / TCVerify
// ---------------------------------------------------------------------

/// Per-participant share bundle for a single LOTS keypair: for each of the
/// `e` secret components, this participant's R Shamir shares (one per
/// BPB-bit block), together with the participant's index `j` (the paper's
/// party identifier, which appears in the verification-key tweaks
/// `index||j||k`) and the field point (`point`) those shares were evaluated
/// at: `omega^j` for `share_ntt`'s root-of-unity domain, not the plain
/// integer `j+1` of a per-point sharing. Carrying `point` with the shares
/// means a caller cannot supply a mismatched point when reconstructing.
pub struct PartyShare {
    /// The participant's index `j` in `0..N` (Algorithm 4's `j`).
    pub id: u32,
    pub point: u64,
    pub blocks: Vec<[u64; shamir::R]>, // length e
}

/// `pi_{ijk} = H_vk(index||j||k, sh_{ij,k})` (Algorithm 4, line 3): the
/// commitment to participant `j`'s share of component `k` of the LOTS key
/// with global index `index`.
pub fn share_commitment(pp: &PubParams, index: LotsTweak, party: u32, k: usize, share: &[u64; shamir::R]) -> Digest32 {
    thash(
        pp.vk(),
        &tweak::vk_share(index, party, k as u32),
        &shamir::share_to_bytes(share),
    )
}

/// `vk_{ij} = H_vk(index||j, (pi_{ijk})_{k in [e]})` (Algorithm 4, line 4):
/// the verification key committing to participant `j`'s whole share vector.
pub fn party_vk(pp: &PubParams, index: LotsTweak, share: &PartyShare) -> Digest32 {
    let mut pis = Vec::with_capacity(share.blocks.len() * 32);
    for (k, blk) in share.blocks.iter().enumerate() {
        pis.extend_from_slice(&share_commitment(pp, index, share.id, k, blk));
    }
    thash(pp.vk(), &tweak::vk_party(index, share.id), &pis)
}

/// The verification-key vector `vk = (vk_{ij})_{j in [N]}` of one LOTS key
/// (Algorithm 4, lines 3-4). Public: every participant receives all of it.
pub fn verification_keys(pp: &PubParams, index: LotsTweak, shares: &[PartyShare]) -> Vec<Digest32> {
    shares.iter().map(|sh| party_vk(pp, index, sh)).collect()
}

/// TCDist(sk, T, N) / Algorithm 4's `IssueShares`: Shamir-share every
/// secret component of `sk` via `shamir::share_ntt` (O(e * M log M) total,
/// M = next_pow2(N), instead of the naive O(e * N * T)) and return the
/// per-participant share bundles together with the verification-key vector
/// `vk = (vk_{ij})_{j in [N]}`, which is public and given to every
/// participant. `index` is the LOTS key's global index (its tweak).
pub fn dist(
    pp: &PubParams,
    index: LotsTweak,
    sk: &[[u8; 32]],
    t: usize,
    n: usize,
    rng: &mut impl RngCore,
) -> (Vec<PartyShare>, Vec<Digest32>) {
    let e = sk.len();
    let domain = shamir::eval_domain(n);
    let mut parties: Vec<PartyShare> = domain
        .iter()
        .enumerate()
        .map(|(j, &point)| PartyShare {
            id: j as u32,
            point,
            blocks: vec![[0u64; shamir::R]; e],
        })
        .collect();
    for i in 0..e {
        let blocks = shamir::secret_to_blocks(&sk[i]);
        for r in 0..shamir::R {
            let shares = shamir::share_ntt(blocks[r], t, n, rng);
            for j in 0..n {
                parties[j].blocks[i][r] = shares[j];
            }
        }
    }
    let vk = verification_keys(pp, index, &parties);
    (parties, vk)
}

/// TCPartSign(msg, j, sh_j): reveal this participant's shares only for the
/// kappa indices in B_msg.
pub fn part_sign(
    params: &LotsParams,
    pp: &PubParams,
    t: LotsTweak,
    msg: &[u8],
    share: &PartyShare,
) -> Vec<[u64; shamir::R]> {
    let block = params.block_for(pp, t, msg);
    block.iter().map(|&i| share.blocks[i]).collect()
}

/// `TLOTS.SigProof` (Algorithm 6): the commitments to the shares this
/// participant did not reveal, `pi = (H_vk(index||j||k, sh_{ij,k}))_{k not in B_msg}`,
/// in ascending order of `k`. Together with a partial signature it reopens
/// the whole of `vk_{ij}`. Algorithm 6 indexes `pi` by all of `[e]` and
/// leaves the `B_msg` positions empty; since only the other `e - kappa`
/// entries carry a value, they are returned here on their own, and
/// `verify_partial` reads them back in the same order.
pub fn sig_proof(
    params: &LotsParams,
    pp: &PubParams,
    index: LotsTweak,
    msg: &[u8],
    share: &PartyShare,
) -> Vec<Digest32> {
    let block = params.block_for(pp, index, msg);
    (0..params.e)
        .filter(|k| block.binary_search(k).is_err())
        .map(|k| share_commitment(pp, index, share.id, k, &share.blocks[k]))
        .collect()
}

/// The check run by `TLOTS.SignFailed` (Algorithm 5, lines 5-9): recompute
/// `vk_{ij}` from the partial signature (for `k in B_msg`) and the proof
/// (for `k not in B_msg`). Returns true iff participant `party`'s partial
/// signature is consistent with the shares it was dealt -- an honest
/// participant always passes, a participant that submitted a malformed
/// partial signature cannot, which is what makes the abort identifiable.
pub fn verify_partial(
    params: &LotsParams,
    pp: &PubParams,
    index: LotsTweak,
    party: u32,
    vk_party: &Digest32,
    msg: &[u8],
    partial: &[[u64; shamir::R]],
    proof: &[Digest32],
) -> bool {
    let block = params.block_for(pp, index, msg);
    if partial.len() != block.len() || proof.len() != params.e - block.len() {
        return false;
    }
    let (mut revealed, mut hidden) = (partial.iter(), proof.iter());
    let mut pis = Vec::with_capacity(params.e * 32);
    for k in 0..params.e {
        let h = if block.binary_search(&k).is_ok() {
            match revealed.next() {
                Some(sh) => share_commitment(pp, index, party, k, sh),
                None => return false,
            }
        } else {
            match hidden.next() {
                Some(h) => *h,
                None => return false,
            }
        };
        pis.extend_from_slice(&h);
    }
    thash(pp.vk(), &tweak::vk_party(index, party), &pis) == *vk_party
}

/// TCCombine(msg, {partial sigs}) / Algorithm 6's
/// `TLOTS.ReconstructSignature`: given >= T partial signatures (each tagged
/// with the participant's field point), reconstruct the ordinary
/// CFF-Lamport signature on `msg` under the key with tweak `tw`.
///
/// The first `t` partials form the quorum S, whose Lagrange coefficients are
/// computed once and reused across all `kappa * R` coordinates, so the cost
/// is Theta(T^2 + kappa*R*T) rather than Theta(kappa*R*T^2) (Section 8).
///
/// Returns `None` if there are too few partials, one has the wrong length,
/// or the quorum repeats a participant's point -- the latter mirroring
/// Algorithm 5's refusal to store two partial signatures from one sender.
pub fn combine(
    params: &LotsParams,
    pp: &PubParams,
    tw: LotsTweak,
    msg: &[u8],
    t: usize,
    partials: &[(u64, Vec<[u64; shamir::R]>)],
) -> Option<Vec<[u8; 32]>> {
    if partials.len() < t {
        return None;
    }
    let kappa = params.block_for(pp, tw, msg).len();
    let quorum = &partials[..t];
    if quorum.iter().any(|(_, p)| p.len() != kappa) {
        return None;
    }
    let xs: Vec<u64> = quorum.iter().map(|&(x, _)| x).collect();
    if !shamir::has_distinct_points(&xs) {
        return None;
    }
    let coeffs = shamir::lagrange_coeffs(&xs);

    let mut column = vec![0u64; t];
    let mut sig = Vec::with_capacity(kappa);
    for pos in 0..kappa {
        let mut blocks_out = [0u64; shamir::R];
        for r in 0..shamir::R {
            for (slot, (_, shares)) in column.iter_mut().zip(quorum) {
                *slot = shares[pos][r];
            }
            blocks_out[r] = shamir::reconstruct_with(&coeffs, &column);
        }
        sig.push(shamir::blocks_to_secret(&blocks_out));
    }
    Some(sig)
}

pub fn tc_verify(
    params: &LotsParams,
    pp: &PubParams,
    t: LotsTweak,
    pk: &[Digest32],
    msg: &[u8],
    sig: &[[u8; 32]],
) -> bool {
    verify(params, pp, t, pk, msg, sig)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    fn pp() -> PubParams {
        PubParams::with_pk(&[5u8; 32], &[3u8; 32])
    }

    #[test]
    fn plain_sign_verify_and_derive_pk() {
        let mut rng = thread_rng();
        let (params, pp) = (LotsParams::new(24, 10), pp());
        let t = 42;
        let kp = key_gen(&params, &pp, t, &mut rng);
        let msg = b"plain lamport";
        let sig = sign(&params, &pp, t, &kp.sk, msg);
        assert!(verify(&params, &pp, t, &kp.pk, msg, &sig));
        assert!(!verify(&params, &pp, t, &kp.pk, b"other", &sig));
        let cmpl = complementary_pk(&params, &pp, t, &kp.pk, msg);
        assert_eq!(cmpl.len(), params.e - params.kappa);
        assert_eq!(derive_pk(&params, &pp, t, msg, &sig, &cmpl), Some(kp.pk.clone()));
    }

    #[test]
    fn tweak_and_parameter_bind_the_key() {
        let mut rng = thread_rng();
        let (params, pp) = (LotsParams::new(24, 10), pp());
        let kp = key_gen(&params, &pp, 7, &mut rng);
        let msg = b"m";
        let sig = sign(&params, &pp, 7, &kp.sk, msg);
        // Same secret, different tweak or parameter: different public key.
        assert!(!verify(&params, &pp, 8, &kp.pk, msg, &sig));
        assert!(!verify(&params, &PubParams::with_pk(&[5u8; 32], &[4u8; 32]), 7, &kp.pk, msg, &sig));
        // The parameters are keyed on the whole public key pk = (Root, id),
        // so a different Root also invalidates the key (Algorithm 4).
        assert!(!verify(&params, &PubParams::with_pk(&[6u8; 32], &[3u8; 32]), 7, &kp.pk, msg, &sig));
    }

    #[test]
    fn threshold_sign_verify() {
        let mut rng = thread_rng();
        let (params, pp) = (LotsParams::new(24, 10), pp());
        let tw = 5;
        let kp = key_gen(&params, &pp, tw, &mut rng);
        let (t, n) = (5, 9);
        let (shares, _vk) = dist(&pp, tw, &kp.sk, t, n, &mut rng);
        let msg = b"hello threshold lamport";
        let partials: Vec<_> = shares
            .iter()
            .take(t)
            .map(|sh| (sh.point, part_sign(&params, &pp, tw, msg, sh)))
            .collect();
        let sig = combine(&params, &pp, tw, msg, t, &partials).unwrap();
        assert!(tc_verify(&params, &pp, tw, &kp.pk, msg, &sig));

        // A different quorum reconstructs the identical signature.
        let partials2: Vec<_> = shares
            .iter()
            .skip(n - t)
            .map(|sh| (sh.point, part_sign(&params, &pp, tw, msg, sh)))
            .collect();
        assert_eq!(combine(&params, &pp, tw, msg, t, &partials2).unwrap(), sig);
    }

    /// Combining shares the quorum's Lagrange coefficients across every
    /// coordinate, so a quorum that repeats one participant must be refused
    /// rather than interpolated through a degenerate coefficient set.
    #[test]
    fn combine_rejects_malformed_quorums() {
        let mut rng = thread_rng();
        let (params, pp) = (LotsParams::new(24, 10), pp());
        let tw = 9;
        let kp = key_gen(&params, &pp, tw, &mut rng);
        let (t, n) = (4, 8);
        let (shares, _vk) = dist(&pp, tw, &kp.sk, t, n, &mut rng);
        let msg = b"quorum hygiene";
        let good: Vec<_> = shares
            .iter()
            .take(t)
            .map(|sh| (sh.point, part_sign(&params, &pp, tw, msg, sh)))
            .collect();
        assert!(combine(&params, &pp, tw, msg, t, &good).is_some());

        // The same participant twice.
        let mut repeated = good.clone();
        repeated[1] = good[0].clone();
        assert!(combine(&params, &pp, tw, msg, t, &repeated).is_none());

        // Too few partials, and a partial of the wrong length.
        assert!(combine(&params, &pp, tw, msg, t, &good[..t - 1]).is_none());
        let mut short = good.clone();
        short[0].1.pop();
        assert!(combine(&params, &pp, tw, msg, t, &short).is_none());
    }

    #[test]
    fn verification_keys_accept_honest_partial_signatures() {
        let mut rng = thread_rng();
        let (params, pp) = (LotsParams::new(24, 10), pp());
        let tw = 11;
        let kp = key_gen(&params, &pp, tw, &mut rng);
        let (t, n) = (4, 7);
        let (shares, vk) = dist(&pp, tw, &kp.sk, t, n, &mut rng);
        assert_eq!(vk.len(), n);
        let msg = b"identifiable abort";
        for sh in &shares {
            let partial = part_sign(&params, &pp, tw, msg, sh);
            let proof = sig_proof(&params, &pp, tw, msg, sh);
            assert_eq!(proof.len(), params.e - params.kappa);
            assert!(verify_partial(&params, &pp, tw, sh.id, &vk[sh.id as usize], msg, &partial, &proof));
            // The verification key binds the participant, the key index and
            // the message's block: none of these may be swapped.
            assert!(!verify_partial(&params, &pp, tw, sh.id ^ 1, &vk[sh.id as usize], msg, &partial, &proof));
            assert!(!verify_partial(&params, &pp, tw + 1, sh.id, &vk[sh.id as usize], msg, &partial, &proof));
        }
    }

    /// The full identifiable-abort flow of Section 3: one participant sends
    /// a malformed partial signature, combination yields an invalid
    /// signature, the culprit is pinned down by its verification key, and a
    /// valid signature is recovered from the remaining honest shares.
    #[test]
    fn malformed_partial_signature_is_identified_and_recovered_from() {
        let mut rng = thread_rng();
        let (params, pp) = (LotsParams::new(24, 10), pp());
        let tw = 3;
        let kp = key_gen(&params, &pp, tw, &mut rng);
        let (t, n) = (4, 8);
        let (shares, vk) = dist(&pp, tw, &kp.sk, t, n, &mut rng);
        let msg = b"one corrupt signer";

        let mut submitted: Vec<(u32, u64, Vec<[u64; shamir::R]>)> = shares
            .iter()
            .map(|sh| (sh.id, sh.point, part_sign(&params, &pp, tw, msg, sh)))
            .collect();
        let culprit = 2usize;
        submitted[culprit].2[0][0] ^= 1; // corrupt one field element

        let quorum: Vec<(u64, Vec<_>)> = submitted[..t].iter().map(|(_, p, s)| (*p, s.clone())).collect();
        let bad = combine(&params, &pp, tw, msg, t, &quorum).unwrap();
        assert!(!tc_verify(&params, &pp, tw, &kp.pk, msg, &bad), "corrupt share must break the signature");

        // Everyone opens their verification key; exactly the culprit fails.
        let blamed: Vec<u32> = submitted
            .iter()
            .zip(&shares)
            .filter(|((id, _, partial), sh)| {
                let proof = sig_proof(&params, &pp, tw, msg, sh);
                !verify_partial(&params, &pp, tw, *id, &vk[*id as usize], msg, partial, &proof)
            })
            .map(|((id, _, _), _)| *id)
            .collect();
        assert_eq!(blamed, vec![culprit as u32]);

        // Dropping the blamed participant, the remaining honest ones sign.
        let honest: Vec<(u64, Vec<_>)> = submitted
            .iter()
            .filter(|(id, _, _)| *id != culprit as u32)
            .take(t)
            .map(|(_, p, s)| (*p, s.clone()))
            .collect();
        let good = combine(&params, &pp, tw, msg, t, &honest).unwrap();
        assert!(tc_verify(&params, &pp, tw, &kp.pk, msg, &good));
    }
}
