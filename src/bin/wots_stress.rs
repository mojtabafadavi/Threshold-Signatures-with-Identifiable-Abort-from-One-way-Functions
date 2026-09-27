//! Regression fuzz check for tweaked WOTS+: random public parameters,
//! seeds, key tweaks and messages; every signature must verify, and must
//! fail under a different tweak.

use rand::{thread_rng, Rng, RngCore};
use tlots::hashing::{tweak::WotsTweak, PubParams, N};
use tlots::wotsplus;

fn main() {
    let rounds: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(20_000);
    let mut rng = thread_rng();
    let mut failures = 0usize;
    for round in 0..rounds {
        let (mut id_mt, mut seed) = ([0u8; N], [0u8; N]);
        rng.fill_bytes(&mut id_mt);
        rng.fill_bytes(&mut seed);
        let pp = PubParams::new(&id_mt);
        let t = WotsTweak { r: rng.gen_range(1..8), idx: rng.gen() };
        let mut msg = vec![0u8; rng.gen_range(0..128)];
        rng.fill_bytes(&mut msg);

        let kp = wotsplus::key_gen(&pp, &seed, &t);
        let sig = wotsplus::sign(&pp, &seed, &t, &msg);
        let ok = wotsplus::verify(&pp, &kp.pk, &t, &msg, &sig);
        let other = WotsTweak { r: t.r, idx: t.idx.wrapping_add(1) };
        let wrong_ok = wotsplus::verify(&pp, &kp.pk, &other, &msg, &sig);
        if !ok || wrong_ok {
            failures += 1;
            eprintln!("round {round}: verify={ok}, verify_with_wrong_tweak={wrong_ok}");
        }
    }
    println!("{rounds} rounds, {failures} failures");
    if failures > 0 {
        std::process::exit(1);
    }
}
