//! Tweaked WOTS+ (Section 2.1 and Appendix A.2, "Winternitz one-time
//! signature scheme"), used for the upper certification layers of the hypertree.
//!
//! With key tweak `t` (`WotsTweak`), public parameters `P_WOTS` and
//! `P_msg` (`PubParams::wots`, `PubParams::wots_msg` -- keyed on the public
//! key's `id`, since these keys are hashed while the long-term root is
//! still being computed; see `hashing`), and the chain
//!
//! ```text
//!   H^0(P, t||i||j, x)   = x
//!   H^{k+1}(P, t||i||j, x) = H(P, t||i||(j+k), H^k(P, t||i||j, x)),
//! ```
//!
//! the algorithms are
//!
//! * `WKeyGen(rho, t)`: `x_i = PRF(rho, t||i)`, `y_i = H^{2^w-1}(P_WOTS, t||i||0, x_i)`;
//! * `WSign(sk, msg, t)`: `D = H(P_msg, t, msg)` in base `2^w` plus checksum,
//!   `sigma_i = H^{D_i}(P_WOTS, t||i||0, x_i)`;
//! * `WDerivePK(msg, sigma, t)`: `y'_i = H^{2^w-1-D_i}(P_WOTS, t||i||D_i, sigma_i)`,
//!   i.e. each chain continues from height `D_i`.
//!
//! Base-`w` encoding and checksum follow XMSS (`wots.c`); the left shift
//! used there does not change the checksum digits.

use crate::hashing::{prf, thash, tweak::WotsTweak, Digest32, PubParams, N};

pub const LOG_W: u32 = 4; // the paper's Winternitz parameter w
pub const W: u32 = 1 << LOG_W; // 2^w = 16
pub const LEN1: usize = (8 * N as u32 / LOG_W) as usize; // l1 = 64
pub const LEN2: usize = 3; // l2 = floor(log2(l1 (2^w - 1)) / w) + 1 = 3
pub const LEN: usize = LEN1 + LEN2; // l = 67

pub struct WotsKeyPair {
    /// Secret seed `rho`; chain `i` starts at `PRF(rho, t||i)`.
    pub sk_seed: [u8; N],
    pub tweak: WotsTweak,
    /// `(y_0, ..., y_{l-1})`.
    pub pk: Vec<Digest32>,
}

fn base_w(msg: &[u8], out_len: usize) -> Vec<u32> {
    let mut out = Vec::with_capacity(out_len);
    let mut bits = 0i32;
    let mut total: u32 = 0;
    let mut in_idx = 0usize;
    while out.len() < out_len {
        if bits <= 0 {
            total = msg[in_idx] as u32;
            in_idx += 1;
            bits += 8;
        }
        bits -= LOG_W as i32;
        out.push((total >> bits) & (W - 1));
    }
    out
}

fn checksum(digits: &[u32]) -> Vec<u32> {
    let csum: u32 = digits.iter().map(|d| W - 1 - d).sum();
    let total_bits = LEN2 as u32 * LOG_W;
    let shift = (8 - (total_bits % 8)) % 8;
    let csum = csum << shift;
    let nbytes = ((total_bits + 7) / 8) as usize;
    let full = csum.to_be_bytes();
    base_w(&full[4 - nbytes..], LEN2)
}

/// `(D_0, ..., D_{l-1})`: base-`2^w` digits of `D = H(P_msg, t, msg)`
/// followed by the checksum digits.
pub fn digits_for(pp: &PubParams, t: &WotsTweak, msg: &[u8]) -> Vec<u32> {
    let digest = thash(pp.wots_msg(), &t.msg(), msg);
    let mut d = base_w(&digest, LEN1);
    d.extend(checksum(&d));
    d
}

/// `H^k(P_WOTS, t||i||j, x)`: `k` chain steps on `x` starting at height `j`.
pub fn chain(pp: &PubParams, t: &WotsTweak, i: u32, j: u32, k: u32, x: &Digest32) -> Digest32 {
    assert!(j + k <= W - 1, "chain would exceed its maximal height");
    let mut v = *x;
    for step in j..j + k {
        v = thash(pp.wots(), &t.chain(i, step), &v);
    }
    v
}

fn chain_start(sk_seed: &[u8; N], t: &WotsTweak, i: u32) -> Digest32 {
    prf(sk_seed, &t.sk(i))
}

/// `WKeyGen(rho, t)`.
pub fn key_gen(pp: &PubParams, sk_seed: &[u8; N], t: &WotsTweak) -> WotsKeyPair {
    let pk = (0..LEN as u32)
        .map(|i| chain(pp, t, i, 0, W - 1, &chain_start(sk_seed, t, i)))
        .collect();
    WotsKeyPair {
        sk_seed: *sk_seed,
        tweak: *t,
        pk,
    }
}

/// `WSign(sk, msg, t)`.
pub fn sign(pp: &PubParams, sk_seed: &[u8; N], t: &WotsTweak, msg: &[u8]) -> Vec<Digest32> {
    let digits = digits_for(pp, t, msg);
    (0..LEN)
        .map(|i| chain(pp, t, i as u32, 0, digits[i], &chain_start(sk_seed, t, i as u32)))
        .collect()
}

/// `WDerivePK(msg, sigma, pk_cmpl = {}, t)`. Returns `None` on a
/// malformed signature (wrong length).
pub fn derive_pk(pp: &PubParams, t: &WotsTweak, msg: &[u8], sig: &[Digest32]) -> Option<Vec<Digest32>> {
    if sig.len() != LEN {
        return None;
    }
    let digits = digits_for(pp, t, msg);
    Some(
        (0..LEN)
            .map(|i| chain(pp, t, i as u32, digits[i], W - 1 - digits[i], &sig[i]))
            .collect(),
    )
}

/// `WVerify(pk, msg, sigma, t)`.
pub fn verify(pp: &PubParams, pk: &[Digest32], t: &WotsTweak, msg: &[u8], sig: &[Digest32]) -> bool {
    derive_pk(pp, t, msg, sig).as_deref() == Some(pk)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, RngCore};

    fn setup() -> (PubParams, [u8; N]) {
        let mut rng = thread_rng();
        let mut id = [0u8; N];
        let mut sk_seed = [0u8; N];
        rng.fill_bytes(&mut id);
        rng.fill_bytes(&mut sk_seed);
        (PubParams::new(&id), sk_seed)
    }

    #[test]
    fn sign_verify() {
        let (pp, sk_seed) = setup();
        let t = WotsTweak { r: 2, idx: 12345 };
        let kp = key_gen(&pp, &sk_seed, &t);
        let msg = b"upper layer message";
        let sig = sign(&pp, &sk_seed, &t, msg);
        assert!(verify(&pp, &kp.pk, &t, msg, &sig));
        assert!(!verify(&pp, &kp.pk, &t, b"other message", &sig));
    }

    #[test]
    fn wrong_tweak_or_parameter_fails() {
        let (pp, sk_seed) = setup();
        let t = WotsTweak { r: 1, idx: 7 };
        let kp = key_gen(&pp, &sk_seed, &t);
        let msg = b"m";
        let sig = sign(&pp, &sk_seed, &t, msg);
        assert!(!verify(&pp, &kp.pk, &WotsTweak { r: 1, idx: 8 }, msg, &sig));
        assert!(!verify(&pp, &kp.pk, &WotsTweak { r: 2, idx: 7 }, msg, &sig));
        let other_pp = PubParams::new(&[1u8; N]);
        assert!(!verify(&other_pp, &kp.pk, &t, msg, &sig));
    }

    /// `BKeyGen` derives the top layer's WOTS+ keys before the long-term
    /// root exists, and `PubParams::bind_root` then fixes `pk = (Root, id)`.
    /// That must not disturb any WOTS+ key already derived, or the root
    /// would no longer be the root of the tree that produced it.
    #[test]
    fn keys_are_independent_of_root_binding() {
        let (pp, sk_seed) = setup();
        let t = WotsTweak { r: 3, idx: 0 };
        let before = key_gen(&pp, &sk_seed, &t);
        let mut bound = pp.clone();
        bound.bind_root(&[0xABu8; N]);
        let after = key_gen(&bound, &sk_seed, &t);
        assert_eq!(before.pk, after.pk);
        let msg = b"m_0";
        assert_eq!(sign(&pp, &sk_seed, &t, msg), sign(&bound, &sk_seed, &t, msg));
    }

    #[test]
    fn chain_composes_across_heights() {
        // H^{a+b}(t||i||0, x) = H^b(t||i||a, H^a(t||i||0, x)): the property
        // that makes WDerivePK recompute y_i.
        let (pp, _) = setup();
        let t = WotsTweak { r: 1, idx: 0 };
        let x = [5u8; N];
        let full = chain(&pp, &t, 3, 0, W - 1, &x);
        for a in 0..W {
            let mid = chain(&pp, &t, 3, 0, a, &x);
            assert_eq!(chain(&pp, &t, 3, a, W - 1 - a, &mid), full);
        }
    }

    #[test]
    fn checksum_digits_match_paper_definition() {
        let digits = vec![0u32; LEN1]; // maximal checksum 64*15 = 960
        let c = checksum(&digits);
        let value = c.iter().fold(0u32, |acc, d| acc * W + d);
        assert_eq!(value, 960);
    }
}
