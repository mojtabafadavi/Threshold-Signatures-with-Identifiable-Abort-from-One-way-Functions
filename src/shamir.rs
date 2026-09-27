//! Shamir (T,N)-threshold secret sharing over the prime field F_p
//! (Section 2.4, Algorithm 2). Each 32-byte CFF-Lamport secret value x_i
//! is parsed into `r = ceil(256 / bpb)` blocks of `bpb` bits each so that
//! every block fits into F_p (p > 2^bpb), exactly as Section 8 does.
//!
//! `share_ntt` uses a number-theoretic transform (NTT) in place of
//! per-point polynomial evaluation, bringing Dist's total cost from
//! O(N*T) to O(N log N): as Section 8 specifies, participant j is
//! assigned the point omega^j for a primitive root of unity omega of
//! order M = next_pow2(N) rather than the plain integer j+1 -- so
//! "evaluate the secret's degree-(T-1) polynomial at every participant's
//! point" becomes exactly "evaluate at every M-th root of unity", which
//! is what a forward NTT computes in O(M log M).
//!
//! Reconstruction interpolates at x = 0 by Lagrange. The coefficients
//! `lambda_{i,S}` depend only on the signer set S and not on which
//! coordinate is being recovered, so `lagrange_coeffs` computes them once
//! per quorum and `reconstruct_with` reuses them across all coordinates of
//! a shared vector, which is the Theta(T^2 + kappa*r*T) cost of Section 8
//! rather than Theta(kappa*r*T^2). The inverse NTT does not help here:
//! an arbitrary quorum of T responders is not a root-of-unity subset.
//! Subproduct-tree interpolation would give O(T log^2 T), and is not
//! implemented.

use rand::RngCore;

/// Section 8's FFT-friendly prime, just above 2^32: p = 4101*2^20 + 1.
/// Every Shamir block is exactly BPB = 32 bits (0 <= block < 2^32), and p
/// exceeds 2^32 by only 5,242,881 (~0.12%), so every element of F_p still
/// needs only ceil(log2(p)/8) = 5 bytes to serialize. Because
/// p-1 = 4101*2^20 has a large power-of-two factor, F_p has a
/// multiplicative subgroup of order 2^20, which `share_ntt` uses.
pub const P: u64 = 4_300_210_177;
pub const BPB: u32 = 32; // bits per block
pub const R: usize = (256 + BPB as usize - 1) / BPB as usize; // blocks per 32-byte secret
/// Minimal number of bytes needed to serialize one F_p element.
pub const ELEM_BYTES: usize = 5;

/// The order of the fixed primitive root `ROOT` below: F_p's multiplicative
/// group has order p-1 = 4101*2^20, so 2^20 is the largest power of two
/// dividing it, and therefore the largest NTT transform size supported --
/// comfortably above any (T,N) this crate benchmarks (N <= 1000).
pub const ROOT_ORDER: u64 = 1 << 20;
/// A primitive `ROOT_ORDER`-th root of unity mod P, i.e. `g^((P-1)/ROOT_ORDER)`
/// for a generator g=5 of F_p^*, precomputed offline (verified: `ROOT^ROOT_ORDER
/// == 1` and `ROOT^(ROOT_ORDER/2) != 1`, so its order is exactly `ROOT_ORDER`,
/// not a smaller divisor).
pub const ROOT: u64 = 307_724_279;

#[inline]
fn mulmod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}
#[inline]
fn addmod(a: u64, b: u64, p: u64) -> u64 {
    let s = a as u128 + b as u128;
    (s % p as u128) as u64
}
#[inline]
fn submod(a: u64, b: u64, p: u64) -> u64 {
    let a = a as i128;
    let b = b as i128;
    let p = p as i128;
    (((a - b) % p + p) % p) as u64
}

pub fn modinv(a: u64, p: u64) -> u64 {
    // extended Euclid, a in [0,p)
    let (mut old_r, mut r) = (a as i128, p as i128);
    let (mut old_s, mut s) = (1i128, 0i128);
    while r != 0 {
        let q = old_r / r;
        let tmp_r = old_r - q * r;
        old_r = r;
        r = tmp_r;
        let tmp_s = old_s - q * s;
        old_s = s;
        s = tmp_s;
    }
    (((old_s % p as i128) + p as i128) % p as i128) as u64
}

#[inline]
fn modexp(mut base: u64, mut exp: u64, p: u64) -> u64 {
    let mut acc = 1u64;
    base %= p;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mulmod(acc, base, p);
        }
        base = mulmod(base, base, p);
        exp >>= 1;
    }
    acc
}

/// A primitive `m`-th root of unity mod P, for any power-of-two `m`
/// dividing `ROOT_ORDER`: `ROOT^(ROOT_ORDER/m) mod P`.
fn root_of_order(m: u64) -> u64 {
    assert!(m.is_power_of_two() && ROOT_ORDER % m == 0, "unsupported NTT size");
    modexp(ROOT, ROOT_ORDER / m, P)
}

/// Split a 32-byte secret into R field-blocks of BPB bits each.
pub fn secret_to_blocks(x: &[u8; 32]) -> [u64; R] {
    let mut blocks = [0u64; R];
    let bits = 256usize;
    for i in 0..R {
        let bit_start = i * BPB as usize;
        let bit_end = (bit_start + BPB as usize).min(bits);
        let mut v: u64 = 0;
        for bit in bit_start..bit_end {
            let byte_idx = bit / 8;
            let bit_idx = 7 - (bit % 8);
            let b = (x[byte_idx] >> bit_idx) & 1;
            v = (v << 1) | (b as u64);
        }
        blocks[i] = v;
    }
    blocks
}

/// Serialized size of one participant's share of one 32-byte secret: `R`
/// field elements at `ELEM_BYTES` bytes each.
pub const SHARE_BYTES: usize = R * ELEM_BYTES;

/// Canonical serialization of one participant's share of one secret:
/// each of the `R` field elements big-endian in `ELEM_BYTES` bytes (every
/// element is `< P < 2^40`, so nothing is truncated). Used as the message
/// of the verification-key commitments `pi_{ijk} = H_vk(i||j||k, sh_{ij,k})`
/// (Algorithm 4), which need a fixed-width encoding of the share.
pub fn share_to_bytes(blocks: &[u64; R]) -> [u8; SHARE_BYTES] {
    let mut out = [0u8; SHARE_BYTES];
    for (i, b) in blocks.iter().enumerate() {
        debug_assert!(*b < P, "share is not a reduced field element");
        let be = b.to_be_bytes();
        out[i * ELEM_BYTES..(i + 1) * ELEM_BYTES].copy_from_slice(&be[8 - ELEM_BYTES..]);
    }
    out
}

pub fn blocks_to_secret(blocks: &[u64; R]) -> [u8; 32] {
    let mut out = [0u8; 32];
    let bits = 256usize;
    for i in 0..R {
        let bit_start = i * BPB as usize;
        let bit_end = (bit_start + BPB as usize).min(bits);
        let width = bit_end - bit_start;
        let v = blocks[i];
        for j in 0..width {
            let bit = bit_start + j;
            let byte_idx = bit / 8;
            let bit_idx = 7 - (bit % 8);
            let b = ((v >> (width - 1 - j)) & 1) as u8;
            out[byte_idx] |= b << bit_idx;
        }
    }
    out
}

/// In-place iterative (Cooley-Tukey, decimation-in-time) NTT of `a`, whose
/// length must be a power of two dividing `ROOT_ORDER`. `root` must be a
/// primitive `a.len()`-th root of unity mod P (e.g. `root_of_order(a.len())`
/// for a forward transform, or its modular inverse for the reverse one).
fn ntt(a: &mut [u64], root: u64) {
    let n = a.len();
    debug_assert!(n.is_power_of_two());

    // bit-reversal permutation
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            a.swap(i, j);
        }
    }

    let mut len = 2usize;
    while len <= n {
        let w_len = modexp(root, (n / len) as u64, P);
        let mut i = 0usize;
        while i < n {
            let mut w = 1u64;
            for k in 0..len / 2 {
                let u = a[i + k];
                let v = mulmod(a[i + k + len / 2], w, P);
                a[i + k] = addmod(u, v, P);
                a[i + k + len / 2] = submod(u, v, P);
                w = mulmod(w, w_len, P);
            }
            i += len;
        }
        len <<= 1;
    }
}

/// The evaluation domain `(omega^0, omega^1, ..., omega^{n-1})` that
/// `share_ntt` assigns to participants `1..=n`, for a primitive
/// `next_pow2(n)`-th root of unity `omega`. Callers needing a given
/// participant's own field point (e.g. to reconstruct with `reconstruct`)
/// should index into this, `domain[j]` for the participant at position `j`.
pub fn eval_domain(n: usize) -> Vec<u64> {
    let m = n.max(1).next_power_of_two() as u64;
    let omega = root_of_order(m);
    let mut domain = Vec::with_capacity(n);
    let mut cur = 1u64;
    for _ in 0..n {
        domain.push(cur);
        cur = mulmod(cur, omega, P);
    }
    domain
}

/// Sample a degree-(T-1) polynomial with constant term `secret`, and
/// evaluate it at every point of `eval_domain(n)` via a forward NTT --
/// O(M log M) for M = next_pow2(n), instead of `share_naive`'s O(n*t).
/// Returns the n shares, in the same order as `eval_domain(n)`.
pub fn share_ntt(secret: u64, t: usize, n: usize, rng: &mut impl RngCore) -> Vec<u64> {
    let m = n.max(1).next_power_of_two();
    let mut coeffs = vec![0u64; m];
    coeffs[0] = secret;
    for c in coeffs.iter_mut().take(t).skip(1) {
        *c = rng.next_u64() % P;
    }
    let omega = root_of_order(m as u64);
    ntt(&mut coeffs, omega);
    coeffs.truncate(n);
    coeffs
}

/// Sample a degree (T-1) polynomial with constant term `secret` and evaluate
/// it at points 1..=N one at a time via Horner's method, returning the N
/// shares f(1), ..., f(N), in O(N*T). `lots_cff::dist` uses `share_ntt`;
/// this is the reference against which `ntt_bench` compares it.
pub fn share_naive(secret: u64, t: usize, n: usize, rng: &mut impl RngCore) -> Vec<u64> {
    let mut coeffs = vec![secret];
    for _ in 1..t {
        coeffs.push(rng.next_u64() % P);
    }
    (1..=n as u64)
        .map(|x| eval_poly(&coeffs, x))
        .collect()
}

fn eval_poly(coeffs: &[u64], x: u64) -> u64 {
    // Horner's method
    let mut acc = 0u64;
    for c in coeffs.iter().rev() {
        acc = addmod(mulmod(acc, x, P), *c, P);
    }
    acc
}

/// The Lagrange coefficients at x = 0 for the quorum `xs`, in the same
/// order: `lambda_{i,S} = prod_{m in S \ {i}} x_m / (x_m - x_i)`
/// (Section 2.4). They are a function of the signer set alone, so a caller
/// recovering many coordinates of one shared vector -- Algorithm 2's
/// `ShamirReconVec` -- computes them once here and passes them to
/// `reconstruct_with` for each coordinate.
///
/// `xs` must be distinct and nonzero; `has_distinct_points` checks that.
pub fn lagrange_coeffs(xs: &[u64]) -> Vec<u64> {
    xs.iter()
        .map(|&xi| {
            let mut num = 1u64;
            let mut den = 1u64;
            for &xj in xs {
                if xj == xi {
                    continue;
                }
                num = mulmod(num, xj, P);
                den = mulmod(den, submod(xj, xi, P), P);
            }
            mulmod(num, modinv(den, P), P)
        })
        .collect()
}

/// Whether `xs` is a usable evaluation set: pairwise distinct and free of
/// the secret's own point 0. A quorum that repeats a participant would
/// otherwise interpolate through a degenerate coefficient set.
pub fn has_distinct_points(xs: &[u64]) -> bool {
    xs.iter().enumerate().all(|(a, &x)| x != 0 && !xs[..a].contains(&x))
}

/// Interpolate one coordinate at x = 0 from coefficients produced by
/// `lagrange_coeffs` and the shares of the same quorum, in the same order:
/// `sum_i lambda_{i,S} y_i`.
pub fn reconstruct_with(coeffs: &[u64], ys: &[u64]) -> u64 {
    assert_eq!(coeffs.len(), ys.len(), "coefficient and share counts differ");
    coeffs
        .iter()
        .zip(ys)
        .fold(0u64, |acc, (&lag, &y)| addmod(acc, mulmod(y, lag, P), P))
}

/// Lagrange-interpolate at x=0 given a set of (point, share) pairs. Generic
/// over whatever x-coordinates it is given -- unaffected by whether those
/// points came from `share_naive` (plain integers) or `share_ntt` (powers
/// of a root of unity). Recovering a whole shared vector one coordinate at
/// a time this way recomputes the coefficients on every call; prefer
/// `lagrange_coeffs` with `reconstruct_with` for that.
pub fn reconstruct(points: &[(u64, u64)]) -> u64 {
    let xs: Vec<u64> = points.iter().map(|&(x, _)| x).collect();
    let ys: Vec<u64> = points.iter().map(|&(_, y)| y).collect();
    reconstruct_with(&lagrange_coeffs(&xs), &ys)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn roundtrip_blocks() {
        let mut x = [0u8; 32];
        for i in 0..32 {
            x[i] = (i as u8).wrapping_mul(37).wrapping_add(5);
        }
        let blocks = secret_to_blocks(&x);
        let back = blocks_to_secret(&blocks);
        assert_eq!(x, back);
    }

    #[test]
    fn share_serialization_is_injective_and_fixed_width() {
        let mut a = [0u64; R];
        a[0] = P - 1;
        a[R - 1] = 1;
        let mut b = a;
        b[R - 1] = 2;
        assert_eq!(share_to_bytes(&a).len(), SHARE_BYTES);
        assert_ne!(share_to_bytes(&a), share_to_bytes(&b));
        // The largest element still fits in ELEM_BYTES bytes.
        assert!((P - 1) < 1u64 << (8 * ELEM_BYTES as u32));
    }

    #[test]
    fn shamir_roundtrip_naive() {
        let mut rng = thread_rng();
        let secret = 123456789u64 % P;
        let (t, n) = (5, 10);
        let shares = share_naive(secret, t, n, &mut rng);
        let pts: Vec<(u64, u64)> = (1..=t as u64).map(|i| (i, shares[i as usize - 1])).collect();
        assert_eq!(reconstruct(&pts), secret);
        // a different quorum of size t should agree
        let pts2: Vec<(u64, u64)> = (n - t + 1..=n)
            .map(|i| (i as u64, shares[i - 1]))
            .collect();
        assert_eq!(reconstruct(&pts2), secret);
    }

    #[test]
    fn root_has_exact_order() {
        assert_eq!(modexp(ROOT, ROOT_ORDER, P), 1, "ROOT^ROOT_ORDER must be 1");
        assert_ne!(
            modexp(ROOT, ROOT_ORDER / 2, P),
            1,
            "ROOT's order must be exactly ROOT_ORDER, not a smaller divisor"
        );
    }

    #[test]
    fn ntt_matches_naive_evaluation() {
        // share_ntt's y-values must equal direct (Horner) evaluation of the
        // same coefficient vector at the same domain points, which is the
        // correctness property required of an NTT.
        let mut rng = thread_rng();
        for &(t, n) in &[(1usize, 1usize), (3, 5), (5, 10), (7, 100), (50, 257)] {
            let secret = (n as u64 * 97 + 13) % P;
            let m = n.max(1).next_power_of_two();
            let mut coeffs = vec![0u64; m];
            coeffs[0] = secret;
            // The identity is checked on a fresh random coefficient vector
            // rather than on share_ntt's internal randomness.
            for c in coeffs.iter_mut().take(t).skip(1) {
                *c = rng.next_u64() % P;
            }
            let domain = eval_domain(n);
            let expected: Vec<u64> = domain.iter().map(|&x| eval_poly(&coeffs, x)).collect();

            let omega = root_of_order(m as u64);
            let mut transformed = coeffs.clone();
            ntt(&mut transformed, omega);
            transformed.truncate(n);

            assert_eq!(transformed, expected, "NTT must match direct evaluation at t={t}, n={n}");
        }
    }

    #[test]
    fn shamir_roundtrip_ntt() {
        let mut rng = thread_rng();
        let secret = 987654321u64 % P;
        let (t, n) = (5, 10);
        let domain = eval_domain(n);
        let shares = share_ntt(secret, t, n, &mut rng);
        assert_eq!(shares.len(), n);

        let pts: Vec<(u64, u64)> = (0..t).map(|j| (domain[j], shares[j])).collect();
        assert_eq!(reconstruct(&pts), secret);

        // a different quorum of size t should also reconstruct correctly
        let pts2: Vec<(u64, u64)> = (n - t..n).map(|j| (domain[j], shares[j])).collect();
        assert_eq!(reconstruct(&pts2), secret);
    }

    /// Reusing one coefficient set across the coordinates of a shared vector
    /// must agree with interpolating each coordinate from scratch, which is
    /// what lets `lots_cff::combine` pay for the coefficients once.
    #[test]
    fn shared_coefficients_match_per_coordinate_interpolation() {
        let mut rng = thread_rng();
        let (t, n) = (6, 11);
        let domain = eval_domain(n);
        let secrets: Vec<u64> = (0..R as u64 * 3).map(|i| (i * 7919 + 11) % P).collect();
        let columns: Vec<Vec<u64>> = secrets
            .iter()
            .map(|&s| share_ntt(s, t, n, &mut rng))
            .collect();

        let xs: Vec<u64> = domain[..t].to_vec();
        assert!(has_distinct_points(&xs));
        let coeffs = lagrange_coeffs(&xs);
        for (coord, shares) in columns.iter().enumerate() {
            let ys: Vec<u64> = shares[..t].to_vec();
            let pts: Vec<(u64, u64)> = xs.iter().copied().zip(ys.iter().copied()).collect();
            assert_eq!(reconstruct_with(&coeffs, &ys), reconstruct(&pts));
            assert_eq!(reconstruct_with(&coeffs, &ys), secrets[coord]);
        }
    }

    #[test]
    fn degenerate_quorums_are_detected() {
        let domain = eval_domain(8);
        assert!(has_distinct_points(&domain[..4]));
        assert!(!has_distinct_points(&[domain[1], domain[2], domain[1]]));
        assert!(!has_distinct_points(&[0, domain[1]]));
    }

    #[test]
    fn ntt_handles_non_power_of_two_n() {
        // n need not itself be a power of two; only the internal transform
        // size next_pow2(n) is.
        let mut rng = thread_rng();
        let secret = 42u64;
        let (t, n) = (4, 7);
        let domain = eval_domain(n);
        let shares = share_ntt(secret, t, n, &mut rng);
        let pts: Vec<(u64, u64)> = (0..t).map(|j| (domain[j], shares[j])).collect();
        assert_eq!(reconstruct(&pts), secret);
    }
}
