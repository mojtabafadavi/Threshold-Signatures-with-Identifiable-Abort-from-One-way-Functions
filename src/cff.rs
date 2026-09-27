//! kappa-uniform 1-cover-free family (e, lambda, kappa): Definition 2 and
//! the CFF-based LOTS of Section 2.1.
//!
//! We realize the bijection g_lambda: {0,1}^lambda -> B (B a family of
//! kappa-subsets of [e] pairwise non-containing) via the standard
//! combinatorial-number-system ("colex") bijection between {0,...,C(e,kappa)-1}
//! and kappa-subsets of [e]. This family is automatically 1-cover-free: two
//! distinct combinations can never contain one another because they are
//! compared elementwise from the largest element down.
//!
//! As long as 2^lambda <= C(e, kappa) the map from lambda-bit digests to
//! subsets is injective, which is all the 1-CFF property requires.

use num_bigint::BigUint;
use num_traits::{One, Zero};

pub struct Cff {
    pub e: usize,
    pub kappa: usize,
    /// binom[n][k] = C(n,k) for n in 0..=e, k in 0..=kappa
    binom: Vec<Vec<BigUint>>,
}

impl Cff {
    pub fn new(e: usize, kappa: usize) -> Self {
        let mut binom = vec![vec![BigUint::zero(); kappa + 1]; e + 1];
        for n in 0..=e {
            binom[n][0] = BigUint::one();
            for k in 1..=kappa.min(n) {
                if k == n {
                    binom[n][k] = BigUint::one();
                } else {
                    binom[n][k] = &binom[n - 1][k - 1] + &binom[n - 1][k];
                }
            }
        }
        Cff { e, kappa, binom }
    }

    pub fn capacity(&self) -> &BigUint {
        &self.binom[self.e][self.kappa]
    }

    fn c(&self, n: usize, k: usize) -> BigUint {
        if k > n {
            BigUint::zero()
        } else {
            self.binom[n][k].clone()
        }
    }

    /// g_lambda(digest): map a digest (interpreted as a big-endian integer,
    /// reduced mod C(e,kappa)) to its kappa-subset of {0,...,e-1}.
    pub fn digest_to_block(&self, digest: &[u8]) -> Vec<usize> {
        let cap = self.capacity();
        let mut idx = BigUint::from_bytes_be(digest) % cap;
        let mut block = Vec::with_capacity(self.kappa);
        let mut k = self.kappa;
        // Standard combinatorial-number-system decoding: find, from the top
        // element down, the largest n such that C(n,k) <= idx.
        let mut n = self.e;
        while k > 0 {
            // find greatest a in [k-1, n-1] with C(a,k) <= idx
            let mut a = n - 1;
            loop {
                let cak = self.c(a, k);
                if cak <= idx {
                    break;
                }
                if a == k - 1 {
                    break;
                }
                a -= 1;
            }
            idx -= self.c(a, k);
            block.push(a);
            n = a;
            k -= 1;
        }
        block.sort_unstable();
        block
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn injective_small() {
        let cff = Cff::new(16, 6);
        let mut seen = std::collections::HashSet::new();
        for i in 0u32..500 {
            let bytes = i.to_be_bytes();
            let block = cff.digest_to_block(&bytes);
            assert_eq!(block.len(), 6);
            assert!(seen.insert(block));
        }
    }

    /// Section 8's parameter choice: `e = 261` is the smallest `e` with
    /// `C(e, floor(e/2)) >= 2^256`, and `kappa = 123` the smallest `kappa`
    /// with `C(261, kappa) >= 2^256`. Since `2^256 <= C(e, kappa)`, the
    /// reduction in `digest_to_block` never reduces a 256-bit
    /// digest, so `g_lambda` is injective as Definition 2 requires.
    #[test]
    fn paper_parameters_are_minimal_and_injective() {
        let two256 = BigUint::one() << 256u32;
        assert!(*Cff::new(261, 123).capacity() >= two256, "2^256 <= C(261,123)");
        assert!(*Cff::new(261, 122).capacity() < two256, "kappa = 123 is the smallest");
        assert!(*Cff::new(260, 130).capacity() < two256, "e = 261 is the smallest");
    }

    #[test]
    fn cover_free() {
        // Any two distinct combination-indices produce sets that don't contain
        // one another, since the combinadic decoding is a bijection onto ALL
        // kappa-subsets and any two distinct kappa-subsets of equal size can
        // never be subsets of each other.
        let cff = Cff::new(10, 4);
        let a = cff.digest_to_block(&7u32.to_be_bytes());
        let b = cff.digest_to_block(&8u32.to_be_bytes());
        assert_ne!(a, b);
    }
}
