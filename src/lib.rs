//! `tlots`: a Rust implementation of the multi-message threshold signature
//! scheme of "Threshold Signatures with Identifiable Abort from One-way
//! Functions", combining
//!   (a) the threshold CFF-based Lamport one-time signature of Sections 2.1
//!       and 3 of that paper, used at the bottom (leaf) layer, together
//!       with the verification keys of Algorithm 4 and the opening check of
//!       Section 5 that make an abort identifiable, and
//!   (b) a WOTS+ certification hypertree (Section 2.3, Algorithm 1) for the
//!       upper certification layers.
//!
//! All hashing goes through one tweakable hash `H(P, t, x)` (Appendix A.1)
//! or the secret-seed PRF: WOTS+, the CFF-Lamport OTS, the verification
//! keys and the Merkle trees each take public parameters and tweaks exactly
//! as in Sections 2.1-2.4, Section 5's Equation (1), Algorithm 4 and
//! Appendix A.3 (see `hashing`).
//!
//! The public parameters of the LOTS layer are keyed on the whole long-term
//! public key `pk = (Root, id)`; those of the hypertree, which is what
//! produces `Root`, are keyed on `id` alone. See `hashing` and `README.md`.

pub mod cff;
pub mod hashing;
pub mod hypertree;
pub mod lots_cff;
pub mod merkle;
pub mod shamir;
pub mod wotsplus;
