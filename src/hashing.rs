//! Hash instantiations, following Appendix A.1 (tweakable hash functions),
//! Section 5's public-parameter convention and Algorithm 4, of "Threshold
//! Signatures with Identifiable Abort from One-way Functions".
//!
//! Every hash in the scheme is either the tweakable hash
//!
//! ```text
//!     H(P, t, x) = SHA-256("TLOTS/TH" || len(P) || P || len(t) || t || x)
//! ```
//!
//! (`thash`, Definition 3) or the secret-seed PRF (`prf`). Length prefixes
//! make the encoding of `(P, t, x)` injective, and every tweak produced by
//! the `tweak` module has fixed-width fields behind a type byte, so two
//! different uses can never share a `(P, t)` pair (the `DIST` condition of
//! Definitions 4 and 5).
//!
//! # Public parameters
//!
//! Section 5, Equation (1) fixes the convention
//!
//! ```text
//!     H_X(t, m) = H(TLOTS.pk || X, t, m)
//! ```
//!
//! for a string literal `X`, and Algorithm 4 (`TLOTS.KeyDist.KeyRefresh`)
//! instantiates the CFF-Lamport OTS with `prm = (pk^B || "LOTS",
//! pk^B || "msg")`. The long-term public key is the pair
//! `pk^B = (Root, id)` (Algorithm 1, `BKeyGen` line 7), serialized here as
//! `Root || id`. `PubParams` therefore keys
//!
//! * `P_LOTS  = pk^B || "LOTS"`   (LOTS key components, Algorithm 4 line 2)
//! * `P_msg   = pk^B || "msg"`    (LOTS message hash,  Algorithm 4 line 2)
//! * `P_vk    = pk^B || "vk"`     (`H_vk`, Algorithm 4 lines 3-4)
//!
//! on the whole public key. The hypertree's own parameters cannot be keyed
//! that way without circularity -- `Root` is by definition the output of
//! the top Merkle tree over WOTS+ public keys, so anything hashed while
//! computing `Root` can only depend on `id`. Algorithm 1 indeed keys those
//! on `id` alone (`MerkleTree(id || r || tau, ...)`), as does this crate:
//!
//! * `P_WOTS  = id || "WOTS"`     (WOTS+ chains)
//! * `P_Wmsg  = id || "WMSG"`     (WOTS+ message hash)
//! * Merkle trees: `id || r || tau` (upper) and `id_msg` (salted bottom).
//!
//! The paper calls both message-hash parameters `P_msg`; keeping the WOTS+
//! one under its own label is what makes them distinct here, since the two
//! are keyed on different values (`id` vs `pk^B`) only once `Root` exists.
//!
//! `PubParams::new` builds the `id`-keyed parameters, which is all that
//! `BKeyGen` needs; `bind_root`/`with_pk` add `Root` and so complete the
//! LOTS and `H_vk` parameters. Before the root is bound, the LOTS
//! parameters are keyed on `id` alone: the degenerate `d = 1` shape, where
//! the LOTS tree is itself the top tree and `Root` is likewise undefined
//! when the LOTS keys are generated.
//!
//! # Tweak namespaces
//!
//! The tweaks fall into three groups that never share a public parameter:
//!
//! * one-time-signature and verification-key tweaks (`0x01..=0x09`), used
//!   with `P_LOTS`, `P_msg`, `P_WOTS`, `P_Wmsg` and `P_vk`;
//! * Merkle tweaks, used only with a tree parameter, and written exactly
//!   as Appendix A.3 writes them: `0 || i` for a leaf and `1 || j || i`
//!   for an internal node at level `j` counted from the root;
//! * PRF tweaks (`0x05`, `0x11..=0x13`), used with the dealer's secret
//!   seed `rho` as the key rather than with any public parameter. WOTS+
//!   secret values and the hypertree's seed/`id_msg`/salt derivations all
//!   use that one key, so their tags are distinct from each other.
//!
//! Since the public parameters above are pairwise distinct byte strings
//! (and `thash` length-prefixes `P`), a Merkle tweak and an OTS tweak can
//! never collide even though `1 || j || i` and the LOTS-key tweak happen
//! to share a leading byte.

use sha2::{Digest, Sha256};

pub const N: usize = 32; // output size in bytes (lambda = 256)

pub type Digest32 = [u8; N];

fn sha256(parts: &[&[u8]]) -> Digest32 {
    let mut h = Sha256::new();
    for p in parts {
        h.update(p);
    }
    let out = h.finalize();
    let mut d = [0u8; N];
    d.copy_from_slice(&out);
    d
}

/// Tweakable hash H(P, t, x) of Definition 3.
pub fn thash(p: &[u8], tweak: &[u8], x: &[u8]) -> Digest32 {
    assert!(p.len() <= u16::MAX as usize, "public parameter too long");
    assert!(tweak.len() <= u8::MAX as usize, "tweak too long");
    sha256(&[
        b"TLOTS/TH",
        &(p.len() as u16).to_be_bytes(),
        p,
        &[tweak.len() as u8],
        tweak,
        x,
    ])
}

/// Keyed PRF used to derive secret values pseudorandomly from a secret
/// seed and a (public) tweak, e.g. `x_i = PRF(rho, t || i)` for WOTS+
/// (Algorithm 1's `WKeyGen(rho, t)`).
pub fn prf(key: &[u8; N], tweak: &[u8]) -> Digest32 {
    sha256(&[b"TLOTS/PRF", key, &[tweak.len() as u8], tweak])
}

const LBL_LOTS: &[u8] = b"LOTS";
const LBL_MSG: &[u8] = b"msg";
const LBL_VK: &[u8] = b"vk";
const LBL_WOTS: &[u8] = b"WOTS";
const LBL_WMSG: &[u8] = b"WMSG";

/// The public parameters of Sections 2.1 and 5: `P_LOTS`, `P_msg`, `P_vk`
/// (keyed on the long-term public key `pk^B = (Root, id)`), and `P_WOTS`,
/// `P_Wmsg` plus the Merkle-tree parameters (keyed on `id`). See the module
/// documentation for why the two groups are keyed differently.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PubParams {
    /// The public key's `id` (Algorithm 1's random public parameter).
    pub id: Digest32,
    /// The public key's `Root`, once it exists.
    root: Option<Digest32>,
    lots: Vec<u8>,
    lots_msg: Vec<u8>,
    vk: Vec<u8>,
    wots: Vec<u8>,
    wots_msg: Vec<u8>,
}

impl PubParams {
    /// The parameters available from `id` alone. Enough for the whole
    /// hypertree (`BKeyGen`, WOTS+, every Merkle tree); the LOTS and `H_vk`
    /// parameters are keyed on `id` until [`PubParams::bind_root`] supplies
    /// the long-term root.
    pub fn new(id: &Digest32) -> Self {
        Self::build(id, None)
    }

    /// The complete parameters of a hypertree whose public key is
    /// `pk^B = (root, id)`, as a verifier holds it.
    pub fn with_pk(root: &Digest32, id: &Digest32) -> Self {
        Self::build(id, Some(*root))
    }

    /// Fix `pk^B = (root, id)` once `BKeyGen` has produced the root,
    /// re-keying `P_LOTS`, `P_msg` and `P_vk` onto it (Algorithm 4 line 2).
    /// The WOTS+ and Merkle parameters are unaffected, so anything already
    /// derived from them stays valid.
    pub fn bind_root(&mut self, root: &Digest32) {
        *self = Self::build(&self.id, Some(*root));
    }

    fn build(id: &Digest32, root: Option<Digest32>) -> Self {
        // `pk^B = (Root, id)`, serialized as `Root || id`; before the root
        // is known the pk-keyed parameters fall back to `id`.
        let key: Vec<u8> = match &root {
            Some(r) => [r.as_slice(), id.as_slice()].concat(),
            None => id.to_vec(),
        };
        let p = |k: &[u8], label: &[u8]| [k, label].concat();
        PubParams {
            id: *id,
            root,
            lots: p(&key, LBL_LOTS),
            lots_msg: p(&key, LBL_MSG),
            vk: p(&key, LBL_VK),
            wots: p(id.as_slice(), LBL_WOTS),
            wots_msg: p(id.as_slice(), LBL_WMSG),
        }
    }

    /// `Root`, if these parameters have been bound to a public key.
    pub fn root(&self) -> Option<&Digest32> {
        self.root.as_ref()
    }
    /// `P_LOTS` (Algorithm 4: `pk^B || "LOTS"`).
    pub fn lots(&self) -> &[u8] {
        &self.lots
    }
    /// `P_msg` for the CFF-Lamport message hash (`pk^B || "msg"`).
    pub fn lots_msg(&self) -> &[u8] {
        &self.lots_msg
    }
    /// `P_vk` for `H_vk` (`pk^B || "vk"`), Algorithm 4 lines 3-4.
    pub fn vk(&self) -> &[u8] {
        &self.vk
    }
    /// `P_WOTS` (`id || "WOTS"`).
    pub fn wots(&self) -> &[u8] {
        &self.wots
    }
    /// The WOTS+ message-hash parameter (`id || "WMSG"`).
    pub fn wots_msg(&self) -> &[u8] {
        &self.wots_msg
    }

    /// Public parameter of an upper-layer (WOTS+) Merkle tree:
    /// `id || r || tau`, with `r` the paper's layer number (1 = the WOTS+
    /// layer directly above the LOTS trees) and `tau` the tree's index
    /// within that layer (Algorithm 1's `id || r || floor(index/L)`).
    pub fn merkle_upper(&self, r: u32, tau: u64) -> Vec<u8> {
        let mut v = Vec::with_capacity(N + 12);
        v.extend_from_slice(&self.id);
        v.extend_from_slice(&r.to_be_bytes());
        v.extend_from_slice(&tau.to_be_bytes());
        v
    }
}

/// Fixed-width tweak encodings (big-endian). Merkle tweaks are Appendix
/// A.3's `0 || i` and `1 || j || i` verbatim; all other tweaks start with a
/// distinct type tag. See the module documentation for the three tweak
/// namespaces.
pub mod tweak {
    // One-time-signature and verification-key tweaks (used with a public
    // parameter of `PubParams`).
    const T_LOTS_KEY: u8 = 0x01;
    const T_LOTS_MSG: u8 = 0x02;
    const T_WOTS_CHAIN: u8 = 0x03;
    const T_WOTS_MSG: u8 = 0x04;
    const T_VK_SHARE: u8 = 0x08;
    const T_VK_PARTY: u8 = 0x09;

    // PRF tweaks (used with the dealer's secret seed as the key).
    const T_WOTS_SK: u8 = 0x05;
    const T_SEED_LEAF: u8 = 0x11;
    const T_SEED_IDMSG: u8 = 0x12;
    const T_SEED_SALT: u8 = 0x13;

    // Merkle tweaks (used with a tree's own public parameter), exactly as
    // Appendix A.3 writes them.
    const T_MERKLE_LEAF: u8 = 0x00;
    const T_MERKLE_NODE: u8 = 0x01;

    /// Tweak of one CFF-Lamport key: `t` = the key's global index.
    /// (A WOTS+ key's tweak is `WotsTweak` below.)
    pub type LotsTweak = u64;

    /// `t || i` for the i-th component of a LOTS key (`y_i = H(P_LOTS, t||i, x_i)`).
    pub fn lots_key(t: LotsTweak, i: u32) -> [u8; 13] {
        let mut b = [0u8; 13];
        b[0] = T_LOTS_KEY;
        b[1..9].copy_from_slice(&t.to_be_bytes());
        b[9..13].copy_from_slice(&i.to_be_bytes());
        b
    }

    /// `t` for the LOTS message hash (`B_msg = g(H(P_msg, t, msg))`).
    pub fn lots_msg(t: LotsTweak) -> [u8; 9] {
        let mut b = [0u8; 9];
        b[0] = T_LOTS_MSG;
        b[1..9].copy_from_slice(&t.to_be_bytes());
        b
    }

    /// `index || j || k`: Algorithm 4's tweak for the commitment
    /// `pi_{ijk} = H_vk(index||j||k, sh_{ij,k})` to participant `j`'s share
    /// of component `k` of the LOTS key with global index `index`.
    pub fn vk_share(t: LotsTweak, party: u32, i: u32) -> [u8; 17] {
        let mut b = [0u8; 17];
        b[0] = T_VK_SHARE;
        b[1..9].copy_from_slice(&t.to_be_bytes());
        b[9..13].copy_from_slice(&party.to_be_bytes());
        b[13..17].copy_from_slice(&i.to_be_bytes());
        b
    }

    /// `index || j`: Algorithm 4's tweak for participant `j`'s verification
    /// key `vk_{ij} = H_vk(index||j, (pi_{ijk})_k)`.
    pub fn vk_party(t: LotsTweak, party: u32) -> [u8; 13] {
        let mut b = [0u8; 13];
        b[0] = T_VK_PARTY;
        b[1..9].copy_from_slice(&t.to_be_bytes());
        b[9..13].copy_from_slice(&party.to_be_bytes());
        b
    }

    /// Tweak of one WOTS+ key: `t = r || index` (Algorithm 1's `t_r`), where
    /// `r` is the layer and `index` the key's position within that layer.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct WotsTweak {
        pub r: u32,
        pub idx: u64,
    }

    impl WotsTweak {
        fn write(&self, tag: u8, b: &mut [u8]) {
            b[0] = tag;
            b[1..5].copy_from_slice(&self.r.to_be_bytes());
            b[5..13].copy_from_slice(&self.idx.to_be_bytes());
        }

        /// `t || i || j`: chain `i`, height `j`.
        pub fn chain(&self, i: u32, j: u32) -> [u8; 21] {
            let mut b = [0u8; 21];
            self.write(T_WOTS_CHAIN, &mut b);
            b[13..17].copy_from_slice(&i.to_be_bytes());
            b[17..21].copy_from_slice(&j.to_be_bytes());
            b
        }

        /// `t` for the WOTS+ message hash (`D = H(P_msg, t, msg)`).
        pub fn msg(&self) -> [u8; 13] {
            let mut b = [0u8; 13];
            self.write(T_WOTS_MSG, &mut b);
            b
        }

        /// `t || i` for deriving the secret chain start `x_i = PRF(rho, t || i)`.
        pub fn sk(&self, i: u32) -> [u8; 17] {
            let mut b = [0u8; 17];
            self.write(T_WOTS_SK, &mut b);
            b[13..17].copy_from_slice(&i.to_be_bytes());
            b
        }
    }

    /// Appendix A.3's leaf tweak `0 || i`.
    pub fn merkle_leaf(i: u32) -> [u8; 5] {
        let mut b = [0u8; 5];
        b[0] = T_MERKLE_LEAF;
        b[1..5].copy_from_slice(&i.to_be_bytes());
        b
    }

    /// Appendix A.3's internal-node tweak `1 || j || i`, where `j` is the
    /// node's level counted from the root (root = level 0).
    pub fn merkle_node(j: u32, i: u32) -> [u8; 9] {
        let mut b = [0u8; 9];
        b[0] = T_MERKLE_NODE;
        b[1..5].copy_from_slice(&j.to_be_bytes());
        b[5..9].copy_from_slice(&i.to_be_bytes());
        b
    }

    /// PRF tweak of a LOTS leaf's secret seed, by hypertree address.
    pub fn seed_leaf(layer: u32, subtree: u64, leaf: u32) -> [u8; 17] {
        let mut b = [0u8; 17];
        b[0] = T_SEED_LEAF;
        b[1..5].copy_from_slice(&layer.to_be_bytes());
        b[5..13].copy_from_slice(&subtree.to_be_bytes());
        b[13..17].copy_from_slice(&leaf.to_be_bytes());
        b
    }

    /// PRF tweak of a bunch's `id_msg`.
    pub fn seed_id_msg(subtree: u64) -> [u8; 9] {
        let mut b = [0u8; 9];
        b[0] = T_SEED_IDMSG;
        b[1..9].copy_from_slice(&subtree.to_be_bytes());
        b
    }

    /// PRF tweak of one salt of a salted bottom tree.
    pub fn seed_salt(subtree: u64, leaf: u32) -> [u8; 13] {
        let mut b = [0u8; 13];
        b[0] = T_SEED_SALT;
        b[1..9].copy_from_slice(&subtree.to_be_bytes());
        b[9..13].copy_from_slice(&leaf.to_be_bytes());
        b
    }
}

/// Deterministically derive a LOTS leaf's secret seed from the master seed
/// `rho` and its hypertree address, so any subtree can be regenerated on
/// demand (Algorithm 1 derives every key from `rho` and its position).
pub fn derive_leaf_seed(master: &[u8; N], layer: u32, subtree: u64, leaf: u32) -> Digest32 {
    prf(master, &tweak::seed_leaf(layer, subtree, leaf))
}

/// Per-bunch `id_msg` (Algorithm 1's `id_msg <- {0,1}^lambda`), derived
/// deterministically so a bunch always regenerates the same value. It is
/// authenticated by the WOTS+ signature on `m_0 = id_msg || MT_0.Root`.
pub fn derive_bunch_id_msg(master: &[u8; N], subtree_idx: u64) -> Digest32 {
    prf(master, &tweak::seed_id_msg(subtree_idx))
}

/// Per-leaf salt of the salted bottom tree (Algorithm 1's `MT_0.Salt(i)`).
pub fn derive_leaf_salt(master: &[u8; N], subtree_idx: u64, leaf_pos: u32) -> Digest32 {
    prf(master, &tweak::seed_salt(subtree_idx, leaf_pos))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thash_separates_p_and_tweak_boundaries() {
        // Moving a byte between P and t must change the output.
        assert_ne!(thash(b"ab", b"c", b"x"), thash(b"a", b"bc", b"x"));
        assert_ne!(thash(b"a", b"bc", b"x"), thash(b"a", b"b", b"cx"));
    }

    #[test]
    fn public_params_are_distinct() {
        let pp = PubParams::with_pk(&[9u8; N], &[7u8; N]);
        let (m10, m20, m11) = (pp.merkle_upper(1, 0), pp.merkle_upper(2, 0), pp.merkle_upper(1, 1));
        let all: Vec<&[u8]> = vec![
            pp.lots(),
            pp.lots_msg(),
            pp.vk(),
            pp.wots(),
            pp.wots_msg(),
            &m10,
            &m20,
            &m11,
        ];
        for a in 0..all.len() {
            for b in a + 1..all.len() {
                assert_ne!(all[a], all[b], "parameters {a} and {b} coincide");
            }
        }
    }

    /// Merkle tweaks reuse Appendix A.3's leading `0`/`1` rather than a
    /// private tag space; that is only safe because no Merkle tree shares a
    /// public parameter with an OTS or vk hash. A Merkle parameter is 32 B
    /// (`id_msg`) or 44 B (`id || r || tau`); every OTS/vk parameter is
    /// 35-36 B (`id || label`) or 66-68 B (`pk || label`).
    #[test]
    fn merkle_parameters_never_coincide_with_ots_parameters() {
        let pp = PubParams::with_pk(&[9u8; N], &[7u8; N]);
        let ots: Vec<&[u8]> = vec![pp.lots(), pp.lots_msg(), pp.vk(), pp.wots(), pp.wots_msg()];
        let merkle: Vec<Vec<u8>> = vec![pp.merkle_upper(1, 0), pp.merkle_upper(3, 77), vec![0u8; N]];
        for m in &merkle {
            for o in &ots {
                assert_ne!(&m[..], *o);
            }
        }
    }

    #[test]
    fn ots_and_prf_tweak_types_are_distinct() {
        let w = tweak::WotsTweak { r: 1, idx: 0 };
        // Tweaks used with a public parameter.
        let ots = [
            tweak::lots_key(0, 0)[0],
            tweak::lots_msg(0)[0],
            w.chain(0, 0)[0],
            w.msg()[0],
            tweak::vk_share(0, 0, 0)[0],
            tweak::vk_party(0, 0)[0],
        ];
        // Tweaks used with the secret seed as a PRF key.
        let prf_tweaks = [
            w.sk(0)[0],
            tweak::seed_leaf(0, 0, 0)[0],
            tweak::seed_id_msg(0)[0],
            tweak::seed_salt(0, 0)[0],
        ];
        for tags in [&ots[..], &prf_tweaks[..]] {
            let set: std::collections::HashSet<_> = tags.iter().collect();
            assert_eq!(set.len(), tags.len());
        }
        // Merkle tweaks are Appendix A.3's literal 0 and 1.
        assert_eq!(tweak::merkle_leaf(3)[0], 0);
        assert_eq!(tweak::merkle_node(1, 3)[0], 1);
    }

    /// Binding the root must re-key only the pk-derived parameters; the
    /// hypertree's own parameters are computed before the root exists.
    #[test]
    fn binding_the_root_leaves_hypertree_parameters_alone() {
        let (id, root) = ([7u8; N], [9u8; N]);
        let unbound = PubParams::new(&id);
        let mut bound = unbound.clone();
        bound.bind_root(&root);
        assert_eq!(bound, PubParams::with_pk(&root, &id));
        assert_eq!(unbound.wots(), bound.wots());
        assert_eq!(unbound.wots_msg(), bound.wots_msg());
        assert_eq!(unbound.merkle_upper(2, 5), bound.merkle_upper(2, 5));
        assert_ne!(unbound.lots(), bound.lots());
        assert_ne!(unbound.lots_msg(), bound.lots_msg());
        assert_ne!(unbound.vk(), bound.vk());
        assert_eq!(bound.root(), Some(&root));
        assert_eq!(unbound.root(), None);
        // A different root gives different LOTS parameters.
        assert_ne!(bound.lots(), PubParams::with_pk(&[8u8; N], &id).lots());
    }
}
