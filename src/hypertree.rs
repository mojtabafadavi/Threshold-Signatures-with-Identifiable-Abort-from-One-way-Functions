//! Threshold multi-message signature from a Lamport-Winternitz hypertree
//! (Section 2.3, Algorithm 1, and Section 3), with every hash tweaked as in
//! Section 5 and Appendix A.1.
//!
//! Layer numbering. The code numbers subtree layers `0` (top) .. `d-1`
//! (bottom). Layer `d-1` holds the Shamir-shared CFF-Lamport keys; layers
//! `0..=d-2` hold WOTS+ keys. Algorithm 1 instead numbers the WOTS+ layers
//! `r = 1` (lowest) .. `r = d-1` (top), so code layer `l` is `r = d-1-l`.
//!
//! Public parameters and tweaks. The LOTS parameters are keyed on the whole
//! long-term public key `pk = (Root, id)` (Section 5's Equation (1) and
//! Algorithm 4); everything hashed while `Root` is still being computed --
//! WOTS+ and every Merkle tree -- is keyed on `id` alone, as Algorithm 1
//! does. See `hashing` for the full argument.
//!
//! | object | public parameter | tweak |
//! |---|---|---|
//! | LOTS key at global index `c` | `P_LOTS`, `P_msg` | `t = c` |
//! | share commitment / verification key | `P_vk` | `c‖j‖k`, `c‖j` |
//! | salted LOTS tree of bunch `beta` | `id_msg` (bunch-specific) | `0‖i`, `1‖j‖i` |
//! | WOTS+ key at position `idx` of layer `r` | `P_WOTS`, `P_Wmsg` | `t_r = r‖idx` (+ `‖i‖j` per chain step) |
//! | WOTS+ tree `tau` of layer `r` | `id‖r‖tau` | `0‖i`, `1‖j‖i` |
//!
//! The WOTS+ key directly above a bunch signs `m_0 = id_msg ‖ MT_0.Root`,
//! which authenticates `id_msg`. Merkle leaf positions are always taken
//! from the leaf index, never from the authentication path.
//!
//! Every subtree is a deterministic function of the master seed and its
//! address, so key generation and signing regenerate only
//! `O(d * 2^(h/d))` leaves.

use crate::hashing::{
    derive_bunch_id_msg, derive_leaf_salt, derive_leaf_seed, tweak::WotsTweak, Digest32, PubParams, N,
};
use crate::lots_cff::{self, LotsKeyPair, LotsParams, PartyShare};
use crate::merkle;
use crate::shamir;
use crate::wotsplus;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::time::{Duration, Instant};

pub struct HyperParams {
    /// Total tree height: capacity is 2^h signatures under one public key.
    pub h: u32,
    /// Number of subtree layers (`d | h`). Algorithm 1's number of WOTS+
    /// layers is `d - 1`.
    pub d: u32,
    pub lots: LotsParams,
    /// The public key's `id` (Algorithm 1's random public parameter); all
    /// public parameters `P_X` are derived from it.
    pub pub_seed: [u8; N],
    /// Secret master seed (Algorithm 1's `rho`), held only by the dealer.
    pub master_seed: [u8; N],
}

impl HyperParams {
    pub fn new(h: u32, d: u32, lots: LotsParams, pub_seed: [u8; N], master_seed: [u8; N]) -> Result<Self, String> {
        if d == 0 {
            return Err("d must be at least 1".into());
        }
        if h % d != 0 {
            return Err(format!("h={h} must be divisible by d={d}"));
        }
        if h > 64 {
            return Err("h > 64 is not supported (leaf indices are u64)".into());
        }
        if h / d > 32 {
            return Err("subtree height h/d > 32 is not supported".into());
        }
        Ok(HyperParams { h, d, lots, pub_seed, master_seed })
    }

    /// Height `h/d` of every subtree.
    pub fn layer_height(&self) -> u32 {
        self.h / self.d
    }
    /// `L = 2^(h/d)`.
    pub fn leaves_per_subtree(&self) -> usize {
        1usize << self.layer_height()
    }
    pub fn capacity(&self) -> u128 {
        1u128 << self.h
    }
    /// Algorithm 1's layer number `r` of code layer `layer` (a WOTS+ layer).
    fn paper_layer(&self, layer: u32) -> u32 {
        self.d - 1 - layer
    }
    /// Tweak of the WOTS+ key at position `idx` of code layer `layer`.
    fn wots_tweak(&self, layer: u32, idx: u64) -> WotsTweak {
        WotsTweak { r: self.paper_layer(layer), idx }
    }
    fn upper_tree_param(&self, pp: &PubParams, layer: u32, subtree_idx: u64) -> Vec<u8> {
        pp.merkle_upper(self.paper_layer(layer), subtree_idx)
    }
}

fn concat_pk(pk: &[Digest32]) -> Vec<u8> {
    pk.concat()
}

/// `m_0 = id_msg || MT_0.Root`: the message the lowest WOTS+ layer signs.
fn bunch_message(id_msg: &Digest32, root: &Digest32) -> Vec<u8> {
    [id_msg.as_slice(), root.as_slice()].concat()
}

/// The salted bottom tree's public parameter `id_msg`. With `d = 1` there
/// is no WOTS+ layer to authenticate a separate `id_msg`, so the public
/// key's own `id` is used instead.
fn bunch_id_msg(params: &HyperParams, pp: &PubParams, subtree_idx: u64) -> Digest32 {
    if params.d == 1 {
        pp.id
    } else {
        derive_bunch_id_msg(&params.master_seed, subtree_idx)
    }
}

/// Regenerate the LOTS key with global index `leaf_index` (its tweak).
fn lots_key(params: &HyperParams, pp: &PubParams, leaf_index: u64) -> LotsKeyPair {
    let l = params.leaves_per_subtree() as u64;
    let seed = derive_leaf_seed(&params.master_seed, params.d - 1, leaf_index / l, (leaf_index % l) as u32);
    let mut rng = StdRng::from_seed(seed);
    lots_cff::key_gen(&params.lots, pp, leaf_index, &mut rng)
}

/// A built subtree: root, per-leaf authentication paths, and (bottom
/// layer only) `id_msg` and salts.
struct Subtree {
    subtree_idx: u64,
    root: Digest32,
    auth_paths: Vec<Vec<Digest32>>,
    id_msg: Digest32,
    salts: Vec<Digest32>,
}

/// Build bunch `subtree_idx` of the bottom layer: `L` LOTS public keys
/// under tweaks `subtree_idx * L + a`, salted and hashed under `id_msg`.
fn build_bottom(params: &HyperParams, pp: &PubParams, subtree_idx: u64) -> (Subtree, Duration, Duration) {
    let l = params.leaves_per_subtree();
    let t0 = Instant::now();
    let id_msg = bunch_id_msg(params, pp, subtree_idx);
    let mut salts = Vec::with_capacity(l);
    let mut leaves = Vec::with_capacity(l);
    for a in 0..l {
        let kp = lots_key(params, pp, subtree_idx * l as u64 + a as u64);
        let salt = derive_leaf_salt(&params.master_seed, subtree_idx, a as u32);
        leaves.push(merkle::salt_leaf(&salt, &concat_pk(&kp.pk)));
        salts.push(salt);
    }
    let keygen = t0.elapsed();
    let t1 = Instant::now();
    let (root, auth_paths) = merkle::build_tree(&id_msg, &leaves);
    let tree = t1.elapsed();
    (Subtree { subtree_idx, root, auth_paths, id_msg, salts }, keygen, tree)
}

/// Build WOTS+ tree `subtree_idx` of code layer `layer`: `L` WOTS+ public
/// keys under tweaks `r || (subtree_idx * L + a)`, hashed under
/// `id || r || subtree_idx`.
fn build_upper(params: &HyperParams, pp: &PubParams, layer: u32, subtree_idx: u64) -> (Subtree, Duration, Duration) {
    let l = params.leaves_per_subtree();
    let t0 = Instant::now();
    let leaves: Vec<Vec<u8>> = (0..l)
        .map(|a| {
            let t = params.wots_tweak(layer, subtree_idx * l as u64 + a as u64);
            concat_pk(&wotsplus::key_gen(pp, &params.master_seed, &t).pk)
        })
        .collect();
    let keygen = t0.elapsed();
    let t1 = Instant::now();
    let (root, auth_paths) = merkle::build_tree(&params.upper_tree_param(pp, layer, subtree_idx), &leaves);
    let tree = t1.elapsed();
    let st = Subtree { subtree_idx, root, auth_paths, id_msg: [0u8; N], salts: Vec::new() };
    (st, keygen, tree)
}

/// Public certificate for one LOTS key: everything in the final signature
/// except the CFF-Lamport signature itself.
#[derive(Clone)]
pub struct CertChain {
    /// The full LOTS public key; the final signature carries only `pk_cmpl`.
    pub pk_lots: Vec<Digest32>,
    pub id_msg: Digest32,
    pub salt: Digest32,
    pub auth_bottom: Vec<Digest32>,
    /// WOTS+ signatures, lowest layer first.
    pub upper_sigs: Vec<Vec<Digest32>>,
    /// Authentication paths of those WOTS+ keys, lowest layer first.
    pub upper_auths: Vec<Vec<Digest32>>,
}

/// Walk from a built bottom subtree up to the root (Algorithm 1's `BSign`
/// loop), signing each child message with the WOTS+ key above it.
/// `upper(layer, subtree_idx)` supplies the subtree at that address.
fn certify<'a>(
    params: &HyperParams,
    pp: &PubParams,
    leaf_index: u64,
    bottom: &Subtree,
    mut upper: impl FnMut(u32, u64) -> &'a Subtree,
) -> (LotsKeyPair, CertChain, Digest32) {
    let l = params.leaves_per_subtree() as u64;
    let kp = lots_key(params, pp, leaf_index);
    let pos = (leaf_index % l) as usize;

    let mut root = bottom.root;
    let mut msg = bunch_message(&bottom.id_msg, &bottom.root);
    let mut idx = bottom.subtree_idx; // position of the signing WOTS+ key in its layer
    let mut upper_sigs = Vec::with_capacity(params.d as usize - 1);
    let mut upper_auths = Vec::with_capacity(params.d as usize - 1);
    for layer in (0..params.d - 1).rev() {
        let t = params.wots_tweak(layer, idx);
        upper_sigs.push(wotsplus::sign(pp, &params.master_seed, &t, &msg));
        let st = upper(layer, idx / l);
        upper_auths.push(st.auth_paths[(idx % l) as usize].clone());
        root = st.root;
        msg = root.to_vec();
        idx /= l;
    }

    let cert = CertChain {
        pk_lots: kp.pk.clone(),
        id_msg: bottom.id_msg,
        salt: bottom.salts[pos],
        auth_bottom: bottom.auth_paths[pos].clone(),
        upper_sigs,
        upper_auths,
    };
    (kp, cert, root)
}

pub struct Dealer {
    pub params: HyperParams,
    /// Public parameters: `id`-keyed for the hypertree, and (once
    /// `key_gen` has the root) `pk`-keyed for LOTS and `H_vk`.
    pub pp: PubParams,
    /// `pk = (root, id)`: this is `id`.
    pub id: Digest32,
    /// `pk = (root, id)`: this is `root`.
    pub pk: Digest32,
    /// Algorithm 1's `counter`.
    next_leaf: u64,
    /// The currently active bottom subtree (public data only).
    bottom_cache: Option<Subtree>,
    /// One slot per upper layer `0..=d-2`.
    upper_cache: Vec<Option<Subtree>>,
}

impl Dealer {
    /// `BKeyGen`: builds (and caches) the single top-layer subtree, whose
    /// root is `pk`.
    pub fn key_gen(params: HyperParams) -> Self {
        let id = params.pub_seed;
        // Only the `id`-keyed parameters exist yet: the top subtree's WOTS+
        // keys and Merkle nodes are precisely what produces `Root`.
        let mut pp = PubParams::new(&id);
        let mut upper_cache: Vec<Option<Subtree>> = (0..params.d - 1).map(|_| None).collect();
        let pk = if params.d == 1 {
            // Degenerate shape: the LOTS tree is itself the top tree, so its
            // keys are also hashed before `Root` exists and the LOTS
            // parameters stay keyed on `id` (see `hashing`).
            build_bottom(&params, &pp, 0).0.root
        } else {
            let top = build_upper(&params, &pp, 0, 0).0;
            let root = top.root;
            upper_cache[0] = Some(top);
            // `pk = (Root, id)` is fixed at this point, so `P_LOTS`, `P_msg` and
            // `P_vk` can be keyed on it (Algorithm 4, line 2). Every LOTS
            // key is generated later, in `resolve_leaf`/`issue_next_leaf`.
            pp.bind_root(&root);
            root
        };
        Dealer { params, pp, id, pk, next_leaf: 0, bottom_cache: None, upper_cache }
    }

    /// Regenerate everything needed for leaf `leaf_index` from scratch,
    /// returning its LOTS keypair (to be Shamir-shared) and certificate.
    pub fn resolve_leaf(&self, leaf_index: u64) -> (LotsKeyPair, CertChain) {
        let (kp, cert, _, _) = self.resolve_leaf_timed(leaf_index);
        (kp, cert)
    }

    /// As `resolve_leaf`, additionally returning per-layer key-generation
    /// and tree-building times, bottom layer first.
    pub fn resolve_leaf_timed(&self, leaf_index: u64) -> (LotsKeyPair, CertChain, Vec<Duration>, Vec<Duration>) {
        assert!((leaf_index as u128) < self.params.capacity(), "leaf index out of range");
        let l = self.params.leaves_per_subtree() as u64;
        let (bottom, kt, tt) = build_bottom(&self.params, &self.pp, leaf_index / l);
        let mut keygen_times = vec![kt];
        let mut tree_times = vec![tt];

        // Build the d-1 ancestor subtrees first (bottom to top), then certify.
        let mut uppers: Vec<Subtree> = Vec::with_capacity(self.params.d as usize - 1);
        let mut sub = leaf_index / l;
        for layer in (0..self.params.d - 1).rev() {
            sub /= l;
            let (st, kt, tt) = build_upper(&self.params, &self.pp, layer, sub);
            keygen_times.push(kt);
            tree_times.push(tt);
            uppers.push(st);
        }
        let d = self.params.d;
        let (kp, cert, root) = certify(&self.params, &self.pp, leaf_index, &bottom, |layer, _| {
            &uppers[(d - 2 - layer) as usize]
        });
        debug_assert_eq!(root, self.pk, "resolved root must match the public key");
        (kp, cert, keygen_times, tree_times)
    }

    /// Issue the next unused leaf (`counter`), or `None` once all `2^h`
    /// leaves are used. Caches one bottom subtree and one subtree per upper
    /// layer (public data only).
    pub fn issue_next_leaf(&mut self) -> Option<(LotsKeyPair, CertChain)> {
        self.issue_next_leaf_verbose().map(|(kp, cert, _)| (kp, cert))
    }

    /// As `issue_next_leaf`, also reporting per layer (bottom first) whether
    /// that layer's subtree was rebuilt on this call.
    pub fn issue_next_leaf_verbose(&mut self) -> Option<(LotsKeyPair, CertChain, Vec<bool>)> {
        if (self.next_leaf as u128) >= self.params.capacity() {
            return None;
        }
        let leaf_index = self.next_leaf;
        self.next_leaf += 1;
        let l = self.params.leaves_per_subtree() as u64;
        let mut rebuilt = Vec::with_capacity(self.params.d as usize);

        let bottom_idx = leaf_index / l;
        let hit = self.bottom_cache.as_ref().map(|c| c.subtree_idx) == Some(bottom_idx);
        if !hit {
            self.bottom_cache = Some(build_bottom(&self.params, &self.pp, bottom_idx).0);
        }
        rebuilt.push(!hit);

        let mut sub = bottom_idx;
        for layer in (0..self.params.d - 1).rev() {
            sub /= l;
            let slot = &mut self.upper_cache[layer as usize];
            let hit = slot.as_ref().map(|c| c.subtree_idx) == Some(sub);
            if !hit {
                *slot = Some(build_upper(&self.params, &self.pp, layer, sub).0);
            }
            rebuilt.push(!hit);
        }

        let caches = &self.upper_cache;
        let (kp, cert, root) = certify(
            &self.params,
            &self.pp,
            leaf_index,
            self.bottom_cache.as_ref().unwrap(),
            |layer, _| caches[layer as usize].as_ref().unwrap(),
        );
        debug_assert_eq!(root, self.pk, "issued leaf's root must match the public key");
        Some((kp, cert, rebuilt))
    }
}

/// Assembled signature `Sigma = (sigma^L, pk_cmpl, index, Sig)`, where
/// `Sig = (counter, id_msg, MT_0.Path(i), MT_0.Salt(i), Sigma^W)`
/// (Algorithm 3's `TLOTS.Verify`, Algorithm 1's `BSign`).
#[derive(Clone)]
pub struct Signature {
    /// The LOTS key's global index, i.e. Algorithm 4's `index` and the LOTS
    /// tweak `t`. Algorithm 1's `counter` (the bunch) and the leaf's
    /// position `i` within its bunch are `index / L` and `index % L`; both
    /// are derived from this one field, never read from a path.
    pub index: u64,
    pub sigma_cff: Vec<[u8; 32]>,
    /// `(y_i)_{i not in B_msg}`.
    pub pk_cmpl: Vec<Digest32>,
    pub id_msg: Digest32,
    pub salt: Digest32,
    pub auth_bottom: Vec<Digest32>,
    pub upper_sigs: Vec<Vec<Digest32>>,
    pub upper_auths: Vec<Vec<Digest32>>,
}

impl Signature {
    /// Serialized size in bytes, including the 8-byte `index` (the LOTS
    /// key's global index, from which the bunch counter and the leaf's
    /// position within it are derived).
    pub fn size_bytes(&self) -> usize {
        8 + self.sigma_cff.len() * 32
            + self.pk_cmpl.len() * N
            + 2 * N
            + self.auth_bottom.len() * N
            + self.upper_sigs.iter().map(|s| s.len() * N).sum::<usize>()
            + self.upper_auths.iter().map(|a| a.len() * N).sum::<usize>()
    }
}

/// Threshold partial signature of the participant holding `share`, for
/// the LOTS key with global index `index`.
pub fn part_sign(
    params: &LotsParams,
    pp: &PubParams,
    index: u64,
    msg: &[u8],
    share: &PartyShare,
) -> (u64, Vec<[u64; shamir::R]>) {
    (share.point, lots_cff::part_sign(params, pp, index, msg, share))
}

/// Combine at least `t` partial signatures into the final signature.
pub fn combine(
    params: &LotsParams,
    pp: &PubParams,
    msg: &[u8],
    t: usize,
    index: u64,
    cert: CertChain,
    partials: &[(u64, Vec<[u64; shamir::R]>)],
) -> Option<Signature> {
    let sigma_cff = lots_cff::combine(params, pp, index, msg, t, partials)?;
    if cert.pk_lots.len() != params.e {
        return None;
    }
    let pk_cmpl = lots_cff::complementary_pk(params, pp, index, &cert.pk_lots, msg);
    Some(Signature {
        index,
        sigma_cff,
        pk_cmpl,
        id_msg: cert.id_msg,
        salt: cert.salt,
        auth_bottom: cert.auth_bottom,
        upper_sigs: cert.upper_sigs,
        upper_auths: cert.upper_auths,
    })
}

/// `Verify` against `pk = (pk_root, pk_id)`: derive the LOTS public key
/// with `LDerivePK`, then the salted bottom root, then walk up the WOTS+
/// layers with `WDerivePK` and `DeriveRoot`, and compare with `pk_root`.
pub fn verify(params: &HyperParams, pk_id: &Digest32, pk_root: &Digest32, msg: &[u8], sig: &Signature) -> bool {
    let d = params.d as usize;
    let l = params.leaves_per_subtree() as u64;
    let hh = params.layer_height() as usize;
    if (sig.index as u128) >= params.capacity()
        || sig.auth_bottom.len() != hh
        || sig.upper_sigs.len() != d - 1
        || sig.upper_auths.len() != d - 1
        || sig.upper_auths.iter().any(|a| a.len() != hh)
    {
        return false;
    }
    // `P_LOTS`, `P_msg` and `P_vk` are keyed on the whole public key
    // `pk = (Root, id)` (Algorithm 4); with `d = 1` the LOTS tree is the top
    // tree, so they are keyed on `id` alone, exactly as the dealer does.
    let pp = if d == 1 {
        PubParams::new(pk_id)
    } else {
        PubParams::with_pk(pk_root, pk_id)
    };
    if d == 1 && sig.id_msg != *pk_id {
        return false;
    }

    let Some(pk_lots) = lots_cff::derive_pk(&params.lots, &pp, sig.index, msg, &sig.sigma_cff, &sig.pk_cmpl) else {
        return false;
    };
    let Some(mut root) =
        merkle::derive_root_salted(&sig.id_msg, sig.index % l, &concat_pk(&pk_lots), &sig.salt, &sig.auth_bottom)
    else {
        return false;
    };

    let mut msg_up = bunch_message(&sig.id_msg, &root);
    let mut idx = sig.index / l;
    for (k, layer) in (0..params.d - 1).rev().enumerate() {
        let t = params.wots_tweak(layer, idx);
        let Some(wpk) = wotsplus::derive_pk(&pp, &t, &msg_up, &sig.upper_sigs[k]) else {
            return false;
        };
        let tree_p = params.upper_tree_param(&pp, layer, idx / l);
        let Some(r) = merkle::derive_root(&tree_p, idx % l, &concat_pk(&wpk), &sig.upper_auths[k]) else {
            return false;
        };
        root = r;
        msg_up = root.to_vec();
        idx /= l;
    }
    root == *pk_root
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{thread_rng, RngCore};

    fn seeds() -> ([u8; N], [u8; N]) {
        let mut rng = thread_rng();
        let (mut a, mut b) = ([0u8; N], [0u8; N]);
        rng.fill_bytes(&mut a);
        rng.fill_bytes(&mut b);
        (a, b)
    }

    fn sign_with(dealer: &Dealer, idx: u64, cert: CertChain, kp: &LotsKeyPair, msg: &[u8]) -> Signature {
        let mut rng = thread_rng();
        let (t, n) = (3, 5);
        let (shares, _vk) = lots_cff::dist(&dealer.pp, idx, &kp.sk, t, n, &mut rng);
        let partials: Vec<_> = shares
            .iter()
            .take(t)
            .map(|sh| part_sign(&dealer.params.lots, &dealer.pp, idx, msg, sh))
            .collect();
        combine(&dealer.params.lots, &dealer.pp, msg, t, idx, cert, &partials).unwrap()
    }

    fn small_dealer(h: u32, d: u32) -> Dealer {
        let (id, master) = seeds();
        Dealer::key_gen(HyperParams::new(h, d, LotsParams::new(20, 8), id, master).unwrap())
    }

    #[test]
    fn end_to_end_small() {
        let dealer = small_dealer(4, 2);
        let idx = 5;
        let (kp, cert) = dealer.resolve_leaf(idx);
        assert_eq!(kp.tweak, idx);
        let msg = b"threshold hypertree message";
        let sig = sign_with(&dealer, idx, cert, &kp, msg);
        assert!(verify(&dealer.params, &dealer.id, &dealer.pk, msg, &sig));
        assert!(!verify(&dealer.params, &dealer.id, &dealer.pk, b"different message", &sig));
        assert_eq!(sig.pk_cmpl.len(), 20 - 8);
    }

    #[test]
    fn tampering_is_rejected() {
        let dealer = small_dealer(6, 3);
        let idx = 9;
        let (kp, cert) = dealer.resolve_leaf(idx);
        let msg = b"m";
        let sig = sign_with(&dealer, idx, cert, &kp, msg);
        assert!(verify(&dealer.params, &dealer.id, &dealer.pk, msg, &sig));

        // Wrong public parameter id.
        let mut other_id = dealer.id;
        other_id[0] ^= 1;
        assert!(!verify(&dealer.params, &other_id, &dealer.pk, msg, &sig));

        // Wrong root: it is both the value the chain must reach and part of
        // the key the LOTS parameters are derived from.
        let mut other_root = dealer.pk;
        other_root[0] ^= 1;
        assert!(!verify(&dealer.params, &dealer.id, &other_root, msg, &sig));

        // Claiming a different index changes every tweak and position.
        for bad in [idx - 1, idx + 1, idx + 4, idx + 16] {
            let mut s = sig.clone();
            s.index = bad;
            assert!(!verify(&dealer.params, &dealer.id, &dealer.pk, msg, &s), "index {bad}");
        }
        let mut s = sig.clone();
        s.index = 1u64 << 6; // out of range
        assert!(!verify(&dealer.params, &dealer.id, &dealer.pk, msg, &s));

        // id_msg is authenticated through m_0 = id_msg || root.
        let mut s = sig.clone();
        s.id_msg[0] ^= 1;
        assert!(!verify(&dealer.params, &dealer.id, &dealer.pk, msg, &s));

        let mut s = sig.clone();
        s.salt[0] ^= 1;
        assert!(!verify(&dealer.params, &dealer.id, &dealer.pk, msg, &s));

        let mut s = sig.clone();
        s.pk_cmpl[0][0] ^= 1;
        assert!(!verify(&dealer.params, &dealer.id, &dealer.pk, msg, &s));

        let mut s = sig.clone();
        s.upper_sigs[1][0][0] ^= 1;
        assert!(!verify(&dealer.params, &dealer.id, &dealer.pk, msg, &s));

        let mut s = sig.clone();
        s.upper_auths.pop();
        assert!(!verify(&dealer.params, &dealer.id, &dealer.pk, msg, &s));
    }

    #[test]
    fn rejects_non_divisible_h_d() {
        assert!(HyperParams::new(10, 3, LotsParams::new(20, 8), [0u8; N], [0u8; N]).is_err());
    }

    #[test]
    fn d_equals_one_works() {
        let dealer = small_dealer(3, 1);
        let (kp, cert) = dealer.resolve_leaf(2);
        assert!(cert.upper_sigs.is_empty());
        assert_eq!(cert.id_msg, dealer.id);
        let msg = b"d=1 case";
        let sig = sign_with(&dealer, 2, cert, &kp, msg);
        assert!(verify(&dealer.params, &dealer.id, &dealer.pk, msg, &sig));
        let mut s = sig.clone();
        s.id_msg[0] ^= 1;
        assert!(!verify(&dealer.params, &dealer.id, &dealer.pk, msg, &s));
    }

    #[test]
    fn four_layer_height_16_matches_spec() {
        let params = HyperParams::new(64, 4, LotsParams::new(20, 8), [0u8; N], [1u8; N]).unwrap();
        assert_eq!(params.layer_height(), 16);
        assert_eq!(params.leaves_per_subtree(), 1 << 16);
        // Code layer 0 is the top: paper layer r = d - 1 = 3.
        assert_eq!(params.paper_layer(0), 3);
        assert_eq!(params.paper_layer(2), 1);
    }

    /// The cached `issue_next_leaf` path must agree with `resolve_leaf` for
    /// every leaf, including across subtree boundaries.
    #[test]
    fn issue_next_leaf_matches_resolve_leaf() {
        let (id, master) = seeds();
        let mk = || HyperParams::new(4, 2, LotsParams::new(20, 8), id, master).unwrap();
        let capacity = mk().capacity() as u64;
        let reference = Dealer::key_gen(mk());
        let mut cached = Dealer::key_gen(mk());
        assert_eq!(cached.pk, reference.pk);
        for i in 0..capacity {
            let (kr, cr) = reference.resolve_leaf(i);
            let (kc, cc) = cached.issue_next_leaf().unwrap();
            assert_eq!((kc.tweak, &kc.sk, &kc.pk), (kr.tweak, &kr.sk, &kr.pk), "leaf {i}");
            assert_eq!(cc.pk_lots, cr.pk_lots, "leaf {i}");
            assert_eq!(cc.id_msg, cr.id_msg, "leaf {i}");
            assert_eq!(cc.salt, cr.salt, "leaf {i}");
            assert_eq!(cc.auth_bottom, cr.auth_bottom, "leaf {i}");
            assert_eq!(cc.upper_sigs, cr.upper_sigs, "leaf {i}");
            assert_eq!(cc.upper_auths, cr.upper_auths, "leaf {i}");
        }
        assert!(cached.issue_next_leaf().is_none());
    }

    #[test]
    fn issue_next_leaf_signs_and_verifies() {
        let mut dealer = small_dealer(6, 3);
        for idx in 0..20u64 {
            let (kp, cert) = dealer.issue_next_leaf().unwrap();
            let msg = format!("message for leaf {idx}");
            let sig = sign_with(&dealer, idx, cert, &kp, msg.as_bytes());
            assert!(verify(&dealer.params, &dealer.id, &dealer.pk, msg.as_bytes(), &sig), "leaf {idx}");
        }
    }

    /// The dealer's LOTS parameters must be exactly the ones a verifier
    /// rebuilds from the published `pk = (Root, id)` (Algorithm 4), while
    /// the hypertree's own parameters stay keyed on `id`.
    #[test]
    fn dealer_parameters_match_the_verifier_s() {
        let dealer = small_dealer(4, 2);
        assert_eq!(dealer.pp, PubParams::with_pk(&dealer.pk, &dealer.id));
        assert_eq!(dealer.pp.root(), Some(&dealer.pk));
        assert_eq!(dealer.pp.wots(), PubParams::new(&dealer.id).wots());
        // d = 1 leaves them id-keyed, on both sides.
        let flat = small_dealer(3, 1);
        assert_eq!(flat.pp, PubParams::new(&flat.id));
    }

    /// Identifiable abort end to end: the dealer issues verification keys
    /// alongside the shares, a corrupt participant is pinned down, and the
    /// honest remainder still produces a signature that verifies under the
    /// long-term public key.
    #[test]
    fn identifiable_abort_through_the_hypertree() {
        let mut rng = thread_rng();
        let dealer = small_dealer(4, 2);
        let idx = 3;
        let (kp, cert) = dealer.resolve_leaf(idx);
        let (t, n) = (3, 6);
        let (shares, vk) = lots_cff::dist(&dealer.pp, idx, &kp.sk, t, n, &mut rng);
        let msg = b"blame the right party";

        let mut partials: Vec<_> = shares
            .iter()
            .map(|sh| part_sign(&dealer.params.lots, &dealer.pp, idx, msg, sh))
            .collect();
        partials[1].1[0][0] ^= 1;

        let blamed: Vec<u32> = shares
            .iter()
            .zip(&partials)
            .filter(|(sh, (_, partial))| {
                let proof = lots_cff::sig_proof(&dealer.params.lots, &dealer.pp, idx, msg, sh);
                !lots_cff::verify_partial(
                    &dealer.params.lots,
                    &dealer.pp,
                    idx,
                    sh.id,
                    &vk[sh.id as usize],
                    msg,
                    partial,
                    &proof,
                )
            })
            .map(|(sh, _)| sh.id)
            .collect();
        assert_eq!(blamed, vec![1]);

        let honest: Vec<_> = partials
            .iter()
            .enumerate()
            .filter(|(j, _)| !blamed.contains(&(*j as u32)))
            .take(t)
            .map(|(_, p)| p.clone())
            .collect();
        let sig = combine(&dealer.params.lots, &dealer.pp, msg, t, idx, cert, &honest).unwrap();
        assert!(verify(&dealer.params, &dealer.id, &dealer.pk, msg, &sig));
    }

    /// Distinct LOTS leaves, and distinct WOTS+ keys, never share a
    /// public key even if a seed were reused, because their tweaks differ.
    #[test]
    fn tweaks_separate_keys() {
        let dealer = small_dealer(4, 2);
        let a = lots_key(&dealer.params, &dealer.pp, 1);
        let b = lots_cff::key_gen(&dealer.params.lots, &dealer.pp, 2, &mut StdRng::from_seed(derive_leaf_seed(
            &dealer.params.master_seed, 1, 0, 1,
        )));
        assert_eq!(a.sk, b.sk);
        assert_ne!(a.pk, b.pk);

        let t1 = dealer.params.wots_tweak(0, 3);
        let t2 = dealer.params.wots_tweak(1, 3);
        let k1 = wotsplus::key_gen(&dealer.pp, &dealer.params.master_seed, &t1);
        let k2 = wotsplus::key_gen(&dealer.pp, &dealer.params.master_seed, &t2);
        assert_ne!(k1.pk, k2.pk);
    }
}
