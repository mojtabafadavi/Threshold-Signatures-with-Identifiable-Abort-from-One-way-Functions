//! Tweaked, domain-separated Merkle trees (Section 2.2 and Appendix A.3).
//!
//! Every hash call is the tweakable hash `H(P, t, m)` (`hashing::thash`):
//! * `P` is the tree's public parameter -- `id || r || tau` for an upper
//!   (WOTS+) tree, or `id_msg` for the salted bottom tree;
//! * leaves use `y_H[i] = H(P, 0 || i, x_i)` (salted: `H(P, 0 || i, r_i || x_i)`);
//! * internal nodes use `y_j[i] = H(P, 1 || j || i, y_{j+1}[2i] || y_{j+1}[2i+1])`,
//!   where `j` is the node's level counted from the root (root = level 0,
//!   leaves = level H), exactly as in Appendix A.3.
//!
//! The leaf position is passed explicitly to `derive_root`, so a verifier
//! takes it from the signature's counter rather than from the path.

use crate::hashing::{thash, tweak, Digest32, N};
use rand::RngCore;

fn height_of(n: usize) -> u32 {
    assert!(n >= 1 && n.is_power_of_two(), "number of leaves must be a power of two");
    let h = n.trailing_zeros();
    assert!(h <= 32, "tree height must fit in u32 leaf indices");
    h
}

fn node(p: &[u8], height: u32, from_leaves: u32, i: usize, left: &Digest32, right: &Digest32) -> Digest32 {
    // A node `from_leaves + 1` levels above the leaves sits at level
    // `height - from_leaves - 1` counted from the root.
    let j = height - from_leaves - 1;
    thash(p, &tweak::merkle_node(j, i as u32), &[left.as_slice(), right.as_slice()].concat())
}

/// `MerkleTree(P, (x_0, ..., x_{n-1}))` with `n` a power of two: returns
/// the root and, for every leaf, its authentication path (bottom to top).
pub fn build_tree(p: &[u8], leaves_data: &[Vec<u8>]) -> (Digest32, Vec<Vec<Digest32>>) {
    let n = leaves_data.len();
    let h = height_of(n);
    let leaves: Vec<Digest32> = leaves_data
        .iter()
        .enumerate()
        .map(|(i, d)| thash(p, &tweak::merkle_leaf(i as u32), d))
        .collect();
    let mut levels: Vec<Vec<Digest32>> = vec![leaves];
    for lvl in 0..h {
        let cur = &levels[lvl as usize];
        let next: Vec<Digest32> = (0..cur.len() / 2)
            .map(|i| node(p, h, lvl, i, &cur[2 * i], &cur[2 * i + 1]))
            .collect();
        levels.push(next);
    }
    let root = levels[h as usize][0];
    let auth_paths = (0..n)
        .map(|leaf_idx| {
            let mut idx = leaf_idx;
            (0..h as usize)
                .map(|lvl| {
                    let s = levels[lvl][idx ^ 1];
                    idx /= 2;
                    s
                })
                .collect()
        })
        .collect();
    (root, auth_paths)
}

/// Prefix a leaf value with its salt: `r_i || x_i`.
pub fn salt_leaf(salt: &[u8; N], data: &[u8]) -> Vec<u8> {
    [salt.as_slice(), data].concat()
}

/// `MerkleTree_salted(P, ...)` with freshly sampled salts; also returns the salts.
pub fn build_tree_salted(
    p: &[u8],
    leaves_data: &[Vec<u8>],
    rng: &mut impl RngCore,
) -> (Digest32, Vec<Vec<Digest32>>, Vec<[u8; N]>) {
    let salts: Vec<[u8; N]> = (0..leaves_data.len())
        .map(|_| {
            let mut s = [0u8; N];
            rng.fill_bytes(&mut s);
            s
        })
        .collect();
    let salted: Vec<Vec<u8>> = leaves_data.iter().zip(&salts).map(|(d, s)| salt_leaf(s, d)).collect();
    let (root, auth_paths) = build_tree(p, &salted);
    (root, auth_paths, salts)
}

/// `MerkleTree.DeriveRoot(P, x, i, path)`: recompute the candidate root
/// from leaf value `x` at position `i`. The tree height is `path.len()`.
/// Returns `None` if the path is too long or `i` lies outside the tree.
pub fn derive_root(p: &[u8], leaf_idx: u64, leaf_data: &[u8], auth: &[Digest32]) -> Option<Digest32> {
    if auth.len() > 32 {
        return None;
    }
    let h = auth.len() as u32;
    if leaf_idx >> h != 0 {
        return None;
    }
    let mut idx = leaf_idx as usize;
    let mut cur = thash(p, &tweak::merkle_leaf(idx as u32), leaf_data);
    for (lvl, sibling) in auth.iter().enumerate() {
        cur = if idx % 2 == 0 {
            node(p, h, lvl as u32, idx / 2, &cur, sibling)
        } else {
            node(p, h, lvl as u32, idx / 2, sibling, &cur)
        };
        idx /= 2;
    }
    Some(cur)
}

/// `MerkleTree_salted.DeriveRoot(P, x, i, path, salt)`.
pub fn derive_root_salted(
    p: &[u8],
    leaf_idx: u64,
    leaf_data: &[u8],
    salt: &[u8; N],
    auth: &[Digest32],
) -> Option<Digest32> {
    derive_root(p, leaf_idx, &salt_leaf(salt, leaf_data), auth)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::thread_rng;

    #[test]
    fn build_and_verify_roundtrip() {
        let p = [7u8; N];
        let leaves: Vec<Vec<u8>> = (0..8u8).map(|i| vec![i; 3]).collect();
        let (root, auths) = build_tree(&p, &leaves);
        for (i, leaf) in leaves.iter().enumerate() {
            assert_eq!(derive_root(&p, i as u64, leaf, &auths[i]), Some(root));
        }
    }

    #[test]
    fn wrong_position_is_rejected() {
        let p = [7u8; N];
        let leaves: Vec<Vec<u8>> = (0..8u8).map(|i| vec![i; 3]).collect();
        let (root, auths) = build_tree(&p, &leaves);
        assert_ne!(derive_root(&p, 3, &leaves[2], &auths[2]), Some(root));
        assert_eq!(derive_root(&p, 10, &leaves[2], &auths[2]), None);
    }

    #[test]
    fn salted_roundtrip() {
        let mut rng = thread_rng();
        let p = [9u8; N];
        let leaves: Vec<Vec<u8>> = (0..4u8).map(|i| vec![i; 5]).collect();
        let (root, auths, salts) = build_tree_salted(&p, &leaves, &mut rng);
        for (i, leaf) in leaves.iter().enumerate() {
            assert_eq!(derive_root_salted(&p, i as u64, leaf, &salts[i], &auths[i]), Some(root));
        }
    }

    #[test]
    fn different_parameters_give_different_roots() {
        let leaves: Vec<Vec<u8>> = (0..4u8).map(|i| vec![i]).collect();
        assert_ne!(build_tree(&[1u8; N], &leaves).0, build_tree(&[2u8; N], &leaves).0);
    }

    #[test]
    fn single_leaf_tree() {
        let p = [3u8; N];
        let leaves = vec![vec![42u8]];
        let (root, auths) = build_tree(&p, &leaves);
        assert!(auths[0].is_empty());
        assert_eq!(derive_root(&p, 0, &leaves[0], &auths[0]), Some(root));
    }
}
