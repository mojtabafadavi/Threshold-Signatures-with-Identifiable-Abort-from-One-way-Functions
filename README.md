# tlots

`tlots` is a `(T,N)`-threshold post-quantum signature scheme with
identifiable abort, implemented purely in Rust and relying on nothing beyond
one-way functions.

Any `T` of `N` participants jointly produce a signature under one long-term
public key; fewer than `T` cannot. When a round fails, the participants that
submitted malformed contributions are named, dropped, and a valid signature
is recovered from the honest remainder — without a second signing round and
without enlarging the signature.

* **Assumption.** One-way functions only: no lattice, code or isogeny
  assumption anywhere in the scheme.
* **Capacity.** `2^h` signatures under one public key, with `h = 64`
  reported and configurable at run time.
* **Signature size.** 16,904 B at the reported parameters.
* **Signing.** One round given consensus on the message and key index; each
  participant's contribution is non-interactive and sub-millisecond.
* **Identifiable abort.** Verification keys distributed once at setup, with
  openings that travel only when a round fails, costing about 10 ms on the
  failure path and nothing on the common one.

## Paper

The scheme, its UC security proof and the parameter choices are from
*Threshold Signatures with Identifiable Abort from One-way Functions*. Module
documentation throughout this crate cites that paper's sections, algorithms
and appendices by number.

## Features

### Core functionality

| Stage | Entry point |
|---|---|
| Hypertree key generation (`BKeyGen`) | `Dealer::key_gen` |
| Leaf issuance, cached (`counter`) | `Dealer::issue_next_leaf` |
| Leaf issuance, arbitrary index | `Dealer::resolve_leaf` |
| Key distribution (`IssueShares`) | `lots_cff::dist` |
| Partial signing (`TLOTS.PartialSignature`) | `hypertree::part_sign` |
| Combining (`TLOTS.ReconstructSignature`) | `hypertree::combine` |
| Verification (`TLOTS.Verify`) | `hypertree::verify` |
| Abort opening (`TLOTS.SigProof`) | `lots_cff::sig_proof` |
| Opening check (`TLOTS.SignFailed`) | `lots_cff::verify_partial` |

### Cryptographic primitives

* One tweakable hash `H(P, t, x)` (Appendix A.1) over SHA-256, plus a
  secret-seed PRF, behind which every hash in the crate goes.
* CFF-based Lamport OTS (§2.1) at the leaves, with the `κ`-uniform
  1-cover-free family realized by combinatorial-number-system decoding.
* Tweaked WOTS+ (Appendix A.2) at the certification layers, `w = 4`,
  base 16.
* Tweaked Merkle trees (§2.2, Appendix A.3), salted at the bottom layer.
* Component-wise Shamir sharing over `F_p` (§2.4, Algorithm 2) for an
  FFT-friendly `p`.

### Security

* Post-quantum: the only assumption is one-wayness, and all parameters are
  set for `λ = 256`.
* Every tweak is a fixed-width encoding behind a type tag, and the public
  parameters are pairwise distinct byte strings, so no two uses in the
  scheme share a `(P, t)` pair — the `DIST` condition the SM-TCR and
  SM-DSPR reductions need. `hashing`'s tests check this.
* Merkle leaf positions and every tweak are derived from the signature's
  own key index, never read from an authentication path; out-of-range
  indices and malformed lengths are rejected.
* Verification keys bind a participant's share vector to its identity, the
  key index and the message's block, so none of the three can be swapped.

### Efficiency

* Key distribution costs `O(N log N)` per secret component rather than
  `O(N·T)`, via an NTT over a root-of-unity evaluation domain, and is
  therefore independent of `T`.
* Combining shares the quorum's Lagrange coefficients across all `κ·r`
  coordinates, for `Θ(T² + κ·r·T)`.
* Every subtree is a deterministic function of the master seed and its
  address, so signing regenerates `O(d · 2^(h/d))` leaves rather than the
  full `2^h`-leaf tree, and `h = 64` is practical.
* Single-threaded throughout.

## Installation

### Prerequisites

* Rust 2021 edition. `let ... else` puts the minimum supported version at
  **1.65.0**; builds and benchmarks here use **1.89.0**.
* `Dealer::issue_next_leaf` caches one bottom subtree and its `d-1`
  ancestors, about **0.9 GiB** at `h = 64, d = 4`. Smaller shapes need
  proportionally less.

### Dependencies

| Crate | Version | Purpose |
|---|---|---|
| `sha2` | 0.10 | SHA-256 behind the tweakable hash and the PRF |
| `rand` | 0.8 | key and salt sampling; `StdRng` for seeded leaf regeneration |
| `num-bigint` | 0.4 | `C(261,123)` is a 257-bit integer, so `g_λ` decodes in big integers |
| `num-traits` | 0.2 | `Zero`/`One` for the binomial table |

### As a library

```toml
[dependencies]
tlots = { path = "../tlots" }
```

## Basic usage

The full version of the following, which compiles and runs, is
`examples/simple.rs` (`cargo run --release --example simple`).

**Dealer setup.** `(e, κ) = (261, 123)` is the paper's CFF-LOTS parameter
pair; `h` and `d` set the hypertree's height and layer count, with `d | h`.

```rust
use tlots::hypertree::{self, Dealer, HyperParams};
use tlots::lots_cff::{self, LotsParams};

let params = HyperParams::new(8, 2, LotsParams::new(261, 123), id, master_seed)?;
let mut dealer = Dealer::key_gen(params);
// dealer.pk, dealer.id are the long-term public key pk = (root, id).
```

**Issue a leaf.** Each call hands out the next unused LOTS keypair with the
certificate chain binding it to the long-term public key.

```rust
let index = 0;
let (keypair, cert) = dealer.issue_next_leaf().unwrap();
```

**Key distribution.** Every secret component is Shamir-shared among `N`
participants. `vk` is the public verification-key vector behind identifiable
abort; every participant receives all of it.

```rust
let (t, n) = (5, 9);
let (shares, vk) = lots_cff::dist(&dealer.pp, index, &keypair.sk, t, n, &mut rng);
```

**Sign.** Each participant's contribution is independent of the others'.

```rust
let partials: Vec<_> = shares.iter().take(t)
    .map(|share| hypertree::part_sign(&dealer.params.lots, &dealer.pp, index, msg, share))
    .collect();

let sig = hypertree::combine(
    &dealer.params.lots, &dealer.pp, msg, t, index, cert, &partials,
).unwrap();
```

**Verify**, against the published public key alone.

```rust
let ok = hypertree::verify(&dealer.params, &dealer.id, &dealer.pk, msg, &sig);
```

**Identifiable abort.** If combining yields a signature that fails
verification, every participant publishes its opening and the check names
the culprits.

```rust
let proof = lots_cff::sig_proof(&dealer.params.lots, &dealer.pp, index, msg, share);
let honest = lots_cff::verify_partial(
    &dealer.params.lots, &dealer.pp, index, share.id,
    &vk[share.id as usize], msg, partial, &proof,
);
```

Dropping the participants that fail and recombining from the remainder gives
a signature that verifies, with no new round of partial signatures.

## Design notes

### Public parameters and tweaks (§5, Algorithm 4)

All hashing is
`thash(P, t, x) = SHA-256("TLOTS/TH" ‖ len(P) ‖ P ‖ len(t) ‖ t ‖ x)`.
The length prefixes make `(P, t, x)` uniquely decodable, and every tweak
(`hashing::tweak`) uses fixed-width big-endian fields, so no two uses in the
scheme share a `(P, t)` pair.

Equation (1) of §5 fixes the convention `H_X(t, m) = H(TLOTS.pk‖X, t, m)`,
and Algorithm 4 instantiates the LOTS with `prm = (pk^B‖"LOTS", pk^B‖"msg")`,
keyed on the long-term public key `pk^B = (Root, id)`, serialized here as
`Root ‖ id`:

| object | public parameter | tweak |
|---|---|---|
| LOTS key with global index `c` | `P_LOTS = pk‖"LOTS"` | `y_i = H(P_LOTS, c‖i, x_i)` |
| LOTS message hash | `P_msg = pk‖"msg"` | `B_msg = g(H(P_msg, c, msg))` |
| share commitment | `P_vk = pk‖"vk"` | `π_cjk = H(P_vk, c‖j‖k, sh_cj,k)` |
| participant verification key | `P_vk = pk‖"vk"` | `vk_cj = H(P_vk, c‖j, (π_cjk)_k)` |
| WOTS+ key at position `idx` of layer `r` | `P_WOTS = id‖"WOTS"` | chain step at height `j` of chain `i`: `r‖idx‖i‖j` |
| WOTS+ message hash | `P_Wmsg = id‖"WMSG"` | `D = H(P_Wmsg, r‖idx, m)` |
| WOTS+ secret values | secret seed `rho` | `x_i = PRF(rho, r‖idx‖i)` |
| WOTS+ tree `tau` of layer `r` | `id‖r‖tau` | leaf `0‖i`, node `1‖j‖i` |
| salted LOTS tree of a bunch | `id_msg` | leaf `0‖i` (on `salt‖pk`), node `1‖j‖i` |

The hypertree's own parameters cannot be keyed on `pk` without circularity,
since `Root` is the output of the top Merkle tree over WOTS+ public keys and
everything hashed while computing `Root` can therefore depend only on `id`.
Algorithm 1 keys those on `id` (`MerkleTree(id‖r‖tau, ...)`), as does this
implementation. `PubParams::new(id)` builds that half of the parameters;
`Dealer::key_gen` calls `bind_root` once `BKeyGen` produces the root, which
re-keys `P_LOTS`, `P_msg` and `P_vk` onto `pk = (Root, id)` before any LOTS
key is generated. A verifier rebuilds the same parameters from the published
`pk` with `PubParams::with_pk`. Under `d = 1` the LOTS tree is itself the top
tree, and the LOTS parameters remain keyed on `id` on both sides.

Merkle tweaks are Appendix A.3's `0‖i` and `1‖j‖i`, with node levels `j`
counted from the root (root = 0). No Merkle tree shares a public parameter
with an OTS or verification-key hash: Merkle parameters are 32 B (`id_msg`)
or 44 B (`id‖r‖tau`), the others 35–36 B (`id‖label`) or 66–68 B
(`pk‖label`). The tests in `hashing` check this. The paper names the LOTS
and WOTS+ message parameters both `P_msg`; since the two are keyed on
different values, the code labels them `"msg"` and `"WMSG"`.

Secret material is derived with the PRF: WOTS+ chain starts, LOTS leaf
seeds, per-bunch `id_msg` and per-leaf salts all use the dealer's master
seed `rho` as the PRF key under distinct tweak tags.

### Subtree regeneration

Every subtree, at every layer, is a function of the secret `master_seed` and
its address `(layer, subtree_index, position)`, via
`hashing::derive_leaf_seed`. Key generation and signing regenerate
`O(d · 2^(h/d))` leaves rather than the full `2^h`-leaf tree, so `h` and `d`
may be chosen at run time subject only to `d | h`, including `h = 64`.

What bounds the cost is that a WOTS+ public key does not depend on the
message it will eventually sign. Building a subtree's Merkle root — needed
for an authentication path, or as the child root that the layer above signs —
therefore takes `O(2^{h/d})` public-key derivations at that subtree, with no
recursion into child subtrees. Only the leaf on the signing path
additionally needs a signature, computed once its child's root is known.

`Dealer::resolve_leaf` rebuilds all `d` subtrees on a leaf's path on every
call, which suits arbitrary one-off access. `Dealer::issue_next_leaf` caches
the active bottom subtree and its chain of `d-1` ancestors across calls,
following Algorithm 1's `counter` pattern with a monotonic `next_leaf` index
that is never reissued, so drawing many leaves from one bunch pays the
`O(2^{h/d})` regeneration cost once per subtree rather than once per leaf.
Only public data is cached: leaf public keys, salts, `id_msg`, the root, and
each leaf's precomputed authentication path. No leaf's secret key is cached;
each is regenerated just in time for the leaf being issued, since building a
subtree's authentication paths does not require any other leaf's secret.
Memory is thereby bounded by one bottom subtree plus `d-1` ancestors,
approximately 0.9 GiB at `h=64, d=4`. A cache slot is replaced, freeing its
contents, when a requested leaf falls outside it.

`Dealer::key_gen` pre-populates the top layer's cache slot, since a
hypertree has exactly one top-layer subtree for its lifetime; without this
the first call to `issue_next_leaf` would rebuild it.
`Dealer::issue_next_leaf_verbose` reports per layer whether that layer's
subtree was rebuilt or reused on a given call.

### Layer numbering

Code layers run `0` (top) to `d-1` (bottom, LOTS). Algorithm 1 numbers the
WOTS+ layers `r = 1` (lowest) to `r = d-1` (top), so code layer `l` is
`r = d-1-l`. The top tree has parameter `id‖(d-1)‖0`, matching `BKeyGen`'s
`id‖h‖0` with `h = d-1`.

### Bottom layer (layer d−1)

Each leaf is a CFF-Lamport key whose tweak is its global index `c`, and
whose secret is Shamir-shared among N participants (`lots_cff::dist`). The
bottom Merkle tree is salted and hashed under a per-bunch `id_msg`. The
WOTS+ key above the bunch signs `m_0 = id_msg ‖ MT_0.Root` (Algorithm 1),
which authenticates `id_msg`. Under `d = 1` there is no WOTS+ layer, so the
bottom tree uses `id` and `verify` requires `id_msg = id`. A signature
carries the global index `c`; Algorithm 1's `counter` and the leaf's
position within its bunch are `c / L` and `c % L`, both derived from `c`
rather than read from an authentication path.

### Upper layers (0..d−2)

Tweaked WOTS+ keys (`wotsplus.rs`) with secret values `PRF(rho, t‖i)` from
the master seed. Each layer's WOTS+ public keys, raw and concatenated, form
a Merkle tree whose root is signed by the WOTS+ key one layer up. Base-`2^w`
encoding and checksum follow XMSS; hashing and addressing follow §2.1,
without XMSS ADRS or bitmasks. Of XMSS only WOTS+'s base-`2^w` encoding and
checksum are used.

### Signing and verification

`hypertree::part_sign` and `hypertree::combine` take the public parameters
and the LOTS key's global index. Combined partial signatures give the LOTS
signature `sigma^L`; `combine` attaches `pk_cmpl = (y_i)_{i ∉ B_msg}`, which
`TLOTS.Verify` accepts in place of the full LOTS public key, together with
the certificate chain, giving `Sigma = (sigma^L, pk_cmpl, index, cert)`.

`hypertree::verify(params, id, root, msg, sig)` derives the public
parameters from `pk = (root, id)`, recomputes the LOTS public key with
`LDerivePK` and the salted bottom root, then per layer the WOTS+ public key
with `WDerivePK` and the tree root with `DeriveRoot`. All tweaks and Merkle
leaf positions come from `sig.index`, never from the authentication paths;
out-of-range indices and malformed lengths are rejected.

### Identifiable abort (§3, §5, Algorithm 4)

Besides the share bundles, `lots_cff::dist` returns the public
verification-key vector `vk = (vk_cj)_{j ∈ [N]}`, where `vk_cj` commits to
each of participant `j`'s `e` shares of LOTS key `c` through the per-share
commitments `π_cjk`. If a round of partial signatures fails to reconstruct,
each participant publishes `lots_cff::sig_proof`, the commitments to the
`e − κ` shares it did not reveal, and `lots_cff::verify_partial` recomputes
`vk_cj` from the revealed shares together with the proof. An honest
participant passes this check; one that submitted a malformed partial
signature does not, and is named, dropped, and the remaining honest shares
recombined. Neither `vk` nor the proof is part of a signature: `vk` is
distributed once with the shares, and a proof travels only when a round
fails. Algorithm 6 indexes `π` by all of `[e]` and leaves the `B_msg`
positions empty; `sig_proof` returns only the `e − κ` entries that carry a
value, and `verify_partial` reads them back in the same ascending order.

### Field size and the cost of Dist

Shamir sharing (§2.4) uses `F_p` with §8's FFT-friendly prime
`p = 4101·2^20 + 1 = 4300210177`, for which `2^20 | p-1` and which exceeds
`2^32` by about 0.12%, so every field element still serializes into 5 bytes.
`lots_cff::dist` assigns participant `j` the field point `omega^j` for a
primitive `next_pow2(N)`-th root of unity `omega` (`shamir::eval_domain`)
rather than the integer `j+1`, so evaluating a secret's degree-`(T-1)`
polynomial at every participant's point is a forward NTT, which
`shamir::share_ntt` computes in `O(N log N)` in place of `O(N*T)` per-point
evaluation. `shamir::share_naive` retains the per-point version for
comparison and as a reference. The cost of `Dist` consequently depends on
`N` alone and not on `T`.

### The cost of Reconstruct

The Lagrange coefficients `lambda_{i,S}` of §2.4 depend only on the signer
set `S`, not on which coordinate of the shared vector is being recovered —
as Algorithm 2's `ShamirReconVec` writes it, `L ← Σ_i λ_{i,S}(x) L_i`, with
one coefficient per participant applied to that participant's whole share
vector. `shamir::lagrange_coeffs` therefore computes them once per quorum
and `lots_cff::combine` reuses them across all `κ·r` coordinates, for the
`Θ(T² + κ·r·T)` cost of §8 rather than `Θ(κ·r·T²)`.
`shamir::reconstruct` keeps the single-coordinate form, which is what the
tests interpolate against and what `bench` times in its last column to show
what the sharing saves.

Reconstruction cannot use the inverse NTT that `Dist` uses in the other
direction: an arbitrary quorum of `T` responders out of `N` is not a
root-of-unity subset, whereas `Dist`'s evaluation domain is fixed and
dealer-chosen. Fast multipoint interpolation via a subproduct tree would
give `O(T log² T)`; it is not implemented here.

Because the coefficients are shared, a quorum that repeats one
participant's field point would interpolate through a degenerate
coefficient set, so `combine` rejects such a quorum
(`shamir::has_distinct_points`), matching Algorithm 5's refusal to store two
partial signatures from one sender.

## Performance

### Parameters

| Parameter | Value | Source |
|---|---|---|
| `λ` (hash output) | 256 bits, SHA-256 | §8 |
| `(e, κ)` | (261, 123) | §8; minimality pinned by tests in `src/cff.rs` |
| `r`, `bpb` (Shamir blocks per preimage) | 8 blocks of 32 bits | §8 |
| `p` (Shamir prime) | `4101·2^20 + 1 = 4300210177` | §8 |
| Evaluation domain | `ω^j`, a `2^⌈log N⌉`-th root of unity | §8 |
| `h`, `d` (hypertree height, layers) | 64, 4 → `L = 2^16`, capacity `2^64` | §8 |
| `(w, ℓ1, ℓ2, ℓ)` (WOTS+) | (4, 64, 3, 67), base `2^w = 16` | §8 |

`(e, κ) = (261, 123)` is the smallest pair with `2^256 ≤ C(e,κ)`, realized in
`src/cff.rs` by combinatorial-number-system decoding of a 256-bit digest
(`num-bigint`, since `C(261,123)` is a ≈257-bit integer). Because
`2^256 ≤ C(e,κ)`, the reduction in `digest_to_block` never reduces a 256-bit
digest, so `g_λ` is injective as Definition 2 requires.

Since `p - 1 = 4101·2^20`, `F_p` has a multiplicative subgroup of order
`2^20`, which `shamir::share_ntt` uses to evaluate each secret's
degree-`(T-1)` polynomial at every participant's point by a radix-2 NTT. `p`
exceeds `2^32` by 0.12%, so a field element serializes into 5 bytes.

`h` and `d` are set at run time, with `d | h`; they are the first two
arguments to `bench`. The figures below are for `h = 64, d = 4`.

### Platform

| | |
|---|---|
| Machine | MacBook Pro `Mac14,7` (Apple M2, 4 performance + 4 efficiency cores) |
| Memory | 16 GiB unified |
| Toolchain | `rustc` 1.89.0, `--release` (`opt-level = 3`) |
| Threading | single-threaded throughout, so one core is used regardless of the machine |

The `sha2` crate is built here without its assembly backend, so SHA-256 does
not use the ARM crypto extensions.

### Signature and key sizes

Grouped as in §8's size table:

| Component | Formula | Size |
|---|---|---|
| `pk^B = (MT_h.Root, id_MT)` | 2·32 | 64 B |
| `σ^L, pk_cmpl` | e·32 | 8,352 B |
| `id_msg, MT_0.Salt, MT_0.Path, index` | (h/d + 2)·32 + h/8 | 584 B |
| `Σ`: d−1 pairs `(σ^W_r, Path_r)` | (d−1)(ℓ + h/d)·32 | 7,968 B |
| **Total** `\|σ\|` | | **16,904 B** |
| Per-participant Shamir share of one leaf (`sh_j`) | e·r·5 | 10,440 B |
| Verification key `vk_ij` (one leaf, one participant) | 32 | 32 B |
| Abort proof `π` (one leaf, one participant) | (e−κ)·32 | 4,416 B |

The second row is the κ revealed preimages (κ·32 = 3,936 B) together with
the e−κ off-block public-key hashes (4,416 B); a signature does not carry
the full LOTS public key, which `LDerivePK` recovers from those two. The
global LOTS index is stored as a `u64`, which is the `h/8` of the third row
at the reported `h = 64`.

The verification key and the abort proof lie off the signing path: `vk` is
public data distributed once at key distribution, and `π` is transmitted
only when a signing round fails. Neither enlarges `|σ|`.

### Round complexity

As in the paper: if all signing requests reach every participant, for
instance through the bulletin board, the scheme is **one round**, each of the
`T`-of-`N` participants producing a partial signature non-interactively for
assembly offline. With the majority-vote bulletin board it signs in **four
rounds**. A misbehaving participant adds one round, the identifiable-abort
opening measured below.

### Dealer costs

One-time and independent of `(T, N)`.

| Stage | Cost |
|---|---|
| `Dealer::key_gen` — top-layer subtree, `2^16` WOTS+ keys | 30.8 s |
| `Dealer::resolve_leaf` — one leaf's chain, uncached | 99.4 s |
| `Dealer::issue_next_leaf` — first call, three layers not cached by `key_gen` | 71.5 s |
| `Dealer::issue_next_leaf` — subsequent calls from the same bunch | 0.95 ms |

Within `resolve_leaf`, key derivation dominates and tree building is a
rounding error:

| Layer | Key derivation | Merkle tree |
|---|---:|---:|
| LOTS (bottom) | 11,448.5 ms | 1,920.1 ms |
| WOTS+ `r=1` | 28,106.8 ms | 491.7 ms |
| WOTS+ `r=2` | 28,531.4 ms | 493.8 ms |
| WOTS+ `r=3` | 27,839.1 ms | 497.4 ms |

The LOTS bunch is cheaper to derive than a WOTS+ layer — `e = 261` hashes
per leaf against `ℓ · (2^w − 1) = 67 · 15 = 1,005` chain steps — but its
Merkle tree costs about four times as much, since its leaves are the `e·32`
= 8,352-byte concatenated LOTS public keys rather than `ℓ·32` = 2,144-byte
WOTS+ ones.

The gap between 71.5 s and 0.95 ms is the bounded cache doing its job.
`issue_next_leaf` reports which layers it rebuilt: the first call rebuilds
the three layers `key_gen` did not cache (`[true, true, true, false]`) and
every later leaf of that bunch reuses all four (`[false, false, false,
false]`).

### Threshold signing

For a single leaf at `e = 261, κ = 123, h = 64, d = 4`. `Dist` is
Algorithm 4's `IssueShares` in full; the next column re-runs just the
verification-key commitments, estimating that portion of `Dist` so that the
Shamir sharing alone is roughly the difference. `PartSign` is the total over
all `T` signers. The last column is the same interpolation with the quorum's
Lagrange coefficients recomputed for every coordinate, which is what sharing
them across the `κ·r` coordinates saves; it is not on the signing path.

| N | T | Dist (ms) | of which vk (ms) | PartSign (ms) | Reconstruct (ms) | Verify (ms) | Total (ms) | per-coord. (ms) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 5 | 3.357 | 2.016 | 0.301 | 0.098 | 0.886 | 4.641 | 1.498 |
| 50 | 25 | 24.575 | 51.864 † | 0.275 | 0.144 | 0.806 | 25.798 | 14.044 |
| 100 | 25 | 43.956 | 18.991 | 0.280 | 0.218 | 0.749 | 45.204 | 9.315 |
| 100 | 50 | 39.158 | 18.952 | 0.583 | 0.233 | 0.748 | 40.721 | 30.020 |
| 400 | 100 | 150.998 | 73.311 | 1.002 | 0.530 | 0.808 | 153.337 | 103.107 |
| 1000 | 200 | 351.488 | 187.669 | 2.090 | 1.105 | 0.748 | 355.432 | 397.004 |

`Total` is `Dist + PartSign + Reconstruct + Verify`; the vk and
per-coordinate columns are diagnostics measured alongside and are not added
in. `|pk| = 64 B` and `|σ| = 16,904 B` for every row, and every signature
verified. Run-to-run variation on this machine is a few percent: three runs
gave 337.0, 350.1 and 351.5 ms for `Dist` at `N = 1000`.

† The vk figure at `N = 50` is a measurement artifact — it exceeds the
`Dist` call it is meant to account for part of, which cannot be right.
Interpolating the neighbouring rows puts the true value near 10 ms.

Reading the table:

* **`Dist` scales with `N` alone.** The NTT is `O(N log N)` and the vk
  commitments are `Θ(e·N)` hashes, both independent of `T`: `N = 100` costs
  39–44 ms at either `T = 25` or `T = 50`, a spread within run-to-run noise,
  while going from `N = 400` to `N = 1000` costs 2.3× more.
* **`Reconstruct` is linear in `T` to within noise**, growing 2.08× from
  `T = 100` to `T = 200`, as `Θ(T² + κ·r·T)` predicts over this range. The
  per-coordinate column, at `Θ(κ·r·T²)`, grows 3.85× over the same step. At
  `T = 200` sharing the coefficients is a **359× saving**, 1.105 ms against
  397.004 ms, and it leaves `Dist` as the dominant `(T,N)`-dependent cost by
  two orders of magnitude.
* **`Verify` stays under 1 ms** and `PartSign` under 2.1 ms, both regardless
  of `N`. `PartSign` is the total over all `T` signers, so one participant's
  own contribution is about 10 µs at `T = 200`: κ lookups into its share
  vector. Verification is a fixed walk up the hypertree.

### The abort path

At `N = 50, T = 25` with participant 7 submitting a corrupted partial
signature: the quorum containing it reconstructs a signature that fails
verification; all 50 participants produce their openings
(`TLOTS.SigProof`) in **4.6 ms** and checking all 50 against their
verification keys (`TLOTS.SignFailed`) takes **5.6 ms**; one participant, the
corrupted one, is blamed; and recombining from the honest remainder gives a
signature that verifies.

So identifiable abort costs one round and about 10 ms on the failure path,
and nothing at all on the common one. The `bench` binary prints this path at
the end of every run, `examples/simple.rs` walks through it, and the same
flow is a unit test in both `lots_cff` and `hypertree`.

## Source code organization

| Module | Paper section | Purpose |
|---|---|---|
| `src/cff.rs` | Definition 2, §2.1 | κ-uniform 1-cover-free family; combinadic bijection `g_λ` |
| `src/hashing.rs` | Appendix A.1, §5 (Eq. 1), Algorithm 4 | the tweakable hash `thash(P, t, x)`, the PRF, `PubParams` (`P_LOTS`, `P_msg`, `P_vk`, `P_WOTS`, Merkle parameters), fixed-width tweak encodings (`hashing::tweak`), and seed/`id_msg`/salt derivation |
| `src/shamir.rs` | §2.4, Algorithm 2 | (T,N) Shamir sharing over §8's FFT-friendly `p`, with an NTT-based `share_ntt` and a per-point `share_naive` for comparison, and Lagrange coefficients (`lagrange_coeffs`) shared across the coordinates of one quorum |
| `src/lots_cff.rs` | §2.1, §3, Algorithm 4 (`dist`), Algorithm 6 (`part_sign`, `combine`, `sig_proof`), §5 (`verification_keys`, `verify_partial`) | plain and threshold CFF-Lamport OTS |
| `src/wotsplus.rs` | §2.1, Appendix A.2 | WOTS+ (w=4, base 16) for the upper hypertree layers; its public key is one raw concatenated Merkle leaf (no L-tree) |
| `src/merkle.rs` | §2.2, Appendix A.3 | tweaked Merkle trees, salted at the bottom layer; `derive_root` takes the leaf position explicitly |
| `src/hypertree.rs` | §2.3 (Algorithm 1), §3 | seed-derived hypertree: `Dealer::key_gen`, `resolve_leaf`, `issue_next_leaf`, `part_sign`, `combine`, `verify` |
| `src/bin/bench.rs` | §8 | benchmark harness; `h` and `d` configurable at the command line |
| `src/bin/ntt_bench.rs` | — | timing comparison of `share_naive` and `share_ntt` |
| `src/bin/wots_stress.rs` | — | WOTS+ fuzz check over random parameters, tweaks and messages |
| `examples/simple.rs` | — | the end-to-end walkthrough above |

## Testing and reproducing the benchmarks

```sh
cargo test --release                      # 43 unit tests per layer, plus
                                          # end-to-end threshold sign/verify
                                          # tests covering identifiable abort,
                                          # parameter keying, h/d validation
                                          # and the d=1 case
cargo run --release --example simple      # the walkthrough above
cargo run --release --bin bench           # every table above (h=64, d=4 default)
cargo run --release --bin bench -- 32 4   # a faster shape
cargo run --release --bin wots_stress     # WOTS+ fuzz check (20,000 rounds)
cargo run --release --bin ntt_bench       # share_naive vs share_ntt timing
```

`bench` takes `h` and `d` as its first two arguments, with `d | h`; an
optional third selects the leaf index to sign, default `0`. One run prints
every figure above: the per-layer dealer costs, the cache behaviour of
`issue_next_leaf` over three successive calls, the size breakdown, the
`(T,N)` table, and the abort path. It asserts as it goes that each signature
verifies and that the serialized size matches the size breakdown, so a run
that completes has checked its own table.

The default shape allocates about 0.9 GiB for the subtree cache and spends
roughly three and a half minutes on the dealer stages before the `(T,N)`
table appears. `Dist`, `PartSign` and both reconstruction columns are
independent of `h` and `d`, so a smaller shape reproduces those rows in
seconds; `Verify` and `Total` are not, since verification walks `d` layers.

## Scope

The bulletin board (§6 and Appendix B, Algorithms 7–10) and the UC
simulation and security proofs are taken from the paper and are not
re-implemented. This crate implements and benchmarks the signing and
verification core of §2–§3 and §5, including the verification keys and
opening protocol behind identifiable abort, but not the message flow that
carries them between parties (`TLOTS.Sign2` and the broadcast logic of
`TLOTS.SignFailed`).

The dealer here plays the trusted key distributor that §5 assumes; no
distributed key generation is implemented, though the paper's use of an
ideal key-generation functionality means one could be substituted.

## Citation

```bibtex
@misc{tlots,
  author = {Fadavi, Mojtaba and Jao, David and Jaques, Sam},
  title  = {Threshold Signatures with Identifiable Abort from One-way Functions},
  year   = {2026},
  note   = {Implementation: \url{https://github.com/}}
}
```

## License

To be chosen before release.

## References

* [AABB+24] Aguilar-Melchor, C., Albrecht, M. R., Bailleux, T., Bindel, N.,
  Howe, J., Hülsing, A., Joseph, D., Manzano, M.: "Batch Signatures,
  Revisited." CT-RSA 2024. https://doi.org/10.1007/978-3-031-58868-6_7
