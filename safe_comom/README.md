# safe_comom - CoMoM that always returns an exact answer

`bin/comom` errors out on degenerate models whose coefficient matrix is
singular (equal-loading / demand relations). `safe_comom` returns the
exact `G(N)`, `X`, `Q` for those models too, using the hybrid MVA/CoMoM
idea of Casale, "CoMoM" (IEEE TSE 2009), Appendix A: when CoMoM is
singular, interleave convolution so the degenerate directions are handled
by an MVA/convolution-style recursion.

```bash
source ./setup-env.sh
make -C safe_comom          # produces bin/safe_comom
./bin/safe_comom models/08_multiclass.qn -e
```

Options are those of `comom`: `-e` (exact `G` num/den), `-g`, `-l`, `-t`,
`-q`.

## Tiers

`safe_comom` orchestrates the exact solvers already in the tree; every
tier is exact, so the result is always the exact answer.

| tier | method | handles |
|---|---|---|
| 1 | plain CoMoM (`bin/comom`) | non-singular models |
| 2 | CoMoM after permuting the class order | singularities that depend only on which class is processed last (the recursion class) |
| 3 | exact convolution (`bin/ca`) | genuine loading degeneracies no class order removes |

- **Tier 2** is the cheap case of the appendix: a singularity removed by
  processing a non-degenerate class last. Permutations are probed with the
  fast `comom -e` path (which errors instantly on a singular system rather
  than running the slow perturbation fallback); only the winning
  permutation is solved with the requested flag. `G(N)` is invariant to
  class order, so `-e`/`-g`/`-l` need no remapping; `-t` and `-q` are
  mapped back to the original class order.
- **Tier 3** is the fully-peeled limit of the appendix's MVA-like
  recursion: every degenerate class swept by convolution. Always exact.

## Validation

Bit-for-bit against the reference solver on every closed model in
`models/` (safe_comom == `comom` when non-singular, == `ca` when
singular):

```
non-singular (Tier 1):  01 02 03 04 09 13 lcfs_2class lcfs_3class
recursion-class (Tier 2): 08_multiclass
genuine degeneracy (Tier 3): 05 07 test_singular test_singular2..5
```

Open/mixed models (`LAMBDA`) are out of CoMoM's closed-network domain and
are reported, not mis-solved by falling to convolution.

## Relation to the appendix, and what is left

This is the exact-but-tiered realisation. The appendix's efficiency
refinement is to peel only the **minimal** set of degenerate classes by
convolution while keeping CoMoM for the rest (reduce the basis by the `d`
linearly dependent constants, add recursion directions `R'`, recover the
dropped constants from `Lambda'(N-1_r)`), rather than falling all the way
back to convolution in Tier 3. That minimal-peel version is a further
optimisation over this one; it requires modifying comom's BTF recursion to
carry a reduced basis and multiple recursion directions, and is left as
future work. The tiered version here already guarantees the exact answer
on every model.
