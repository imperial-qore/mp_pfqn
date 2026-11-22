# safe_mom - MoM that always returns an exact answer

`bin/mom` degrades silently on models whose coefficient matrix is singular:
it warns and re-solves a model perturbed at digit 20, so `-g`, `-l`, `-t`,
`-q` return an APPROXIMATE result (with `-e` it refuses outright).
`safe_mom` returns the exact `G(N)`, `X`, `Q` on those models too, applying
the tiered hybrid MVA/MoM construction of Casale, "CoMoM" (IEEE TSE 2009),
Appendix A, to the MoM basis. It is the MoM-basis counterpart of
`safe_comom` and shares its tier driver (`util/safe_tiers.c`).

```bash
source ./setup-env.sh
make -C safe_mom            # produces bin/safe_mom
./bin/safe_mom models/08_multiclass.qn -e
```

Options are those of `mom`: `-e` (exact `G` num/den), `-g`, `-l`, `-t`, `-q`.

## Tiers

Every tier is exact, so the result is always the exact answer.

| tier | method | handles |
|---|---|---|
| 1 | plain MoM (`bin/mom -X`) | non-singular models |
| 2 | MoM with a different recursion class (O(R) search) | singularities that depend only on which class is processed last |
| 3 | exact convolution (`bin/ca`) | genuine loading degeneracies no class order removes |

- **`mom -X` (`--no-perturb`)** is added by this module: on a singular
  system `mom` reports the singularity and exits non-zero instead of
  falling back to the perturbed re-solve. Without it, Tier 1 and Tier 2
  would pay for a full approximate solve before being rejected. An explicit
  `-p` is still honoured.
- **Tier 2** is the cheap case of the appendix: the recursion class's
  loadings do not enter the coefficient matrix, so a singularity removed by
  processing a non-degenerate class last is found by trying each class in
  that position -- O(R) orders, not O(R!). Each order is probed with
  `mom -e -X` (fails at once when singular); only the winning order is
  re-solved with the requested flag. `G(N)` is class-order invariant, so
  `-e`/`-g`/`-l` need no remapping; `-t` and `-q` are mapped back to the
  original class order.
- **Tier 3** is the fully-peeled limit of the appendix's MVA-like
  recursion. Unlike `safe_gmom` would be, this keeps mom's full domain:
  `bin/ca` handles multiserver stations exactly, by station expansion.

The loading-only oracle in `util/pfqn_sing.h` is exact for gmom's b=1
matrix, not for mom's basis, so Tier 2 probes class orders by running `mom`
rather than screening them with the oracle.

## Validation

On every model in `models/`, for each of `-e`, `-l`, `-t`, `-q`:
`safe_mom` == `mom` where `mom` is exact, == `ca` where `mom` is singular.
No mismatches, no perturbed output.

```
Tier 1 (13): 01 02 03 04 09 13_gld_small 14_ld_multi 15_repairman
             16_mixed 17_mixed_ld lcfs_1class lcfs_2class lcfs_3class
Tier 2 (2):  08_multiclass 11_swapped
Tier 3 (13): 05 06 07 10 12 test_singular test_singular2..8
```

Open/mixed models (`LAMBDA`) are reported, not handed to the closed-network
convolution fallback.

## What is left

Same refinement as for `safe_comom`: the appendix peels only the MINIMAL
set of degenerate classes by convolution and keeps MoM for the rest
(reduced basis, extra recursion directions), instead of falling all the way
back to convolution in Tier 3. That requires modifying mom's recursion to
carry a reduced basis, and is future work; the tiered version here already
guarantees exactness on every model.
