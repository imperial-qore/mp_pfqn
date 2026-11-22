# safe_gmom - generalized MoM that always returns an exact answer

`bin/gmom` degrades silently on a model that is singular under every class
order: it warns and re-solves a model perturbed at digit 20, so `-g`, `-l`,
`-t`, `-q` return an APPROXIMATE result (with `-e` it refuses). `safe_gmom`
returns the exact `G(N)`, `X`, `Q` on those models too, applying the tiered
hybrid MVA/MoM construction of Casale, "CoMoM" (IEEE TSE 2009), Appendix A,
to gmom's b=1 basis. It is the generalized-MoM counterpart of `safe_mom` and
`safe_comom` and shares their tier driver (`util/safe_tiers.c`).

```bash
source ./setup-env.sh
make -C safe_gmom           # produces bin/safe_gmom
./bin/safe_gmom models/test_singular7.qn -e
```

Options are those of `gmom`: `-e` (exact `G` num/den), `-g`, `-l`, `-t`, `-q`.

## Tiers

Every tier is exact, so the result is always the exact answer.

| tier | method | handles |
|---|---|---|
| 1 | plain gmom (`bin/gmom -X`) | non-singular models inside gmom's domain |
| 2 | gmom with a different recursion class (O(R) search) | singularities that depend only on which class is processed last |
| 3 | exact convolution (`bin/ca`) | genuine loading degeneracies, and every model outside gmom's domain |

- **`gmom -X` (`--no-perturb`)** is added by this module: on a system that is
  singular under every class order, `gmom` reports it and exits non-zero
  instead of falling back to the perturbed re-solve. Without it, Tier 1 would
  pay for a full approximate solve before being rejected. An explicit `-p` is
  still honoured, so `bin/gmom -p` remains the way to get a perturbed answer
  (and `bin/gmom -p -b` to bracket it).
- **Tier 2 is expected never to fire.** `gmom` already performs the class
  search internally: its init-time oracle (`util/pfqn_sing.h`) is exact for
  the b=1 coefficient matrix, so `gmom` itself moves a degenerate class to the
  recursion position and un-permutes `X`/`Q` on the way out. Tier 2 is kept
  for uniformity with `safe_mom`/`safe_comom` and as a backstop if the two
  class searches ever disagree. This is the difference from `safe_mom`, where
  Tier 2 does real work (`08_multiclass`, `11_swapped`), because mom's basis
  has no exact loading-only oracle.
- **Tier 3 also covers gmom's domain restrictions**, not only its
  singularities. `gmom` is limited to closed models with `R>=2` and distinct
  single-server queues; `bin/ca` handles multiserver stations exactly by
  station expansion and single-class models directly, so `safe_gmom` answers
  the whole `models/` corpus while `gmom` alone rejects part of it.

## Validation

`safe_gmom -e` on all 28 models in `models/`: 26 exact and equal to `bin/ca
-e`, 2 reported as out of scope.

```
Tier 1 (11): 02 03 08 09 13_gld_small lcfs_2class lcfs_3class
             test_singular test_singular2
Tier 3 (15): 01_single 04 05 06 07 10 11 12 14_ld_multi 15_repairman
             lcfs_1class test_singular3..8
Tier 2 ( 0): never fires, as expected (gmom picks the class order itself)
out of scope (2): 16_mixed, 17_mixed_ld
```

Tier 3 carries most of the corpus because it absorbs gmom's domain
restrictions as well as its singularities: `01_single` and `lcfs_1class` are
`R=1`, `04`/`11`/`12` have multiserver stations, `14`/`15` are
load-dependent.

**Open/mixed models are reported, not answered.** `gmom` rejects them and
`bin/ca` is a closed-network convolution, so no tier is exact; `safe_gmom`
prints `open/mixed model is out of scope` and exits non-zero rather than
emit a bogus number. This differs from `safe_mom`, whose Tier 1 (`mom`)
handles open classes directly and never reaches that check. Use `bin/mom` or
`bin/safe_mom` for mixed models.

## Relation to the other wrappers

| wrapper | Tier 1/2 solver | Tier 2 does real work? | reason |
|---|---|---|---|
| `safe_mom` | `mom -X` | yes | mom's basis has no exact loading-only oracle, so orders must be probed |
| `safe_comom` | `comom` | yes | same, on the class-oriented basis |
| `safe_gmom` | `gmom -X` | no | gmom's b=1 oracle already picks the class order internally |

The payoff of `safe_gmom` over `gmom` is therefore Tier 3 alone: exactness in
place of perturbation on genuinely degenerate models, plus coverage of the
models outside gmom's domain. It is built for uniformity across the family,
not for speed.

## What is left

Same refinement as for `safe_mom` and `safe_comom`: the appendix peels only
the MINIMAL set of degenerate classes by convolution and keeps the moment
recursion for the rest, instead of falling all the way back to convolution in
Tier 3. That requires carrying a reduced basis through the recursion and is
future work; the tiered version here already guarantees exactness on every
model.
