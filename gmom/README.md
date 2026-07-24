# gmom - Generalized Method of Moments (divide-and-conquer, b=1)

`gmom` implements the b=1 case of the generalized MoM (Casale, "A
Generalized Method of Moments for Closed Queueing Networks", Perf. Eval.
2011; the divide-and-conquer / multi-branch MoM). It is a faithful C port
of the reference MATLAB `mbmom1.m` / `setup1.m`, reusing the `gmpla`
exact-rational linear algebra and `util` combinatorics that back `bin/mom`
and `bin/comom`.

```bash
source ./setup-env.sh
make -C gmom              # produces bin/gmom
./bin/gmom models/09_asymmetric.qn --validate
./bin/gmom models/02_bottleneck.qn -e
```

## What it computes

The original MoM (`bin/mom`) recurses on the population alone. The
generalized MoM adds the **Convolution recursion**, which recurses on the
number of queues: the model on the prefix `{1..m}` is built from the model
on `{1..m-1}` plus one queue-removal (subtractive-convolution) branch.
`gmom` is the b=1 case - one queue removed at a time along the prefix
chain `1 -> 2 -> ... -> M`.

Per level the CE+PC system is overdetermined and is solved in exact
rational arithmetic by normal equations `(A^T A) x = A^T b` (eq 8 of the
paper), matching how the reference and `bin/procomom` solve it. The result
is the moment **basis** `V{M,l}`, from which `gmom` recovers the plain
normalizing constant `G(N)`, throughputs `X_r` and mean queue lengths
`Q_kr` by the mdecrease replica-descent (the `mom/mdecrease.c` post-step,
started one level lower because the generalized basis is more compact).

## Options

| flag | meaning |
|---|---|
| `-e` | exact `G(N)` numerator/denominator |
| `-g` | `G(N)` as a double |
| `-l` | log `G(N)` |
| `-t` | throughputs `X_r`, one per line |
| `-q` | mean queue lengths `Q_kr`, one row per station |
| `--validate` | check the whole basis, at every level, against exact convolution |

`G(N)`, `X` and `Q` are bit-for-bit identical to `bin/ca` on every
in-scope non-singular model (02, 03, 09, 13, lcfs_2class, lcfs_3class).

## Correctness

`--validate` recomputes every basis entry `V{m}` (all levels `m=1..M`, all
replica combinations, all class decrements) by independent Buzen
convolution of the corresponding replica-augmented model, and compares
bit-for-bit. It passes on every closed, non-singular model in `models/`
with `R>=2`:

```
02_bottleneck R=2   all 20 basis entries match
03_think      R=2   all 20 basis entries match   (Z != 0)
04_replicated R=2   all 12 basis entries match
09_asymmetric R=3   all 168 basis entries match  (M=6, Z != 0)
13_gld_small  R=2   all  6 basis entries match
14_ld_multi   R=2   all  6 basis entries match
15_repairman  R=2   all  6 basis entries match
lcfs_2class   R=2   all  6 basis entries match
lcfs_3class   R=3   all 12 basis entries match
```

## Scope and two findings from the port

- **Distinct single-server closed queues only.** Open/mixed models (a
  class with `N<0` / `LAMBDA`), load-dependent stations (`MU`), and
  multiserver/replicated stations (`mi>1`) are rejected with a clear
  message rather than mis-solved. The b=1 basis assumes `m=(1,...,1)`,
  the reference's stated domain; the base convolution is single-server.
- **Singular models.** Models whose MoM coefficient matrix is singular
  (equal-loading / demand relations: `05_sparse`, `06_large`,
  `08_multiclass`, the `test_singular*` family, ...) make the normal
  equations singular, exactly as they make `bin/mom` auto-perturb. `gmom`
  reports the singularity. This is the demand-related singularity of the
  MoM system, not a defect of the port. (Confirmed by construction: 09's
  demands validate; only specific demand relations are singular.)
- **The reference's new-class shift is Z=0-only.** Between classes the
  reference initialises the new class by a cheap component shift, and the
  MATLAB explicitly `error('Z<>0')` on that path. That shift zeroes the
  decrement component of the class just completed, which is wrong once
  that class has a think time. `gmom` instead initialises each new class
  exactly by convolution at the class boundary (`fill_level_all`): correct
  for all `Z`, and costing only `R-1` convolutions off the population
  sweep, so the recursion still carries the `O(N)` cost. With the shift,
  `09_asymmetric` (Z != 0) fails; with the exact init it validates.

## Performance measures (mdecrease)

`measures.c` recovers `G(N)`, `X_r` and `Q_kr` from the basis by the
replica-descent of the original MoM (a port of `mom/mdecrease.c` +
`perfindices`). Two adaptations: the descent starts at level `R-1` rather
than `R` (gmom's basis is one level more compact), and the class-`R`
decrement `G(N-1_R)` that the descent needs is supplied by capturing the
top basis one population step back (`grM`, the analogue of mom's `gr`).
Validated bit-for-bit against `bin/ca`.

## Files

- `base.c` - exact Buzen convolution `gmva` (the reference's base solver),
  valid for `Z=0` and `Z!=0`.
- `setup1.c` - the six per-level matrices `A,B,C,D,E,F` (port of
  `setup1.m`), with the overdetermined row count computed.
- `main.c` - the `mbmom1` recursion, normal-equations solve, exact
  new-class init, `grM` capture, CLI and validation harness.
- `measures.c` - mdecrease replica-descent and perfindices for
  `G(N)`/`X`/`Q`.
