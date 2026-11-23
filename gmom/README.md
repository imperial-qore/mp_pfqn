# gmom - Generalized Method of Moments (divide-and-conquer, MoM(Bmax))

`gmom` implements the generalized MoM (Casale, "A Generalized Method of
Moments for Closed Queueing Networks", Perf. Eval. 2011, the
divide-and-conquer / multi-branch MoM), MoM(Bmax) for any maximum branching
factor (`--bmax`, default 1). Bmax = 1 began as a C port of the reference
MATLAB `mbmom1.m` / `setup1.m`; Bmax > 1 follows Sec. 4.1.2-4.2 of the
paper (no reference code). It reuses the `gmpla`
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
With the default `--bmax 1` one queue is removed at a time along the prefix
chain `1 -> 2 -> ... -> M`; larger Bmax branches into several sub-models
per step (see "Branching factor" below).

Per level the CE+PC system is overdetermined and is solved in exact
rational arithmetic by normal equations `(A^T A) x = A^T b` (eq 8 of the
paper), matching how the reference and `bin/procomom` solve it. The result
is the moment **basis** `V{M,l}`, from which `gmom` recovers the plain
normalizing constant `G(N)`, throughputs `X_r` and mean queue lengths
`Q_kr` from replica levels 0 and 1 of the full basis, which are carried
through the recursion (see "Performance measures" below).

## Options

| flag | meaning |
|---|---|
| `-e` | exact `G(N)` numerator/denominator |
| `-g` | `G(N)` as a double |
| `-l` | log `G(N)` |
| `-t` | throughputs `X_r`, one per line |
| `-q` | mean queue lengths `Q_kr`, one row per station |
| `--bmax B` | maximum branching factor (default 1, the original gmom); see "Branching factor" below |
| `-j n` | threads (OpenMP): sub-models with the same queue count in parallel (default 1) |
| `--validate` | check both basis slots of every sub-model against exact convolution |
| `-p digit` | perturb demands/think-times at that decimal digit (approximate solution for degenerate models) |
| `-s seed` | RNG seed for the perturbation (default 23000) |
| `-b` | bracket a perturbed run by re-solving with `-eps` as well as `+eps` (requires `-p`) |
| `--force-perturb` | allow the perturbed solve under `-e` |

`G(N)`, `X` and `Q` are bit-for-bit identical to `bin/ca` on every in-scope
non-singular model (02, 03, 09, 13, lcfs_2class, lcfs_3class), and on the
auto-perturbing ones too (05_sparse, 07, 08, test_singular2/3/5/7/8). Two
auto-perturbing models differ in the last printed digit only
(`test_singular4` -t/-q, `test_singular6` -q): the digit-20 perturbation is
below double precision, so the difference is one ulp. `06_large` (auto-perturbed)
now finishes in 76 s at B = 1 and 1.8 s at B = M (was over 300 s), matching
`bin/ca` to the printed digit; `12_expanded` (M = 3, R = 10) is singular under every class order at every
Bmax (`-e` refuses), and its perturbed `-t` still exceeds 300 s at B = 1, 2, 3 = M.

Extra GCEs can also restore rank: `07_mixed` (zero demands) is singular at
B = 1 under every class order, but at B >= 2 it is solved exactly, `-e`, `-t`
and `-q` bit-for-bit equal to `bin/ca`, and validates. `safe_gmom` still
calls gmom at the default B = 1.

## Branching factor (`--bmax`)

`--bmax B` runs the paper's MoM(Bmax) with variable branching (Sec. 4.1.2
and 4.2). A sub-model with `d` queues adds `B(d) = min{d, B}` general
convolution expressions (GCEs, eq. 20), one for each of its last `B(d)`
queues; each GCE ties the sub-model to the one with that queue removed.
Every extra GCE lowers the basis level carried as unknowns, which is
`l(d) = max{1, r - B(d)}` while `r` classes are processed (Theorem 2). The
default `B = 1` is the fixed-branching prefix chain the original gmom
implemented (`l = R-1`); `B >= R-1` brings the full model to level 1, i.e.
`M*R` unknowns. `R = 2` has `l = 1` for every `B`, so `--bmax` only changes
anything from `R = 3` on.

The price is the number of sub-models: those with `d` queues number
`C(B(d+1)+M-d-1, M-d)` (the paper's `card(m(d))`), about `M^2/2` in all for
`B = 2` and every non-empty subset of the queues (`2^M - 1`) for `B = M`. Each
is solved at every population step, children before parents.

Measured (exact `G`, identical to `bin/ca` bit for bit; `Z = 0`, distinct
random demands, one core; `orig` is gmom before this change; the models are
`bench/m4r4.qn`, `m5r5.qn`, `m6r3.qn`, `m8r3.qn`, next to `bench/mixed_r4*.qn`
and `bench/nr0_last.qn`, which cover mixed think times at R=4 and empty
classes):

| model | orig | B=1 | B=2 | B=3 | B=M |
|---|---|---|---|---|---|
| M=4 R=4 N=12 each | 8.5 s | 0.15 s | 0.04 s | 0.03 s | 0.03 s |
| M=5 R=5 N=8 each | >300 s | 58 s | 10.2 s | 0.61 s | 0.20 s |
| M=6 R=3 N=(60,30,40) | 15.5 s | 0.34 s | 0.16 s | 0.25 s | 0.94 s |
| M=8 R=3 N=30 each | 13.1 s | 0.84 s | 0.12 s | 0.42 s | 1.86 s |

As the paper finds, B = M wins at M = R, where the basis shrinks most. With
M > R, levels stop falling once `B >= R-1` while the sub-model count keeps
growing, so B = R-1 is about optimal there.

Implementation (`setupb.c`, the model tree in `main.c`):

- `setupb()` builds one sub-model's system: GCE rows (for each tagged
  queue `j` and replica multiset `delta` with `delta_j = 0`), then CE and PC
  rows as in `setup1`. A GCE row's convolution term is component 0 of the
  child's level-`l` basis, which lives in the child's `vl` slot when the
  child carries level `l` and in its `vl1` slot when it carries `l+1`
  (Theorem 2 guarantees one of the two). With `B = 1` the rows are exactly
  `setup1`'s.
- Single-queue sub-models are advanced by their own level system: with no
  GCEs it has `r` rows and `r` unknowns at any level, square, as Sec. 4.2
  suggests. They used to be refilled by a lattice convolution per entry at
  every population step, which was 88% of the run time at B = 1. A queue
  with a zero demand in a class `s < r` makes that system singular; such a
  sub-model falls back to convolution.
- New classes are carried over exactly instead of by convolution (except
  class 2, whose lattice is one-dimensional): each basis keeps its `r-1`
  entries per combination and gains the new decrement `G(n-1_{r-1})` from a
  snapshot taken one step before the end of the previous sweep. When the
  level rises, the old level becomes `vl1` and the new `vl` solves the level
  system with `G(n-1_r) = 0`, the paper's reduced system. The paper's
  pseudocode has two index errors here (see the review of the paper), and
  the reference's shift initialisation is Z=0-only.
- The measures read levels 0 and 1 of the full model: directly from its two
  slots when `l = 1`, else from the lower-level carry (mixed think times) or
  the downward pass started at level `l`.
- The singularity screen tests every sub-model's matrix for every class
  count (`gb_singular`); `pfqn_sing.c` builds only the `B = 1` matrices, so it
  is still used for `B = 1`.

Every level matrix depends on the demands only, so its `A^T A` is formed
and LU-factored once per (sub-model, class) and reused along the class
sweep, instead of at every population step as before.

**Threads** (`-j n`). Setup and factorization of a class read no other
sub-model and run for all of them at once; a level solve reads the
children (one queue fewer), so the class start and each population step
run the sub-models of equal queue count in parallel, bottom-up. `nck()`'s
lazily filled table is warmed first, since it is not thread safe. Output is
byte-identical for every `n`. With `-j 8` (shared host): M=12 R=3 B=M
38.9 s -> 10.3 s, M=10 R=4 B=M 10.7 s -> 5.1 s, M=8 R=5 B=M 4.2 s -> 1.7 s.
`B = 1` is a chain, one sub-model per queue count, and does not gain.

## Perturbation

A perturbation replaces each demand `L_kr` by `L_kr*10^d + eps_kr` with
`eps_kr` a small integer, so the perturbed model is solved exactly and the
answer is descaled analytically: `G` is homogeneous of degree `sum(N)` in the
demands and is divided by `10^(d*sum(N))`, `X_r = G(N-1_r)/G(N)` picks up one
factor `10^d`, and in `Q_kr` the two powers cancel (`measures.c`). The
demands are therefore *not* restored before the measures are computed - they
must not be, since the basis was built from the perturbed model.

- **Default output** reports the model as read together with the `eps=`
  actually applied per class, via `printmodel_with_perturbation`, so a
  perturbed run is self-describing rather than only flagged on stderr. This
  covers the automatic digit-20 fallback as well as an explicit `-p`.
- **`-e`** normally refuses on a genuinely singular model, since a perturbed
  rational is not the exact `G`. `--force-perturb` prints it anyway, with a
  stderr note that the rational is exact for the *perturbed* model; useful for
  measuring the perturbation error against `bin/ca`.
- **`-b`** solves twice at the same digit and seed, once with `+eps` and once
  with `-eps`, and prints `midpoint [lower, upper] relhw`. `G` is a polynomial
  with non-negative coefficients in the demands, hence monotone in them, so
  `[G(-eps), G(+eps)]` **encloses** the exact `G(N)` (same for `log G`). `X`
  and `Q` are ratios of such polynomials and are not monotone in general, so
  their intervals are labelled as error indicators, not proven enclosures. A
  demand of `0` is left at `0` under `-eps`, since a negative demand is not a
  model. This is a strict improvement on `mom -b`, which is still the stub it
  announces ("bounds computation using dual perturbation is in development").

```
$ ./bin/gmom models/08_multiclass.qn -p 3 -b
G            4.766175608676675e+34  [4.763085381152656e+34, 4.769265836200694e+34]  relhw=6.48e-04
```
exact `G = 4.766087716788053e+34`, inside the bracket as expected.

## Correctness

`--validate` recomputes both slots of every sub-model at `N` (all replica
combinations, all class decrements) by independent Buzen convolution of the
corresponding replica-augmented model, and compares bit-for-bit. It passes
on every in-scope non-singular model in `models/`, at B = 1, 2 and 3:

```
                      entries (sub-models)
                      B=1        B=2        B=3
02_bottleneck R=2      28 (4)     60 (10)    86 (14)
03_think      R=2      28 (4)     60 (10)    86 (14)   (Z != 0)
09_asymmetric R=3     231 (6)    231 (21)   456 (41)   (M=6, Z != 0)
18_zero_first R=2      10 (2)     14 (3)     14 (3)
19_zero_last  R=2      10 (2)     14 (3)     14 (3)
20_zero_middle R=3     21 (2)     21 (3)     21 (3)
21_mixed_think R=2     10 (2)     14 (3)     14 (3)    (mixed Z: lower-level carry)
lcfs_3class   R=3      21 (2)     21 (3)     21 (3)
```

## Scope and two findings from the port

- **Distinct single-server closed queues only.** Open/mixed models (a
  class with `N<0` / `LAMBDA`), load-dependent stations (`MU`), and
  multiserver/replicated stations (`mi>1`) are rejected with a clear
  message rather than mis-solved. The basis assumes `m=(1,...,1)`, the
  reference's stated domain.
- **Init-time singularity check and auto-reorder.** `util/pfqn_sing.c` provides
  a shared oracle `pfqn_recursion_singular` / `pfqn_nonsingular_recclass`:
  the coefficient matrix depends only on the loadings (not the population),
  so singularity is decided once at startup by factorizing each level's
  `A^T A`. If the recursion would go singular, gmom moves the offending
  class to the recursion position (its loadings then leave the coefficient
  matrix) and solves in that order, un-permuting `X`/`Q` back; if no single
  class reorder removes it, gmom falls back to perturbation (see below).
  This replaces the old mid-recursion abort. NB the
  oracle is exact for gmom's Bmax = 1 matrix (Bmax > 1 screens its own
  matrices, `gb_singular`); mom/comom use a different basis
  whose degeneracy is stronger, so they must probe their own matrix.
- **Singular models.** Models whose MoM coefficient matrix is singular
  (equal-loading / demand relations: `05_sparse`, `06_large`,
  `08_multiclass`, the `test_singular*` family, ...) make the normal
  equations singular, exactly as they make the `bin/mom` system singular.
  `gmom` reports the singularity; `mom`, `momf`, `mommod` and `comommod`
  solve such models through the tiered fallback of `util/singular.h`,
  which `gmom` does not use. This is the demand-related singularity of the
  MoM system, not a defect of the port. (Confirmed by construction: 09's
  demands validate; only specific demand relations are singular.)
- **The reference's new-class shift is Z=0-only.** Between classes the
  reference initialises the new class by a cheap component shift, and the
  MATLAB explicitly `error('Z<>0')` on that path. That shift zeroes the
  decrement component of the class just completed, which is wrong once
  that class has a think time. `gmom` instead carries every basis over
  exactly, taking the decrement of the class just completed from a
  snapshot one step before its end (see "Branching factor"); correct for
  all `Z`. With the shift, `09_asymmetric` (Z != 0) fails; with the exact
  carry-over it validates.

## Performance measures (full basis)

The full model's two slots hold replica levels `l` (`vl`) and `l-1` (`vl1`),
`l = max{1, r - min{M, Bmax}}` (`r-1` under `--bmax 1`). When think times are
mixed, `main.c` also carries levels `0..l-2` of the full model (Casale,
Perf. Eval. 2011, Def. 1 and App. C; tables in `util/lvtab.h`). They are
carried over between classes like the two slots, and advanced along the class-`r` sweep by the class-`r` population constraint,
which holds for any think time. The class-`R` decrement `G(N-1_R)` comes
from snapshots of all levels one population step back (`pvl`, `pvl1`,
`lowr`). Otherwise the downward pass (`gmom_descend`) starts at level `l`. `measures.c` then reads `G(N)`, `X_r` and `Q_kr` directly from
levels 0 and 1.

This replaces the earlier `mdecrease` replica-descent, whose `Z_s=0` branch
used an identity that holds only when no class has a positive think time,
so it was wrong for models mixing zero and positive think times. Validated
bit-for-bit against `bin/ca`.

## Files

- `base.c` - exact Buzen convolution `gmva` (the reference's base solver),
  valid for `Z=0` and `Z!=0`.
- `setup1.c` - the six per-level matrices `A,B,C,D,E,F` (port of
  `setup1.m`), with the overdetermined row count computed.
- `setupb.c` - one sub-model's level system for `--bmax` (GCE, CE, PC rows).
- `main.c` - the MoM(Bmax) recursion over the sub-model tree, cached
  normal-equations factorization, exact new-class init, lower-level carry
  and `N-1_R` snapshots, singularity screen, CLI and validation harness.
- `measures.c` - `G(N)`/`X`/`Q` output from levels 0 and 1.
