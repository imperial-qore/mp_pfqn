# gcomom - CoMoM basis with multi-branch queue removal

`gcomom` is an exact solver for closed product-form networks. It keeps
CoMoM's population recursion on classes and adds a generalized-MoM-style
queue-removal recursion, but removes **r-1 stations per level** instead of
gmom's one. That brings the CoMoM decrement depth from M down to **1**: each
sub-model is solved by small independent exact systems instead of comom's
single dense per-class system.

```bash
source ./setup-env.sh
make -C gcomom            # produces bin/gcomom (plain C, no C++)
./bin/gcomom models/09_asymmetric.qn
./bin/gcomom -e models/04_replicated.qn
```

## Options

| flag | meaning |
|---|---|
| `-l` | log G |
| `-e` | exact G, numerator and denominator |
| `-g` | G as a double |
| `-t` | throughputs X_r, one per line |
| `-q` | mean queue lengths Q_kr, one row per station (total over its `mi` copies) |
| `-d` | metrics as exact rationals |
| `-X`, `--no-perturb` | fail on a singular system instead of running the exact fallback |
| `-j n` | threads (OpenMP): sub-models of one level in parallel (default 1) |
| `--stats` | sub-model count, depth, unknowns and timing on stderr |

## Algorithm

While class c is swept, classes s < c are at full population and class c
at n. A sub-model is a copy vector v (v_k <= mi_k, so multiplicity is
native: a station of multiplicity m is m identical single-server queues).
A fractional multiplicity e (`mom`-style rational `mi`, see the top-level
CLAUDE.md) splits into a base frac(e) that is never removed plus floor(e)
removable copies: the support counts stations with a base, depth 1 needs
R-1 stations with a removable copy, and the PC coefficients are frac(e_k) +
v_k. Sub-models with a fractional base never switch to the integer
convolution fallback (it adds whole copies); a singular one goes to
`util/singular.h`'s exact tiers instead.
For decrements d over the classes s < c, it carries `P[d] = G_v(n-d)` and
`Q[k,d] = G_{v+e_k}(n-d)` and moves them one step with:

| relation | equation | role |
|---|---|---|
| class-c PC | `n P[d] = Z_c P'[d] + sum_k v_k L_kc Q'[k,d]` | plain constants, explicit |
| CE | `Q[k,d] - sum_{s<c} L_ks Q[k,d+1_s] = P[d] + L_kc Q'[k,d]` | |
| PC, s < c | `sum_k v_k L_ks Q[k,d+1_s] = (N_s-d_s) P[d] - Z_s P[d+1_s]` | |
| removal of one copy of j | `Q[k,d] - sum_{s<c} L_js Q[k,d+1_s] = Q^{v-e_j}[k,d] + L_jc Q'[k,d]` | |

(' = previous step). Removal minus CE gives
`sum_s (L_ks - L_js) Q[k,d+1_s]`, so r-1 removals determine each station's
depth-1 values, and the PC rows supply the one row each removed station
lacks. Why fewer removals do not work: with b < r-1 removals the depth
must stay at m, CoMoM's own depth (rank argument, checked numerically for
r <= 5), so gmom's b=1 scheme gives nothing on the CoMoM basis.

- Sub-models with >= R-1 support stations run at depth 1 and remove one copy
  of each of their top R-1 support stations; the closure of that rule is
  the family. Smaller supports run at CoMoM depth |support|, with the
  removal rows as extra rows; every child then has depth >= parent - 1,
  which is what its removal rows read.
- A class idle at every support station adds the exact rows of its
  decoupled population constraint (`Z_s Q[k,d+1_s] = (N_s-d_s) Q[k,d]`, or
  `Q = 0` while the class has jobs if also `Z_s = 0`).
- The coefficient matrix depends on (L, v, Z, N) but not on the step. Once
  per class a full set of independent rows is chosen (mod a 61-bit prime, a
  second prime, then exactly over Q if both are inconclusive), the square
  system is split into its independent blocks (one per station outside the
  removal set, one coupled block for the removed stations), and each block
  is LU-factored exactly. A step is one back-substitution per block.

## Singular models

Equal demands of two stations in one class, proportional demand rows, zero
demands, or other loading relations can make a sub-model's system rank
deficient. The singular sub-models are almost always small (two or three
stations), so gcomom first handles them locally:

1. **Convolution sub-models.** If the population lattice prod_r (N_r+1) times
   (support+1) is at most 4e6, the singular sub-model's tables come from
   exact Buzen convolution over the lattice from that class on, in
   integers (`h(n) = G(n) prod n_s!`, so each queue copy adds
   `sum_s L_ks n_s h(n-1_s)`). Every other sub-model keeps the recursion.
   On `10_diverse` (stations 4 and 6 proportional) 14 of 677 sub-models go
   this way and the model takes 5.7 s; without it, it went to interpolation
   and ran for more than 7 minutes.
2. **`util/singular.h`**, as in comom, when the lattice is too large: other
   class orders, identical stations folded into one of summed `mi`,
   proportional classes aggregated, then the exact limit of a perturbation
   by interpolation. The interpolation runs over Q (comom uses comommod's
   modular engine for it), so on large singular models it is the slow path.

`-X` disables both and fails instead.

A class with jobs but no think time and no demand anywhere has G = 0;
gcomom rejects such models with a message instead of interpolating to 0.

## Validation

- `models/`: G exact and X, Q to 1e-12 against `bin/ca` on every closed
  non-LD model (see the table below).
- Random models (R <= 4, M <= 8, N_r <= 4, `mi` in {1,2,3}, think times,
  demands drawn from a small pool so that zero and equal demands are
  common): 293/293 well-posed models exact against `bin/ca` (G, and Q to
  1e-12), 123 of them through the singular fallback.

## Performance

Total unknowns per step equal comom's basis `M*C(M+R-1,R-1)` for R=3 and are
1-2x larger for R >= 4; the gain is structural. comom factors one dense
rational system of that size per class, gcomom only small blocks.

R=3, N=(2,2,N_R), random demands in 1..999, Z=0, wall-clock seconds on
picard09 (shared, load ~13 of 20 cores):

| M | N_R | comom | gcomom | gcomom route |
|---|---|---|---|---|
| 4 | 20 | 0.04 | 0.017 | direct |
| 8 | 20 | 0.69 | 0.090 | direct |
| 12 | 20 | 19.5 | 0.19 | direct |
| 12 | 40 | 21.3 | 0.25 | direct |
| 16 | 20 | 421 | 0.28 | class order (tier 1) |
| 16 | 40 | 509 | 0.60 | class order (tier 1) |
| 24 | 20 | not run | 2.1 | direct |

**Threads.** Setting up a class (row selection, block LU) reads no other
sub-model, so `-j n` does it for all of them at once; a step reads the
children's new tables, which have one copy fewer, so it runs the
sub-models of equal total copies in parallel, one level after the other.
Convolution sub-models are built in parallel too. Output is byte-identical
for every `n`. Wall clock with `-j 8`, same host under load: M=24 N_R=40
1.24 s -> 0.33 s, `10_diverse` 5.6 s -> 1.7 s; `06_large` (15 sub-models,
its time in one convolution) does not gain.

At M=16 two stations share a class-2 demand (387); the sub-model that removes
both has the single form `(L_3 - L_15)` with a zero class-2 component, so it
is singular, and tier 1 (another class order) solves it exactly at no
noticeable cost. Choosing the removed stations per sub-model to avoid such
pairs would keep these models on the direct path.
