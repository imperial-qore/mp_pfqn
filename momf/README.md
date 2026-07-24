# momf - fixed-precision MoM, with Wilkinson/Moler iterative refinement

`momf` is a floating-point twin of `mom`. It runs the same Method of
Moments recursion on the same linear systems, but carries every quantity
as an MPFR float at a run-time settable working precision instead of an
exact rational. Its purpose is to answer one question: **is the exact
rational arithmetic of the SIGMETRICS'06 algorithm actually necessary, or
would double precision plus iterative refinement do?**

The answer measured here is that the exact arithmetic is necessary, and
that iterative refinement does not substitute for it.

## Build and use

```bash
source ./setup-env.sh
make -C fpla && make -C momf     # produces bin/momf
```

`momf` accepts every `mom` option (`-l -g -t -q -d -b -p -s`) plus:

| flag | meaning | default |
|---|---|---|
| `-P bits` | working precision | 53 (IEEE double mantissa) |
| `-W bits` | residual precision used by refinement | `2*P` |
| `-I n` | max Wilkinson/Moler sweeps per block solve | 0 (refinement off) |
| `-T tol` | stop refining below this relative correction | 0 (run full budget) |
| `-D` | per-step recursion diagnostics and block condition numbers, on stderr | off |

`-e` (exact numerator/denominator) has no meaning in floating point and
falls back to printing `G` as a double.

MPFR is used rather than `double` or GMP `mpf` on purpose. At `-P 53`
MPFR reproduces IEEE double *precision* exactly (round-to-nearest, 53-bit
significand) while leaving the exponent range effectively unbounded. That
separates the two things the paper's Section 5.2 conflates: the *range*
problem, which is that `G` overflows a double, and the *precision*
problem, which is that roundoff accumulates along the population
recursion. Everything below is a measurement of precision alone.

## Reproducing the experiment

`sweep.py` scales the population of a base model, runs `mom` and `momf`
on each instance, and reports the floating-point error against the exact
answer:

```bash
python3 momf/sweep.py --base <model.qn> --scale 1 2 3 4 5 --momf-args="-P 53"
```

Metrics: `dlogG` = |log G_fp - log G_exact| (the relative error of `G`),
`eX` = max relative error over class throughputs, `eQ` = max relative
error over per-station-per-class queue lengths.

The model used below is Table 1 of the paper: `R=6` classes, `q=3` queue
types with multiplicities 12/5/8, integer loadings, scaled from 1 to 5
jobs per class.

## Result 1: double precision is unstable, as the paper claims

| jobs | eX at 53 bits | eX at 106 bits | eX at 212 bits |
|---|---|---|---|
| 6 | 2.2e-15 | 3.5e-16 | 3.5e-16 |
| 12 | 1.0e-10 | 3.7e-16 | 3.7e-16 |
| 18 | 2.2e-03 | 2.0e-16 | 2.0e-16 |
| 24 | 1.4e+01 | 2.3e-14 | 1.9e-16 |
| 30 | 1.2e+01 | 2.8e-08 | 8.6e-16 |

At 53 bits the throughputs are meaningless past 18 jobs and the queue
lengths are wrong by three orders of magnitude. This reproduces the
paper's statement that the floating-point technique "typically fails if
the network contains more than a few tens of jobs".

The failure is roundoff and nothing else: doubling the working precision
buys roughly a doubling of the reachable population. Bits lost are
therefore *linear* in the population, about 2.8 bits per job on this
model. That is the signature of a relative error amplified by a constant
factor at every population step, not of a fixed per-step loss.

## Result 2: Wilkinson/Moler refinement does not fix it

`-I 8 -T 1e-15` refines each diagonal-block solve with the residual
evaluated at 106 bits, and converges: the loop stops after 1.2 sweeps per
solve on average, i.e. the relative correction is already below 1e-15.
Raising the residual precision to 424 bits and the tolerance to 1e-30
changes nothing.

| jobs | 53 bits, no refinement | 53 bits + refinement |
|---|---|---|
| 6 | 2.2e-15 | 3.5e-16 |
| 12 | 1.0e-10 | 1.1e-11 |
| 18 | 2.2e-03 | 1.4e-04 |
| 24 | 1.4e+01 | 1.6e+01 |
| 30 | 1.2e+01 | 1.2e+01 |

Refinement buys a uniform factor of roughly ten, one decimal digit, and
leaves the growth rate untouched. The breakdown point moves by less than
one job per class.

## Why refinement cannot help

Run with `-D`, the solver reports the infinity-norm condition number of
each diagonal block and the composition of each right-hand side. On the
Table 1 model:

```
diag class r=3 block h=2 order=6   cond_inf=5.6e+04
diag class r=4 block h=2 order=12  cond_inf=2.6e+06
diag class r=5 block h=3 order=30  cond_inf=1.0e+07
diag class r=6 block h=3 order=60  cond_inf=4.5e+07
```

and, at every population step,

```
diag n=(1,1,1,1,1,1) |g|max=1.045e+13 |g|min=1.555e+10 range=6.7e+02 cancel=9.4e-01
```

Two facts follow. First, there is no catastrophic cancellation in forming
`b = B1r*g - A12*G` (the cancellation ratio stays below 1) and no large
dynamic range in the unknown vector (about 1e3), so the loss is not a
scaling or accumulation defect that extended-precision accumulation could
repair. Second, the diagonal blocks are ill conditioned, up to 4.5e7 for
the order-60 block of the last class.

Iterative refinement makes each solve backward stable: the computed `x`
solves a nearby system exactly, and its forward error relative to the
true solution *of the given right-hand side* falls to the working unit
roundoff. What it cannot do is undo the amplification of error that is
*already present in that right-hand side*. The right-hand side at
population `n` is built from the solution at population `n-1`, so an
inherited relative error `eps` emerges as up to `cond(C) * eps`. The
recursion applies `C^{-1}` once per population step, and the amplification
compounds over all `N` steps.

Refinement removes the fresh roundoff of the current solve, which is a
one-off gain of `cond(C)*u -> u`, i.e. the observed single digit. It
cannot change the compounding factor. The clearest evidence is the
converged run itself: every linear system is solved to full working
precision, and the answer is still wrong.

```
$ ./bin/momf tab1_n3.qn -P 53 -I 5 -T 1e-15   # refinement converged
X[1] = 6.228451200152182e-03
$ ./bin/mom  tab1_n3.qn -t                    # exact
       6.228391782095788e-03
```

## Consequences

1. The exact rational arithmetic in `mom` is not defensive engineering.
   It is required by the forward conditioning of the population
   recursion, which needs a number of bits growing linearly in `N`.
2. Mixed-precision iterative refinement is the wrong remedy, because the
   instability is not backward instability of the linear solves.
3. What remains open is whether a *reformulation* of the recursion has a
   better-conditioned forward map, for instance by carrying normalized
   unknowns whose ratios are the quantities actually consumed, rather
   than the normalizing constants themselves. That is a change to the
   algorithm, not to its arithmetic.
4. Where exactness is kept, the cost can still be reduced by changing the
   exact solver rather than the precision: Dixon p-adic lifting with
   rational reconstruction, or modular arithmetic with CRT, instead of
   rational Gaussian elimination whose operands grow throughout the
   elimination. The coefficient matrix does not depend on the population
   of the recursion class, so a single factorization or lifting setup
   serves all of its populations.

## Note on an output-path defect in `mom`

`mom -l` converts `G` to a double before taking its logarithm, so it
prints `inf` for any model beyond a few hundred jobs even though the
rational `G` it computed is exact. `momf` takes the logarithm inside
MPFR and is unaffected. `sweep.py` therefore reads `mom -e` and takes the
logarithm on the integers. This is incidentally a small illustration of
the range/precision distinction above: the exact solver loses the
answer's magnitude at the very last step, in the one place where it
stopped being exact.
