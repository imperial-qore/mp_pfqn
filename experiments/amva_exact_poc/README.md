# AMVA-assisted exact reconstruction: proof of concept

**Exact recovery works on the six test models, but AMVA initialization did not
reduce the correction counts compared with a zero initial guess.** Closing the
equations with an exact single-class boundary was necessary. Binary64 correction
solves stalled on the ill-conditioned case; an 80-decimal-digit correction solver
then recovered a certified rational answer in two further corrections.

This is a standalone experiment, with no changes to production solvers and no
third-party Python dependencies. It implements actual two-class CoMoM CE/PC
equations, not an unrelated linear system.

```bash
python3 experiments/amva_exact_poc/reproduce.py --check-ca \
  --output experiments/amva_exact_poc/results.json
```

Use `--case accurate_amva` or another case name for one model. `--check-ca`
requires the existing `bin/ca`; omitting it still performs independent exact
convolution in Python. Exit status is nonzero if any run remains uncertified
after the prescribed correction attempts and higher-precision retry.

## The exact system

Scope: two queues, two classes, positive demands, distinct class-1 demands,
single-server queues, nonnegative think times, and populations `(p,q)` with
`p >= 2`, `q >= 1`. In this experiment all input parameters are integers.

At each class-2 population `n`, carry the nine CoMoM basis entries

```
G(p-t,n), G^(+1)(p-t,n), G^(+2)(p-t,n),  t = 0,1,2.
```

All entries are divided by the SAME unknown target constant `G(p,q)`. Writing
the normalized entries as `g_t(n), h_kt(n)`, the equations are:

```
CE, t=0,1; k=1,2:
  h_kt(n) - g_t(n) - D_k1 h_k,t+1(n) - D_k2 h_kt(n-1) = 0

PC for class 1, t=0,1:
  (p-t) g_t(n) - Z1 g_t+1(n) - sum_k D_k1 h_k,t+1(n) = 0

PC for class 2, t=0,1,2:
  n g_t(n) - Z2 g_t(n-1) - sum_k D_k2 h_kt(n-1) = 0.
```

These are the two-class specialization of `comom/setupls.c`. The six mean
metrics are linear functions of normalized entries:

```
X1 = g_1(q),             X2 = g_0(q-1)
Q_k1 = D_k1 h_k1(q),     Q_k2 = D_k2 h_k0(q-1).
```

The first test leaves the preceding population's basis unknown. For the five
test models with positive think times, exact rank calculations find:

| Equations available | Unknowns | Rank | Nullity | Undetermined metric directions |
|---|---:|---:|---:|---:|
| Last population step + target normalization | 18 | 10 | 8 | 4 |
| Last three steps + target normalization | 36 | 28 | 8 | 4 |

In the zero-think-time model, the coefficient ranks and nullities are the same,
but there are three undetermined metric directions with one step and two with
three steps. Extra equations remove some ambiguity here, but do not close it.

The report includes explicit exact nullspace witnesses: `H v = 0` but
`C v != 0`. Thus nearby equations do not determine the metrics, even when the
target normalizing constant is set to one. AMVA observations are NOT inserted
as purportedly exact constraints.

For closure, exact one-class convolution computes the nine boundary entries
`b0` at `n=0`. Set `v0 = s*b0`, with `s = 1/G(p,q)` an additional unknown, and
impose `g_0(q)=1`. The resulting global system has `9q+1` unknowns and equations.
Its rank is checked exactly and is full in all six cases. This system includes
every class-2 step back to the known boundary; it is not a compact final-step
shortcut.

## Approximate initialization and exact certification

Three initializations are compared with exactly the same coefficient matrix,
equilibration, correction factorization, and reconstruction schedule:

1. **Zero:** no approximate solution information.
2. **Target AMVA only:** target normalization and the six target metrics are
   initialized from Bard-Schweitzer AMVA; all other unknowns start at zero.
3. **Nearby AMVA:** initialize the full basis using AMVA at neighboring
   populations, products of throughput estimates, and
   `G^(+k)/G = 1 + sum_r Q_kr`. The target metric entries are overwritten with
   the same target estimates used in initialization 2.

Initialization 3 needs additional AMVA model evaluations (9 to 36 here); it
does not assume that means for the original network determine higher moments.
AMVA uses the class-specific Bard-Schweitzer arrival approximation
`Q_kr(N-e_r) ~= (1-1/N_r) Q_kr(N)`, with other classes unchanged, damping 0.5,
and convergence tolerance `1e-13`.

The correction loop:

- Stores the current estimate as exact rationals.
- Computes `r = d-H*y` in exact arithmetic.
- Solves `H*delta ~= r` with a reused, equilibrated floating LU factorization.
- Adds the returned finite-precision correction exactly to the rational estimate.
- Tries continued-fraction reconstruction with denominator caps `2^8, 2^16,
  ..., 2^512`. A cap is only a candidate-search parameter, not an assumed bound
  on the answer.
- Accepts a candidate ONLY when `H*y=d` holds exactly. The exact full-rank
  certificate then proves uniqueness.

After 16 unsuccessful binary64 corrections, the experiment retries from the
current iterate with a predetermined 80-decimal-digit LU factorization. It
still requires the same exact certificate; a small residual is never sufficient.

The final multiclass reference is not computed until all reconstruction attempts
are over. Neither exact final metrics, exact final basis values, nor reference
denominator sizes influence initialization, precision choice, or stopping.
Only the legitimate single-class boundary is exact solver input.

This follows the general numeric-symbolic strategy of approximate correction,
exact residuals, and certified reconstruction; it is not a claim that AMVA can
solve arbitrary signed residual equations. See [Saunders, Wood and Youse,
*Numeric-Symbolic Exact Rational Linear System Solver*, ISSAC 2011](https://www.eecis.udel.edu/~youse/itersolve.pdf).
The AMVA formula is also documented in the [LINE Bard-Schweitzer implementation](https://line-solver.sourceforge.net/doxygen-cpp/pfqn__bs_8h.html).

## Observed results

The AMVA error is the largest relative error among the two throughputs and four
queue lengths. Correction counts exclude candidate reconstruction attempts.

| Case | Populations | System size | Initial AMVA error | Binary64 corrections: zero / target / nearby | 80-digit retry |
|---|---|---:|---:|---|---|
| small | (3,3) | 28 | 11.08% | 1 / 1 / 1 | unnecessary |
| medium | (6,6) | 55 | 12.29% | 3 / 3 / 3 | unnecessary |
| larger | (12,12) | 109 | 9.50% | 7 / 7 / 7 | unnecessary |
| accurate_amva | (6,6) | 55 | **0.816%** | 6 / 6 / 6 | unnecessary |
| low_delay | (8,8) | 73 | 9.87% | 3 / 3 / 3 | unnecessary |
| close_class1 | (3,6) | 55 | 16.75% | all uncertified after 16 | **2 / 2 / 2** |

The first four cases have demands `[[1,2],[3,1]]`. Think times are `[3,2]`
except for `accurate_amva`, which uses `[100,80]`. `low_delay` has demands
`[[10,5],[5,9]]` and zero think times. `close_class1` has demands
`[[1000001,5],[1000000,9]]` and think times `[1,1]`.

The stressed case matters: with nearby AMVA initialization, the final binary64
residual was about `1.55e-31`, yet the largest relative metric error was about
**30.3**. It was rejected. After two 80-digit corrections, exact reconstruction
and verification succeeded. An accurate approximate initial solution does not
make an ill-conditioned correction operator reliable.

For example, the small case recovers the exact throughput
`X1 = 13877/68084`; every metric and basis entry is recorded or checked exactly.
The largest reduced denominator in the normalized bases ranges from 17 to 99
bits across these tests.

## Validation and limits

- All 18 initialization/model combinations eventually produced certified exact
  solutions; all reconstructed basis entries matched independent positive
  convolution in Python.
- `--check-ca` verified seven normalizing constants per model against `bin/ca`
  (42 exact comparisons), covering the target and every throughput/queue-length
  numerator. The binary's decimal metric output is not used as an exact oracle.
- Rank calculations and explicit nullspace witnesses distinguish local
  underdetermination from globally unique reconstruction.

This establishes feasibility of certified correction for this restricted
CoMoM family. It does not establish a speedup, a general AMVA-to-exact conversion,
or convergence with fixed correction precision. The implementation assembles
the global system densely for floating factorization and pays for exact rank
audits. Its cost should not be compared directly with an optimized modular
solver. No MoM or GMoM implementation was changed or benchmarked.

The useful outcome is that exact reconstruction can be certified without using
the exact answer as a seed, while the tested AMVA estimates did not remove the
need for full closure or reduce the observed correction counts.
