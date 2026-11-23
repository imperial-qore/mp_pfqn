# automom - run the fastest exact MoM-family solver for a model

```bash
make -C automom                           # plain C
./bin/automom -t models/09_asymmetric.qn  # same output flags as every solver
./bin/automom --explain -q model.qn       # also show the candidates and predictions
./bin/automom --dry-run model.qn          # only show the choice
```

## Options

| flag | meaning |
|---|---|
| `-l -e -g -t -q` | output, passed to the chosen solver (all candidates accept them) |
| `-d` | passed through; the exact-rational format is the chosen solver's own |
| `-j n` | threads, given to the solvers that take `-j` |
| `--explain` | candidates, predicted times and choice on stderr |
| `--dry-run` | the same, without running anything |
| `--solver NAME` | force one candidate |
| `--features` | the cost features of the model (used by `calibrate.py`) |

## How it chooses

Candidates: `ca`, `mom`, `comom`, `gcomom`, `gmom`, `binmom`, `bincomom`,
`mommod`, `comommod`. For each one the work W is computed in closed form
from the solver's own basis sizes, which are exact (checked against
`binmom -v`/`bincomom -v`):

| solver | W |
|---|---|
| ca | prod_r (N_r+1) R sum_k mi_k |
| mom, mommod | sum_r d_m(r)^2 M + N_r d_m(r) (M+1), d_m(r) = r [C(M+r-1,r) + C(M+r-2,r-1)] |
| comom, comommod | sum_r d_c(r)^3 + N_r d_c(r)^2, d_c(r) = (M+1) C(M+r-1,r-1) |
| binmom, bincomom | as mom/comom, with d^3 log N_r for a class long enough to be split |
| gcomom | sum_r U_r (N_r + R^2) R, U_r from its sub-model family (built in full) |
| gmom | sum_r N_r sum_m colA(m,r)^3, colA = r C(m+r-2,r-1) |

The modular solvers multiply W by the prime count (about 2b/61). The
predicted time is `exp(c0) W^c1 b^c2`, where b estimates the bit size of G
(`Ntot log2(Lmax sum mi + Zmax)`) and (c0, c1, c2) per solver come from
`calibrate.py`: the formulas give how the work scales, while the constants
convert it to seconds (rational vs word arithmetic, coefficient growth in
exact elimination, per-solver overheads).

Domain filters: `gmom` only for R >= 2 and all `mi = 1`, and always run with
`-X` so that it fails rather than returning an approximate perturbation;
`binmom` not with mixed zero/non-zero think times; `mommod` not with `-d`.
The candidates are run in order of predicted time until one succeeds; a
refusal (binmom or gmom on a singular system) moves on to the next, and the
output of a failed attempt is discarded.

Outside the MoM family: closed load-dependent models go to `comomld` then
`gld` for `-l/-e/-g`, and to `mvaldmx` for `-t/-q`; mixed models to `mvamx`,
or `mvaldmx` if they also have a MU section.

## Calibration

```bash
python3 automom/calibrate.py --workdir /scratch/bitbucket/automom-cal -p 4
make -C automom
```

Random non-degenerate models (R = 1..5, M = 1..24, four population
profiles, half with think times, 30% with `mi` in {1,2}); each solver
single-threaded with a 20 s cap, a chain stopping at its first timeout.
The fit is least squares on log t, written to `calib.h`. Calibrated on
picard09 (20 cores, shared): 1333 runs, median prediction error 1.5x to 2.3x
by solver. On the grid itself the predicted-fastest solver was the fastest
on 72/164 models, with picked/best time 1.06 at the median and 1.95 at the
90th percentile.

The predictions assume a non-singular model: a singular one adds the cost
of the chosen solver's exact fallback, which the model does not price.
