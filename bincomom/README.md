# bincomom - CoMoM by binary splitting

`bincomom` computes exactly what `bin/comom` computes, with the same output
byte for byte, but multiplies the population steps of each class together
by binary splitting instead of applying them one at a time. It is the CoMoM
counterpart of `binmom` (see `binmom/README.md` for the method).

```bash
source ./setup-env.sh
make -C bincomom               # produces bin/bincomom
./bin/bincomom -e models/03_think.qn
./bin/bincomom -v -q model.qn  # per-class route, basis size, timing on stderr
```

Options: every `comom` option (`-l -e -g -t -q -d -X -p -s -j`); `-j n` also
threads the matrix products; `-v`; `--cross c` (split a class when
`N_r - M - 1 >= c*d^2`, default 2); `--leaf k` (steps multiplied directly at
a leaf, default 8); `--seq` (never split).

## How it is built

`bincomom` is `bin/comom` with one hook. Its Makefile compiles
`comom/main.c` with `main` renamed to `comom_cli_main` and links comom's
other objects and its `ccengine.o`. `bincomom`'s own `main` sets
`COMOM_ADVANCE` and calls `comom_cli_main`. Once the system of class `r`
is set up, comom calls the hook, which may move the basis from `n_r = 0`
to `n_r = N_r - (M+1)` in one go. The last `M+1` steps stay with comom,
because they feed the next class's initial basis and the final marginals.
Every report, flag and the exact singular-model fallback (tiers of
`util/singular.h`, whose derived models also go through the hook) are
therefore comom's own.

comom's step was moved, unchanged, into `comom_step()` (`comom/main.c`) so
that both comom's loop and the hook call it:

```
G(n)  = B2 g(n-1) / n
Gk(n) = A11^{-1} (B1 g(n-1) - A12 G(n))
```

Here A11 is solved by BTF or LU, with comom's handling of impossible
combinations and of non-filtered combinations. Every operation is linear
in `g`, and every pivot is chosen from the coefficients alone, so

```
n g(n) = (n P + Q) g(n-1)
```

with `P` and `Q` fixed within the class. `bincomom` reads them off the step
itself: column `j` of `P + Q` is `1*step(e_j, 1)`, and column `j` of
`2P + Q` is `2*step(e_j, 2)`, so `2d` calls of `comom_step` in all.
Building them through the step, instead of rebuilding them from the blocks,
keeps every special path of comom inside the matrices. A guard compares
`(P + Q) g(0)` with comom's own first step; on a mismatch the class falls
back to comom's steps with a warning (never seen in validation). The
product is then formed by `util/binsplit.h`, the code shared with `binmom`.

The step matrices are dense `d x d`, because `A11^{-1}` fills them in,
while comom steps on sparse blocks. So, as in `binmom`, a class is
multiplied out only when the population left for the hook is at least
`cross*d^2`.

Unlike `binmom`, `bincomom` has no refusals of its own. Mixed think times,
singular systems and every other case comom handles go through comom's
code.

## Validation

Byte-identical to `bin/comom` (stdout and exit code) on every run, both
with the default route and with every class forced through binary splitting
(`--cross 0`):

- 979 models (`models/` plus the generated test models), `-e -l -t -q`:
  3912 runs, identical in both modes, 4 timeouts (in comom too); 287 of
  these models reach the split path;
- 400 random models aimed at the special paths (impossible combinations,
  zero demands, replicated stations (`mi > 1`), zero and mixed think times,
  empty classes), with `--cross 0`:
  1543 runs identical, 57 timeouts; all 400 split, 185 went through the
  singular fallback.

The refactored `bin/comom` is byte-identical to the build before the
refactoring on the same 3912 runs.

## Benchmark

Exact `G` (`-e`), wall seconds, the 23 models of the `binmom` benchmark
(`RxMy_Nz`: `x` classes of equal population, `y` single-server stations,
total population `z`; random demands 1..1000), 300 s limit. Every output
equals comom's. Selected rows:

| model | bincomom | bincomom -j 8 | comom | binmom | mom | mommod | comommod |
|---|---|---|---|---|---|---|---|
| R1M3_N1000 | 0.020 | 0.018 | 0.248 | 0.020 | 0.302 | 0.574 | 0.428 |
| R1M3_N3000 | 0.049 | 0.085 | 4.18 | 0.037 | 5.08 | 8.27 | 7.36 |
| R1M3_N10000 | 0.396 | 0.347 | 131.0 | 0.154 | 96.3 | 93.6 | 82.6 |
| R2M3_N1000 | 1.34 | 1.49 | 1.41 | 1.32 | 1.31 | 1.21 | 1.72 |
| R2M3_N3000 | 1.15 | 0.421 | 24.7 | 1.24 | 20.9 | 11.5 | 15.0 |
| R3M3_N600 | 0.920 | 0.976 | 1.06 | 1.04 | 1.09 | 0.971 | 1.33 |
| R2M6_N1000 | 6.28 | 6.28 | 6.44 | 4.54 | 4.58 | 2.80 | 6.58 |
| R2M12_N120 | 0.268 | 0.263 | 0.266 | 0.132 | 0.209 | 0.176 | 0.525 |
| R4M8_N12 | 42.8 | 42.8 | 43.4 | 0.054 | 0.148 | 0.125 | 1.19 |
| R4M8_N30 | timeout | timeout | timeout | 0.097 | 0.178 | 0.370 | 28.0 |

Reading:

- small basis, large population: 330x faster than comom at `N = 10000`
  (`d = 4`), 21x at `R2M3_N3000` (59x with `-j 8`), and ahead of every
  solver but `binmom`;
- elsewhere it never splits and costs what comom costs: the steps it takes
  are comom's own. So it dominates comom (never slower beyond noise), but
  inherits comom's weak spots, such as the large CoMoM setup on `R4M8`
  (`R4M8_N12` 300x slower than the MoM solvers, `R4M8_N30` past the 300 s
  limit, as comom);
- against `binmom` it is about 2.5x slower at `R1M3_N10000`: its basis is
  the same size there, but its step matrices carry wider entries after
  comom's eliminations.
