# binmom - MoM by binary splitting

`binmom` computes exactly what `bin/mom` computes, with the same output
byte for byte, but multiplies the population steps of each class together
by binary splitting instead of applying them one at a time.

```bash
source ./setup-env.sh
make -C binmom               # produces bin/binmom
./bin/binmom -e models/03_think.qn
./bin/binmom -v -q model.qn  # per-class route, basis size, timing on stderr
```

Options: `-l -e -g -t -q -d` as `mom`; `-v`; `-j n` (threads for the matrix
products); `--cross c` (split a class when `N_r >= c*d^2`, default 2);
`--leaf k` (steps multiplied directly at a leaf, default 8); `--seq` (never
split).

## Why it can be faster

Within class `r` the MoM blocks `A11`, `A12`, `B1r`, `B2r` do not depend
on `n_r` (`mom` factors them once per class). One step is

```
G(n)  = B2r g(n-1) / n
Gk(n) = A11^{-1} (B1r g(n-1) - A12 G(n))
```

so with `g = [Gk; G]`

```
n g(n) = (n P + Q) g(n-1),   P = [A11^{-1} B1r; 0],  Q = [-A11^{-1} A12 B2r; B2r]
```

and, clearing the common denominator `D` of `P` and `Q`,

```
g(N_r) = (N_r P' + Q') ... (P' + Q') g(0) / (D^N_r N_r!),   P' = D P, Q' = D Q integer.
```

`mom` applies the `N_r` steps to a vector whose entries grow by
`Theta(log n)` bits per step: about `d^2 N^2 log N` bit operations per class
(`d` the basis size). The modular solvers do not change that order (`d^2 N`
per prime, `Theta(N log N)` primes). `binmom` forms the product of the
integer matrices in a balanced tree, so the large multiplications are
between operands of equal size and GMP's subquadratic multiplication
applies: about `d^omega M(N log N)` bit operations, quasi-linear in `N`. The
output itself (`G` has `Theta(N log N)` bits) bounds what any method can do
in `N`, so in `N` this is near optimal.

It loses when `d` is large against `N_r`: the products are dense `d x d`
with entries widened by `D` (40 to 140 bits here), while `mom` steps on
sparse blocks of about `d*k` entries, so the break-even population grows
like `d^2`. Measured (split vs sequential on one class, single thread):

| d  | N_r  | split   | sequential | |
|----|------|---------|------------|---|
| 18 | 1500 | 0.94 s  | 4.59 s     | split 4.9x faster |
| 28 | 1500 | 5.9 s   | 7.3 s      | about even (N_r = 1.9 d^2) |
| 48 | 1200 | 44.2 s  | 10.0 s     | sequential 4.4x faster |
| 48 | 2400 | 118.9 s | 54.8 s     | sequential 2.2x faster |
| 54 | 1500 | 60.1 s  | 19.8 s     | sequential 3.0x faster |
| 54 | 3000 | 208.1 s | 117.2 s    | sequential 1.8x faster |

Sequential cost grows like `N^2` (x5.5 per doubling), the split like
`N^1.4` (x2.7), so the crossover lies near `2 d^2`. A class is multiplied
out only when `N_r >= cross*d^2` (default 2); shorter classes take `mom`'s
own sequential steps, so in that regime `binmom` costs what `mom` costs.
The products parallelize where the steps cannot: `-j 8` brings the
`d = 48, N_r = 1200` split from 44.5 s to 6.9 s (below the sequential
10 s), so with threads a lower `--cross` pays.

## Limits of the prototype

Refused with a message, never mis-solved:

- mixed zero/non-zero think times: `mom` then carries lower replica levels
  (`util/lvtab.h`); they are linear in the previous state with the same
  `1/n` factor, so they can join the product state, but that is not done;
- singular systems, including those `mom` resolves by another class order:
  use `bin/mom`, whose `util/singular.h` fallback solves them exactly;
- open classes and load-dependent stations, as `mom`.

The binary-splitting branch builds `P` and `Q` densely (`2d` block solves)
and holds `d x d` matrices of `O(N log N)`-bit integers at the top of the
tree, so its memory is `O(d^2 N log N)` bits.

## Validation

Byte-identical to `bin/mom` on every run it accepts, both with the default
route and with every class forced through binary splitting (`--cross 0`):
904 runs each (`models/` plus the generated small models, `-e -l -t -q`),
and every benchmark run below.

## Benchmark

Exact `G` (`-e`), wall seconds, single thread unless noted, random integer
demands 1..1000 and think times (non-singular in the given class order),
300 s limit. `RxMy_Nz`: `x` classes of equal population, `y` single-server
stations, total population `z`.

| model | binmom | binmom -j 8 | mom | comom | mommod | comommod | binmom vs best other |
|---|---|---|---|---|---|---|---|
| R1M3_N300 | 0.004 | 0.012 | 0.050 | 0.021 | 0.060 | 0.046 | 5.2x |
| R1M3_N1000 | 0.009 | 0.019 | 0.323 | 0.250 | 0.595 | 0.452 | 27.8x |
| R1M3_N3000 | 0.031 | 0.042 | 4.377 | 4.413 | 6.713 | 5.428 | 141x |
| R1M3_N10000 | 0.147 | 0.163 | 95.4 | 96.4 | 89.1 | 74.4 | 506x |
| R2M3_N1000 | 1.139 | 1.140 | 1.231 | 1.345 | 1.125 | 1.605 | 1.0x |
| R2M3_N3000 | 1.160 | 0.371 | 19.19 | 22.04 | 10.88 | 14.27 | 9.4x (29x with -j 8) |
| R3M3_N600 | 0.882 | 0.885 | 0.931 | 0.883 | 0.828 | 1.200 | 0.9x |
| R2M6_N300 | 0.247 | 0.246 | 0.316 | 0.390 | 0.237 | 0.516 | 1.0x |
| R2M6_N1000 | 4.509 | 4.539 | 4.488 | 6.292 | 2.759 | 6.424 | 0.6x |
| R2M12_N120 | 0.136 | 0.133 | 0.215 | 0.262 | 0.179 | 0.547 | 1.3x |
| R4M8_N12 | 0.062 | 0.053 | 0.150 | 42.06 | 0.119 | 1.177 | 1.9x |
| R4M8_N30 | 0.094 | 0.100 | 0.185 | timeout | 0.324 | 27.06 | 2.0x |

(Full table of 23 models: every model under 0.03 s is a tie within noise.)

Reading:

- small basis, large population: the regime it was built for, and the gain
  grows with `N` as the analysis predicts (5x at `N = 300`, 506x at
  `N = 10000` for `d = 4`);
- moderate basis (`d` about 50): below `2 d^2` it takes `mom`'s steps and
  matches `mom`; the multi-modular `mommod` is then faster (2.8 s vs 4.5 s on
  `R2M6_N1000`), because that regime is bound by the length of the numbers,
  not by the step count, and the sequential steps over `Q` pay for it;
- large basis, small population (`R4M8`, `d = 1800`): `binmom` never
  splits and takes `mom`'s steps (here somewhat faster than `bin/mom`
  itself, whose run includes its own setup and progress handling).

So `binmom` dominates `mom` (never slower beyond noise, up to 500x faster)
but not the whole family: `mommod` keeps the moderate-basis regime. The
direct way to close that gap is to take `mommod`'s modular steps
(`mommod/mmengine.h`, already linked into `bin/mom`) for the classes that do
not split, instead of `mom`'s steps over `Q`.
