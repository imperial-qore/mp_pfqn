# comommod - multi-modular CoMoM

`comommod` runs the entire CoMoM class-oriented recursion independently
in `Z/p` for a sequence of word-size primes and recovers the rational
answers at the end by incremental CRT plus rational reconstruction. It is
to `comom` what `mommod` is to `mom`: the parallelisation of Makaronidis
and Casale, "Efficient Parallelization of the Method of Moments for
Queueing Networks Using Multi-Modular Algebra" (VALUETOOLS 2011), applied
to the class-oriented recursion instead of the population one.

```bash
source ./setup-env.sh
make -C zpla && make -C comommod    # produces bin/comommod
./bin/comommod models/03_think.qn -t -j 8 -v
```

## Options

Those of `comom` (`-e -l -g -t -q -d`), plus:

| flag | meaning | default |
|---|---|---|
| `-j n` | primes solved concurrently (OpenMP) | 1 |
| `-b bits` | bit size of the primes drawn | 61 |
| `-B`, `--bounded` | bound the denominator during rational reconstruction | **off** |
| `-W n`, `--witness n` | primes solved but held out of the CRT, used only to verify bounded reconstructions | 1 when `-B`, else 0 |
| `--denbits n` | denominator bound for `-B`, as a power of two | a third of the modulus so far |
| `-v` | prime count, modulus size and timing on stderr | off |

`comom`'s `-p`/`-s` perturbation flags are deliberately absent: see
**Singular models** below.

## Design

Every source here is a transliteration of its `comom/` counterpart with
`mpq_t` replaced by a residue in a machine word, so the two solvers walk
the same recursion, build the same rows in the same order and must agree.
What changed:

- **`zpla/`** supplies the `Z/p` arithmetic, unchanged from `mommod`, so
  `setupls.c`, `btf_decompose.c` and `btf_solve.c` are transliterations
  rather than reimplementations. Residues are in Montgomery form.
- **The basis is built once, not once per prime.** `Dn` depends only on
  `(r,M)`, never on `p`, so `makeDn` runs in `main` and every image shares
  it read-only. `comom` rebuilds it per class and never frees it, which
  costs nothing in a single solve but would leak once per prime here.
- **`hash` is no longer a linear scan.** `comom` finds a population by
  scanning the whole basis (`int_matmatchrow`), which the impossible-combo
  fixup calls `O(card*M)` times per population step. Since the basis is
  now immutable and shared, `makeDn` also builds a sorted mixed-radix
  encoding of it and `hash` binary-searches that, falling back to the
  linear scan when the encoding would overflow 64 bits. The returned index
  is the same either way.
- **Unlucky primes.** A BTF diagonal block that is singular mod `p` need
  not be singular over `Q`, so the solver first falls back to a direct LU
  of `A11` -- which solves the same system -- and only discards the prime
  when `A11` itself is singular mod `p`. A matrix singular over `Q` is
  singular modulo every prime and is reported after 40 consecutive
  failures.
- **The rank-revealing sub-solve is guarded.** On models where some basis
  population is unreachable (`N_s < M` for a class `s` processed before
  the last), `comom` re-solves on the reachable sub-system with a
  Gaussian elimination that simply stops when a pivot column is empty and
  sets the remaining unknowns to zero. Modulo `p` that is dangerous: an
  entry that is non-zero over `Q` can vanish, the rank drops, and a
  different particular solution comes out with nothing to signal it. The
  image therefore discards the prime when a pivot goes missing. Since
  rank mod `p` is at most the rank over `Q`, this only rejects primes that
  are genuinely unlucky, and a system that is rank-deficient over `Q`
  fails for every prime and is reported as singular.
- **Raw constants are exported, not the measures.** `X_r = G_r/G` is a
  rational whose numerator and denominator are each about as long as `G`,
  so reconstructing it needs roughly four times the modulus the integers
  need. `comommod` reconstructs `G(N)`, `G(N-e_r)` and `G^k(N-e_r)` and
  forms the ratios once, exactly, afterwards.
- **Thread safety**: the modulus is `__thread`-local, and the first batch
  of primes runs sequentially to warm the lazily-memoized `nck` cache,
  which is read-only thereafter.
- **Termination** is not by an a priori height bound. Primes are added
  until `G` reconstructs to the same rational twice, then the full set of
  constants is reconstructed and must reproduce itself at a later prime
  count.

## `-B`: bounded reconstruction, and why its default bound differs from mommod's

Plain Wang reconstruction splits the modulus evenly between numerator and
denominator, so it needs `2B` bits of modulus to recover a `B`-bit value
however narrow the denominator is. In `mommod` every exported quantity
has a denominator of at most a few dozen bits, because `mdecrease`
divides only by populations, station multiplicities and think times, so a
fixed `2^64` bound recovers the value from about `B` bits and halves the
prime count.

CoMoM has no such bound. Every population step divides by `n_r` and every
basis solve divides by the pivots of `A11`, so the denominator grows with
the population along with the numerator; measured on the model below it
settles at about a fifth of the length of `G`:

| population scale | digits of num(G) | digits of den(G) |
|---|---|---|
| x1 | 255 | 46 |
| x2 | 537 | 122 |
| x4 | 1141 | 318 |
| x8 | 2424 | 790 |

A fixed 64-bit bound therefore never fires and `-B` degenerates into
plain Wang. Unless `--denbits` says otherwise the bound is instead taken
as a third of the modulus accumulated so far, which leaves the numerator
two thirds of the modulus rather than the half Wang leaves it, and cuts
the prime count by about 30%:

| population scale | primes, unbounded | primes, `-B` |
|---|---|---|
| x1 | 39 | 30 |
| x2 | 81 | 57 |
| x4 | 163 | 129 |
| x8 | 363 | 258 |

That split is a heuristic, not a theorem, which is exactly what the
witness primes are there to check: `-W n` solves `n` further modular
images and deliberately keeps them out of the CRT accumulator, so their
residues carry information the reconstruction cannot have used. A
candidate `num/den` is accepted only if `num == c_w * den (mod p_w)` for
every witness `p_w`, which a wrong candidate survives with probability
about `2^-61` per witness. If the check fails, the solver falls back to
unbounded reconstruction, whose size condition `2|num|den < M` needs no
auxiliary evidence.

`-B` is off by default for the same reason as in `mommod`: unbounded
reconstruction is correct by a size argument with no probabilistic step,
and for a solver whose entire purpose is exactness that is the right
default.

## Singular models

`comom` auto-perturbs when its coefficient matrix is singular, multiplying
the demands by `10^20` and adding a permutation, and returns an
approximation with a warning on stderr; with `-e` it refuses instead.
`comommod` always refuses, reporting

```
Error: the linear system is singular modulo 41 consecutive primes;
the coefficient matrix is singular over the rationals.
```

This is the same stance `mommod` takes against `mom`, and it agrees with
`comom -e` on every model in `models/`: the two classify exactly the same
models as singular (`05 06 07 08 10 11 12`, `test_singular*`). For an
exact answer on those models the tiered `safe_comom` is the right tool;
teaching it to drive `comommod` instead of `comom` on its first two tiers
would be the natural next step, and is not done here.

## Validation

Bit-for-bit identical to `bin/comom` on `-e`, `-l`, `-t` and `-q` for
every model in `models/` that `comom` can solve exactly, under every
combination of `-B`, `-W 2`, `-j 8` and `-b 40`:

```
01_single  02_bottleneck  03_think  04_replicated  09_asymmetric
13_gld_small  15_repairman  lcfs_1class  lcfs_2class  lcfs_3class
```

`09_asymmetric` is the model that exercises the unreachable-population
sub-solve (`N_2 = 5 < M = 6`), so that path is covered.

Two notes on `-l`. `mommod` computes `log(G)` by splitting off the binary
exponent, which rounds twice and lands one ulp away from every other
solver in the tree on small models; `comommod` takes `comom`'s route
(`mpq -> mpf -> double`, then `log`) instead and keeps the exponent split
only as a fallback. That fallback is the one place the two solvers'
printed output differs: once `G` exceeds what a double can hold, `comom`
prints `inf` and `comommod` prints the true logarithm. On the scaled
model below that happens from about 310 digits of `G` upwards, so `-e`,
`-t` and `-q` still agree there and only `-l` does not.

## Performance

`09_asymmetric` (`R=3`, `M=6`) with the populations scaled, best of three
runs, 16-core machine:

| population scale | digits(G) | `comom` | `comommod` | `comommod -j8` | `comommod -j8 -B` | primes (plain / `-B`) |
|---|---|---|---|---|---|---|
| x1 | 255 | 0.39 s | 0.51 s | 0.18 s | 0.18 s | 39 / 30 |
| x2 | 537 | 0.24 s | 0.23 s | 0.07 s | 0.06 s | 81 / 57 |
| x4 | 1141 | 1.06 s | 0.76 s | 0.17 s | 0.16 s | 163 / 129 |
| x8 | 2424 | 5.22 s | 3.00 s | 0.58 s | 0.44 s | 363 / 258 |

Both solvers are slower at x1 than at x2, which is not a measurement
artefact: at x1 the second class has `N_2 = 5 < M = 6`, so some basis
populations are unreachable and every population step pays the cubic
sub-system elimination described above. From x2 on that path is not taken.

Memory at x8, and parallel scaling of the same run:

| | maxRSS |
|---|---|
| `comom` | 8.5 MB |
| `comommod -j1` | 3.4 MB |
| `comommod -j8` | 5.6 MB |

| threads | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| wall | 2.99 s | 1.64 s | 0.92 s | 0.57 s | 0.53 s |

5.2x on 8 threads, 66% efficiency, against `mommod`'s 3.8x and 47% on the
same machine: CoMoM's images are longer and more uniform, so the batch
barrier costs proportionally less. It is still a barrier: primes are
drawn and joined in batches of `-j`, and the CRT and the reconstruction
probe are serial, so a work queue with per-thread accumulators merged at
the end would recover the rest.

The shape is the same as `mommod`'s. Multi-modular is not asymptotically
cheaper: the number of primes grows linearly with `log G`, and each prime
replays the whole recursion, so the total work is of the same order as
the rational solver, whose operands also grow linearly. What changes is
the constant, the memory, and above all the fact that the work becomes
embarrassingly parallel, whereas the recursion itself is a strictly
sequential chain of dependent steps.

Two things are more favourable here than in `mommod`. The per-prime
constant is better, because CoMoM's inner loop is dominated by dense
block solves rather than by GMP's small-by-large multiply, so `comommod`
overtakes `comom` on one core once `G` passes about five hundred digits,
rather than staying several times slower as `mommod` does. And the
memory is `O(cardG*M)` words per image and independent of `N`, against
`O(cardG*M * N log N)` bits for the rational solver.
