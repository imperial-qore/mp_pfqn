# mommod - multi-modular Method of Moments

`mommod` runs the entire MoM population recursion independently in `Z/p`
for a sequence of word-size primes and recovers the rational answers at
the end by incremental CRT plus rational reconstruction. It is the
"parallelisation via the Chinese Remainder Theorem" option, built as a
working solver rather than a sketch.

```bash
source ./setup-env.sh
make -C zpla && make -C mommod       # produces bin/mommod
./bin/mommod models/03_think.qn -t -j 8 -v
```

## Options

Those of `mom` (`-l -e -g -t -q`), plus:

| flag | meaning | default |
|---|---|---|
| `-j n` | primes solved concurrently (OpenMP) | 1 |
| `-b bits` | bit size of the primes drawn | 61 |
| `-B`, `--bounded` | bound the denominator during rational reconstruction | **off** |
| `-W n`, `--witness n` | primes solved but held out of the CRT, used only to verify bounded reconstructions | 1 when `-B`, else 0 |
| `--denbits n` | denominator bound for `-B`, as a power of two | 64 |
| `-v` | prime count, modulus size and timing on stderr | off |

### `-B`: bounded reconstruction, and why it is off by default

Plain Wang reconstruction splits the modulus evenly between numerator and
denominator, so it needs `2B` bits of modulus to recover a `B`-bit value
even when the denominator is only a few bits wide. Every quantity
exported here is of that kind: `mdecrease` divides only by populations,
station multiplicities and think times. Bounding the denominator recovers
the value from about `B` bits and therefore from **half as many primes**,
which is the largest single speedup available to this solver.

The catch is that a bounded reconstruction which succeeds is not
necessarily correct. While the modulus is still too small for the true
value to satisfy the bounds, the Euclid nevertheless returns *some* pair
that does; and for a quantity whose denominator exceeds the bound it can
keep returning different spurious pairs indefinitely, so the run never
terminates. This was observed directly: an early version using a `2^256`
bound returned garbage on `models/01_single.qn` and failed to terminate
on `models/09_asymmetric.qn`.

`-B` is therefore sound only in combination with an independent check,
which the **witness primes** supply. `-W n` solves `n` further modular
images and deliberately keeps them out of the CRT accumulator, so their
residues carry information the reconstruction cannot have used. A
candidate `num/den` is accepted only if `num == c_w * den (mod p_w)` for
every witness `p_w`. A wrong candidate survives one 61-bit witness with
probability about `2^-61` per quantity; `-W 2` squares that. If the check
fails, the solver falls back to unbounded reconstruction, whose size
condition `2|num|den < M` guarantees uniqueness with no auxiliary
evidence at all.

Even so the default is off. Unbounded reconstruction is correct by a
size argument that needs no probabilistic step, and for a solver whose
entire purpose is exactness that is the right default; `-B` trades a
`2^-61` failure probability for a factor of two, and that trade should be
made deliberately by the caller rather than silently.

## Design

- **`zpla/`** is a linear-algebra layer over `Z/p` mirroring the `gmpla`
  `mpq_*` API, so `setupls.c`, `blocksolve.c` and `mdecrease.c` are
  transliterations of the `mom` sources rather than reimplementations.
  Residues are held in **Montgomery form**: the inner-loop multiply is
  three 64x64 products instead of a 128-by-64 hardware division.
- **Pivoting** in `zp_ludcmp` is by first nonzero. Over a finite field
  there is no magnitude and no growth to control, so any nonzero pivot
  serves.
- **Unlucky primes** (some diagonal block singular mod `p`) are detected
  and discarded. A matrix that is singular over `Q` is singular modulo
  every prime, and is reported after 40 consecutive failures. This is
  strictly better behaviour than `mom`, which silently substitutes a
  perturbed model.
- **Termination** is not by an a priori height bound. Primes are added
  until `G` reconstructs to the same rational twice, then the full set of
  constants is reconstructed and must reproduce itself at a later prime
  count.
- **Raw constants are exported, not the measures.** `X_r = G_r/G` is a
  rational whose numerator and denominator are each about as long as `G`,
  so reconstructing it needs roughly four times the modulus the integers
  need. `mommod` reconstructs `G`, `G_r` and `G^{+k}_r` and forms the
  ratios once, exactly, afterwards.
- **Thread safety**: the modulus is `__thread`-local, and the first batch
  of primes runs sequentially to warm the lazily-memoized `nck` cache,
  which is read-only thereafter.

### Prime supply and exhaustion

Primes are drawn descending from `2^b` (`-b`, default 61) and must stay
above a floor of `max(3, max_r N_r + 1)`: below the largest class
population the division by `n_r` is not invertible, and `p = 2` would
break the Montgomery reduction, which inverts `p` modulo `2^64`.

At `-b 61` exhaustion is unreachable. The prime count scales with the
number of *digits* of `G`, not with its magnitude: about `0.055` primes
per decimal digit with `-B`, `0.11` without.

| digits of G | primes with -B | primes unbounded |
|---|---|---|
| 2 924 | 159 | 318 |
| 9 939 | 541 | 1 083 |
| 11 328 (the paper's Table 1 result, `G = 1.21e11328`) | 617 | 1 234 |
| 100 000 | 5 446 | 10 891 |

There are about `2^60/ln(2^61) = 2.7e16` primes in `[2^60, 2^61]`, so
exhausting them would take a `G` of roughly `5e17` decimal digits. A
model of that size is unreachable for every other reason first.

Measured at `G = 1e9939` (Table 1 model, `Ntot = 5262`):

| | wall | maxRSS |
|---|---|---|
| `mom` | 4.76 s | 29.4 MB |
| `mommod -j8` | 6.85 s | 56.9 MB |
| `mommod -j8 -B` | 3.65 s | 30.8 MB |

720 primes, 43920-bit modulus, `G` identical to `mom` in all 9939
digits. At this size the binding constraint is time, not prime supply:
720 images of a 5262-step recursion. Memory is dominated by the CRT
accumulators rather than by the images, which is why `-B` halves it.

With a small `-b` it is reachable, and the solver now reports it:

```
$ ./bin/mommod exh.qn -e -b 8
Error: exhausted the primes of 8 bits after 26 of them (184 bits of modulus).
The reconstruction needs a larger modulus than this prime size can supply;
rerun with a larger -b
```

Before the floor was enforced this was silent and fatal in three
different ways, all worth recording because none of them announces
itself:

1. **Wrap-around and prime reuse.** The descending search ran past 2 into
   negative candidates; `mpz_probab_prime_p` tests `|n|` and
   `mpz_get_ui` returns the low limb of the absolute value, so primes
   already consumed re-entered the sequence. A repeated prime makes
   `M mod p` zero, so `zp_crt_step` finds `inv_plain(0,p) = 0`, adds
   nothing to the accumulator, and yet `M` is multiplied by `p` again.
   The residue silently stops being a valid CRT residue for the claimed
   modulus, reconstruction never converges, and the solver **hangs with
   no output**. This was the observed behaviour at `-b 8`.
2. **Misdiagnosis as a singular system.** Once the primes fall below
   `max_r N_r`, every image fails on a non-invertible `n_r` and is
   counted as unlucky. After 40 consecutive failures the solver reports
   that the coefficient matrix is singular over the rationals, which
   would be a false accusation against a perfectly well-conditioned
   model.
3. **Silently wrong arithmetic at `p = 2`.** The Newton iteration in
   `zp_setmod` computes `p^{-1} mod 2^64` and requires `p` odd; at
   `p = 2` it does not converge and every subsequent `redc` is wrong,
   with no error raised.

## Validation

Bit-for-bit identical to `bin/mom` on `-e`, `-t` and `-q` for every
non-singular model in `models/`, plus the paper's Table 1 model at
several populations. On the singular models the two differ by design:
`mom` auto-perturbs and returns an approximation, `mommod` reports the
singularity.

## A latent bug in `mom` that this exposed

`mom` was **not deterministic**: on `models/03_think.qn -q`, 5 runs in 80
returned different answers, including negative queue lengths, with
nothing on stderr. AddressSanitizer localised it to `setupls.c:217`
reading four bytes from the one-byte `USE_LINBOX`.

The cause was in `util/util.h`, which defined `bool` as `unsigned int`
whenever `<stdbool.h>` had not already been included in that translation
unit. `main.c` includes `<stdbool.h>` and so allocated a 1-byte `_Bool`;
`setupls.c` and `blocksolve.c` do not, and read the same global as a
4-byte `unsigned int`, picking up three adjacent bytes of unrelated data.
`USE_LINBOX` therefore read as true at random, making `setupls` skip the
LU factorisation while `blocksolve` still used it, or the reverse.

Fixed by including `<stdbool.h>` unconditionally. `mom` is now
deterministic over 200 runs and agrees with `mommod`. The defect affected
every solver in the repository that mixes translation units this way, not
only `mom`.

## Performance

Paper Table 1 model (`R=6`, `q=3`, multiplicities 12/5/8), scaled
population, 16-core machine:

Best of three runs each:

| Ntot | digits(G) | mom | mommod -j8 | mommod -j8 -B | primes with -B |
|---|---|---|---|---|---|
| 192 | 384 | 0.080 s | 0.030 s | 0.023 s | 30 |
| 384 | 749 | 0.116 s | 0.069 s | 0.047 s | 57 |
| 768 | 1475 | 0.206 s | 0.185 s | 0.122 s | 115 |
| 1536 | 2924 | 0.515 s | 0.608 s | 0.345 s | 205 |

Bounded reconstruction halves the prime count (407 to 205 at
`Ntot=1536`) and with it the wall clock, taking the solver from 1.2x
slower than `mom` to 1.5x faster. Single-threaded, `mommod` remains
several times slower than `mom`: the gain comes from parallelism, not
from the arithmetic.

Memory at `Ntot=1536`: 23.8 MB for `mommod` regardless of thread count,
15.7 MB for `mom`. Parallel scaling: 3.8x on 8 threads (47% efficiency);
16 threads is slower than 8.

**The honest reading: multi-modular is not asymptotically cheaper.** The
number of primes grows linearly with `log G = Theta(N log N)`, and each
prime replays the whole `N`-step recursion, so total work is
`Theta(theta^2 N^2 log N)` word operations, the same order as the
rational solver, whose operands also grow linearly in `N`. What
multi-modular changes is not the exponent but the constant, the memory,
and above all the fact that the work becomes embarrassingly parallel,
whereas the recursion itself is a strictly sequential chain of `N`
dependent steps. GMP's small-by-large multiply is efficient enough that
on one core the constant favours `mom` by several times; the speedup
seen above is bought entirely with cores, and with the halved prime count
that `-B` provides.

The remaining inefficiency is **parallel efficiency**: primes are drawn
and joined in batches of `-j`, so every batch is a barrier and the CRT
and reconstruction probe are serial. Scaling is 3.8x on 8 threads (47%),
and 16 threads is slower than 8. A work queue with per-thread CRT
accumulators merged at the end would remove both barriers.

## Where multi-modular is unambiguously better

- Memory is `O(theta)` words per image and independent of `N`, against
  `O(theta N log N)` bits for the rational solver. On models where `mom`
  runs out of memory this is the difference between an answer and none.
- The work is embarrassingly parallel across primes, with no
  communication until the final CRT. The rational recursion has no such
  axis: it is a strictly sequential chain of `N` dependent steps.
- No gcd, no allocation churn, no operand growth in the inner loop, so
  performance is flat and predictable rather than degrading as the
  recursion proceeds.
