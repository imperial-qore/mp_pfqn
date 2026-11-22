# camarg - exact marginal queue-length distributions by convolution

`bin/promom` and `bin/procomom` compute marginal queue-length distributions
from a moment basis, and both inherit that basis's singularities: `promom`
declines a degenerate model, `procomom` perturbs it at digit 20 and returns an
approximation. Convolution has no coefficient matrix, so it never goes
singular. `camarg` is therefore

- an **exact oracle** for `promom`/`procomom` on every closed model, including
  the degenerate ones where neither can be trusted, and
- the exact **Tier 3** of `safe_promom`, which had none before it existed.

It is to distributions what `bin/ca` is to means: the same exponential-in-`R`
lattice cost, meant as a reference rather than a fast path.

```bash
source ./setup-env.sh
make -C camarg                          # produces bin/camarg
./bin/camarg models/05_sparse.qn -P     # distribution, one row per station
./bin/camarg models/05_sparse.qn -q     # mean queue length per station
```

## What it computes

Conditioning a single-server station `k` on its contents:

```
q_k(N, n) = sum_{|j|=n} (n! / prod_r j_r!) prod_r L_{k,r}^{j_r} G_{-k}(N - j)
```

`G_{-k}` is the constant of the network with ONE copy of station `k` removed,
and the multinomial counts the orderings of `j` jobs in that station's queue.
Then `P(n_k = n) = q_k(N,n) / G(N)`, and `sum_n q_k(N,n) = G(N)` -- the
normalisation used internally.

Per reference station: one Buzen convolution over the population lattice
omitting one copy of that station, then one pass over the lattice accumulating
the sum above. Cost `O(M * Meff * prod_r (N_r+1))`.

`mi_k > 1` means `mi_k` replicated single-server queues (as `bin/ca` expands
them), so `-P` is the marginal at ONE copy and `-q` multiplies its mean by
`mi_k` to give the station total -- the same convention as `promom` and
`procomom`.

## Validation

`camarg -q` equals `bin/ca -q` summed over classes on ALL 24 closed models in
`models/`, max relative error `3.8e-16`, and `camarg -P` is bit-for-bit
identical to `promom -P` on each of the 9 that `promom` solves.

The 15 remaining are exactly the ones no moment-basis solver answers exactly:

```
05 06 07 08 10 11 12 and the 8 test_singular*   promom declines (singular),
                                                procomom perturbs at digit 20,
                                                camarg exact (<4e-16 vs ca)
```

The 4 not covered are `14_ld_multi` and `15_repairman` (load-dependent) and
`16_mixed`, `17_mixed_ld` (open/mixed), all rejected by design.

## Two traps worth recording

**`factorial()` overflows.** `util`'s `factorial` returns `long int`, so it is
silent nonsense past 20!. The multinomial here needs `n!` with `n` up to
`sum(N)` -- 50+ on `02_bottleneck` -- and the delay station needs `1/n_r!` on
the same scale. Both use `mpz_fac_ui` instead. The symptom was subtle: station
1 of `02_bottleneck` was exact while station 2, the bottleneck with the long
tail, read 19.16 against the true 49.45, because only the high-`n` terms were
corrupted.

**Sweep order matters.** The Buzen update `g_m(pop) = g_{m-1}(pop) + sum_r
L_mr g_m(pop - 1_r)` reads the *updated* value at `pop - 1_r`, so the lattice
must be swept in increasing population. `util`'s `nextpop` does that, which is
why the in-place update is correct.

## Limits

Closed models only; open/mixed (`LAMBDA`) and load-dependent (`MU`) stations
are rejected rather than mis-solved. The lattice is `prod_r (N_r + 1)` entries
of exact rationals, so memory and time grow exponentially in the number of
classes -- expected for a convolution reference, and the reason `safe_promom`
reaches for it only after the moment-basis tiers have failed.
