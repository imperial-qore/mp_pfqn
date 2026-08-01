# promom - marginal queue-length probabilities on the MoM basis

`promom` is the MoM-basis counterpart of `procomom`: it computes the exact
marginal queue-length distribution `P(n_k = n)` of every station, and the mean
queue lengths that follow from it, in exact rational arithmetic.

```bash
source ./setup-env.sh
make -C promom                        # produces bin/promom
./bin/promom models/02_bottleneck.qn -P    # distribution, one row per station
./bin/promom models/02_bottleneck.qn -q    # mean queue lengths
```

| flag | meaning |
|---|---|
| `-P`, `--prob` | marginal distribution per station, one row of `sumN+1` probabilities per station |
| `-q`, `--qlen` | mean queue length per station (station total, comparable with `ca -q` summed over classes) |
| (none) | both, in readable form |

## The quantity

Fix a reference station, rotated to be station M. For the model augmented
with the replica vector `d` (station `j` has `mi_j + d_j` servers), at
population `P`,

```
q_d(P, n) = sum_{|k|=n} (n! / prod_r k_r!) prod_r L_{M,r}^{k_r} g_{d,-M}(P-k)
```

is the unnormalized probability that the reference COPY of station M holds
exactly `n` jobs; `g_{d,-M}` is the constant of the same model with that copy
removed. Summing over `n` gives the ordinary constant, and
`P(n_M = n) = q(N,n) / G(N)`.

## Three identities

Everything follows from the multinomial form above.

**n-recursion.** Splitting off the class of one of the `n` jobs, via the
multinomial Pascal rule `n!/prod k_r! = sum_r (n-1)!/(k_1!..(k_r-1)!..k_R!)`:

```
q_d(P, n) = sum_r L_{M,r} q_d(P - 1_r, n-1),     n >= 1
q_d(P, 0) = g_{d,-M}(P)
```

**CE.** Removing a replica of station `j` from `W - M` and applying mom's
convolution equation under the sum:

```
q_d(P, n) = q_{d-1_j}(P, n) + sum_r L_{j,r} q_d(P - 1_r, n)
```

identical to mom's, since a replica of any station leaves the count at the
reference copy alone.

**PC.** Applying mom's population constraint to `g_{d,-M}` at `P-k` and
converting back. The population there is `P - k`, not `P`, and the `k_s` term
that this produces is again resolved by the multinomial Pascal rule, which is
where the `n` factor comes from:

```
P_s q_d(P,n) - n L_{M,s} q_d(P-1_s, n-1)
    = Z_s q_d(P-1_s, n) + sum_j m'_j L_{j,s} q_{d+1_j}(P-1_s, n)
m'_j = mi_j + d_j - [j = M]
```

The `-[j = M]` matters: the station sum runs over `W - M`, which has one copy
of station M fewer. Without it the system reduces to mom's at `n = 0` and
returns `G(P)` where `q(P,0) = g_{-M}(P)` is wanted. That was the one real
error in the first implementation, and it is visible immediately on a
two-station single-class model.

Summing the PC over all classes with `Z = 0` and substituting the n-recursion
gives the conditioned form of mom's server-count identity, used by the
descent:

```
sum_j m'_j q(m+1_j, Q, n) = (|Q| - n + |m| - 1) q(m, Q, n)
```

mom's divisor `|Q| + |m|`, less the `n` jobs pinned at the reference copy and
less that copy itself.

## Structure

The per-step system therefore has exactly procomom's shape,

```
A q(n) = B q_prev(n) + n [ DC q(n-1) + DD q_prev(n-1) ]
```

only assembled over the MoM basis rather than the CoMoM one. Columns are

```
[0, cardk*r)                 level-r,     index w*r + s
[cardk*r, (cardk+cardi)*r)   level-(r-1), index (cardk+i)*r + s
```

with `w` over `multichoose(M, r)`, `i` over `multichoose(M, r-1)`, and `s = 0`
meaning `P`, `s >= 1` meaning `P - 1_s`. Rows are mom's CE rows, its class-`s`
population constraints, and its class-`r` population constraints (the `B2r`
block), all with the `n`-coupled terms added. The system is square: the row
count identity `ce + cardi*(r-1) = cardk*r` is mom's `A12` row count.

- `genqmatrix.c` builds the four blocks densely, in the manner of
  `procomom/genpmatrix.c`.
- `main.c` sweeps the population as `mom/main.c` does, solving for
  `n = 0..sum(Ncur)` at each step, with the class transition handled as mom
  handles it (the old top level becomes the new level-`(r-1)` block).
- `qdecrease.c` is `mom/mdecrease.c` carrying the `n` index, descending to
  level 0 where the base component is `q(N,n)`.

Dropping every `n`-coupled term recovers mom's own recursion, which is the
structural check that the port is faithful.

## Validation

Two independent oracles, since the obvious one turned out to be unreliable.

**Against exact convolution.** `promom -q` equals `bin/ca -q` summed over
classes to full double precision on all 9 models it solves:

```
01_single 0.0e+00   02_bottleneck 4.6e-16   03_think     3.1e-16
04_replicated 1.8e-16   09_asymmetric 2.3e-16   13_gld_small 0.0e+00
lcfs_1class 0.0e+00   lcfs_2class 4.3e-16     lcfs_3class  3.5e-16
```

**Against procomom.** `promom -P` is bit-for-bit identical to
`bin/procomom -P` on every model with all `mi_k = 1`.

**Where they disagree, procomom is wrong.** On `04_replicated` (station 1 has
`mi = 3`), against the exact per-station queue lengths from `ca`:

| station | mi | exact (ca) | procomom | promom |
|---|---|---|---|---|
| 1 | 3 | 2.534588116251616 | 0.8609086881853183 | 2.534588116251616 |
| 2 | 1 | 7.025688729838590 | 7.711274616873007 | 7.025688729838591 |
| 3 | 1 | 10.43972315390979 | 11.42781669494167 | 10.43972315390979 |

`procomom` is off at every station once any `mi > 1` is present, not only at
the replicated one. This is a pre-existing `procomom` defect, unrelated to
this module.

**Replicated stations.** `mi_k > 1` means `mi_k` replicated single-server
queues (this is how `bin/ca` expands them), so the distribution `-P` prints is
the marginal at ONE copy; `-q` multiplies its mean by `mi_k` to report the
station total, which is what makes the `ca` comparison above line up.

## Limits

- Closed models only; open/mixed (`LAMBDA`) and load-dependent (`MU`) stations
  are rejected with a message rather than mis-solved.
- The step matrix is mom's, so `promom` is singular on exactly the models that
  `bin/mom -X` rejects. Of the 28 models in `models/`, `promom` solves 9 and
  declines 19: 15 mom-singular (`05`, `06`, `07`, `08`, `10`, `11`, `12`, the
  8 `test_singular*`), 2 load-dependent (`14_ld_multi`, `15_repairman`) and 2
  open/mixed (`16_mixed`, `17_mixed_ld`). The singular set is exactly
  `safe_mom`'s Tier 2 + Tier 3 set, which is the check that the shared
  coefficient matrix behaves as it should.

  It reports the singularity and exits non-zero rather than perturbing;
  `procomom` perturbs at digit 20 there and returns an approximation. Two
  cheap ways to close the gap, in order of value:

  1. a `safe_promom` wrapper over `util/safe_tiers` -- the class-permutation
     tier alone recovers `08_multiclass` and `11_swapped`, which are exactly
     the two models `safe_mom` answers at Tier 2. There is no Tier-3 exact
     fallback for marginals, though: `bin/ca` returns queue lengths, not
     distributions, so a genuinely degenerate model would still be declined.
  2. a `-p digit` flag matching procomom's, for an approximate answer
     everywhere.
