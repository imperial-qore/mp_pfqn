# procomom - ProCoMoM, marginal queue-length probabilities on the CoMoM basis

`procomom` computes the marginal queue-length distribution `P(n_k = n)` of
every station and the mean queue lengths that follow, by carrying a
probability vector alongside CoMoM's population recursion:

```
A p(n) = B p_prev(n) + n [ DC p(n-1) + DD p_prev(n-1) ]
```

`bin/promom` is the same construction on the MoM basis; see `promom/README.md`
for the derivation of the three identities behind that recursion, which apply
here too.

```bash
source ./setup-env.sh
make -C procomom
./bin/procomom models/02_bottleneck.qn -P    # distribution, one row per station
./bin/procomom models/02_bottleneck.qn -q    # mean queue lengths
./bin/procomom models/02_bottleneck.qn -p 20 # perturbed (approximate) solve
```

Unlike `promom`, `procomom` auto-perturbs at digit 20 on a singular system
rather than declining, so it answers the degenerate models too -- but
approximately.

## Multiserver stations: three defects, all fixed

Found by parity testing against `bin/promom` and exact `bin/ca` on
`04_replicated`. Before the fixes, `procomom -q` was wrong at EVERY station of
any model containing a replicated station, not only the replicated one:

| station | mi | exact (ca) | procomom before | procomom now |
|---|---|---|---|---|
| 1 | 3 | 2.534588116251616 | 0.8609086881853183 | 2.534588116251616 |
| 2 | 1 | 7.025688729838590 | 7.711274616873007 | 7.025688729838591 |
| 3 | 1 | 10.43972315390979 | 11.42781669494167 | 10.43972315390979 |

`procomom -P` is now bit-for-bit identical to `promom -P` on this model, the
two bases agreeing exactly as they should.

**1. Missing multiplicities.** `genpmatrix.c` never referenced `qnm->mi`,
while the CoMoM system it mirrors carries `mi[k-1]*L[k-1][s-1]` in its
population constraints (`comom/setupls.c:103`) and `mi[k-1]*L[k-1][r-1]` in
the class-`r` block (`comom/setupls.c:135`). Both are now weighted. Any `mi`
factor is the identity when all `mi = 1`, so this cannot change a
single-server result, and none changed.

**2. No column for the reference station's own replica.** The station being
measured is rotated to position M, and its replica term carries multiplicity
`mi_M - 1` (see the PC identity in `promom/README.md`). The basis holds the
base component plus stations `1..M-1`, so station M had nowhere to live.

The stride is now set per rotation (`set_stride` in `main.c`): `M+1`
components when the reference station is replicated, `M` otherwise. Widening
it unconditionally does NOT work -- with a single-server reference the weight
`mi_M - 1` is zero everywhere, so the extra column is unconstrained and
`A^T A` goes singular. That is presumably why station M was dropped in the
first place; the omission is only correct while `mi_M = 1`.

**3. Per-copy versus per-station means.** `mi_k > 1` means `mi_k` replicated
single-server queues, so the distribution is the marginal at ONE copy and the
station total is `mi_k` times its mean. `-q` now applies that factor, which is
what makes it comparable with `ca -q` summed over classes. `-P` remains the
per-copy distribution, the same convention as `promom -P`.

## Fixed: duplicated rows on the perturbation retry

`main.c` printed each station's result inside the rotation loop, but a
singular solve `goto`s back to `solve_attempt` and recomputes ALL stations
with the perturbed model. Rows already emitted were therefore printed twice.
On `11_swapped` (M=2) `-q` produced three rows, the first two identical, so
any downstream parse silently mis-associated stations. Both `-q` and `-P` now
buffer their results and print only after all `M` stations have succeeded.

## Relation to promom

| | procomom | promom |
|---|---|---|
| basis | CoMoM: population shift x station | MoM: replica combination x class decrement |
| singular models | auto-perturbs at digit 20 (approximate) | declines; see `safe_promom` |
| replicated station elsewhere | correct (after the fixes above) | correct |
| replicated REFERENCE station | correct (after the fixes above) | correct |
| `-P` convention | marginal at one copy | marginal at one copy |
| `-q` convention | station total (`mi` x per-copy mean) | station total |
