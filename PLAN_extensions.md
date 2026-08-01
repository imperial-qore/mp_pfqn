# Plan: filling the gaps in the MoM/CoMoM solver matrix

Baseline: worktree `mp_pfqn-dev`, branch `dev`, forked from `main` at `e068567`.

Current matrix (X = exists, . = missing):

| basis            | plain | generalized (d&c) | probabilities | singular-safe |
|------------------|-------|-------------------|---------------|---------------|
| MoM (`mom`)      | X     | `gmom` X          | `promom` X    | `safe_mom` X  |
| CoMoM (`comom`)  | X     | `gcomom` .        | `procomom` X  | `safe_comom` X|

`gmom` also has `safe_gmom` X (Step 2).

Four extensions below, ordered by value/effort. Each is independent after
step 0 and can be committed separately.

Status: Steps 0, 1, 1b, 2 and 3 are IMPLEMENTED on this branch (see the notes
at the end of each section). Step 4 is GATED OFF by its own spike; the spike
result and the alternative it suggests are recorded there.

## Step 0 - shared plumbing (prerequisite, small) [DONE]

- Move `gmom/pfqn_sing.c` and the two prototypes (`pfqn_recursion_singular`,
  `pfqn_nonsingular_recclass`, currently `gmom/gmom.h:64-65`) into `util/`
  as `util/pfqn_sing.{c,h}`. `gmom` keeps working via the util include; the
  oracle becomes callable from `mom` and from the safe orchestrators.
- Factor the tier driver out of `safe_comom/safe_comom.c` into
  `util/safe_tiers.{c,h}`: argv rewriting, tier probing, class-order
  remapping of `-t`/`-q`, exec of a child solver. `safe_comom.c` becomes a
  thin table of (solver, tiers). This is what makes steps 1 and 2 cheap.
- Caveat to encode in the header comment: the `A^T A` oracle is exact for
  gmom's b=1 matrix only (`gmom/README.md:86`); `mom`/`comom` use a
  different basis, so for those the oracle is a *hint* and each candidate
  class order must still be probed by running the solver.

Implemented: `setup1.c` and `pfqn_sing.c` moved to `util/` behind
`util/pfqn_sing.h` (`gmom.h` now includes it); tier driver extracted to
`util/safe_tiers.{c,h}` with a `safe_cfg` table, `safe_comom.c` reduced to
that table. `gmom` and `safe_comom` outputs verified bit-for-bit against
the pre-refactor binaries on all 28 models x `-e`/`-l`/`-t`/`-q`; only the
two reworded stderr diagnostics differ.

## Step 1 - `safe_mom` (highest value) [DONE]

Today `mom` silently degrades: on a singular coefficient matrix it
auto-perturbs at digit 20 (`mom/main.c:435-484`) and returns an
*approximate* answer. `safe_mom` makes the answer exact on every model.

1. Add `-X` (no-perturb) to `mom/main.c`: on singularity, print the
   existing diagnostic and exit with a distinct status instead of entering
   the perturbation branch. Fast-fail probe, same role `comom -e` plays for
   `safe_comom`.
2. `safe_mom/safe_mom.c` over `util/safe_tiers`:
   - Tier 1: `bin/mom -X` (non-singular models).
   - Tier 2: O(R) recursion-class permutations, each probed with `mom -X`;
     only the winning order re-run with the requested flag. `G(N)` is class-
     order invariant, so `-e`/`-g`/`-l` need no remapping; `-t`/`-q` are
     permuted back.
   - Tier 3: `bin/ca` (exact convolution). Covers genuine loading
     degeneracies, and `ca` handles multiserver via station expansion, so
     `safe_mom` keeps mom's full domain (unlike `safe_gmom`, below).
3. Validation: bit-for-bit `safe_mom` == `mom` on non-singular models
   (01-04, 09, 13, lcfs_*), == `ca` on `05`, `07`, `test_singular*`.
   Open/mixed (`LAMBDA`) models reported, not mis-solved.

Implemented as specified: `mom -X`/`--no-perturb`, `safe_mom/` over
`util/safe_tiers.c`, registered in `DIRS`/`EXE` and in the top-level build
recipe (which previously listed `bin/gmom` and `bin/safe_comom` in `EXE`
but never built them; those `cd` lines were added too). Validation on all
28 models x `-e`/`-l`/`-t`/`-q`: `safe_mom` == `mom` where mom is exact,
== `ca` where mom is singular, no mismatches. Tier 1: 13 models, Tier 2:
`08_multiclass`, `11_swapped`, Tier 3: 13 models. `mom`'s own output is
bit-identical to before the change.

## Step 1b - numerical perturbation parity for gmom

Requested: give `gmom` the same numerical-perturbation escape that `mom`
and `comom` take when they detect a singularity.

Already present as of `e068567` (`gmom/main.c:150,191,217,226-238`):

- manual `-p digit` / `-s seed`, same semantics and same
  `apply_perturbation_to_model` construction as `mom`;
- automatic perturbation at digit 20 when the init-time oracle reports a
  degeneracy that no class reorder removes, matching mom's auto-fallback;
- suppressed under `-e`, where a perturbed rational is not the exact `G`;
  the message points at `bin/ca` / `bin/safe_comom` instead.

Residual sub-items [ALL DONE, see the note at the end of this section]:

1. `-b`/`--bounds` parity: `mom` computes +/- perturbation ranges under
   `-b`; `gmom` has no equivalent, so a perturbed `gmom` run gives a point
   estimate with no error bracket.
2. Report the perturbation in the output header the way
   `printmodel_with_perturbation` does for `mom`, so a perturbed `gmom`
   result is self-describing rather than only flagged on stderr.
3. ~~Restore the unperturbed `L`/`Z` before computing `X` and `Q`.~~ STALE:
   this was a misreading of `measures.c`. `gmom` does not restore the
   demands and must not: the basis was built from the perturbed model, so
   the descent has to use the same demands. Instead `measures.c:107-120`
   descales analytically -- `G` is homogeneous of degree `sum(N)` and is
   divided by `scale^sum(N)`, `X_r = G(N-1_r)/G(N)` is multiplied by
   `scale`, and in `Q_kr` the two powers cancel. Verified: `gmom -p 20`
   reproduces `bin/ca` to 15-16 digits on `-g`/`-t`/`-q` for 02, 08, 09,
   and the auto-perturbed runs match `ca` on 05_sparse and
   test_singular6/7/8.
4. Optional `--force-perturb` to allow the perturbed solve under `-e`, for
   benchmarking the perturbation error against the exact `ca` answer.

Implemented:

- (1) `-b`/`--bounds` is NOT a port of `mom -b`, which is a stub that says
  so in its own output ("bounds computation using dual perturbation is in
  development", `mom/perfindices.c:154-157`). `gmom -b` solves twice at the
  same digit and seed, once with `+eps` and once with `-eps`, and prints
  `midpoint [lower, upper] relhw`. `G` is a polynomial with non-negative
  coefficients in the demands, hence monotone in them, so the interval
  ENCLOSES the exact `G(N)` (and `log G`); `X` and `Q` are ratios and are
  labelled as error indicators, not proven enclosures. A demand of 0 is left
  at 0 under `-eps`. Checked on 08_multiclass at digit 3: exact
  `G = 4.766087716788053e+34` lies inside `[4.7631e34, 4.7693e34]`.
  Required extracting the solve into `gmom_compute()` (main.c) and a capture
  mode in `gmom_measures()` (`gmom_out` in gmom.h) so it can run twice.
- (2) the default output now prints the model as read plus the per-class
  `eps=` actually applied, via `printmodel_with_perturbation`, off a view of
  the model in the ORIGINAL class order (the singularity screen may permute
  it). This covers the automatic digit-20 fallback, not only `-p`.
- (4) `--force-perturb` prints the perturbed rational under `-e`, with a
  stderr note that it is exact for the perturbed model. On test_singular7 it
  gives 9.049316350000000e+11, matching `ca` to 16 digits.
- also added: `-X`/`--no-perturb`, the exact-only probe that Step 2 needs.

## Step 2 - `safe_gmom` (small, low payoff)

`gmom` already implements the Tier-2 analogue internally (init-time
singularity check plus auto-reorder of the recursion class,
`gmom/pfqn_sing.c`), and falls back to `-p` perturbation only when no
reorder helps. So the only missing tier is the exact last resort.

- Replace the perturbation fallback with a Tier-3 hand-off to `bin/ca`,
  behind a flag so the current `-p` behaviour stays reachable for
  benchmarking.
- Either a `safe_gmom` wrapper (consistent with `safe_comom`) or an
  in-process `--exact-fallback` in `gmom`. Prefer the wrapper for symmetry.
- Note the domain mismatch: gmom is restricted to distinct single-server
  closed queues, `ca` is not, so `safe_gmom` == `gmom`-or-`ca` adds little
  beyond what the reorder already gives. Do it for uniformity, not speed.

Effort: half a day.

Implemented as the wrapper: `safe_gmom/` over `util/safe_tiers.c` with
`solver="gmom"`, `solver_flags="-X"`, `probe_flags="-e -X"`,
`fallback="ca"`; registered in `DIRS`/`EXE` and the top-level recipe.
Validation on all 28 models: 26 equal to `bin/ca -e`, Tier 1 on 11 and
Tier 3 on 15, Tier 2 never fires as predicted. The domain mismatch turned
out to be the main payoff rather than a caveat -- Tier 3 absorbs `R=1`,
multiserver and load-dependent models that `gmom` alone rejects. The two
open/mixed models are reported as out of scope (no tier is exact for them:
`gmom` rejects open classes and `ca` is a closed-network convolution), so
`safe_gmom` exits non-zero instead of printing a bogus number; `safe_mom`
covers those because its Tier 1 `mom` handles open classes directly.

## Step 3 - `promom` (moderate)

Marginal queue-length probabilities on the MoM basis, mirroring
`procomom/` (`genpmatrix.c`, `pexact.c`, `phash.c`, ~24K of C).

- Port `genpmatrix.c`'s P-matrix construction from the CoMoM basis to
  mom's basis: same `A`/`B` block recursion, different basis vector and
  different population-index map.
- Reuse `phash.c` verbatim (basis-independent memoisation) and `pexact.c`
  for the final marginal extraction.
- Validation is free and strict: `promom` must agree bit-for-bit with
  `bin/procomom` on every closed model in `models/`, since both compute the
  same marginals from different bases.

Effort: 2-3 days. Main risk is the index map between bases, which the
bit-for-bit oracle catches immediately.

Implemented as `promom/` (`promom.h`, `genqmatrix.c`, `qdecrease.c`,
`main.c`), registered in `DIRS`/`EXE` and the top-level recipe. Two
corrections to the plan as written:

1. There is no MATLAB `dmomprob.m` in `gmom-reference/`, so the equations
   were derived rather than ported. They are in `promom/README.md`; the
   quantity is

     q_d(P,n) = sum_{|k|=n} (n!/prod k_r!) prod_r L_{M,r}^{k_r} g_{d,-M}(P-k)

   the unnormalized probability that the reference COPY of station M holds n
   jobs. The n-recursion, CE and PC all follow from the multinomial Pascal
   rule, the PC picking up the n factor that is procomom's DC/DD. The system
   then has procomom's exact shape over the MoM basis, and dropping the
   n-coupled terms recovers mom's own recursion.

   The one real trap: the PC's station sum runs over the model with the
   reference copy REMOVED, so station M carries multiplicity mi_M + d_M - 1.
   Without that, the n = 0 slice is mom's system and returns G(P) instead of
   g_{-M}(P). It shows up immediately on a 2-station single-class model.

2. `phash.c`/`pexact.c` were NOT reusable: the MoM basis has no station
   index, so the hash is the (replica combination, class decrement) map and
   the initial condition is "every constant of the empty model is 1, at
   n = 0 only". Both are a few lines inline.

Validation, and a finding about the oracle. `promom -q` equals `bin/ca -q`
summed over classes to < 5e-16 relative on every model it solves, and
`promom -P` is bit-for-bit `bin/procomom -P` on every model with all mi_k=1.
On models with a replicated station the two disagree and `procomom` is the
one that is wrong: on `04_replicated` it misses the exact per-station queue
lengths at ALL THREE stations, not only the replicated one (0.861 / 7.711 /
11.428 against the exact 2.535 / 7.026 / 10.440), while `promom` reproduces
all three. Worth a separate look at `procomom`; `promom` is unaffected.

Coverage: `promom` solves 9 of the 28 models and declines 19 -- 15
mom-singular, 2 load-dependent, 2 open/mixed. The singular set is exactly
`safe_mom`'s Tier 2 + Tier 3 set, i.e. `promom` is singular precisely where
`mom -X` is, which is the check that the shared coefficient matrix behaves as
it should. It reports the singularity instead of perturbing, where
`procomom` perturbs at digit 20.

Follow-on (small, not done): a `safe_promom` over `util/safe_tiers` would
recover `08_multiclass` and `11_swapped` from the class-permutation tier
alone. Note there is no exact Tier-3 for marginals -- `bin/ca` gives queue
lengths, not distributions -- so a genuinely degenerate model would still be
declined; that tier would need a convolution-based marginal, which does not
exist in the tree yet.

## Step 4 - `gcomom` (largest, do last)

Generalized CoMoM: apply gmom's divide-and-conquer (prefix-chain queue
removal, one queue per step, b=1) to the class-oriented basis.

- Level construction: replace gmom's `setup1.c` per-level system with the
  BTF-structured class-oriented basis from `comom/setupls.c` +
  `comom/btf_decompose.c`.
- Per-level solve: gmom's overdetermined normal equations (`A^T A`) applied
  block-wise, so the BTF sparsity is not destroyed. This is the open
  design question - naive normal equations fill in the BTF blocks and
  throw away CoMoM's entire complexity advantage. Prototype this on one
  level before committing to the rest.
- Recovery: port `gmom/measures.c` `mdecrease` to the class-oriented basis
  to get `G(N)`, `X`, `Q`.
- Validation: against `bin/comom` on non-singular models and `bin/gmom` on
  the b=1 intersection of their domains.

Effort: 1-2 weeks, and it may not pay off if the fill-in problem above has
no clean answer. Gate it on a spike: measure BTF fill-in under the normal
equations on `08_multiclass` first, and abandon if the block structure is
lost.

SPIKE RUN, GATE CLOSED for the normal-equations route. Building comom's own
per-class `A11` (`comom/setupls.c`) and its block partition
(`comom/btf_decompose.c`), then forming `A11^T A11`:

```
08_multiclass  r=2  size=12   nnz 27 -> 48    (18.8% -> 33.3% dense)
08_multiclass  r=3  size=30   nnz 90 -> 192   (10.0% -> 21.3% dense)
09_asymmetric  r=2  size=42   nnz 108 -> 294  ( 6.1% -> 16.7% dense)
09_asymmetric  r=3  size=168  nnz 630 -> 1734 ( 2.2% ->  6.1% dense)
```

a 2.1x-2.8x nnz growth, with new entries in the region the block ordering
keeps empty. This is what the structure predicts: for `A = [A11 0; A21 A22]`
the Gram matrix has the block `A21^T A22`, generically nonzero, so the block
triangular solve is gone and with it CoMoM's advantage. Do not implement
`gcomom` this way.

The spike also suggests the way around it. gmom uses normal equations only
because its per-level system is overdetermined -- but that system is
CONSISTENT by construction (its rows are true identities about the same
constants), so any square nonsingular row subset has the same solution.
Selecting rows that preserve the block ordering, instead of squaring the
system, would keep the BTF structure intact. That turns Step 4 from "solve
an overdetermined sparse system exactly" into "choose a square subsystem",
which is a much smaller problem and is worth prototyping before the rest.
The spike program is in the session scratchpad as `gcomom_spike.c`.

## Suggested commit sequence

1. `refactor: move pfqn_sing to util, extract safe tier driver`   [done]
2. `feat: add -X no-perturb probe mode to mom`                    [done]
3. `feat: add safe_mom, exact MoM on singular models`             [done]
4. `fix: restore unperturbed demands before gmom measures`        (Step 1b.3)
5. `feat: add safe_gmom exact fallback to convolution`
6. `feat: add promom, marginal probabilities on the MoM basis`        [done]
7. (gated) `feat: add gcomom, generalized class-oriented MoM`         [gate
   closed for the normal-equations route; see Step 4]

Each step updates the module table in `CLAUDE.md` and adds a module
`README.md` in the style of `safe_comom/README.md`.
