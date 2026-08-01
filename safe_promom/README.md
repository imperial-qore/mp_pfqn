# safe_promom - marginal probabilities, exact wherever a class order makes them so

`bin/promom` inherits mom's coefficient matrix, so it is singular on exactly
the models `bin/mom -X` rejects, and it declines them rather than perturbing
(`bin/procomom` perturbs at digit 20 and returns an approximation).
`safe_promom` recovers the models whose singularity depends only on which
class is processed last, using the class-permutation tier of Casale, "CoMoM"
(IEEE TSE 2009), Appendix A, and the tier driver shared with `safe_mom`,
`safe_comom` and `safe_gmom` (`util/safe_tiers.c`).

```bash
source ./setup-env.sh
make -C safe_promom          # produces bin/safe_promom
./bin/safe_promom models/08_multiclass.qn -q
```

Options are those of `promom`: `-P` (distribution per station), `-q` (mean
queue length per station).

## Tiers

| tier | method | handles |
|---|---|---|
| 1 | plain `bin/promom` | non-singular models |
| 2 | `promom` with a different recursion class (O(R) search) | singularities that depend only on which class is processed last |
| 3 | none | see below |

**There is no Tier 3.** The other wrappers fall back to `bin/ca`, exact for
`G`, `X` and `Q`, but convolution in this tree returns mean queue lengths, not
marginal DISTRIBUTIONS. A model degenerate under every class order is
therefore declined with a message rather than answered approximately. Closing
that would need a convolution-based marginal, which does not exist here yet.
This is the one wrapper in the family with a genuinely open bottom tier.

## Output is not un-permuted

The other wrappers map a Tier-2 `-t`/`-q` result back to the original class
order, because those results are indexed by class. `promom`'s `-P` and `-q`
are indexed by STATION and carry no class index at all: the marginal at a
station does not depend on the order the classes are processed in. So this
wrapper declares `permuted_flags = ""` and prints a Tier-2 result exactly as
it comes back. That field was added to `safe_cfg` for this module; the three
existing wrappers declare `"-t -q"` and are unaffected.

## Validation

`safe_promom -q` equals `bin/ca -q` summed over classes on all 11 of the 28
models in `models/` that it answers, `promom` alone answering 9:

```
answered (11): 01_single 02_bottleneck 03_think 04_replicated 08_multiclass
               09_asymmetric 11_swapped 13_gld_small lcfs_1class lcfs_2class
               lcfs_3class                      max relative error 4.6e-16
declined (17): 13 singular under every class order (05 06 07 10 12 and the 8
               test_singular*), 2 load-dependent (14_ld_multi 15_repairman),
               2 open/mixed (16_mixed 17_mixed_ld)
```

Tier 2 recovers the two models `promom` alone declines for a
class-order-dependent singularity, the same two `safe_mom` answers at its
Tier 2:

```
08_multiclass  3.970840712887463e+00  7.058318574225073e+00  3.970840712887463e+00
ca (summed)    3.970840712887463e+00  7.058318574225074e+00  3.970840712887463e+00

11_swapped     2.695487897605068e-03  7.996452677998317e+00
ca (summed)    2.695487897605068e-03  7.996452677998317e+00
```

`safe_mom`, `safe_comom` and `safe_gmom` are bit-for-bit unchanged by the
`safe_cfg` extension.

A load-dependent model never reaches the class search (the driver skips it
for `isLD`), so with no Tier 3 it is reported as out of domain rather than as
singular -- otherwise `14_ld_multi` and `15_repairman` would be blamed on a
degeneracy they do not have.

## Relation to the other wrappers

| wrapper | Tier 1/2 solver | Tier 3 | Tier 2 does real work? |
|---|---|---|---|
| `safe_mom` | `mom -X` | `ca` | yes |
| `safe_comom` | `comom` | `ca` | yes |
| `safe_gmom` | `gmom -X` | `ca` | no, gmom reorders internally |
| `safe_promom` | `promom` | none | yes |
