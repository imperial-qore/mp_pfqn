# CoMoM boundary-recurrence falsification

Status: **counterexample found; investigation stopped**. This is an isolated
reproducer of the proposed two-class scalar reformulation, not a production
solver or a general CoMoM implementation.

Run from the repository root, using only Python's standard library:

```bash
python3 experiments/comom_boundary_poc/reproduce.py
```

For an independent check with the repository's exact convolution solver:

```bash
./bin/ca experiments/comom_boundary_poc/counterexample.qn -e
```

The model file multiplies all four demands by `2^48` to accommodate the
installed binary's integer-only parser. Its exact answer is therefore the
normalized answer multiplied by `(2^48)^3`, reported by the reproducer as
`reference_G_integer_model`.

The experiment grants the proposal an exact terminal value, exact parameters,
and exact arithmetic after one deliberate perturbation. The reference is a
positive two-station convolution, independent of the proposed recurrence.

## Counterexample

Take no delays, populations `(1, 2)`, and demand matrix

```
    1 + epsilon    1
    1              1 + epsilon
```

with `epsilon = 2^-48`. These parameters satisfy the original assumptions
`a > b > 0`, `d > c > 0`; there is no structural degeneracy or zero population.

In the proposed recurrence, `lambda = 2 + epsilon` and

```
g_n = lambda * g_(n-1) + f_n
f_n = binom(n+p,p) * (a^(p+1) c^n - b^(p+1) d^n) / (a-b).
```

The sign switch is exactly at `n = 2`. Recover `g_2` backward from an exact
`g_3`, using `f_3 < 0`. The control with exact forcing recovers the convolution
answer exactly. A relative perturbation to the propagated `g_3` has gain about
`0.8333333333`, consistent with the scalar nonamplification proof.

Now perturb just the first auxiliary term in `f_3` by `u = 2^-53` relatively.
Writing `a = 1 + epsilon`, the terms are `A = 4a^2` and `B = 4a^3`:

```
f_3          = (A - B) / epsilon
perturbed_f3 = ((1+u) A - B) / epsilon
g_2          = (g_3 - f_3) / lambda.
```

The perturbation does not change the forcing sign. Nevertheless, the output
relative error divided by `u` is exactly

```
K = 4(1+epsilon)^2
    / [epsilon (2+epsilon) (12+18epsilon+8epsilon^2+epsilon^3)],
```

so `K ~ 1/(6 epsilon)` as the queues approach equality. At `epsilon = 2^-48`,
the gain is about `4.69125e13` and the output relative error is **0.520833%**.
This is sensitivity to an injected auxiliary error, not a claim that an
ordinary floating-point evaluation necessarily incurs that exact error.

The physical answer tends smoothly to 12. It is a positive degree-three
polynomial in the demands, so relative componentwise perturbations of the
model itself have first-order gain at most three. The large gain above is
introduced by the auxiliary representation, not by the queueing model.

## Scope of the negative result

The conditional proof for nonnegative scalar propagation remains valid. The
counterexample shows that the full proposed representation still has unbounded
amplification when its auxiliary sequences suffer independent small errors.
It does **not** demonstrate exponential growth with population under exact
forcing, or rule out a different theoretical revision.

In this particular symmetric example, symbolic cancellation gives `f_3 = -4a^2`
and avoids this failure. Thus this is not an impossibility theorem: it shows
that the proposed direction change alone is insufficient. Per the requested
stop condition, no replacement representation or broader implementation was
pursued after finding the counterexample.
