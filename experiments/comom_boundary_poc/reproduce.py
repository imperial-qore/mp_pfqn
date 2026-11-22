#!/usr/bin/env python3
"""Falsify the complete scalar boundary-recurrence representation.

All operations use Fraction. One auxiliary forcing term receives a deliberate
relative perturbation of 2**-53. There is no floating-point evaluation, tail
truncation, or inaccurate terminal condition in the experiment itself.

This is a reproducer, not a parameter search: stop after the known counterexample.
"""

import json
from fractions import Fraction
from math import comb
from pathlib import Path


def coefficient(p, n, a, b, c, d):
    """Positive two-station convolution, independent of the proposed recurrence."""
    return sum(
        (
            comb(i + j, i) * a**i * c**j
            * comb(p - i + n - j, p - i) * b**(p - i) * d**(n - j)
            for i in range(p + 1)
            for j in range(n + 1)
        ),
        Fraction(0),
    )


def main():
    epsilon = Fraction(1, 2**48)
    perturbation = Fraction(1, 2**53)
    a = d = 1 + epsilon
    b = c = Fraction(1)
    p, n = 1, 2
    j = n + 1
    rate = (a * d - b * c) / (a - b)
    reference = coefficient(p, n, a, b, c, d)
    terminal = coefficient(p, j, a, b, c, d)

    # Each term includes the common binomial coefficient. Perturbing it as
    # a whole is equivalent to perturbing one of the two power products.
    term_a = comb(p + j, p) * a ** (p + 1) * c**j
    term_b = comb(p + j, p) * b ** (p + 1) * d**j
    forcing = (term_a - term_b) / (a - b)
    perturbed_forcing = ((1 + perturbation) * term_a - term_b) / (a - b)

    # Establish that the model and the chosen backward branch obey the stated
    # assumptions. This example's switch is exactly n = p + 1 = 2.
    assert a > b > 0 and d > c > 0
    assert rate == 2 + epsilon
    assert forcing == -4 * a**2
    assert forcing < 0 and perturbed_forcing < 0
    assert reference == 12 + 18 * epsilon + 8 * epsilon**2 + epsilon**3

    recovered = (terminal - forcing) / rate
    assert recovered == reference

    # Control: perturb the propagated scalar, keeping the forcing exact.
    scalar_result = ((1 + perturbation) * terminal - forcing) / rate
    scalar_gain = abs(scalar_result / reference - 1) / perturbation
    assert scalar_gain == terminal / (rate * reference)
    assert 0 <= scalar_gain <= 1

    # Falsification: perturb an auxiliary term, with an exact scalar boundary.
    result = (terminal - perturbed_forcing) / rate
    relative_error = abs(result / reference - 1)
    gain = relative_error / perturbation
    predicted_gain = term_a / ((a - b) * rate * reference)
    assert gain == predicted_gain
    assert gain > 10**13
    assert relative_error > Fraction(1, 200)  # More than 0.5 percent.

    root = Path(__file__).resolve().parent
    # The installed CA binary accepts integer demands. Uniformly scale them;
    # G is homogeneous of degree p+n, so its exact reference scales likewise.
    demand_scale = 2**48
    model = (
        f"2\n{p} {n}\n0 0\n2\n"
        f"1 {a * demand_scale} {c * demand_scale}\n"
        f"1 {b * demand_scale} {d * demand_scale}\n"
    )
    assert (root / "counterexample.qn").read_text() == model
    output = {
        "status": "COUNTEREXAMPLE_FOUND_STOPPED",
        "arithmetic": "exact rational, with one deliberate relative perturbation",
        "epsilon": str(epsilon),
        "relative_auxiliary_perturbation": str(perturbation),
        "populations": [p, n],
        "reference_G": str(reference),
        "integer_model_demand_scale": demand_scale,
        "reference_G_integer_model": str(reference * demand_scale ** (p + n)),
        "perturbed_G": str(result),
        "scalar_propagation_gain": float(scalar_gain),
        "auxiliary_to_output_gain": float(gain),
        "relative_output_error": float(relative_error),
        "relative_output_error_percent": float(100 * relative_error),
        "failure": "auxiliary forcing construction, not scalar backward propagation",
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
