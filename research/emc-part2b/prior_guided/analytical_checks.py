"""Exact, small arithmetic witnesses for the prior-guided EMC design.

These functions evaluate one assay vector; they do not fit a mineral model.
"""

from __future__ import annotations

from fractions import Fraction
from math import isfinite

MAX_DIMENSION = 4
ZERO = Fraction(0)
ONE = Fraction(1)
HALF = Fraction(1, 2)


def _number(value: object, name: str) -> Fraction:
    if type(value) is int or isinstance(value, Fraction):
        return Fraction(value)
    if type(value) is float and isfinite(value):
        return Fraction(str(value))
    raise ValueError(f"{name} must be a finite int, float, or Fraction")


def _vector(value: object, length: int, name: str) -> list[Fraction]:
    if not isinstance(value, (list, tuple)) or len(value) != length:
        raise ValueError(f"{name} must be a vector of length {length}")
    return [_number(item, name) for item in value]


def _chemistry(value: object) -> list[list[Fraction]]:
    if not isinstance(value, (list, tuple)) or not 1 <= len(value) <= MAX_DIMENSION:
        raise ValueError("C must have one to four rows")
    if not isinstance(value[0], (list, tuple)) or not 1 <= len(value[0]) <= MAX_DIMENSION:
        raise ValueError("C must have one to four columns")
    columns = len(value[0])
    matrix = [_vector(row, columns, "C row") for row in value]
    if any(entry < ZERO for row in matrix for entry in row):
        raise ValueError("C must be nonnegative")
    if any(all(row[j] == ZERO for row in matrix) for j in range(columns)):
        raise ValueError("every C column must have a positive entry")
    return matrix


def _order(value: object, columns: int) -> list[int]:
    if (
        not isinstance(value, (list, tuple))
        or len(value) != columns
        or any(type(item) is not int for item in value)
        or set(value) != set(range(columns))
    ):
        raise ValueError("order must be a permutation of phase indices")
    return list(value)


def _nonnegative(values: list[Fraction], name: str) -> None:
    if any(value < ZERO for value in values):
        raise ValueError(f"{name} must be nonnegative")


def _reconstruct(matrix: list[list[Fraction]], mass: list[Fraction]) -> list[Fraction]:
    return [sum((coefficient * amount for coefficient, amount in zip(row, mass)), ZERO)
            for row in matrix]


def _capacity(
    matrix: list[list[Fraction]], residual: list[Fraction], remaining: Fraction, phase: int
) -> Fraction:
    return min(
        remaining,
        *(residual[i] / row[phase] for i, row in enumerate(matrix) if row[phase] > ZERO),
    )


def decode(
    C: object, a: object, order: object, fractions: object, budget: object = 1
) -> dict[str, object]:
    """Allocate canonical phase fractions in ``order`` within assay and mass limits.

    A negative assay returns an explicit unsupported status after all inputs pass
    validation. No negative channel is clipped or fitted.
    """
    matrix = _chemistry(C)
    channels, phases = len(matrix), len(matrix[0])
    assay = _vector(a, channels, "a")
    sequence = _order(order, phases)
    shares = _vector(fractions, phases, "fractions")
    total_budget = _number(budget, "budget")
    if total_budget < ZERO or any(not ZERO <= share <= ONE for share in shares):
        raise ValueError("budget must be nonnegative and fractions must be in [0, 1]")
    if any(value < ZERO for value in assay):
        return {
            "status": "unsupported_negative_assay",
            "m": None,
            "reconstructed": None,
            "residual": None,
            "remaining": None,
        }

    mass = [ZERO] * phases
    residual = assay.copy()
    remaining = total_budget
    for phase in sequence:
        amount = shares[phase] * _capacity(matrix, residual, remaining, phase)
        mass[phase] = amount
        residual = [value - row[phase] * amount for value, row in zip(residual, matrix)]
        remaining -= amount
    reconstructed = _reconstruct(matrix, mass)
    return {
        "status": "ok",
        "m": mass,
        "reconstructed": reconstructed,
        "residual": [value - fit for value, fit in zip(assay, reconstructed)],
        "remaining": remaining,
    }


def encode(
    C: object, a: object, order: object, target: object, budget: object = 1
) -> list[Fraction]:
    """Find canonical fractions that decode to a feasible target mass vector."""
    matrix = _chemistry(C)
    channels, phases = len(matrix), len(matrix[0])
    assay = _vector(a, channels, "a")
    sequence = _order(order, phases)
    mass = _vector(target, phases, "target")
    total_budget = _number(budget, "budget")
    _nonnegative(assay, "a")
    _nonnegative(mass, "target")
    if total_budget < ZERO:
        raise ValueError("budget must be nonnegative")
    if sum(mass, ZERO) > total_budget or any(
        fit > value for fit, value in zip(_reconstruct(matrix, mass), assay)
    ):
        raise ValueError("target exceeds assay or budget")

    shares = [ZERO] * phases
    residual = assay.copy()
    remaining = total_budget
    for phase in sequence:
        capacity = _capacity(matrix, residual, remaining, phase)
        shares[phase] = mass[phase] / capacity if capacity > ZERO else ZERO
        residual = [value - row[phase] * mass[phase] for value, row in zip(residual, matrix)]
        remaining -= mass[phase]
    return shares


def quadratic_score(
    C: object,
    a: object,
    mass: object,
    weights: object,
    prior_mean: object,
    prior_precision: object,
    strength: object = 1,
) -> dict[str, Fraction]:
    """Evaluate exact diagonal assay and prior quadratic penalties.

    This objective accepts signed assays and masses and imposes no feasibility
    constraint. Precisions are diagonal entries, not standard deviations.
    """
    matrix = _chemistry(C)
    channels, phases = len(matrix), len(matrix[0])
    assay = _vector(a, channels, "a")
    amounts = _vector(mass, phases, "mass")
    assay_weights = _vector(weights, channels, "weights")
    center = _vector(prior_mean, phases, "prior_mean")
    precision = _vector(prior_precision, phases, "prior_precision")
    multiplier = _number(strength, "strength")
    _nonnegative(assay_weights, "weights")
    _nonnegative(precision, "prior_precision")
    if multiplier < ZERO:
        raise ValueError("strength must be nonnegative")

    fitted = _reconstruct(matrix, amounts)
    assay_term = HALF * sum(
        (weight * (value - fit) ** 2
         for weight, value, fit in zip(assay_weights, assay, fitted)), ZERO
    )
    prior_term = multiplier * HALF * sum(
        (weight * (amount - mean) ** 2
         for weight, amount, mean in zip(precision, amounts, center)), ZERO
    )
    return {"assay": assay_term, "prior": prior_term, "total": assay_term + prior_term}
