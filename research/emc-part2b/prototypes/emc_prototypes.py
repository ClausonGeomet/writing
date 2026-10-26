"""Small proposed assay allocation and order-search heuristics for EMC feasibility.

Assays and chemistry use the same unnormalised analyte mass-fraction basis.
These heuristics do not infer unique mineralogy or return uncertainty intervals.
"""

from __future__ import annotations

from itertools import permutations

import numpy as np

FEASIBILITY_TOL = 1e-12
TIE_TOL = 1e-12
SCALE_FLOOR = 0.01
HIST_EDGES = np.r_[-np.inf, np.linspace(0.0, 1.0, 9), np.inf]


def _array(value: object, shape: tuple[int | None, ...], name: str) -> np.ndarray:
    try:
        result = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if result.ndim != len(shape) or any(
        size is not None and result.shape[i] != size for i, size in enumerate(shape)
    ) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} has invalid shape or nonfinite entries")
    return result


def _chemistry(C: object) -> np.ndarray:
    matrix = _array(C, (None, None), "C")
    if matrix.shape[0] < 1 or not 1 <= matrix.shape[1] <= 5:
        raise ValueError("C must have p>=1 and 1<=k<=5")
    if np.any(matrix < 0) or np.any(~np.any(matrix > 0, axis=0)):
        raise ValueError("C must be nonnegative with positive support in each column")
    return matrix


def _order(order: object, k: int) -> tuple[int, ...]:
    try:
        values = tuple(order)
    except TypeError as exc:
        raise ValueError("order must be a permutation") from exc
    if len(values) != k or any(type(i) not in (int, np.int64, np.int32) for i in values):
        raise ValueError("order must be a permutation of integer phase indices")
    if set(values) != set(range(k)):
        raise ValueError("order must be a permutation of integer phase indices")
    return tuple(int(i) for i in values)


def _positive_vector(value: object, p: int, name: str) -> np.ndarray:
    vector = _array(value, (p,), name)
    if np.any(vector <= 0):
        raise ValueError(f"{name} must be positive")
    return vector


def _matrix(value: object, p: int, name: str, minimum_rows: int) -> np.ndarray:
    matrix = _array(value, (None, p), name)
    if matrix.shape[0] < minimum_rows:
        raise ValueError(f"{name} needs at least {minimum_rows} rows")
    return matrix


def _positive_integer(value: object, name: str, minimum: int = 1) -> int:
    if type(value) not in (int, np.int64, np.int32) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def allocate(C: object, a: object, order: object) -> dict:
    """Allocate one assay (p,) to k whole-sample phase fractions in order.

    Return mass fractions, reconstructed assay, original-basis residual, remaining
    unallocated mass, status, and largest arithmetic correction. Assays below
    -1e-12 return ``unsupported_negative_assay`` with NaN-valued results.
    """
    matrix = _chemistry(C)
    assay = _array(a, (matrix.shape[0],), "a")
    sequence = _order(order, matrix.shape[1])
    if np.any(assay < -FEASIBILITY_TOL):
        return {
            "m": np.full(matrix.shape[1], np.nan),
            "reconstructed": np.full(matrix.shape[0], np.nan),
            "residual": np.full(matrix.shape[0], np.nan),
            "unallocated": float("nan"),
            "status": "unsupported_negative_assay",
            "correction_max": 0.0,
        }
    residual = assay.copy()
    correction = float(max(0.0, -np.min(residual)))
    residual = np.maximum(residual, 0.0)
    m = np.zeros(matrix.shape[1], dtype=np.float64)
    remaining = 1.0
    for j in sequence:
        positive = matrix[:, j] > 0
        delta = min(remaining, float(np.min(residual[positive] / matrix[positive, j])))
        m[j] = delta
        residual -= delta * matrix[:, j]
        remaining -= delta
        minimum = min(float(np.min(residual)), remaining)
        if minimum < -FEASIBILITY_TOL:
            raise ArithmeticError("allocation exceeded available assay or mass")
        correction = max(correction, -minimum, 0.0)
        residual = np.maximum(residual, 0.0)
        remaining = max(remaining, 0.0)
    reconstructed = matrix @ m
    return {
        "m": m,
        "reconstructed": reconstructed,
        "residual": assay - reconstructed,
        "unallocated": float(remaining),
        "status": "ok",
        "correction_max": correction,
    }


def predict_order(C: object, A: object, order: object) -> dict:
    """Apply one phase order to assay rows (n,p), preserving unsupported status."""
    matrix = _chemistry(C)
    assays = _matrix(A, matrix.shape[0], "A", 1)
    sequence = _order(order, matrix.shape[1])
    rows = [allocate(matrix, assay, sequence) for assay in assays]
    return {
        "M": np.stack([row["m"] for row in rows]),
        "A_hat": np.stack([row["reconstructed"] for row in rows]),
        "residual": np.stack([row["residual"] for row in rows]),
        "unallocated": np.array([row["unallocated"] for row in rows]),
        "status": [row["status"] for row in rows],
    }


def score_assays(
    A: object, A_hat: object, *, sigma: object = None,
    objective: str = "row", scales: object = None,
) -> dict:
    """Score observed and reconstructed (n,p) mass fractions, with n>=2.

    ``row`` is mean squared assay error scaled by sigma. ``distribution``
    additionally compares training-scale means, covariances and histogram JS
    divergence. Supplied scales permit validation without refitting on test rows.
    """
    observed = _array(A, (None, None), "A")
    if observed.shape[0] < 2 or observed.shape[1] < 1:
        raise ValueError("A needs n>=2 and p>=1")
    fitted = _matrix(A_hat, observed.shape[1], "A_hat", 2)
    if fitted.shape != observed.shape:
        raise ValueError("A_hat must match A")
    if objective not in ("row", "distribution"):
        raise ValueError("unsupported objective")
    weights = np.ones(observed.shape[1]) if sigma is None else _positive_vector(sigma, observed.shape[1], "sigma")
    training_scales = (
        np.maximum(np.std(observed, axis=0, ddof=1), SCALE_FLOOR)
        if scales is None else _positive_vector(scales, observed.shape[1], "scales")
    )
    row = float(np.mean(((observed - fitted) / weights) ** 2))
    mean = float(np.mean(((observed.mean(axis=0) - fitted.mean(axis=0)) / training_scales) ** 2))
    obs_cov = np.atleast_2d(np.cov(observed, rowvar=False, ddof=1))
    fit_cov = np.atleast_2d(np.cov(fitted, rowvar=False, ddof=1))
    covariance = float(np.mean(((obs_cov - fit_cov) / np.outer(training_scales, training_scales)) ** 2))
    divergences = []
    for i in range(observed.shape[1]):
        p = np.histogram(observed[:, i], bins=HIST_EDGES)[0] / observed.shape[0]
        q = np.histogram(fitted[:, i], bins=HIST_EDGES)[0] / fitted.shape[0]
        mix = (p + q) / 2
        divergences.append(0.5 * (np.sum(p[p > 0] * np.log(p[p > 0] / mix[p > 0])) +
                                  np.sum(q[q > 0] * np.log(q[q > 0] / mix[q > 0]))))
    js = float(np.mean(divergences))
    return {
        "row": row, "mean": mean, "covariance": covariance, "js": js,
        "total": row if objective == "row" else row + mean + covariance + js,
        "scales": training_scales.copy(),
    }


def _better(score: float, order: tuple[int, ...], best_score: float, best_order: tuple[int, ...] | None) -> bool:
    return best_order is None or score < best_score - TIE_TOL or (
        abs(score - best_score) <= TIE_TOL and order < best_order
    )


def fit_order(
    C: object, A_train: object, *, method: str = "exhaustive",
    objective: str = "row", sigma: object = None, seed: int = 0,
    budget: int = 24, beam_width: int = 2,
) -> dict:
    """Fit one proposed global order using training assays only, without labels.

    Search methods are fixed, exhaustive, greedy, beam, random, and ga. The
    reported trace counts proposals, including cached duplicates, and actual
    score evaluations. Budget constrains proposals except exhaustive search.
    """
    matrix = _chemistry(C)
    assays = _matrix(A_train, matrix.shape[0], "A_train", 2)
    if method not in ("fixed", "exhaustive", "greedy", "beam", "random", "ga"):
        raise ValueError("unsupported method")
    if objective not in ("row", "distribution"):
        raise ValueError("unsupported objective")
    weights = np.ones(matrix.shape[0]) if sigma is None else _positive_vector(sigma, matrix.shape[0], "sigma")
    seed = _positive_integer(seed, "seed", 0)
    budget = _positive_integer(budget, "budget")
    beam_width = _positive_integer(beam_width, "beam_width")
    scales = np.maximum(np.std(assays, axis=0, ddof=1), SCALE_FLOOR)
    if np.any(assays < -FEASIBILITY_TOL):
        return {
            "order": None, "train_score": float("nan"), "components": {},
            "trace": [], "evaluations": 0, "proposals": 0,
            "status": "unsupported_negative_assay", "seed": seed,
            "scales": scales.copy(),
        }
    cache: dict[tuple[int, ...], dict] = {}
    trace: list[dict] = []
    best_order: tuple[int, ...] | None = None
    best_score = float("inf")

    def propose(candidate: tuple[int, ...]) -> float:
        nonlocal best_order, best_score
        if candidate not in cache:
            fitted = predict_order(matrix, assays, candidate)["A_hat"]
            cache[candidate] = score_assays(
                assays, fitted, sigma=weights, objective=objective, scales=scales
            )
        score = cache[candidate]["total"]
        if _better(score, candidate, best_score, best_order):
            best_order, best_score = candidate, score
        trace.append({
            "proposal": len(trace) + 1, "evaluations": len(cache),
            "best_score": best_score, "order": best_order,
        })
        return score

    lexical = tuple(range(matrix.shape[1]))
    rng = np.random.Generator(np.random.PCG64(seed))
    if method == "fixed":
        propose(lexical)
    elif method == "exhaustive":
        for candidate in permutations(lexical):
            propose(candidate)
    elif method == "random":
        for _ in range(budget):
            propose(tuple(int(i) for i in rng.permutation(lexical)))
    elif method == "greedy":
        current = lexical
        current_score = propose(current)
        while len(trace) < budget:
            next_order, next_score = current, current_score
            for i in range(len(lexical) - 1):
                if len(trace) >= budget:
                    break
                swapped = list(current)
                swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
                candidate = tuple(swapped)
                score = propose(candidate)
                if score < next_score - TIE_TOL or (
                    abs(score - next_score) <= TIE_TOL and score < current_score - TIE_TOL
                    and candidate < next_order
                ):
                    next_order, next_score = candidate, score
            if next_order == current:
                break
            current, current_score = next_order, next_score
    elif method == "beam":
        prefixes: list[tuple[int, ...]] = [()]
        for _ in lexical:
            ranked = []
            for prefix in prefixes:
                for phase in lexical:
                    if phase in prefix:
                        continue
                    child = prefix + (phase,)
                    completed = child + tuple(i for i in lexical if i not in child)
                    if len(trace) >= budget:
                        break
                    ranked.append((propose(completed), child))
                if len(trace) >= budget:
                    break
            if not ranked:
                break
            # Scores closer than tolerance retain lexical priority.
            ranked.sort(key=lambda item: item[1])
            ordered = []
            for item in ranked:
                pos = next((i for i, old in enumerate(ordered)
                            if item[0] < old[0] - TIE_TOL), len(ordered))
                ordered.insert(pos, item)
            prefixes = [child for _, child in ordered[:beam_width]]
    else:  # ga
        population = [lexical, lexical[::-1]]
        population.extend(tuple(int(i) for i in rng.permutation(lexical)) for _ in range(4))
        for generation in range(4):
            scored = []
            for candidate in population:
                if len(trace) >= budget:
                    break
                scored.append((propose(candidate), candidate))
            if len(trace) >= budget or generation == 3:
                break
            scored.sort(key=lambda item: item[1])
            ranked = []
            for item in scored:
                pos = next((i for i, old in enumerate(ranked)
                            if item[0] < old[0] - TIE_TOL), len(ranked))
                ranked.insert(pos, item)
            elites = [candidate for _, candidate in ranked[:2]]
            def parent() -> tuple[int, ...]:
                draws = rng.integers(len(ranked), size=2)
                return ranked[min(draws)][1]
            children = []
            for _ in range(4):
                first, second = parent(), parent()
                if len(lexical) == 1:
                    child = first
                else:
                    left, right = sorted(int(x) for x in rng.choice(len(lexical) + 1, size=2, replace=False))
                    segment = first[left:right]
                    unused = iter(i for i in second if i not in segment)
                    child = tuple(segment[i - left] if left <= i < right else next(unused)
                                  for i in range(len(lexical)))
                    if rng.random() < 0.3:
                        i, j = rng.choice(len(lexical), size=2, replace=False)
                        values = list(child)
                        values[i], values[j] = values[j], values[i]
                        child = tuple(values)
                children.append(child)
            population = elites + children
    assert best_order is not None
    return {
        "order": best_order, "train_score": best_score,
        "components": cache[best_order].copy(), "trace": trace,
        "evaluations": len(cache), "proposals": len(trace), "status": "ok",
        "seed": seed, "scales": scales.copy(),
    }
