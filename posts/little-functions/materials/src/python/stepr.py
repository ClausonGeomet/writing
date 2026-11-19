"""Scoped Gaussian translation of easyCODA::STEPR method 1.

Method attribution: Michael Greenacre, easyCODA 0.40.2, pinned source commit
3a9ea3953d8a6a7753c2c54298d93ef54558b3d2; Coenders and Greenacre
(2021). The upstream source and this translated selection flow are GPL.
"""
from __future__ import annotations

import math
from numbers import Integral

from lr_ml import _composition, _ratios

__all__ = ["stepr_py", "predict_stepr"]

_BOUNDARY_TOL = 1e-10
_RSS_TOL = 1e-12


def _ols(y, columns):
    import numpy as np

    design = np.column_stack((np.ones(len(y)), *columns))
    coefficients, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    residuals = y - design @ coefficients
    rss = float(residuals @ residuals)
    return coefficients, rank, rss


def stepr_py(X, y, method=1, family="gaussian", criterion="AIC", nsteps=3,
             top=1, previous=None, previousparts=None, denom=None):
    """Select Gaussian method-1 STEPR ratios, with explicit parity boundaries.

    Exact later-step ties differ from R STEPR, which randomly chooses among them.
    """
    import numpy as np

    if (type(method) is not int or method != 1 or family != "gaussian"
            or criterion not in ("AIC", "BIC") or type(top) is not int or top != 1
            or previous is not None or previousparts is not None or denom is not None):
        raise NotImplementedError("UNSUPPORTED STEPR combination")
    values, parts = _composition(X)
    n, d = values.shape
    if not isinstance(nsteps, Integral) or isinstance(nsteps, bool) or not 1 <= nsteps <= min(3, d - 1):
        raise NotImplementedError("UNSUPPORTED nsteps for STEPR scope")
    if n <= nsteps + 2:
        raise ValueError("STEPR requires n > nsteps + 2")
    try:
        response = np.asarray(y, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("y must be finite numeric") from exc
    if response.shape != (n,) or not np.isfinite(response).all():
        raise ValueError("y must contain one finite value per row")
    pairs = [(parts[i], parts[j]) for j in range(1, d) for i in range(j)]
    logs = np.log(values)
    ratio_values = {(a, b): logs[:, parts.index(a)] - logs[:, parts.index(b)] for a, b in pairs}
    total_ss = float(np.sum((response - response.mean()) ** 2))
    selected, path = [], []
    previous_score = None
    tie_sensitive = False
    for step in range(1, nsteps + 1):
        candidates = []
        for pair in pairs:  # STEPR method 1 rescans even previously chosen pairs.
            columns = [ratio_values[p] for p in selected + [pair]]
            _, _, rss = _ols(response, columns)
            if rss <= _RSS_TOL * max(1.0, total_ss):
                raise ValueError("UNSUPPORTED numerical boundary: effectively zero RSS")
            neg2loglik = n * (math.log(2 * math.pi) + 1 + math.log(rss / n))
            candidates.append((neg2loglik, pair, rss))
        candidates.sort(key=lambda item: item[0])  # stable sort keeps pair enumeration.
        best, pair, rss = candidates[0]
        if len(candidates) > 1 and candidates[1][0] - best <= _BOUNDARY_TOL:
            tie_sensitive = True
        aic = best + 2 * (step + 2)
        bic = best + math.log(n) * (step + 2)
        score = aic if criterion == "AIC" else bic
        if previous_score is not None and abs(score - previous_score) <= _BOUNDARY_TOL:
            raise ValueError("UNSUPPORTED numerical boundary: criterion difference near zero")
        accepted = previous_score is None or score < previous_score
        path.append({"step": step, "pair": pair, "rss": rss, "neg2loglik": best,
                     "aic": aic, "bic": bic, "accepted": accepted})
        if not accepted:
            break
        coefficients, rank, retained_rss = _ols(response, [ratio_values[p] for p in selected + [pair]])
        if rank != step + 1:
            raise ValueError("UNSUPPORTED rank-deficient accepted STEPR fit")
        if retained_rss <= _RSS_TOL * max(1.0, total_ss):
            raise ValueError("UNSUPPORTED numerical boundary: effectively zero RSS")
        selected.append(pair)
        previous_score = score
    coefficients, rank, _ = _ols(response, [ratio_values[p] for p in selected])
    if rank != len(selected) + 1:
        raise ValueError("UNSUPPORTED rank-deficient retained STEPR fit")
    return {"parts": parts, "selected_pairs": selected, "path": path,
            "coefficients": coefficients.tolist(), "criterion": criterion,
            "parity_scope": "tie-sensitive" if tie_sensitive else "gaussian-method1-no-ties"}


def predict_stepr(fit, new_composition):
    """Predict from complete compositions using retained linear coefficients."""
    import numpy as np

    values, _ = _composition(new_composition, fit["parts"])
    design = np.column_stack((np.ones(len(values)),
                              *[_ratios(values, fit["parts"], [pair]).ravel()
                                for pair in fit["selected_pairs"]]))
    return design @ np.asarray(fit["coefficients"], dtype=float)
