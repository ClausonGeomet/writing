"""Grouped out-of-fold selection of pairwise compositional log ratios."""
from __future__ import annotations

import math
from numbers import Integral, Real

__all__ = ["select_lr_ml", "predict_lr_ml", "paired_contrast"]

_CLOSURE_TOL = 1e-10
_TIE_TOL = 1e-12


def _composition(X, parts=None):
    import numpy as np
    import pandas as pd

    if not isinstance(X, pd.DataFrame) or X.empty:
        raise ValueError("X must be a nonempty named DataFrame")
    if any(not isinstance(c, str) or not c for c in X.columns) or not X.columns.is_unique:
        raise ValueError("composition columns must have unique nonempty names")
    if parts is not None:
        if set(X.columns) != set(parts) or len(X.columns) != len(parts):
            raise ValueError("composition part names differ from training parts")
        X = X.loc[:, parts]
    if any(not pd.api.types.is_numeric_dtype(dtype) or pd.api.types.is_bool_dtype(dtype)
           for dtype in X.dtypes):
        raise ValueError("composition must be numeric")
    try:
        values = X.to_numpy(dtype=float, copy=True)
    except (TypeError, ValueError) as exc:
        raise ValueError("composition must be numeric") from exc
    if values.shape[1] < 2 or not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("composition needs at least two finite, strictly positive parts")
    if (np.abs(values.sum(axis=1) - 1) > _CLOSURE_TOL).any():
        raise ValueError("composition rows must sum to one")
    return values, list(X.columns)


def _ratios(values, parts, pairs):
    import numpy as np

    indices = [(parts.index(a), parts.index(b)) for a, b in pairs]
    return np.column_stack([np.log(values[:, i]) - np.log(values[:, j]) for i, j in indices])


def _forest(n_trees, n_features, seed):
    from sklearn.ensemble import RandomForestRegressor

    return RandomForestRegressor(
        n_estimators=n_trees, max_features=1.0, min_samples_leaf=5,
        bootstrap=True, n_jobs=1, random_state=seed,
    )


def select_lr_ml(X, y, groups, fold_id, max_ratios=3, min_improvement=0.01,
                 n_trees=50, seed=104729):
    """Select at most three ratios with grouped, pooled out-of-fold MSE."""
    import numpy as np
    import pandas as pd
    import sklearn

    values, parts = _composition(X)
    n, d = values.shape
    if n < 6:
        raise ValueError("at least six rows are required")
    if not isinstance(max_ratios, Integral) or isinstance(max_ratios, bool) or not 0 <= max_ratios <= min(3, d * (d - 1) // 2):
        raise ValueError("max_ratios is outside the allowed integer range")
    if not isinstance(n_trees, Integral) or isinstance(n_trees, bool) or not 1 <= n_trees <= 100:
        raise ValueError("n_trees must be an integer in [1, 100]")
    if not isinstance(seed, Integral) or isinstance(seed, bool) or seed < 0 or seed > 2**31 - 1:
        raise ValueError("seed must be a nonnegative 32-bit integer")
    if not isinstance(min_improvement, Real) or isinstance(min_improvement, bool) or not math.isfinite(min_improvement) or not 0 <= min_improvement <= 1:
        raise ValueError("min_improvement must be finite and in [0, 1]")
    try:
        response = np.asarray(y, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("y must be finite numeric") from exc
    if response.shape != (n,) or not np.isfinite(response).all():
        raise ValueError("y must contain one finite value per row")
    if len(groups) != n or len(fold_id) != n:
        raise ValueError("groups and fold_id must contain one value per row")
    g = pd.Series(groups).reset_index(drop=True)
    f = pd.Series(fold_id).reset_index(drop=True)
    if g.isna().any() or f.isna().any() or any(not isinstance(v, Integral) or isinstance(v, bool) for v in f):
        raise ValueError("groups must be nonmissing and fold_id must be nonmissing integers")
    if (pd.DataFrame({"g": g, "f": f}).groupby("g", sort=False)["f"].nunique() != 1).any():
        raise ValueError("a group spans multiple folds")
    folds = list(dict.fromkeys(f.tolist()))
    if any(seed + int(fold) < 0 or seed + int(fold) > 2**31 - 1 for fold in folds):
        raise ValueError("seed plus fold_id must be a valid forest seed")
    if len(folds) < 2 or any((f == fold).sum() < 2 or (f != fold).sum() < 2 for fold in folds):
        raise ValueError("at least two folds with two training and validation rows each are required")
    pairs = [(parts[i], parts[j]) for j in range(1, d) for i in range(j)]
    log_values = np.log(values)
    ratio_values = {(a, b): log_values[:, parts.index(a)] - log_values[:, parts.index(b)] for a, b in pairs}
    def oof(candidate_pairs):
        pred = np.empty(n)
        for fold in folds:
            train = (f != fold).to_numpy()
            valid = ~train
            if candidate_pairs:
                design = np.column_stack([ratio_values[p] for p in candidate_pairs])
                model = _forest(n_trees, len(candidate_pairs), seed + int(fold))
                model.fit(design[train], response[train])
                pred[valid] = model.predict(design[valid])
            else:
                pred[valid] = response[train].mean()
        return float(np.mean((response - pred) ** 2))
    old_loss = oof([])
    selected = []
    path = []
    for step in range(1, max_ratios + 1):
        if old_loss == 0:
            break
        best, best_loss = None, float("inf")
        for pair in pairs:
            if pair in selected:
                continue
            loss = oof(selected + [pair])
            if loss < best_loss - _TIE_TOL:
                best, best_loss = pair, loss
        gain = (old_loss - best_loss) / old_loss
        accepted = gain > 0 and gain >= min_improvement
        path.append({"step": step, "pair": best, "loss": best_loss,
                     "improvement": gain, "accepted": accepted})
        if not accepted:
            break
        selected.append(best)
        old_loss = best_loss
    model = None
    if selected:
        model = _forest(n_trees, len(selected), seed)
        model.fit(np.column_stack([ratio_values[p] for p in selected]), response)
    return {"parts": parts, "selected_pairs": selected, "path": path,
            "model": model, "mean_y": float(response.mean()),
            "folds": pd.DataFrame({"group": g, "fold_id": f}).drop_duplicates().to_dict("records"),
            "settings": {"max_ratios": max_ratios, "min_improvement": min_improvement,
                         "n_trees": n_trees, "seed": seed, "sklearn_version": sklearn.__version__,
                         "learner": "RandomForestRegressor(max_features=1.0,min_samples_leaf=5,bootstrap=True,n_jobs=1)"}}


def predict_lr_ml(fit, new_composition):
    """Predict from complete, named fractions using the selected pair definitions."""
    import numpy as np

    values, _ = _composition(new_composition, fit["parts"])
    if not fit["selected_pairs"]:
        return np.full(len(values), fit["mean_y"])
    return np.asarray(fit["model"].predict(_ratios(values, fit["parts"], fit["selected_pairs"])), dtype=float)


def paired_contrast(fit, X, delta=0.05, donor="Quartz", recipient="Kaolinite"):
    """Return finite predicted changes after a valid donor-to-recipient exchange."""
    import pandas as pd

    values, parts = _composition(X, fit["parts"])
    if not isinstance(delta, Real) or isinstance(delta, bool) or not math.isfinite(delta) or delta < 0:
        raise ValueError("delta must be finite and nonnegative")
    if donor == recipient or donor not in parts or recipient not in parts:
        raise ValueError("donor and recipient must be distinct training parts")
    if (values[:, parts.index(donor)] <= delta).any():
        raise ValueError("donor must remain strictly positive")
    after = X.loc[:, parts].copy()
    after[donor] -= delta
    after[recipient] += delta
    before_pred = predict_lr_ml(fit, X)
    after_pred = predict_lr_ml(fit, after)
    return pd.DataFrame({"before": before_pred, "after": after_pred,
                         "contrast": after_pred - before_pred}, index=X.index)
