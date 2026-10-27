"""Finite calibration and globally solved QAR quadratic estimators.

Masses are whole-sample fractions in Quartz, Anatase, Rutile order. The
learned quadratic is a plug-in regularizer; it is not a Bayesian posterior.
"""

from __future__ import annotations

import numbers

import numpy as np


C = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 1.0]])
U = np.array([[2 / np.sqrt(6), 0.0],
              [-1 / np.sqrt(6), 1 / np.sqrt(2)],
              [-1 / np.sqrt(6), -1 / np.sqrt(2)]])
PSD_TOL = 1e-12
FEAS_TOL = 1e-9
KKT_TOL = 1e-8
ACTIVE_TOL = 1e-9


def _array(value: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    raw = np.asarray(value)
    original = np.asarray(value, dtype=object)
    if any(not isinstance(x, numbers.Real) or isinstance(x, (bool, np.bool_))
           for x in original.flat):
        raise ValueError(f'{name} must contain finite real numbers, not booleans')
    if raw.shape != shape or raw.dtype.kind not in 'iuf':
        raise ValueError(f'{name} must be a finite real array of shape {shape}')
    result = raw.astype(np.float64)
    if not np.isfinite(result).all():
        raise ValueError(f'{name} must be finite')
    return result


def _scalar(value: object, name: str, lower: float | None = None,
            upper: float | None = None, strict_lower: bool = False) -> float:
    raw = np.asarray(value)
    if not isinstance(value, numbers.Real) or isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{name} must be a finite real scalar')
    if raw.shape != () or raw.dtype.kind not in 'iuf':
        raise ValueError(f'{name} must be a finite real scalar')
    result = float(raw)
    if not np.isfinite(result) or (lower is not None and
            (result <= lower if strict_lower else result < lower)) or (
            upper is not None and result > upper):
        raise ValueError(f'{name} is outside its domain')
    return result


def _cov(value: object, size: int, name: str, positive: bool = False) -> np.ndarray:
    matrix = _array(value, (size, size), name)
    if np.max(np.abs(matrix - matrix.T)) > PSD_TOL:
        raise ValueError(f'{name} must be symmetric')
    matrix = (matrix + matrix.T) / 2
    eig = np.linalg.eigvalsh(matrix)
    if eig[0] < (-PSD_TOL if not positive else PSD_TOL):
        raise ValueError(f'{name} must be {"positive definite" if positive else "PSD"}')
    if positive:
        return matrix
    # Only numerical roundoff is removed; materially negative eigenvalues fail.
    if eig[0] < 0:
        d, v = np.linalg.eigh(matrix)
        matrix = (v * np.maximum(d, 0)) @ v.T
    return matrix


def _constraints(support: object, fixed_quartz: object) -> tuple[np.ndarray, float | None, np.ndarray, np.ndarray]:
    if not isinstance(support, (tuple, list, np.ndarray)) or np.asarray(support).shape != (3,) or not all(
            isinstance(x, (bool, np.bool_)) for x in support):
        raise ValueError('support must contain three booleans')
    allowed = np.asarray(support, dtype=bool)
    if not allowed.any():
        raise ValueError('support must be nonempty')
    q = None if fixed_quartz is None else _scalar(fixed_quartz, 'fixed_quartz', 0, 1)
    if q is not None and ((not allowed[0] and q != 0) or
                          (not allowed[1:3].any() and q != 1)):
        raise ValueError('fixed_quartz is infeasible with support')
    rows = [np.ones(3)]
    rhs = [1.0]
    if q is not None:
        rows.append(np.array([1.0, 0.0, 0.0]))
        rhs.append(q)
    for j in np.flatnonzero(~allowed):
        row = np.zeros(3)
        row[j] = 1.0
        rows.append(row)
        rhs.append(0.0)
    return allowed, q, np.array(rows), np.array(rhs)


def _violation(m: np.ndarray, allowed: np.ndarray, E: np.ndarray,
               b: np.ndarray) -> float:
    return float(max(np.max(np.abs(E @ m - b)), max(0.0, -float(np.min(m))),
                     np.max(np.abs(m[~allowed])) if (~allowed).any() else 0.0))


def _independent(E: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rows: list[np.ndarray] = []
    rhs: list[float] = []
    for row, value in zip(E, b):
        trial = np.array(rows + [row])
        if np.linalg.matrix_rank(trial, tol=PSD_TOL) > len(rows):
            rows.append(row)
            rhs.append(float(value))
    return np.array(rows), np.array(rhs)


def _face_solution(H: np.ndarray, h: np.ndarray, E: np.ndarray,
                   b: np.ndarray) -> np.ndarray | None:
    A, rhs = _independent(E, b)
    block = np.block([[H, A.T], [A, np.zeros((len(A), len(A)))]])
    solution = np.linalg.lstsq(block, np.concatenate([h, rhs]), rcond=None)[0]
    m = solution[:3]
    if np.max(np.abs(A @ m - rhs)) > FEAS_TOL or np.max(
            np.abs(H @ m + A.T @ solution[3:] - h)) > KKT_TOL * max(
            1.0, np.linalg.norm(H, np.inf) + np.linalg.norm(h, np.inf)):
        return None
    return m


def _project(y: np.ndarray, allowed: np.ndarray, E: np.ndarray,
             b: np.ndarray) -> np.ndarray | None:
    best = None
    best_loss = np.inf
    indices = np.flatnonzero(allowed)
    for mask in range(1, 1 << len(indices)):
        zero = [j for bit, j in enumerate(indices) if not mask & (1 << bit)]
        extra = np.eye(3)[zero]
        candidate = _face_solution(np.eye(3), y, np.vstack((E, extra)),
                                   np.concatenate((b, np.zeros(len(zero)))))
        if candidate is None or _violation(candidate, allowed, E, b) > FEAS_TOL:
            continue
        loss = float(np.sum((candidate - y) ** 2))
        if loss < best_loss:
            best, best_loss = candidate, loss
    return best


def _failure_prior(diagnostics: dict) -> dict:
    return dict(status='numerical_failure', center=None, tangent_cov=None,
                mean_cov=None, precision=None, diagnostics=diagnostics)


def fit_prior(y_cal: object, measurement_cov: object, *,
              support: object = (True, True, True), fixed_quartz: object = None,
              shrinkage: object = 0.25, variance_floor: object = 0.0001) -> dict:
    """Fit a projected center and measurement-error-adjusted tangent penalty."""
    raw = np.asarray(y_cal)
    if raw.ndim != 2 or raw.shape[1] != 3 or raw.shape[0] < 2:
        raise ValueError('y_cal must be n by 3 with n >= 2')
    y = _array(y_cal, raw.shape, 'y_cal')
    R = _cov(measurement_cov, 3, 'measurement_cov')
    allowed, _, E, b = _constraints(support, fixed_quartz)
    alpha = _scalar(shrinkage, 'shrinkage', 0, 1)
    floor = _scalar(variance_floor, 'variance_floor', 0, strict_lower=True)
    n = len(y)
    ybar = y.mean(axis=0)
    Sz = np.cov(y @ U, rowvar=False, ddof=1)
    Rz = U.T @ R @ U
    d, V = np.linalg.eigh((Sz - Rz + (Sz - Rz).T) / 2)
    df = np.maximum(d, floor)
    G = (V * df) @ V.T
    T = (1 - alpha) * G + alpha * np.trace(G) * np.eye(2) / 2
    M = (T + Rz) / n
    B = T + M
    center = _project(ybar, allowed, E, b)
    diagnostics = dict(n=n, raw_mean=ybar.tolist(), raw_tangent_cov=Sz.tolist(),
                       deconvolved_eigenvalues=d.tolist(),
                       floored_eigenvalues=df.tolist(),
                       projection_distance=None if center is None else float(
                           np.linalg.norm(center - ybar)),
                       shrinkage=alpha, variance_floor=floor)
    if center is None or not all(np.isfinite(x).all() for x in (T, M, B)):
        return _failure_prior(diagnostics)
    try:
        K = U @ np.linalg.solve(B, U.T)
    except np.linalg.LinAlgError:
        return _failure_prior(diagnostics)
    if not np.isfinite(K).all():
        return _failure_prior(diagnostics)
    return dict(status='ok', center=center.tolist(), tangent_cov=T.tolist(),
                mean_cov=M.tolist(), precision=K.tolist(), diagnostics=diagnostics)


def _kkt(m: np.ndarray, H: np.ndarray, h: np.ndarray, allowed: np.ndarray,
         E: np.ndarray, b: np.ndarray) -> float:
    A, _ = _independent(E, b)
    active = [j for j in np.flatnonzero(allowed) if m[j] <= ACTIVE_TOL]
    g = H @ m - h
    scale = max(1.0, np.linalg.norm(H, np.inf) * max(1.0, np.linalg.norm(m, np.inf))
                + np.linalg.norm(h, np.inf))
    best = np.inf
    # Every possible positive-multiplier subset is checked. This tiny exact
    # enumeration avoids relying on an iterative bounded least-squares status.
    for mask in range(1 << len(active)):
        selected = [j for bit, j in enumerate(active) if mask & (1 << bit)]
        system = np.hstack((A.T, -np.eye(3)[:, selected]))
        multipliers = np.linalg.lstsq(system, -g, rcond=None)[0]
        z = np.zeros(3)
        z[selected] = multipliers[len(A):]
        stationarity = np.max(np.abs(g + A.T @ multipliers[:len(A)] - z))
        dual = max(0.0, -float(np.min(z)))
        complementary = np.max(np.abs(z * m))
        best = min(best, max(stationarity, dual, complementary) / scale)
    return float(best)


def _failure_estimate() -> dict:
    return dict(status='numerical_failure', mass=None, reconstructed=None,
                residual=None, assay_loss=None, prior_loss=None, objective=None,
                constraint_violation=None, kkt_residual=None, active_support=[])


def estimate(a: object, assay_cov: object, prior: object, *, method: str,
             support: object = (True, True, True), fixed_quartz: object = None,
             strength: object = 1.0) -> dict:
    """Solve A0, center-only A1, or global constrained quadratic A2."""
    if not isinstance(method, str) or method not in ('A0', 'A1', 'A2'):
        raise ValueError('method must be A0, A1, or A2')
    allowed, q_fixed, E, b = _constraints(support, fixed_quartz)
    lam = _scalar(strength, 'strength', 0)
    assay = None if a is None else _array(a, (2,), 'a')
    cov = None if assay_cov is None else _cov(assay_cov, 2, 'assay_cov', positive=True)
    if method != 'A1' and (assay is None or cov is None):
        raise ValueError('A0 and A2 require assay and assay_cov')
    if method == 'A0' and prior is not None:
        raise ValueError('A0 requires prior=None')
    center = K = None
    if method != 'A0':
        if not isinstance(prior, dict) or prior.get('status') != 'ok':
            raise ValueError('A1 and A2 require a successful fitted prior')
        center = _array(prior.get('center'), (3,), 'prior center')
        K = _cov(prior.get('precision'), 3, 'prior precision')
        if _violation(center, allowed, E, b) > FEAS_TOL:
            raise ValueError('prior center is infeasible')
    W = np.linalg.inv(cov) if cov is not None else None
    if method == 'A1':
        m = center.copy()
        H = h = None
    elif method == 'A0' or lam == 0:
        # The assay observes q and total Ti only. Equal splitting gives the
        # unique least-Euclidean-norm element of the full minimizer set.
        if q_fixed is not None:
            q = q_fixed
        elif not allowed[0]:
            q = 0.0
        elif not allowed[1:3].any():
            q = 1.0
        else:
            direction = np.array([1.0, -1.0])
            q_unconstrained = float(direction @ W @ (assay - np.array([0.0, 1.0]))
                                    / (direction @ W @ direction))
            q = min(1.0, max(0.0, q_unconstrained))
        m = np.zeros(3)
        m[0] = q
        ti_indices = np.flatnonzero(allowed[1:3]) + 1
        if len(ti_indices):
            m[ti_indices] = (1 - q) / len(ti_indices)
        H = C.T @ W @ C
        h = C.T @ W @ assay
    else:
        H = C.T @ W @ C + lam * K
        h = C.T @ W @ assay + lam * K @ center
        m = None
        best_loss = np.inf
        indices = np.flatnonzero(allowed)
        for mask in range(1, 1 << len(indices)):
            zero = [j for bit, j in enumerate(indices) if not mask & (1 << bit)]
            candidate = _face_solution(H, h, np.vstack((E, np.eye(3)[zero])),
                                       np.concatenate((b, np.zeros(len(zero)))))
            if candidate is None or _violation(candidate, allowed, E, b) > FEAS_TOL:
                continue
            residual = assay - C @ candidate
            delta = candidate - center
            loss = float(0.5 * residual @ W @ residual + 0.5 * lam * delta @ K @ delta)
            if loss < best_loss:
                m, best_loss = candidate, loss
        if m is None:
            return _failure_estimate()
    violation = _violation(m, allowed, E, b)
    kkt = 0.0 if method == 'A1' and violation <= FEAS_TOL else _kkt(
        m, H, h, allowed, E, b)
    if violation > FEAS_TOL or kkt > KKT_TOL or not np.isfinite(kkt):
        return _failure_estimate()
    reconstructed = C @ m
    residual = None if assay is None else assay - reconstructed
    assay_loss = None if method == 'A1' else float(0.5 * residual @ W @ residual)
    prior_loss = 0.0 if method == 'A0' else float(
        0.5 * lam * (m - center) @ K @ (m - center))
    objective = prior_loss if method == 'A1' else assay_loss + prior_loss
    return dict(status='ok', mass=m.tolist(), reconstructed=reconstructed.tolist(),
                residual=None if residual is None else residual.tolist(),
                assay_loss=assay_loss, prior_loss=prior_loss, objective=objective,
                constraint_violation=violation, kkt_residual=kkt,
                active_support=[int(j) for j in np.flatnonzero(m > ACTIVE_TOL)])
