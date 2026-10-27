"""Three-phase, full-covariance convex estimators for the bounded QAR study.

All masses are whole-sample fractions in Quartz, Anatase, Rutile order.
The fitted prior is a working quadratic penalty, not a calibrated posterior.
"""
from __future__ import annotations

import numbers
from collections.abc import Mapping, Sequence

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog, nnls

B = np.full(3, 1.0 / 3.0)
U = np.array([[2 / np.sqrt(6), 0], [-1 / np.sqrt(6), 1 / np.sqrt(2)],
              [-1 / np.sqrt(6), -1 / np.sqrt(2)]])
SYMMETRY_TOL = 1e-12
FEAS_TOL = 1e-9
KKT_TOL = 1e-8
RANK_TOL = 1e-10
REDUCED_RANK_TOL = 1e-12
TIE_TOL = 1e-12
ACTIVE_TOL = 1e-9


def _array(value: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    try:
        raw = np.asarray(value)
        objects = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError(f'{name} must be a rectangular real array') from exc
    if raw.shape != shape or raw.dtype.kind not in 'iuf' or any(
        not isinstance(x, numbers.Real) or isinstance(x, (bool, np.bool_))
        for x in objects.flat
    ):
        raise ValueError(f'{name} must be a real array of shape {shape}')
    try:
        result = raw.astype(np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'{name} cannot be represented as float64') from exc
    if not np.isfinite(result).all():
        raise ValueError(f'{name} must be finite')
    return result


def _scalar(value: object, name: str, low: float | None = None,
            high: float | None = None, strict: bool = False) -> float:
    if not isinstance(value, numbers.Real) or isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{name} must be a real scalar')
    result = float(value)
    if not np.isfinite(result) or (low is not None and
        (result <= low if strict else result < low)) or (high is not None and result > high):
        raise ValueError(f'{name} is outside its domain')
    return result


def _symmetric(value: object, n: int, name: str, *, positive: bool = False,
               semidefinite: bool = False) -> np.ndarray:
    x = _array(value, (n, n), name)
    if np.max(np.abs(x - x.T)) > SYMMETRY_TOL:
        raise ValueError(f'{name} is nonsymmetric')
    x = (x + x.T) / 2
    least = np.linalg.eigvalsh(x)[0]
    if positive and least <= SYMMETRY_TOL:
        raise ValueError(f'{name} must be positive definite')
    if semidefinite and least < -SYMMETRY_TOL:
        raise ValueError(f'{name} must be positive semidefinite')
    return x


def _support(value: object) -> np.ndarray:
    if not isinstance(value, (tuple, list, np.ndarray)):
        raise ValueError('support must contain three booleans')
    try:
        if np.asarray(value).shape != (3,) or not all(
            isinstance(x, (bool, np.bool_)) for x in value
        ):
            raise ValueError('support must contain three booleans')
    except (TypeError, ValueError) as exc:
        raise ValueError('support must contain three booleans') from exc
    result = np.asarray(value, dtype=bool)
    if not result.any():
        raise ValueError('support must be nonempty')
    return result


def _constraints(support: object, fixed_quartz: object = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    allowed = _support(support)
    q = None if fixed_quartz is None else _scalar(fixed_quartz, 'fixed_quartz', 0, 1)
    if q is not None and ((not allowed[0] and q != 0) or
                          (not allowed[1:].any() and q != 1)):
        raise ValueError('fixed_quartz is infeasible')
    E = [np.ones(3)]
    rhs = [1.0]
    if q is not None:
        E.append(np.array([1.0, 0, 0]))
        rhs.append(q)
    for j in np.flatnonzero(~allowed):
        E.append(np.eye(3)[j])
        rhs.append(0.0)
    return allowed, np.array(E), np.array(rhs)


def _independent(E: np.ndarray, rhs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rows: list[np.ndarray] = []
    vals: list[float] = []
    for row, val in zip(E, rhs):
        trial = np.array(rows + [row])
        if np.linalg.matrix_rank(trial, tol=SYMMETRY_TOL) > len(rows):
            rows.append(row)
            vals.append(float(val))
    return np.array(rows), np.array(vals)


def _violation(m: np.ndarray, allowed: np.ndarray, E: np.ndarray, rhs: np.ndarray) -> float:
    return float(max(np.max(np.abs(E @ m - rhs)), max(0.0, -np.min(m)),
                     np.max(np.abs(m[~allowed])) if (~allowed).any() else 0.0))


def _face(H: np.ndarray, f: np.ndarray, E: np.ndarray, rhs: np.ndarray) -> np.ndarray | None:
    A, v = _independent(E, rhs)
    m0 = np.linalg.lstsq(A, v, rcond=None)[0]
    if np.max(np.abs(E @ m0 - rhs)) > FEAS_TOL:
        return None
    N = null_space(A, rcond=SYMMETRY_TOL)
    if N.shape[1] == 0:
        m = m0
    else:
        Q = N.T @ H @ N
        Q = (Q + Q.T) / 2
        eigen, V = np.linalg.eigh(Q)
        # A restriction of a singular Hessian can acquire a tiny eigenvalue
        # from forming H and rotating into the face. Compare with the full
        # Hessian scale: a relative cutoff based only on Q would call a
        # one-dimensional roundoff eigenvalue identified.
        cutoff = REDUCED_RANK_TOL * np.linalg.norm(H, ord=2)
        if eigen[0] < -cutoff:
            return None
        identified = eigen > cutoff
        gradient = N.T @ (f - H @ m0)
        z = V[:, identified] @ ((V[:, identified].T @ gradient) / eigen[identified])
        m = m0 + N @ z
        # m0 is orthogonal to ker(A), so the least-norm reduced solution
        # gives the minimum Euclidean norm on this stationary affine face.
        if np.max(np.abs(Q @ z - gradient)) > KKT_TOL * max(
            1., np.linalg.norm(H, np.inf) + np.linalg.norm(f, np.inf)
        ):
            return None
    return m if np.isfinite(m).all() else None


def _kkt(m: np.ndarray, H: np.ndarray, f: np.ndarray, allowed: np.ndarray,
         E: np.ndarray, rhs: np.ndarray) -> float:
    A, _ = _independent(E, rhs)
    N = null_space(A, rcond=SYMMETRY_TOL)
    active = [int(j) for j in np.flatnonzero(allowed) if m[j] <= ACTIVE_TOL]
    g = H @ m - f
    if active and N.shape[1]:
        Z = N.T @ np.eye(3)[:, active]
        z_active = nnls(Z, N.T @ g)[0]
    else:
        z_active = np.zeros(len(active))
    z = np.zeros(3)
    z[active] = z_active
    nu = np.linalg.lstsq(A.T, z - g, rcond=None)[0]
    scale = max(1., np.linalg.norm(H, np.inf) * max(1., np.linalg.norm(m, np.inf))
                + np.linalg.norm(f, np.inf))
    return float(max(np.max(np.abs(g + A.T @ nu - z)),
                     np.max(np.abs(z * m))) / scale)


def _solve(H: np.ndarray, f: np.ndarray, allowed: np.ndarray, E: np.ndarray,
           rhs: np.ndarray) -> tuple[np.ndarray, float, float] | None:
    best: tuple[np.ndarray, float, float] | None = None
    best_obj = np.inf
    indices = np.flatnonzero(allowed)
    for mask in range(1, 1 << len(indices)):
        zeros = [j for bit, j in enumerate(indices) if not (mask & (1 << bit))]
        F = np.vstack((E, np.eye(3)[zeros]))
        v = np.concatenate((rhs, np.zeros(len(zeros))))
        m = _face(H, f, F, v)
        if m is None:
            continue
        violation = _violation(m, allowed, E, rhs)
        if violation > FEAS_TOL:
            continue
        cert = _kkt(m, H, f, allowed, E, rhs)
        if cert > KKT_TOL or not np.isfinite(cert):
            continue
        obj = float(0.5 * m @ H @ m - f @ m)
        if best is None or obj < best_obj - TIE_TOL * (1 + abs(best_obj)) or (
            abs(obj - best_obj) <= TIE_TOL * (1 + abs(best_obj)) and
            np.linalg.norm(m) < np.linalg.norm(best[0])
        ):
            best = (m, violation, cert)
            best_obj = obj
    return best


def _failure_prior(reason: str, diagnostics: dict) -> dict:
    return dict(status='numerical_failure', reason=reason, center=None,
                tangent_cov=None, mean_cov=None, precision=None, diagnostics=diagnostics)


def _prior(center: np.ndarray, T: np.ndarray, M: np.ndarray, diagnostics: dict) -> dict:
    try:
        K = U @ np.linalg.solve(T + M, U.T)
    except np.linalg.LinAlgError:
        return _failure_prior('prior precision solve failed', diagnostics)
    if not np.isfinite(K).all():
        return _failure_prior('nonfinite prior precision', diagnostics)
    return dict(status='ok', reason='', center=center.tolist(), tangent_cov=T.tolist(),
                mean_cov=M.tolist(), precision=((K + K.T) / 2).tolist(),
                diagnostics=diagnostics)


def _validate_prior(prior: object, allowed: np.ndarray, E: np.ndarray,
                    rhs: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not isinstance(prior, Mapping) or prior.get('status') != 'ok' or not isinstance(
        prior.get('diagnostics'), Mapping
    ):
        raise ValueError('successful prior with diagnostics required')
    center = _array(prior.get('center'), (3,), 'prior.center')
    T = _symmetric(prior.get('tangent_cov'), 2, 'prior.tangent_cov', positive=True)
    M = _symmetric(prior.get('mean_cov'), 2, 'prior.mean_cov', positive=True)
    K = _symmetric(prior.get('precision'), 3, 'prior.precision', semidefinite=True)
    expected = U @ np.linalg.solve(T + M, U.T)
    if np.max(np.abs(K - expected)) > SYMMETRY_TOL:
        raise ValueError('prior precision is inconsistent with T and M')
    diagnostics = prior['diagnostics']
    required = ('n_cal', 'raw_mean', 'raw_tangent_cov', 'deconvolved_eigenvalues',
                'floored_eigenvalues', 'projection_distance', 'included_unit_ids',
                'excluded_unit_id', 'n_update', 'mean_information', 'mean_score',
                'mean_constraint_violation', 'mean_kkt_residual')
    if any(key not in diagnostics for key in required):
        raise ValueError('prior diagnostics are incomplete')
    info = _symmetric(diagnostics['mean_information'], 2, 'prior.mean_information', positive=True)
    if np.max(np.abs(info @ M - np.eye(2))) > SYMMETRY_TOL:
        raise ValueError('prior mean information is inconsistent with M')
    _array(diagnostics['mean_score'], (2,), 'prior.mean_score')
    if not isinstance(diagnostics['included_unit_ids'], list) or not all(
        isinstance(x, str) and x for x in diagnostics['included_unit_ids']
    ) or diagnostics['excluded_unit_id'] is not None and not isinstance(
        diagnostics['excluded_unit_id'], str
    ):
        raise ValueError('prior unit provenance is invalid')
    if _violation(center, allowed, E, rhs) > FEAS_TOL:
        raise ValueError('prior center is infeasible on requested support')
    return center, T, M


def fit_initial(Y: object, R_cal: object, *, center_support: object = (True, True, True),
                shrinkage: object = .25, variance_floor: object = .0001) -> dict:
    """Fit the finite-calibration covariance and projected initial center."""
    try:
        shape = np.asarray(Y).shape
    except (TypeError, ValueError) as exc:
        raise ValueError('Y must be n by 3') from exc
    if len(shape) != 2 or shape[0] < 2 or shape[1] != 3:
        raise ValueError('Y must be n by 3 with n >= 2')
    y = _array(Y, shape, 'Y')
    R = _symmetric(R_cal, 3, 'R_cal', positive=True)
    allowed, E, rhs = _constraints(center_support)
    alpha = _scalar(shrinkage, 'shrinkage', 0, 1)
    floor = _scalar(variance_floor, 'variance_floor', 0, strict=True)
    n = len(y)
    raw_mean = y.mean(axis=0)
    S = np.cov(y @ U, rowvar=False, ddof=1)
    Rz = U.T @ R @ U
    eigen, V = np.linalg.eigh((S - Rz + (S - Rz).T) / 2)
    floored = np.maximum(eigen, floor)
    G = (V * floored) @ V.T
    T = (1 - alpha) * G + alpha * np.trace(G) * np.eye(2) / 2
    M = (T + Rz) / n
    fitted = _solve(np.eye(3), raw_mean, allowed, E, rhs)
    diagnostics = dict(n_cal=n, raw_mean=raw_mean.tolist(), raw_tangent_cov=S.tolist(),
        deconvolved_eigenvalues=eigen.tolist(), floored_eigenvalues=floored.tolist(),
        projection_distance=None if fitted is None else float(np.linalg.norm(fitted[0] - raw_mean)),
        included_unit_ids=[], excluded_unit_id=None, n_update=0,
        mean_information=None, mean_score=None, mean_constraint_violation=None,
        mean_kkt_residual=None)
    if fitted is None or not np.isfinite(T).all() or not np.isfinite(M).all():
        return _failure_prior('initial center or covariance solve failed', diagnostics)
    center, violation, kkt = fitted
    P = np.linalg.inv(M)
    diagnostics.update(mean_information=P.tolist(),
        mean_score=(P @ (U.T @ (center - B))).tolist(),
        mean_constraint_violation=violation, mean_kkt_residual=kkt)
    return _prior(center, T, M, diagnostics)


def _observations(y: object, D: object, cov: object) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    try:
        shape = np.asarray(D).shape
    except (TypeError, ValueError) as exc:
        raise ValueError('D must be p by 3') from exc
    if len(shape) != 2 or shape[0] < 1 or shape[1] != 3:
        raise ValueError('D must be p by 3 with p >= 1')
    design = _array(D, shape, 'D')
    return _array(y, (shape[0],), 'y'), design, _symmetric(cov, shape[0], 'cov', positive=True)


def update_mean(prior: object, units: object, *, support: object = (True, True, True),
                exclude_unit: str | None = None) -> dict:
    """Update a working domain mean, excluding a complete specimen stack first."""
    allowed, E, rhs = _constraints(support)
    center, T, M0 = _validate_prior(prior, allowed, E, rhs)
    if not isinstance(units, (list, tuple)):
        raise ValueError('units must be a list of unit mappings')
    if exclude_unit is not None and (not isinstance(exclude_unit, str) or not exclude_unit):
        raise ValueError('exclude_unit must be a nonempty unit ID')
    seen: set[str] = set()
    selected: list[Mapping] = []
    for unit in units:
        if not isinstance(unit, Mapping) or not isinstance(unit.get('unit_id'), str) or not unit['unit_id']:
            raise ValueError('each unit requires a nonempty unit_id')
        uid = unit['unit_id']
        if uid in seen:
            raise ValueError('duplicate unit_id')
        seen.add(uid)
        if uid != exclude_unit:
            selected.append(unit)
    if exclude_unit is not None and exclude_unit not in seen:
        raise ValueError('exclude_unit is unknown')
    diagnostics = dict(prior['diagnostics'])
    diagnostics.update(included_unit_ids=[u['unit_id'] for u in selected],
        excluded_unit_id=exclude_unit, n_update=len(selected))
    if not selected:
        return dict(status='ok', reason='', center=list(prior['center']),
            tangent_cov=np.array(T).tolist(), mean_cov=np.array(M0).tolist(),
            precision=np.array(prior['precision']).tolist(), diagnostics=diagnostics)
    try:
        P = np.linalg.inv(M0)
        v = P @ (U.T @ (center - B))
        for unit in selected:
            y, D, R = _observations(unit.get('y'), unit.get('D'), unit.get('cov'))
            A = D @ U
            V = A @ T @ A.T + R
            np.linalg.cholesky(V)
            P = P + A.T @ np.linalg.solve(V, A)
            v = v + A.T @ np.linalg.solve(V, y - D @ B)
        M = np.linalg.inv(P)
        H = U @ P @ U.T
        f = U @ v + H @ B
        fitted = _solve(H, f, allowed, E, rhs)
    except np.linalg.LinAlgError:
        return _failure_prior('mean information solve failed', diagnostics)
    diagnostics.update(mean_information=P.tolist(), mean_score=v.tolist())
    if fitted is None or not np.isfinite(M).all():
        return _failure_prior('updated center solve failed', diagnostics)
    updated, violation, kkt = fitted
    diagnostics.update(mean_constraint_violation=violation, mean_kkt_residual=kkt,
        update_projection_distance=float(np.linalg.norm(updated - (B + U @ (M @ v)))))
    return _prior(updated, T, M, diagnostics)


def _failure_estimate(reason: str) -> dict:
    return dict(status='numerical_failure', reason=reason, mass=None, reconstructed=None,
        residual=None, assay_loss=None, prior_loss=None, objective=None,
        constraint_violation=None, kkt_residual=None, active_support=None)


def estimate(y: object, D: object, cov: object, prior: object = None, *, method: str,
             strength: object = 1., support: object = (True, True, True),
             fixed_quartz: object = None) -> dict:
    """Fit A0/A1/A2 on the declared simplex with full known assay covariance."""
    if not isinstance(method, str) or method not in ('A0', 'A1', 'A2'):
        raise ValueError('method must be A0, A1, or A2')
    allowed, E, rhs = _constraints(support, fixed_quartz)
    lam = _scalar(strength, 'strength', 0)
    if method == 'A0' and prior is not None:
        raise ValueError('A0 requires prior=None')
    if method == 'A1':
        if y is not None or D is not None or cov is not None:
            raise ValueError('A1 requires y=D=cov=None')
    else:
        obs, design, R = _observations(y, D, cov)
    if method != 'A0':
        center, T, M = _validate_prior(prior, allowed, E, rhs)
        K = _symmetric(prior['precision'], 3, 'prior.precision', semidefinite=True)
    if method == 'A1':
        m = center.copy()
        return dict(status='ok', reason='', mass=m.tolist(), reconstructed=None,
            residual=None, assay_loss=None, prior_loss=0., objective=0.,
            constraint_violation=_violation(m, allowed, E, rhs), kkt_residual=None,
            active_support=[int(j) for j in np.flatnonzero(m > ACTIVE_TOL)])
    Wd = np.linalg.solve(R, design)
    Wy = np.linalg.solve(R, obs)
    H = design.T @ Wd
    f = design.T @ Wy
    if method == 'A2' and lam > 0:
        H = H + lam * K
        f = f + lam * K @ center
    fitted = _solve(H, f, allowed, E, rhs)
    if fitted is None:
        return _failure_estimate('global convex solve or KKT certificate failed')
    m, violation, cert = fitted
    reconstructed = design @ m
    residual = obs - reconstructed
    assay_loss = float(.5 * residual @ np.linalg.solve(R, residual))
    prior_loss = float(.5 * lam * (m - center) @ K @ (m - center)) if method == 'A2' else 0.
    return dict(status='ok', reason='', mass=m.tolist(), reconstructed=reconstructed.tolist(),
        residual=residual.tolist(), assay_loss=assay_loss, prior_loss=prior_loss,
        objective=assay_loss + prior_loss, constraint_violation=violation,
        kkt_residual=cert, active_support=[int(j) for j in np.flatnonzero(m > ACTIVE_TOL)])


def _spectrum(X: np.ndarray) -> tuple[int, list[float]]:
    singular = np.linalg.svd(X, compute_uv=False)
    threshold = RANK_TOL * max(1., float(singular[0]) if singular.size else 0.)
    return int(np.sum(singular > threshold)), singular.tolist()


def information(D: object, cov: object, *, support: object = (True, True, True)) -> dict:
    """Report chemical, closure-augmented, and whitened tangent ranks."""
    try:
        shape = np.asarray(D).shape
    except (TypeError, ValueError) as exc:
        raise ValueError('D must be p by 3') from exc
    if len(shape) != 2 or shape[0] < 1 or shape[1] != 3:
        raise ValueError('D must be p by 3 with p >= 1')
    _, design, R = _observations(np.zeros(shape[0]), D, cov)
    allowed = _support(support)
    support_tangent = null_space(np.vstack((np.ones(3), np.eye(3)[~allowed])), rcond=SYMMETRY_TOL)
    W = np.linalg.solve(np.linalg.cholesky(R), design)
    matrices = (design, np.vstack((design, np.ones(3))), W @ U, W @ support_tangent)
    spectra = [_spectrum(x) for x in matrices]
    return dict(rank_raw=spectra[0][0], rank_augmented=spectra[1][0],
        rank_tangent=spectra[2][0], rank_support_tangent=spectra[3][0],
        singular_raw=spectra[0][1], singular_augmented=spectra[1][1],
        singular_tangent=spectra[2][1], singular_support_tangent=spectra[3][1])


def abundance_bounds(y: object, D: object, cov: object, *,
                     support: object = (True, True, True), fixed_quartz: object = None,
                     multiplier: object = 3.) -> dict:
    """Find marginal extrema in a measurement-consistent box, not intervals."""
    obs, design, R = _observations(y, D, cov)
    allowed, E, rhs = _constraints(support, fixed_quartz)
    width = _scalar(multiplier, 'multiplier', 0, strict=True) * np.sqrt(np.diag(R))
    A_ub = np.vstack((design, -design))
    b_ub = np.concatenate((obs + width, width - obs))
    kwargs = dict(A_ub=A_ub, b_ub=b_ub, A_eq=E, b_eq=rhs,
                  bounds=[(0., None) if x else (0., 0.) for x in allowed], method='highs')
    lower: list[float] = []
    upper: list[float] = []
    max_violation = 0.
    for j in range(3):
        for sign, target in ((1., lower), (-1., upper)):
            objective = sign * np.eye(3)[j]
            result = linprog(objective, **kwargs)
            if result.status == 2:
                return dict(status='infeasible_envelope', reason='measurement box is empty',
                    lower=None, upper=None, width=None, max_constraint_violation=None)
            if result.status != 0 or result.x is None:
                return dict(status='numerical_failure', reason=f'LP failed: {result.message}',
                    lower=None, upper=None, width=None, max_constraint_violation=None)
            m = result.x
            violation = max(_violation(m, allowed, E, rhs),
                            float(np.max(np.maximum(np.abs(design @ m - obs) - width, 0.))))
            if violation > FEAS_TOL:
                return dict(status='numerical_failure', reason='LP certificate failed',
                    lower=None, upper=None, width=None, max_constraint_violation=violation)
            max_violation = max(max_violation, violation)
            target.append(float(m[j]))
    return dict(status='ok', reason='', lower=lower, upper=upper,
                width=(np.array(upper) - lower).tolist(), max_constraint_violation=max_violation)
