"""Frozen least-squares comparators for EMC feasibility version 1."""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize, nnls

METHODS = ('pinv', 'nnls', 'wnnls', 'cls_eq', 'cls_le')

def solve(C, a, *, method, sigma=None):
    C = np.asarray(C, dtype=np.float64)
    a = np.asarray(a, dtype=np.float64)
    if C.ndim != 2 or C.shape[0] < 1 or C.shape[1] < 1 or C.shape[1] > 5 or a.shape != (C.shape[0],) or not np.isfinite(C).all() or not np.isfinite(a).all() or np.any(C < 0) or np.any(np.sum(C,axis=0) <= 0):
        raise ValueError('invalid C or a')
    if method not in METHODS:
        raise ValueError('invalid method')
    s = np.ones(C.shape[0]) if sigma is None else np.asarray(sigma,dtype=np.float64)
    if s.shape != (C.shape[0],) or not np.isfinite(s).all() or np.any(s <= 0):
        raise ValueError('invalid sigma')
    iterations = 0
    status = 'ok'
    if method == 'pinv':
        m = np.linalg.pinv(C) @ a
    elif method == 'nnls':
        m, _ = nnls(C, a)
    elif method == 'wnnls':
        m, _ = nnls(C / s[:,None], a / s)
    else:
        B, y = C / s[:,None], a / s
        def objective(x):
            z = B @ x - y
            return float(z @ z)
        def jac(x):
            return 2 * B.T @ (B @ x - y)
        k = C.shape[1]
        constraint = {'type':'eq' if method == 'cls_eq' else 'ineq',
                      'fun':(lambda x: np.sum(x)-1) if method == 'cls_eq' else (lambda x: 1-np.sum(x)),
                      'jac':(lambda x: np.ones(k)) if method == 'cls_eq' else (lambda x: -np.ones(k))}
        result = minimize(objective,np.full(k,1/k if method=='cls_eq' else 1/(2*k)),jac=jac,bounds=[(0,None)]*k,constraints=[constraint],method='SLSQP',options={'ftol':1e-12,'maxiter':1000})
        m = result.x
        iterations = int(result.nit)
        if not result.success or np.min(m)<-1e-8 or (abs(np.sum(m)-1)>1e-8 if method=='cls_eq' else np.sum(m)>1+1e-8):
            status = 'solver_failure'
    reconstructed = C @ m
    residual = reconstructed-a
    return dict(m=m,reconstructed=reconstructed,residual=residual,status=status,objective=float(np.sum((residual/s)**2)),iterations=iterations,unallocated=float(1-np.sum(m)))
