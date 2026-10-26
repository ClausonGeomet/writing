"""Descriptive summaries of the frozen rowwise Part A output."""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
b=pd.read_csv(ROOT/'results/baselines.csv'); p=pd.read_csv(ROOT/'results/predictions.csv')
rows=[]; phases=[]
for (case,r,split,method),g in b.groupby(['scenario','replicate','split','method'],sort=True):
    q=p[(p.scenario==case)&(p.replicate==r)&(p.split==split)&(p.method==method)]
    piv=q.pivot(index='sample_id',columns='phase',values=['true_mass','estimated_mass'])
    truth=piv['true_mass'].to_numpy(); est=piv['estimated_mass'].to_numpy()
    eligible=(np.min(truth,axis=1)>0)&(np.min(est,axis=1)>0)&(abs(truth.sum(axis=1)-1)<=1e-8)&(abs(est.sum(axis=1)-1)<=1e-8)
    if eligible.any():
        # Orthogonal projection onto zero-sum log space gives the ILR norm.
        d=np.log(est[eligible])-np.log(truth[eligible]); d-=d.mean(axis=1,keepdims=True)
        ilr=float(np.sqrt(np.mean(np.sum(d*d,axis=1))))
    else: ilr=np.nan
    rows.append(dict(scenario=case,replicate=r,split=split,method=method,n=len(g),failure_count=int(g.failure_count.sum()),unsupported_count=int(g.unsupported_count.sum()),assay_rmse=float(np.sqrt(np.mean(g.assay_rmse**2))),standardized_assay_rmse=float(np.sqrt(np.mean(g.standardized_assay_rmse**2))),mineral_rmse=float(np.sqrt(np.mean(g.mineral_rmse**2))),mean_closure_error=float(g.closure_error.mean()),mean_unallocated=float(g.unallocated.mean()),negative_mass=float(g.negative_mass.sum()),max_exact_residual=float(g.exact_residual.max()),wall_seconds=float(g.wall_seconds.sum()),cpu_seconds=float(g.cpu_seconds.sum()),ilr_rms=ilr,ilr_eligible=int(eligible.sum())))
    for phase,h in q.groupby('phase'):
        e=h.error.to_numpy()
        phases.append(dict(scenario=case,replicate=r,split=split,method=method,phase=phase,n=len(e),bias=float(e.mean()),mae=float(abs(e).mean()),rmse=float(np.sqrt(np.mean(e*e)))))
pd.DataFrame(rows).to_csv(ROOT/'results/summary.csv',index=False,float_format='%.17g')
pd.DataFrame(phases).to_csv(ROOT/'results/phase_summary.csv',index=False,float_format='%.17g')
