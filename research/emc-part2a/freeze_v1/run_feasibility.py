#!/usr/bin/env python3
"""Run fixed, paired EMC baseline feasibility cohorts; freeze before B."""
from __future__ import annotations
import csv, hashlib, json, math, platform, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linprog
import scipy
from baselines import METHODS, solve

ROOT=Path(__file__).resolve().parent
FIX=ROOT/'fixtures'; OUT=ROOT/'results'
CASES=('easy','ambiguous','noisy','missing_phase','wrong_chemistry','infeasible')
BASE=('Quartz','Hematite','Kaolinite','Forsterite')
AMB=('Quartz','Anatase','Rutile')
SIGMA=np.array([.010,.008,.003,.006])

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def write_csv(path,rows,columns):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore'); w.writeheader(); w.writerows(rows)
def svd_record(C,case,r):
    s=np.linalg.svd(C,compute_uv=False)
    aug=np.vstack([C,np.ones(C.shape[1])]); t=np.linalg.svd(aug,compute_uv=False)
    def rank(x,shape,relative): return int(np.sum(x>max(shape)*np.finfo(float).eps*x[0] if not relative else x>1e-8*x[0]))
    return dict(scenario=case,replicate=r,singular_values=json.dumps(s.tolist()),augmented_singular_values=json.dumps(t.tolist()),rank=rank(s,C.shape,False),rank_rel1e8=rank(s,C.shape,True),augmented_rank=rank(t,aug.shape,False),augmented_rank_rel1e8=rank(t,aug.shape,True))

def main():
    OUT.mkdir(exist_ok=True)
    cbase=pd.read_csv(FIX/'C_base.csv'); camb=pd.read_csv(FIX/'C_ambiguous.csv')
    bases=[]; preds=[]; ident=[]; splits=[]; counter=[]; case_meta={}
    start=time.perf_counter(); cpu=time.process_time()
    for si,case in enumerate(CASES):
      for r in range(3):
        df=pd.read_csv(FIX/f'{case}_r{r}.csv').sort_values('sample_id').reset_index(drop=True)
        original=pd.read_csv(FIX/f'{case}_r{r}.csv')['sample_id'].tolist()
        trueparts=AMB if case=='ambiguous' else BASE
        fitparts=AMB if case=='ambiguous' else BASE[:3] if case=='missing_phase' else BASE
        cdf=camb if case=='ambiguous' else cbase
        channels=cdf['analyte'].tolist()
        Ct=cdf[list(trueparts)].to_numpy(dtype=float)
        C=Ct[:,[trueparts.index(p) for p in fitparts]].copy()
        if case=='wrong_chemistry': C[channels.index('Al2O3'),fitparts.index('Kaolinite')]*=.8
        M=df[[f'm_{p}' for p in trueparts]].to_numpy(dtype=float)
        A0=df[[f'a0_{c}' for c in channels]].to_numpy(dtype=float)
        A=A0.copy()
        seedinfo=None
        if case=='noisy':
            ss=np.random.SeedSequence([20261026,si,r,0]); rng=np.random.Generator(np.random.PCG64(ss)); A += rng.normal(size=A.shape)*SIGMA
            seedinfo=dict(entropy=[20261026,si,r,0],bit_generator='PCG64')
        if case=='infeasible': A=1.2*A
        sigma=SIGMA if case=='noisy' else np.ones(len(channels))
        assert np.max(abs(A0-M@Ct.T))<=1e-10
        ident.append(svd_record(C,case,r))
        case_meta[f'{case}_r{r}']=dict(C_true=Ct.tolist(),C_fit=C.tolist(),channels=channels,true_phases=trueparts,fit_phases=fitparts,sigma=sigma.tolist(),noise_seed=seedinfo,source_file=f'fixtures/{case}_r{r}.csv')
        for i,row in df.iterrows():
          split='train' if i<30 else 'validation' if i<40 else 'test'
          sid=row['sample_id']; splits.append(dict(scenario=case,replicate=r,sample_id=sid,original_index=original.index(sid),sorted_index=i,split=split))
          for method in METHODS:
            tick=time.perf_counter(); ctick=time.process_time()
            result=solve(C,A[i],method=method,sigma=sigma)
            m=np.zeros(len(trueparts)); m[[trueparts.index(p) for p in fitparts]]=result['m']
            ah=result['reconstructed']
            bases.append(dict(scenario=case,replicate=r,sample_id=sid,split=split,method=method,objective='weighted_least_squares' if method in ('wnnls','cls_eq','cls_le') else 'ordinary_least_squares',search_seed='',status=result['status'],iterations=result['iterations'],solver_objective=result['objective'],assay_rmse=float(np.sqrt(np.mean((ah-A[i])**2))),standardized_assay_rmse=float(np.sqrt(np.mean(((ah-A[i])/sigma)**2))),mineral_rmse=float(np.sqrt(np.mean((m-M[i])**2))),closure_error=float(abs(np.sum(m)-1)),unallocated=result['unallocated'],negative_mass=float(np.sum(np.maximum(-m,0))),exact_residual=float(np.max(abs(ah-A[i]))),wall_seconds=time.perf_counter()-tick,cpu_seconds=time.process_time()-ctick,unsupported_count=0,failure_count=int(result['status']!='ok')))
            for j,p in enumerate(trueparts):
              preds.append(dict(scenario=case,replicate=r,sample_id=sid,split=split,method=method,objective='',search_seed='',phase=p,true_mass=M[i,j],estimated_mass=m[j],error=m[j]-M[i,j]))
          if case in ('easy','ambiguous','missing_phase','wrong_chemistry'):
            for closure in ('eq','le'):
              lp=[]; feasible=True
              for j in range(C.shape[1]):
                v=np.eye(C.shape[1])[j]
                kwargs=dict(A_eq=C,b_eq=A[i],bounds=[(0,None)]*C.shape[1],method='highs')
                if closure=='eq': kwargs['A_eq']=np.vstack([C,np.ones(C.shape[1])]); kwargs['b_eq']=np.r_[A[i],1.]
                else: kwargs['A_ub']=np.ones((1,C.shape[1])); kwargs['b_ub']=[1.]
                lo=linprog(v,**kwargs); hi=linprog(-v,**kwargs)
                if not lo.success or not hi.success: feasible=False; break
                lp.append((float(lo.fun),float(-hi.fun)))
              for j,p in enumerate(fitparts):
                ident.append(dict(scenario=case,replicate=r,sample_id=sid,split=split,closure=closure,phase=p,feasible=feasible,lower=lp[j][0] if feasible else '',upper=lp[j][1] if feasible else ''))
        if case=='ambiguous':
          swapped=M[:,[0,2,1]]
          assert np.max(abs(M@Ct.T-swapped@Ct.T))<=1e-12
          for i,sid in enumerate(df.sample_id):
            counter.append(dict(scenario=case,replicate=r,sample_id=sid,split=splits[-60+i]['split'],source='simCoDaR paired swap',assay_max_diff=float(np.max(abs(M[i]@Ct.T-swapped[i]@Ct.T))),anatase_original=M[i,1],anatase_swapped=swapped[i,1],rutile_original=M[i,2],rutile_swapped=swapped[i,2],phase_specific_diff=float(abs(M[i,1]-swapped[i,1]))))
    write_csv(OUT/'baselines.csv',bases,list(bases[0]))
    write_csv(OUT/'predictions.csv',preds,list(preds[0]))
    write_csv(OUT/'identifiability.csv',ident,sorted(set().union(*(x.keys() for x in ident))))
    write_csv(FIX/'split_manifest.csv',splits,list(splits[0]))
    write_csv(OUT/'observational_equivalence.csv',counter,list(counter[0]))
    # Independent analytic illustration, distinct from simCoDaR cohorts.
    analytic=[]
    for r in range(3):
      rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([20261026,1,r,1])))
      Q=rng.uniform(.4,.8,60); U=rng.beta(2,8,60)
      for i in range(60):
        T=1-Q[i]; x=np.array([Q[i],T*U[i],T*(1-U[i])]); y=x[[0,2,1]]
        analytic.append(dict(replicate=r,sample_id=f'analytic_{i+1}',Q=Q[i],U=U[i],original_anatase=x[1],swapped_anatase=y[1],total_tio2_diff=float((x[1]+x[2])-(y[1]+y[2])),phase_specific_diff=float(x[1]-y[1])))
    write_csv(OUT/'analytic_counterexample.csv',analytic,list(analytic[0]))
    ds=[('Fortescue CSIRO collection6753','whole-rock XRD/chemistry; prior 23 modeled, 14/5/4 split','proxy chemistry; independent holdout and basis audit absent'),('Mary-Wildman CSIRO','nine prior bulk XRD/XRF pairs','dry basis and censoring audit needed'),('Cannington AuScope Cloncurry METAL','TIMA vol% with surface pXRF','volume to mass and specimen basis mismatch'),('Mantovani 10.3389/fenvs.2023.1179272','ash/grain-bin groups','amorphous material and group join'),('Qiong local','historical local data','provenance incomplete'),('2026 lithium 10.3390/min16020139','article','raw data confidential'),('Koashva 10.1038/s41598-026-52558-9','article accessible','paired raw matrix unverified'),('simCoDaR IODP','bundled assays','no mineral truth')]
    write_csv(OUT/'dataset_screen.csv',[dict(dataset=a,scope=b,blocker=c,admitted=False) for a,b,c in ds],['dataset','scope','blocker','admitted'])
    manifest=dict(schema='emc-feasibility-1',created_utc=datetime.now(timezone.utc).isoformat(),dimensions=dict(scenarios=6,replicates=3,rows_per_cohort=60,analytes='varies by case',units='whole-sample mass fraction'),cases=case_meta,settings=dict(methods=METHODS,split=[30,10,20],base_seed=20261026,RNG='R Mersenne-Twister/Inversion/Rejection; Python PCG64 SeedSequence'),runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__),time=dict(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu))
    files=[ROOT/'generation.R',ROOT/'baselines.py',ROOT/'run_feasibility.py',*FIX.glob('*.csv'),*FIX.glob('*.json'),*OUT.glob('*.csv')]
    manifest['sha256']={str(p.relative_to(ROOT)):sha(p) for p in sorted(files)}
    (FIX/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    frozen=dict(schema='emc-feasibility-1',frozen_utc=datetime.now(timezone.utc).isoformat(),dependency_commit='d9a0b27676a0cb3a4acb300a3f6b2ebc845d205d',settings=manifest['settings'],runtime=manifest['runtime'],hashes=manifest['sha256'])
    (ROOT/'baseline_freeze.json').write_text(json.dumps(frozen,indent=2)+'\n')
    print(json.dumps(dict(baseline_rows=len(bases),prediction_rows=len(preds),wall_seconds=manifest['time']['wall_seconds'],freeze_sha256=sha(ROOT/'baseline_freeze.json'))))
if __name__=='__main__': main()
