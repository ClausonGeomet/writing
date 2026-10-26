#!/usr/bin/env python3
"""Evaluate the builder black box after the Part A freeze and worktree merge."""
from __future__ import annotations
import csv, hashlib, importlib.util, json, math, platform, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import scipy

ROOT=Path(__file__).resolve().parent
BROOT=ROOT.parent/'emc-part2b'
FIX=ROOT/'fixtures'; OUT=BROOT/'results'
CASES=('easy','ambiguous','noisy','missing_phase','wrong_chemistry','infeasible')
BASE=('Quartz','Hematite','Kaolinite','Forsterite'); AMB=('Quartz','Anatase','Rutile')
SIGMA=np.array([.010,.008,.003,.006])

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def write_csv(path,rows,fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def main():
    frozen=ROOT/'baseline_freeze.json'
    if not frozen.exists(): raise RuntimeError('freeze Part A first')
    proto=BROOT/'prototypes/emc_prototypes.py'
    if not proto.exists():
        print('PENDING: prototype absent; run after builder merge')
        return
    spec=importlib.util.spec_from_file_location('emc_prototypes',proto)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    OUT.mkdir(exist_ok=True)
    cbase=pd.read_csv(FIX/'C_base.csv'); camb=pd.read_csv(FIX/'C_ambiguous.csv')
    records=[]; predictions=[]; traces=[]; times={}
    start=time.perf_counter(); cpu=time.process_time()
    configs=[(m,o,seed) for o in ('row','distribution') for m in ('fixed','exhaustive','greedy','beam','random','ga') for seed in ((101,102,103) if m in ('random','ga') else (0,))]
    for si,case in enumerate(CASES):
      for r in range(3):
        df=pd.read_csv(FIX/f'{case}_r{r}.csv').sort_values('sample_id').reset_index(drop=True)
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
        if case=='noisy':
            rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([20261026,si,r,0])))
            A+=rng.normal(size=A.shape)*SIGMA
        if case=='infeasible': A*=1.2
        sigma=SIGMA if case=='noisy' else np.ones(len(channels))
        for method,obj,seed in configs:
          if time.perf_counter()-start>600:
            status='time_limit_partial'; break
          tick=time.perf_counter(); ctick=time.process_time()
          fit=mod.fit_order(C,A[:30],method=method,objective=obj,sigma=sigma,seed=seed,budget=24,beam_width=2)
          fit_wall=time.perf_counter()-tick; fit_cpu=time.process_time()-ctick
          key=(case,r,method,obj,seed)
          times[str(key)]=dict(fit_wall_seconds=fit_wall,fit_cpu_seconds=fit_cpu)
          order=fit.get('order')
          for ti,tr in enumerate(fit.get('trace',[])):
            traces.append(dict(scenario=case,replicate=r,method=method,objective=obj,search_seed=seed,trace_index=ti,proposal=tr.get('proposal'),evaluations=tr.get('evaluations'),best_score=tr.get('best_score'),order=json.dumps(tr.get('order'))))
          for split,lo,hi in (('train',0,30),('validation',30,40),('test',40,60)):
            At=A[lo:hi]; Mt=M[lo:hi]
            if order is None:
              mh=np.full((hi-lo,len(fitparts)),np.nan); ah=np.full_like(At,np.nan); statuses=['unsupported_negative_assay']*(hi-lo); ua=np.full(hi-lo,np.nan)
            else:
              pred=mod.predict_order(C,At,order)
              mh=np.asarray(pred['M'],dtype=float); ah=np.asarray(pred['A_hat'],dtype=float); statuses=pred['status']; ua=np.asarray(pred['unallocated'],dtype=float)
            full=np.zeros((hi-lo,len(trueparts)))
            full[:,[trueparts.index(p) for p in fitparts]]=mh
            good=np.isfinite(ah).all(axis=1)&np.isfinite(full).all(axis=1)
            score=mod.score_assays(At[good],ah[good],sigma=sigma,objective=obj,scales=fit.get('scales')) if good.sum()>=2 else {}
            assay_err=ah[good]-At[good]; mineral_err=full[good]-Mt[good]
            rec=dict(scenario=case,replicate=r,split=split,method=method,objective=obj,search_seed=seed,order=json.dumps(order),status=fit.get('status'),n_total=hi-lo,n_supported=int(good.sum()),unsupported_count=int((~good).sum()),failure_count=sum(x not in ('ok','unsupported_negative_assay') for x in statuses),assay_rmse=float(np.sqrt(np.mean(assay_err**2))) if good.any() else math.nan,standardized_assay_rmse=float(np.sqrt(np.mean((assay_err/sigma)**2))) if good.any() else math.nan,mineral_rmse=float(np.sqrt(np.mean(mineral_err**2))) if good.any() else math.nan,closure_error=float(np.mean(np.abs(full[good].sum(axis=1)-1))) if good.any() else math.nan,unallocated=float(np.mean(ua[good])) if good.any() else math.nan,negative_mass=float(np.sum(np.maximum(-full[good],0))) if good.any() else math.nan,exact_residual=float(np.max(np.abs(assay_err))) if good.any() else math.nan,fit_wall_seconds=fit_wall,fit_cpu_seconds=fit_cpu,evaluations=fit.get('evaluations'),proposals=fit.get('proposals'),train_score=fit.get('train_score'),row_loss=score.get('row'),mean_loss=score.get('mean'),covariance_loss=score.get('covariance'),js_loss=score.get('js'),total_loss=score.get('total'),scales=json.dumps(np.asarray(fit.get('scales')).tolist()) if fit.get('scales') is not None else '')
            records.append(rec)
            for i in range(hi-lo):
              for j,p in enumerate(trueparts):
                predictions.append(dict(scenario=case,replicate=r,sample_id=df.iloc[lo+i]['sample_id'],split=split,method=method,objective=obj,search_seed=seed,phase=p,true_mass=Mt[i,j],estimated_mass=full[i,j],error=full[i,j]-Mt[i,j],status=statuses[i]))
        else: continue
        break
      else: continue
      break
    oracle={(x['scenario'],x['replicate'],x['objective']):x['train_score'] for x in records if x['method']=='exhaustive' and x['split']=='train'}
    for x in records:
      ref=oracle.get((x['scenario'],x['replicate'],x['objective']))
      x['objective_gap']=x['train_score']-ref if ref is not None and x['train_score'] is not None else math.nan
    write_csv(OUT/'metrics.csv',records,list(records[0]) if records else ['scenario','status'])
    write_csv(OUT/'predictions.csv',predictions,list(predictions[0]) if predictions else ['scenario','status'])
    write_csv(OUT/'search_trace.csv',traces,list(traces[0]) if traces else ['scenario','status'])
    source=ROOT/'results/observational_equivalence.csv'
    (OUT/'observational_equivalence.csv').write_bytes(source.read_bytes())
    files=[proto,frozen,source,*OUT.glob('*.csv')]
    manifest=dict(schema='emc-feasibility-1',created_utc=datetime.now(timezone.utc).isoformat(),baseline_freeze_sha256=sha(frozen),status='complete' if len(records)==18*20*3 else 'time_limit_partial',configurations_per_cohort=20,fit_calls=len(records)//3,split=[30,10,20],budget=24,beam_width=2,search_seeds=[101,102,103],runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__),time=dict(wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu),fit_times=times,sha256={str(p.relative_to(ROOT.parent)):sha(p) for p in files if p.is_relative_to(ROOT.parent)})
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=True)+'\n')
    print(manifest['status'],manifest['fit_calls'])
if __name__=='__main__': main()
