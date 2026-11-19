#!/usr/bin/env python3
"""Compare bounded Python STEPR calls with actual pinned R references."""
import os,sys,json
from pathlib import Path
import numpy as np,pandas as pd
base=Path('examples/copper'); d=pd.read_csv(base/'copper.csv'); parts=['Chalcopyrite','Chalcocite','Pyrite','Quartz','Kaolinite']; tr=d.split.eq('train'); te=~tr
sys.path.insert(0,str(Path(os.environ.get('LF_COPPER_SRC','src'))/'python'))
from stepr import stepr_py,predict_stepr
refs=json.loads((base/'evidence/r_stepr_reference.json').read_text()); out=[]
for response in ('linear','nonlinear'):
 for criterion in ('AIC','BIC'):
  for nsteps in (1,3):
   key=f'{response}_{criterion}_{nsteps}'; r=refs[key]
   row=dict(response=response,criterion=criterion,nsteps=nsteps,r_status=r['status'],python_status='success',selected_pairs_match=False,prediction_max_abs_diff=None,parity_status='UNVERIFIED',error='')
   try:
    if r['status']!='success':raise ValueError('R reference failed: '+r.get('error',''))
    f=stepr_py(d.loc[tr,parts],d.loc[tr,f'y_{response}'],method=1,family='gaussian',criterion=criterion,nsteps=nsteps,top=1)
    pairs=f['selected_pairs'] if isinstance(f,dict) else f.selected_pairs
    names=['/'.join(p) for p in pairs];row['r_names']=';'.join(r['names'] if isinstance(r['names'],list) else [r['names']]);row['python_names']=';'.join(names)
    row['selected_pairs_match']=row['r_names']==row['python_names']
    p=np.asarray(predict_stepr(f,d.loc[te,parts])); rr=np.asarray(r['test_prediction']);diff=np.max(np.abs(p-rr));row['prediction_max_abs_diff']=float(diff)
    row['parity_status']='PASS' if row['selected_pairs_match'] and diff<=1e-8 else 'FAIL'
   except Exception as e:row.update(python_status='failed',parity_status='FAIL',error=f'{type(e).__name__}: {e}')
   out.append(row)
pd.DataFrame(out).to_csv(base/'evidence/stepr_parity.csv',index=False)
print(pd.DataFrame(out).to_string(index=False))

# Three-step comparator rows remain visible as diagnostics when pair parity fails.
m=base/'evidence/python_metrics.csv'
if m.exists():
    z=pd.read_csv(m); z['selected_ratios']=z['selected_ratios'].fillna(''); z['error']=z['error'].fillna('')
    for response in ('linear','nonlinear'):
        ok=all(q['parity_status']=='PASS' for q in out if q['response']==response and q['criterion']=='AIC' and q['nsteps']==3)
        if not ok:
            mask=z.response.eq(response)&z.method.isin(['stepr_py_ols','stepr_py_rf'])
            z.loc[mask,'status']='diagnostic_unverified'
            z.loc[mask,'error']='Actual R STEPR selected third ratio differs; no method-equivalence claim'
    z.to_csv(m,index=False)
