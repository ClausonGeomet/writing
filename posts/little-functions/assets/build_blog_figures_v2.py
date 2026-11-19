"""Three editorial figures from pinned frozen rows; never imports or runs model code."""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import polars as pl
import plotnine as pn
import pandas as pd
import matplotlib

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'assets'
F = ROOT / 'materials/examples/copper_v2/copper_v2.csv'
E = ROOT / 'materials/examples/publication_sc020/artifacts/SC020-v1-reproduce-002/root/evidence/SC020-v1'
L, D = E / 'levels.csv', E / 'lf03_effects.csv'
SOURCES = {F: 'a6454551e639edcf472a5a1794ddd0c367cc278d3b09a4601c5c39cccf27b02c', L: 'cb468f0e63bed471d4c185874f8d2ff6c2291606895ca8b67f493bbf06a28a1d', D: '0cf7ea4b1a8eff08c4ef04b50b8080d6079cc8b5f6cea1a78cf3c19b08ff1c0a'}
PARTS = ['Chalcopyrite', 'Chalcocite', 'Pyrite', 'Quartz', 'Kaolinite']
METHODS = ['stepr_ols', 'stepr_rf', 'select_lr_ml_rf', 'ilr_rf', 'raw_rf']
LABELS = {'stepr_ols':'STEPR ratios + linear regression', 'stepr_rf':'STEPR ratios + forest', 'select_lr_ml_rf':'Forest-selected ratios + forest', 'ilr_rf':'Four ILR coordinates + forest', 'raw_rf':'Four raw mineral fractions + forest'}
COLORS = {'Chalcopyrite':'#0072B2','Chalcocite':'#56B4E9','Pyrite':'#E69F00','Quartz':'#999999','Kaolinite':'#009E73'}

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def check_source(path):
    if sha(path) != SOURCES[path]: raise ValueError(f'Pinned source changed: {path}')
def unique(df, keys):
    if df.group_by(keys).len().select(pl.col('len').max()).item() != 1: raise ValueError(f'Duplicate key: {keys}')
def save(df, stem, plot, width, height):
    inp = OUT / f'{stem}-input.csv'; png = OUT / f'{stem}.png'
    df.write_csv(inp)
    plot.save(png, width=width, height=height, dpi=160, verbose=False)
    return {'input':str(inp.relative_to(ROOT)), 'input_sha256':sha(inp), 'image':str(png.relative_to(ROOT)), 'image_sha256':sha(png), 'rows':df.height}

for p in SOURCES: check_source(p)
base = pl.read_csv(F).filter(pl.col('sample_id') == 'W01_001')
if base.height != 1 or base['split'][0] != 'train': raise ValueError('Expected one training W01_001')
values = [float(base[p][0]) for p in PARTS]
if not math.isclose(sum(values), 1, abs_tol=1e-10) or values[3] <= .05: raise ValueError('Invalid starting composition')
feed = pl.DataFrame([{'sample_id':'W01_001','stage':stage,'part':part,'fraction':value + (delta if stage == 'After transfer' else 0), 'mass_percent':100*(value + (delta if stage == 'After transfer' else 0))} for stage in ['Before transfer','After transfer'] for part,value,delta in zip(PARTS,values,[0,0,0,-.05,.05])])
unique(feed,['sample_id','stage','part'])
feed_plot = (pn.ggplot(feed.to_pandas(),pn.aes(x='stage',y='mass_percent',fill='part')) + pn.geom_col(width=.62,position=pn.position_stack(reverse=True),color='white',size=.4) + pn.scale_fill_manual(values=COLORS,breaks=PARTS) + pn.scale_x_discrete(limits=['Before transfer','After transfer'],labels=['Before','After']) + pn.scale_y_continuous(limits=(0,100),breaks=[0,25,50,75,100]) + pn.labs(x='',y='Mineral mass (% of feed)',fill='',title='Five minerals, one feed') + pn.guides(fill=pn.guide_legend(ncol=1)) + pn.theme_minimal(base_size=12) + pn.theme(legend_position='bottom',plot_title=pn.element_text(weight='bold',size=12),legend_text=pn.element_text(size=11),axis_text=pn.element_text(size=11),axis_title_y=pn.element_text(size=11)))
artifacts = {'feed':save(feed,'blog-v2-feed-transfer',feed_plot,3.6,5.3)}

levels = pl.read_csv(L).filter((pl.col('study')=='LF03') & (pl.col('response')=='observed_recovery_pct') & pl.col('language').is_in(['R','Python']) & pl.col('method').is_in(METHODS))
unique(levels,['study','language','response','method','model_id','sample_id'])
counts = levels.group_by(['language','method']).agg(pl.len().alias('n'),pl.col('model_id').n_unique().alias('models'),pl.col('sample_id').n_unique().alias('ids'))
if counts.height != 10 or counts.filter((pl.col('n')!=256)|(pl.col('models')!=1)|(pl.col('ids')!=256)).height: raise ValueError('Unexpected level cell')
rmse = (levels.group_by(['language','method','model_id']).agg(pl.len().alias('n_heldout'),((pl.col('pred')-pl.col('observed')).pow(2).mean().sqrt()).alias('rmse_pp')).with_columns(pl.col('method').replace(LABELS).alias('method_label')).sort(['language','method']))
order = [LABELS[m] for m in METHODS]
rpdf = rmse.to_pandas()
short = {'STEPR ratios + linear regression':'STEPR + linear', 'STEPR ratios + forest':'STEPR + forest', 'Forest-selected ratios + forest':'Selected + forest', 'Four ILR coordinates + forest':'ILR + forest', 'Four raw mineral fractions + forest':'Raw parts + forest'}
rpdf['plot_label'] = pd.Categorical(rpdf['method_label'].map(short),categories=[short[x] for x in order[::-1]],ordered=True)
rmse_plot = (pn.ggplot(rpdf,pn.aes(x='rmse_pp',y='plot_label',color='language',shape='language')) + pn.geom_point(size=3.6) + pn.scale_color_manual(values={'R':'#0072B2','Python':'#D55E00'}) + pn.scale_x_continuous(limits=(0,5),breaks=[0,1,2,3,4,5]) + pn.labs(x='Recovery RMSE (pp)',y='',color='',shape='',title='Heldout error') + pn.theme_minimal(base_size=12) + pn.theme(legend_position='bottom',plot_title=pn.element_text(weight='bold',size=12),legend_text=pn.element_text(size=11),axis_text=pn.element_text(size=11),axis_title_x=pn.element_text(size=11)))
artifacts['level'] = save(rmse,'blog-v2-level-rmse',rmse_plot,3.6,4.3)

raw = pl.read_csv(D, schema_overrides={'dose_token':pl.String}).filter((pl.col('dose_token')=='0.05') & pl.col('language').is_in(['R','Python']) & pl.col('method').is_in(['select_lr_ml_rf','stepr_ols']))
unique(raw,['language','method','model_id','sample_id','dose_token'])
if raw.height != 1024: raise ValueError('Expected 4x256 dose rows')
status = raw.group_by(['language','method']).agg((pl.col('status')=='SUCCESS').sum().alias('success'),(pl.col('status')=='OFF_SUPPORT').sum().alias('excluded'))
if status.height != 4 or status.filter((pl.col('success')!=235)|(pl.col('excluded')!=21)).height: raise ValueError('Unexpected support counts')
good = raw.filter(pl.col('status')=='SUCCESS')
if good.filter(pl.col('Delta_oracle')>=-1e-8).height: raise ValueError('Oracle sign changed')
classified = good.with_columns(pl.when(pl.col('Delta_pred').abs()<=1e-8).then(pl.lit('Near zero')).when(pl.col('Delta_pred')<0).then(pl.lit('Decrease')).otherwise(pl.lit('Increase')).alias('direction'))
direction = classified.group_by(['language','method','model_id','dose_token','direction']).agg(pl.len().alias('count')).with_columns(pl.lit(235).alias('n_supported'),pl.lit(21).alias('n_excluded'),pl.lit('pp').alias('effect_unit')).sort(['language','method','direction'])
expected={('R','select_lr_ml_rf'):(191,0,44),('Python','select_lr_ml_rf'):(194,4,37),('R','stepr_ols'):(235,0,0),('Python','stepr_ols'):(235,0,0)}
for (lang,method),(neg,zero,pos) in expected.items():
    ss=direction.filter((pl.col('language')==lang)&(pl.col('method')==method))
    actual={r['direction']:r['count'] for r in ss.to_dicts()}
    if (actual.get('Decrease',0),actual.get('Near zero',0),actual.get('Increase',0)) != (neg,zero,pos): raise ValueError('Direction count changed')
# Add zero-count categories only for plotting; CSV retains observed categories and identities.
dpdf=direction.to_pandas(); dpdf['fit']=dpdf['language']+' '+dpdf['method'].map({'stepr_ols':'STEPR linear','select_lr_ml_rf':'selected forest'})
dpdf['fit']=pd.Categorical(dpdf['fit'],categories=['R STEPR linear','R selected forest','Python STEPR linear','Python selected forest'],ordered=True)
dpdf['direction']=pd.Categorical(dpdf['direction'],categories=['Decrease','Near zero','Increase'],ordered=True)
direction_plot=(pn.ggplot(dpdf,pn.aes(x='fit',y='count',fill='direction')) + pn.geom_col(width=.66,position=pn.position_stack(reverse=True)) + pn.scale_fill_manual(values={'Decrease':'#0072B2','Near zero':'#E69F00','Increase':'#D55E00'},breaks=['Decrease','Near zero','Increase']) + pn.scale_y_continuous(limits=(0,250),breaks=[0,50,100,150,200]) + pn.labs(x='',y='Supported rows (of 235)',fill='',title='Transfer direction') + pn.coord_flip() + pn.guides(fill=pn.guide_legend(ncol=1)) + pn.theme_minimal(base_size=12) + pn.theme(legend_position='bottom',plot_title=pn.element_text(weight='bold',size=12),legend_text=pn.element_text(size=11),axis_text=pn.element_text(size=10.5),axis_title_x=pn.element_text(size=11)))
artifacts['direction']=save(direction,'blog-v2-effect-direction',direction_plot,3.6,5.2)
script=Path(__file__)
record={'purpose':'Plot frozen values and a closed arithmetic transfer only; no scientific model calls','sources':{str(p.relative_to(ROOT)):sha(p) for p in SOURCES},'filters':{'feed':'W01_001 training row by sample_id','level':'LF03 observed_recovery_pct; R/Python; five methods; 256 unique heldouts per cell','direction':'dose_token 0.05; R/Python; selected forest/STEPR OLS; SUCCESS only; 235 supported and 21 OFF_SUPPORT per cell'},'keys':{'feed':['sample_id','stage','part'],'level':['study','language','response','method','model_id','sample_id'],'direction':['language','method','model_id','sample_id','dose_token']},'derivation':{'feed':'Quartz fraction minus 0.05 and Kaolinite plus 0.05; mass percent = 100*fraction','level':'sqrt(mean((pred-observed)^2)) per model/language; units percentage points','direction':'count Delta_pred < -1e-8, within ±1e-8, > 1e-8 pp among SUCCESS; oracle negative in all 235'},'presentation':'full-precision CSV arithmetic; no numeric rounding before plotting; printed article values rounded separately','plotnine_boundary':'Polars to pandas once per figure after filtering and aggregation','software':{'polars':pl.__version__,'plotnine':pn.__version__,'pandas':pd.__version__,'matplotlib':matplotlib.__version__},'artifacts':artifacts,'script':str(script.relative_to(ROOT))}
record['script_sha256']=sha(script)
(ROOT/'blog-v2-provenance.json').write_text(json.dumps(record,indent=2)+'\n')
