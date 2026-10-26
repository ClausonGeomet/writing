"""Small exportable diagnostics backed by CSV data."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
x=pd.read_csv(ROOT/'results/summary.csv').query("split=='test'")
fig,ax=plt.subplots(figsize=(6,4))
for method,g in x.groupby('method'):
    ax.scatter(g.assay_rmse,g.mineral_rmse,label=method,alpha=.8,s=24)
ax.set(xlabel='Held-out assay RMSE (mass fraction)',ylabel='Held-out mineral RMSE (mass fraction)')
ax.legend(fontsize=7,ncol=2);fig.tight_layout();fig.savefig(ROOT/'results/heldout_errors.png',dpi=180);plt.close(fig)
y=pd.read_csv(ROOT/'results/observational_equivalence.csv').query('replicate==0')
fig,ax=plt.subplots(figsize=(6,4))
ax.scatter(y.anatase_original,y.anatase_swapped,s=18,alpha=.8)
ax.plot([0,1],[0,1],color='gray',linewidth=1)
ax.set(xlabel='Original Anatase fraction',ylabel='Swapped Anatase fraction',title='Identical SiO2 and TiO2 assays for every pair')
fig.tight_layout();fig.savefig(ROOT/'results/ambiguous_swap.png',dpi=180);plt.close(fig)
