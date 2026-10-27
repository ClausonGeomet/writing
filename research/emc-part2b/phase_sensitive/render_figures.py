"""Redraw the two report figures from a saved correction-002 paired summary."""
import argparse
import csv
import pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--paired-summary', type=pathlib.Path, required=True,
                    help='restored correction-002 results/paired_summary.csv')
parser.add_argument('--output-dir', type=pathlib.Path, required=True,
                    help='directory for the two report PNGs')
args = parser.parse_args()
dest = args.output_dir
dest.mkdir(parents=True, exist_ok=True)
with args.paired_summary.open(newline='') as saved_csv:
    rows = list(csv.DictReader(saved_csv))
labels = ['1 accurate', '2 weak', '3 mismatch', '4 soft', '4 hard']
keys = [('1','full'),('2','full'),('3','full'),('4','soft'),('4','hard')]
def get(key, stage, stratum, contrast):
    return next(r for r in rows if (r['cell_id'],r['support_group'],r['stage'],r['budget'],r['stratum'],r['contrast'],r['metric']) == (*key,stage,'5' if stage=='2' else '-1',stratum,contrast,'mass_RMSE'))
fig, ax = plt.subplots(figsize=(8,4.7))
for i,key in enumerate(keys):
    r=get(key,'2','unmeasured','A2_P-A2_R')
    ax.errorbar(float(r['mean_difference']),i,xerr=float(r['descriptive_se']),fmt='o',color='#17639c',capsize=3)
ax.axvline(0,color='0.4',lw=.8)
ax.set_yticks(range(5),labels); ax.invert_yaxis()
ax.set_xlabel('A2 phase purchase − A2 bulk repeat: mass RMSE')
ax.set_title('Independent unmeasured targets, k=5; bars = descriptive SE')
fig.tight_layout(); fig.savefig(dest/'transfer_summary.png',dpi=180); plt.close(fig)
fig, ax = plt.subplots(figsize=(8,4.7))
for i,key in enumerate(keys):
    d=get(key,'1','direct','A0_H-A0_T')
    t=get(key,'2','unmeasured','A2_P-A2_R')
    ax.scatter(float(d['mean_difference']),i-.13,color='#b04d37',marker='s',label='Direct A0_H−A0_T' if i==0 else None)
    ax.scatter(float(t['mean_difference']),i+.13,color='#17639c',marker='o',label='Unmeasured A2_P−A2_R' if i==0 else None)
ax.axvline(0,color='0.4',lw=.8)
ax.set_yticks(range(5),labels); ax.invert_yaxis()
ax.set_xlabel('Mean paired mass-RMSE difference; negative favors phase channel')
ax.set_title('Direct measurement and independent transfer are different estimands')
ax.legend(loc='lower right',frameon=False); fig.tight_layout(); fig.savefig(dest/'direct_vs_transfer.png',dpi=180); plt.close(fig)
