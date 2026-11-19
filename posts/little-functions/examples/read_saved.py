"""Read frozen compositions, predict with one admitted saved fit, evaluate frozen rows."""
import csv
import hashlib
import json
import math
import pickle
import sys
from pathlib import Path

if len(sys.argv) != 4:
    raise SystemExit('usage: read_saved.py REPO ROOT OUTPUT')
if not all(Path(a).is_absolute() for a in sys.argv[1:]):
    raise SystemExit('REPO, ROOT, OUTPUT must be absolute')
repo, root, out = (Path(a).resolve(strict=True) if i < 2 else Path(a).resolve() for i, a in enumerate(sys.argv[1:]))
if out == repo or out == root or repo in out.parents or root in out.parents:
    raise SystemExit('REPO, ROOT, OUTPUT must be absolute; OUTPUT must be outside inputs')
out.mkdir(parents=True, exist_ok=True)
parts = ['Chalcopyrite', 'Chalcocite', 'Pyrite', 'Quartz', 'Kaolinite']
E = repo / 'examples/publication_sc020/artifacts/SC020-v1-reproduce-002/root'
registry = json.loads((root / 'registry/SC020-v1.json').read_text())
model_id = 'SC020-v1-LF01-Python-y_nonlinear-select_lr_ml_rf'
entry = next(m for m in registry['models'] if m['model_id'] == model_id)
assert entry['response'] == 'y_nonlinear' and entry['root_object_path'].startswith('models/SC020-v1/')
model_file = repo / 'examples/publication_sc020' / entry['root_object_path']
assert hashlib.sha256(model_file.read_bytes()).hexdigest() == entry['object_sha256']
sys.path.insert(0, str(repo / 'src/python'))
from lr_ml import predict_lr_ml, paired_contrast
import pandas as pd
with model_file.open('rb') as stream:
    wrapper = pickle.load(stream)
fit = wrapper['fit']
assert [list(p) for p in fit['selected_pairs']] == entry['selected_pairs']
fixture = pd.read_csv(repo / 'examples/copper/copper.csv')
assert len(fixture) == 240 and list(fixture.loc[:, parts].columns) == parts
assert (fixture.loc[:, parts] > 0).all().all()
assert ((fixture.loc[:, parts].sum(axis=1) - 1).abs() <= 1e-10).all()
assert fixture['split'].value_counts().to_dict() == {'train': 180, 'test': 60}
ids = ['W07_sample_1', 'W07_sample_2', 'W07_sample_3', 'W07_sample_4', 'W07_sample_9', 'W07_sample_10']
rows = fixture.set_index('sample_id').loc[ids]
X = rows.loc[:, parts].copy()
levels_pct = predict_lr_ml(fit, X)
assert len(levels_pct) == 6
print('Python selected pairs:', fit['selected_pairs'])
print('Python first complete composition:', X.iloc[0].to_dict())
print('Python first recovery prediction (%):', f'{levels_pct[0]:.6f}')
D = pd.read_csv(E / 'evidence/SC020-v1/lf01_effects.csv', dtype={'dose_token': str})
records = []
for dose in (0, .01, .025, .05):
    result = paired_contrast(fit, X, delta=dose, donor='Quartz', recipient='Kaolinite')
    frozen = D[(D.language == 'Python') & (D.dose_token == ('0' if dose == 0 else str(dose)))].set_index('sample_id').loc[ids]
    for sid, row in result.iterrows():
        expected = frozen.loc[sid]
        assert max(abs(row.before-float(expected.before)), abs(row.after-float(expected.after)), abs(row.contrast-float(expected.effect))) <= 1e-10
        assert abs(row.contrast - (row.after - row.before)) <= (1e-12 if dose == 0 else 1e-10)
        if dose == 0:
            assert abs(row.contrast) <= 1e-12 and abs(row.after-row.before) <= 1e-12
        records.append({'sample_id': sid, 'dose_token': '0' if dose == 0 else str(dose), 'dose': dose, 'before_percent': row.before, 'after_percent': row.after, 'effect_pp': row.contrast})
pd.DataFrame(records).to_csv(out / 'lf01_predictions.csv', index=False)
print('Python LF01 paired rows:', len(records))

def read_rows(name):
    with (E / 'evidence/SC020-v1' / name).open(newline='') as stream:
        return list(csv.DictReader(stream))

levels, effects = read_rows('levels.csv'), read_rows('lf03_effects.csv')
assert len({r['sample_id'] for r in levels if r['study'] == 'LF03' and r['language'] == 'Python' and r['method'] == 'stepr_ols'}) == 256
summary = []
for language in ('R', 'Python'):
    for method in ('stepr_ols', 'stepr_rf', 'select_lr_ml_rf', 'ilr_rf', 'raw_rf'):
        L = [r for r in levels if r['study'] == 'LF03' and r['language'] == language and r['method'] == method]
        F = [r for r in effects if r['language'] == language and r['method'] == method and r['dose_token'] == '0.05']
        assert len(L) == len(F) == 256 and len({r['sample_id'] for r in L}) == len({r['sample_id'] for r in F}) == 256
        assert {r['sample_id'] for r in L} == {r['sample_id'] for r in F}
        assert len({r['model_id'] for r in L+F}) == 1
        good = [r for r in F if r['status'] == 'SUCCESS']
        excluded = [r for r in F if r['status'] == 'OFF_SUPPORT']
        assert len(good) == 235 and len(excluded) == 21
        assert all(r['error'] == '' and r['Delta_pred'] == '' for r in excluded)
        obs = math.sqrt(sum((float(r['pred'])-float(r['observed']))**2 for r in L)/256)
        quiet = math.sqrt(sum((float(r['pred'])-float(r['noiseless']))**2 for r in L)/256)
        err = math.sqrt(sum((float(r['Delta_pred'])-float(r['Delta_oracle']))**2 for r in good)/235)
        signs = [sum(('zero' if abs(float(r['Delta_pred'])) <= 1e-8 else ('negative' if float(r['Delta_pred']) < 0 else 'positive')) == label for r in good) for label in ('negative', 'zero', 'positive')]
        assert all(float(r['Delta_oracle']) < -1e-8 for r in good)
        summary.append({'language': language, 'method': method, 'n_level': 256, 'n_effect': 235, 'excluded': 21, 'observed_rmse_pp': obs, 'noiseless_rmse_pp': quiet, 'effect_rmse_pp': err, 'negative': signs[0], 'zero': signs[1], 'positive': signs[2]})
with (out / 'lf03_table.csv').open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(summary[0])); writer.writeheader(); writer.writerows(summary)
print('LF03 frozen-row arithmetic:')
for row in summary:
    print(row['language'], row['method'], f"{row['observed_rmse_pp']:.6f}", f"{row['noiseless_rmse_pp']:.6f}", f"{row['effect_rmse_pp']:.6f}", f"{row['negative']}/{row['zero']}/{row['positive']}")
