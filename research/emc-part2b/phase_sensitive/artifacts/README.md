# Phase-sensitive run archive

These six deterministic gzip tar archives preserve every regular file in the completed `run/` tree, including the root blocked attempt, `correction-001` and its saved replay, and the final `correction-002` attempt. Archive member names begin with `run/`. No raw run files are staged in Git.

The complete source index is `source-files.sha256` (116 files, 475,212,777 bytes, index SHA256 `2876fedc6ebc72a3c9f83874ab5f717f51c4259eea6d1824e03b93a9bd9973dd`). `archive-files.sha256` authenticates the six compressed archives. Each archive is under 90 MiB. Archive ordering and member ordering are lexical; gzip has zero modification time and no embedded filename; tar members have zero modification time, uid/gid 0, blank owner/group names, and mode 0644. The source file bytes are unchanged.

## Restore and check file bytes

Run from the repository root. Restore into a checkout that does not already have `research/emc-part2b/phase_sensitive/run/`:

```bash
set -euo pipefail
study=research/emc-part2b/phase_sensitive
(cd "$study/artifacts" && sha256sum -c archive-files.sha256)
for archive in "$study"/artifacts/*.tar.gz; do
  tar -xzf "$archive" -C "$study"
done
(cd "$study" && sha256sum -c artifacts/source-files.sha256)
```

All 116 `sha256sum` checks must report `OK`. The archive index itself can be checked with `sha256sum "$study/artifacts/source-files.sha256"` against the digest above. Archive entries were also independently streamed and checked against every source file's SHA256 and byte count during packaging.

## Study verifier and legacy prerequisites

The frozen verifier includes absolute paths in its historical seed inventory and requires the original checkout at `/home/matt/Documents/GitHub/writing`. It also requires the existing accepted PR4/PR5 files, `/tmp/emc-part2-workspace/writing/runs/` historical specifications, pinned simCoDaR source at `/home/matt/Documents/GitHub/simCoDaR`, its isolated installation at `/tmp/emc-simcodar-d9/lib`, and the Python/R environment recorded in `run/attempts/correction-002/manifests/environment.json`. A relocated checkout can pass the file-byte check above while the verifier reports `seed inventory changed after freeze`; relocation is not evidence of altered archived files. Do not edit the frozen seed or source manifests to make relocation pass.

At the original path, with those prerequisites available, run this **read-only** check from the repository root. The named Python environment is the executable recorded by the final manifest; another environment must reproduce that manifest exactly.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/tmp/bif-v3-world1-validation-env-20261023/bin/python \
research/emc-part2b/phase_sensitive/measurement_experiment.py \
  --config research/emc-part2b/phase_sensitive/preregistration.json \
  --output research/emc-part2b/phase_sensitive/run \
  --review research/emc-part2b/phase_sensitive/PREEXECUTION_REVIEW.json \
  --attempt correction-002 --verify
```

The verifier returned `PASS` with 350,250 predictions, 21,250 cohort metric rows, and zero prediction failures on the unchanged original tree. Since the restored tree has byte-identical files, restoration at the same path with the recorded prerequisites presents the same inputs to `--verify`. This is an artifact check; it does not constitute final scientific or publication clearance.

## Reproduce the archives

This command creates the same archive bytes from the original absolute run path. It refuses an oversized archive and checks every member against its source before finishing. Run it only with the original source tree present; it does not remove or edit source files.

```bash
python3 - <<'PY'
from pathlib import Path
import gzip, hashlib, tarfile
src = Path('/home/matt/Documents/GitHub/writing/research/emc-part2b/phase_sensitive/run')
dst = Path('research/emc-part2b/phase_sensitive/artifacts')
dst.mkdir(exist_ok=True)
def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''): h.update(b)
    return h.hexdigest()
files = sorted(p for p in src.rglob('*') if p.is_file())
assert all(not p.is_symlink() for p in files)
source = {'run/' + p.relative_to(src).as_posix(): (sha(p), p.stat().st_size) for p in files}
index = ''.join(f'{h}  {rel}\n' for rel,(h,n) in source.items()).encode()
(dst/'source-files.sha256').write_bytes(index)
def group(rel):
    parts = rel.split('/')
    if len(parts)<3 or parts[1]!='attempts': return 'run-root'
    a = parts[2]
    if len(parts)>=4 and parts[3]=='results': return f'run-{a}-results'
    if len(parts)>=5 and parts[3]=='audit' and parts[4]=='replay': return f'run-{a}-replay'
    return f'run-{a}-evidence'
groups = {}
for rel in source: groups.setdefault(group(rel), []).append(rel)
archives = []
for name,members in sorted(groups.items()):
    out = dst/(name+'.tar.gz')
    with out.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',filename='',mtime=0,compresslevel=9) as gz:
            with tarfile.open(fileobj=gz,mode='w',format=tarfile.GNU_FORMAT) as tar:
                for rel in members:
                    p = src/rel.removeprefix('run/')
                    t = tarfile.TarInfo(rel)
                    t.size=p.stat().st_size; t.mode=0o644
                    t.uid=t.gid=0; t.uname=t.gname=''; t.mtime=0
                    with p.open('rb') as f: tar.addfile(t,f)
    assert out.stat().st_size < 90*1024*1024
    archives.append(out)
seen = {}
for out in archives:
    with tarfile.open(out,'r:gz') as tar:
        for m in tar:
            assert m.isfile() and m.name in source and m.name not in seen
            assert not m.name.startswith('/') and '..' not in Path(m.name).parts
            h=hashlib.sha256(); n=0
            f=tar.extractfile(m)
            for b in iter(lambda:f.read(1024*1024),b''): h.update(b); n+=len(b)
            assert (h.hexdigest(),n)==source[m.name]
            seen[m.name]=(h.hexdigest(),n)
assert seen==source
assert all((sha(src/rel.removeprefix('run/')),(src/rel.removeprefix('run/')).stat().st_size)==value for rel,value in source.items())
(dst/'archive-files.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in archives))
PY
```
