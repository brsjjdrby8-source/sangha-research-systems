"""Deterministic archive entries; package target must run tests before this script."""
import hashlib
import json
import zipfile
from pathlib import Path
root=Path(__file__).resolve().parents[1]; dest=root/'dist'; dest.mkdir(exist_ok=True)
name='SANGHA_Field_Geometry_v0.1.0'; archive=dest/(name+'.zip')
files=[p for p in sorted(root.rglob('*')) if p.is_file() and not any(x in p.relative_to(root).parts for x in ('dist','__pycache__','.git')) and '.pyc'!=p.suffix and not any(x.endswith('.egg-info') for x in p.parts)]
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p in files:
        info=zipfile.ZipInfo(name+'/'+str(p.relative_to(root)),date_time=(2026,1,1,0,0,0))
        info.compress_type=zipfile.ZIP_DEFLATED; info.external_attr=0o100644<<16
        z.writestr(info,p.read_bytes())
sha=hashlib.sha256(archive.read_bytes()).hexdigest()
archive.with_suffix('.zip.sha256').write_text(sha+'  '+archive.name+'\n')
print(json.dumps({'archive':str(archive),'sha256':sha,'files':len(files)},indent=2))
