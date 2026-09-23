import argparse
import csv
import hashlib
import html
import io
import json
import platform
import sys
from pathlib import Path
import numpy as np
import scipy
from . import __version__
from .core import Config, Refusal, analyze, require_keys, positive, ingest


def canonical(obj):
    return json.dumps(obj,sort_keys=True,indent=2,allow_nan=False)+'\n'

def digest(b): return hashlib.sha256(b).hexdigest()

def write_json(path,obj): path.write_text(canonical(obj))

def table(path,rows):
    if not rows: return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def source_hash():
    return digest(b''.join(p.name.encode()+p.read_bytes() for p in sorted(Path(__file__).parent.glob('*.py'))))


def run(dataset,config_path,out):
    data=json.loads(dataset.read_text()); config=Config(**json.loads(config_path.read_text()))
    result,objects,fields,widths=analyze(data,config)
    # Write only after successful input validation and inference. Never overwrite a run.
    out.mkdir(parents=True,exist_ok=False)
    write_json(out/'dataset.json',data); write_json(out/'config.json',result['config'])
    write_json(out/'results.json',result); table(out/'statistics.csv',result['statistics'])
    table(out/'objects.csv',objects)
    field_rows=[]
    for (s,f),meta in fields.items():
        rr=[r for r in objects if r['specimen_id']==s and r['field_id']==f]
        row={'specimen_id':s,'field_id':f,'area':meta['area'],'n_objects':len(rr),'object_density':len(rr)/meta['area']}
        for k in range(1,config.max_order+1):
            val=sum(r['feret_mean']**k for r in rr)
            row[f'feret_mean_raw_{k}']=val/len(rr) if rr else None
            row[f'feret_mean_power_{k}_per_area']=val/meta['area']
        field_rows.append(row)
    table(out/'fields.csv',field_rows)
    np.save(out/'directional_widths.npy',widths,allow_pickle=False)
    bad=sum(r['status']!='OK' for r in result['statistics'])
    write_json(out/'gates.json',{'input':'PASS','scope':'OBSERVED_2D_ONLY','mode':'EXACT_POLYGON_KERNEL',
                              'intervals':'REVIEW' if bad else 'PASS_WITH_ASSUMPTIONS',
                              'non_ok_statistics':bad,'legacy_kernel_parity':'NOT_CLAIMED'})
    esc=lambda x:html.escape(str(x))
    trs=''.join('<tr>'+''.join('<td>'+esc(r[k])+'</td>' for k in ('feature','statistic','estimate','ci_low','ci_high','unit','status'))+'</tr>' for r in result['statistics'])
    report='''<!doctype html><meta charset="utf-8"><title>SANGHA | Field Geometry</title>
<style>body{background:#f5f2e9;color:#202824;font:15px system-ui;margin:40px auto;max-width:1200px;padding:0 24px}h1{font-size:38px}small{letter-spacing:2px}table{border-collapse:collapse;width:100%;background:#fffdf6}td,th{padding:9px;border-bottom:1px solid #d6dace;text-align:left}th{position:sticky;top:0;background:#dfe7dd}li{margin:8px 0}.note{border-left:4px solid #355b47;padding:16px;background:#e6ecdf}a{color:#355b47}</style>
<small>SANGHA RESEARCH SYSTEMS / OBSERVATION → ESTIMATE → WITNESS</small>
<h1>Field geometry</h1>'''
    report+=f"<p>{result['n_specimens']} specimens · {result['n_fields']} fields · {result['n_objects']} objects · {esc(data['unit'])}</p>"
    report+=f"<div class='note'>{esc(config.weighting)} · {config.confidence:.0%} pointwise percentile intervals · {config.bootstrap} specimen-cluster replicates. Observed 2D geometry; synthetic data if using the bundled demo.</div>"
    report+='<ul>'+''.join('<li>'+esc(w)+'</li>' for w in result['warnings'])+'</ul>'
    report+='<p><a href="statistics.csv">Statistics CSV</a> · <a href="objects.csv">Objects CSV</a> · <a href="fields.csv">Fields CSV</a></p>'
    report+='<table><thead><tr><th>Feature</th><th>Statistic</th><th>Estimate</th><th>CI low</th><th>CI high</th><th>Unit</th><th>Status</th></tr></thead><tbody>'+trs+'</tbody></table>'
    (out/'report.html').write_text(report)
    files={p.name:digest(p.read_bytes()) for p in sorted(out.iterdir()) if p.is_file()}
    write_json(out/'receipt.json',{'schema':'sangha.geometry.receipt.v1','version':__version__,
         'source_sha256':source_hash(),'input_sha256':digest(dataset.read_bytes()),
         'config_input_sha256':digest(config_path.read_bytes()),
         'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
         'outputs_sha256':files,'rng':'numpy.default_rng/PCG64','claim':'computational reproducibility; no signature or coverage certification'})
    return result


def verify(out):
    receipt=json.loads((out/'receipt.json').read_text())
    for name,expected in receipt['outputs_sha256'].items():
        if Path(name).name!=name or digest((out/name).read_bytes())!=expected:
            raise Refusal('receipt mismatch: '+name)
    return {'status':'PASS','scope':'listed output byte integrity; receipt itself is unsigned'}


def from_labels(manifest_path,out):
    """Raster labels → pixel-cell boundaries selected by centroid reference frame."""
    m=json.loads(manifest_path.read_text())
    require_keys(m,['schema','unit','observation','fields'])
    if m['schema']!='sangha.geometry.labels.v1': raise Refusal('unknown labels schema')
    data={'schema':'sangha.geometry.dataset.v1','unit':m['unit'],'observation':m['observation'],'fields':[],'objects':[]}
    for f in m['fields']:
        require_keys(f,['specimen_id','field_id','image','pixel_size_x','pixel_size_y','frame'])
        sx=positive(f['pixel_size_x']); sy=positive(f['pixel_size_y'])
        path=(manifest_path.parent/f['image']).resolve()
        if path.suffix.lower()=='.npy': a=np.load(path,allow_pickle=False)
        else:
            try: import tifffile
            except ImportError as e: raise Refusal('TIFF requires pip install tifffile') from e
            a=tifffile.imread(path)
        if a.ndim!=2 or a.dtype.kind not in 'ui' or np.any(a<0):
            raise Refusal('image must be a 2D nonnegative integer label image; background=0')
        frame=np.asarray(f['frame'],float)
        if frame.shape!=(4,) or not np.isfinite(frame).all(): raise Refusal('frame must be [x0,y0,x1,y1]')
        x0,y0,x1,y1=frame
        if not (0<=x0<x1<=a.shape[1] and 0<=y0<y1<=a.shape[0]): raise Refusal('invalid reference frame')
        data['fields'].append({'specimen_id':f['specimen_id'],'field_id':f['field_id'],
           'area':float((x1-x0)*(y1-y0)*sx*sy),'pixel_size_x':sx,'pixel_size_y':sy,
           'sampling':'complete_objects_by_reference_point'})
        for label in np.unique(a):
            if label==0: continue
            yy,xx=np.where(a==label); cx=xx.mean()+.5; cy=yy.mean()+.5
            if not (x0<=cx<x1 and y0<=cy<y1): continue
            if xx.min()==0 or yy.min()==0 or xx.max()==a.shape[1]-1 or yy.max()==a.shape[0]-1:
                raise Refusal(f'label {label} touches image edge: enlarge guard region')
            # Pixel cells, not centers: a single-pixel object has nonzero area.
            p=np.column_stack([xx,yy]); corners=np.concatenate([p+[dx,dy] for dx,dy in ((0,0),(1,0),(1,1),(0,1))])
            from scipy.spatial import ConvexHull
            p=corners[ConvexHull(corners).vertices]
            data['objects'].append({'specimen_id':f['specimen_id'],'field_id':f['field_id'],
               'object_id':str(label),'points':p.tolist(),'complete':True})
    # Validate converted data before writing.
    from .core import ingest
    ingest(data,Config())
    if out.exists() or out.with_suffix('.import.json').exists(): raise Refusal('output already exists')
    write_json(out,data)
    provenance={'schema':'sangha.geometry.import.v1','manifest_sha256':digest(manifest_path.read_bytes()),
        'images':{f['image']:digest((manifest_path.parent/f['image']).read_bytes()) for f in m['fields']},
        'output_sha256':digest(out.read_bytes()),'boundary':'convex hull of foreground pixel cells',
        'selection':'centroid in half-open reference frame; selected edge-touching labels refused'}
    write_json(out.with_suffix('.import.json'),provenance)
    return {'status':'PASS','objects':len(data['objects'])}


def main():
    p=argparse.ArgumentParser(description='SANGHA Field Geometry v'+__version__)
    sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('run'); r.add_argument('dataset',type=Path); r.add_argument('--config',required=True,type=Path); r.add_argument('--out',required=True,type=Path)
    v=sub.add_parser('verify'); v.add_argument('out',type=Path)
    for verb in ('validate','measure'):
        c=sub.add_parser(verb); c.add_argument('dataset',type=Path); c.add_argument('--config',required=True,type=Path)
        if verb=='measure': c.add_argument('--out',required=True,type=Path)
    sub.add_parser('doctor')
    i=sub.add_parser('from-labels'); i.add_argument('manifest',type=Path); i.add_argument('--out',required=True,type=Path)
    args=p.parse_args()
    try:
        if args.command=='run':
            result=run(args.dataset,args.config,args.out)
            print(canonical({'status':'COMPLETE','out':str(args.out),'n_statistics':len(result['statistics'])}))
        elif args.command=='verify': print(canonical(verify(args.out)))
        elif args.command=='doctor':
            print(canonical({'status':'PASS','python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'source_sha256':source_hash()}))
        elif args.command in ('validate','measure'):
            data=json.loads(args.dataset.read_text()); cfg=Config(**json.loads(args.config.read_text()))
            fields,rows,widths=ingest(data,cfg)
            if args.command=='measure':
                args.out.mkdir(parents=True,exist_ok=False)
                table(args.out/'objects.csv',rows); np.save(args.out/'directional_widths.npy',widths,allow_pickle=False)
                write_json(args.out/'receipt.json',{'schema':'sangha.geometry.receipt.v1','source_sha256':source_hash(),
                    'input_sha256':digest(args.dataset.read_bytes()),'config_input_sha256':digest(args.config.read_bytes()),
                    'outputs_sha256':{p.name:digest(p.read_bytes()) for p in sorted(args.out.iterdir())}})
            print(canonical({'status':'PASS','operator':args.command,'n_objects':len(rows),'n_fields':len(fields)}))
        else: print(canonical(from_labels(args.manifest,args.out)))
    except (Refusal,ValueError,TypeError,KeyError,OSError) as e:
        print(canonical({'status':'REFUSE','reason':str(e)}),file=sys.stderr); return 2
    return 0

if __name__=='__main__': sys.exit(main())
