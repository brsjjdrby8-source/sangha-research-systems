from dataclasses import dataclass, asdict
import numpy as np
from scipy.spatial import ConvexHull, QhullError
from scipy.spatial.distance import pdist

class Refusal(ValueError):
    """An input or inferential contract was not satisfied."""

@dataclass(frozen=True)
class Config:
    schema: str = 'sangha.geometry.config.v1'
    mode: str = 'EXACT'
    max_order: int = 4
    angles: int = 180
    bootstrap: int = 2000
    seed: int = 1729
    confidence: float = 0.95
    weighting: str = 'specimen_equal'
    resampling: str = 'specimen_cluster'
    target: str = 'observed_2d'

    def __post_init__(self):
        if self.schema != 'sangha.geometry.config.v1' or self.mode != 'EXACT':
            raise Refusal('unsupported schema or mode')
        for name, lo, hi in [('max_order',4,12),('angles',16,4096),('bootstrap',200,100000),('seed',0,2**32-1)]:
            v = getattr(self,name)
            if type(v) is not int or not lo <= v <= hi:
                raise Refusal(f'{name} must be an integer in [{lo}, {hi}]')
        if type(self.confidence) not in (int,float) or not .5 < self.confidence < 1:
            raise Refusal('confidence must lie between 0.5 and 1')
        if self.weighting not in ('specimen_equal','object_pooled'):
            raise Refusal('unsupported weighting')
        if self.resampling != 'specimen_cluster' or self.target != 'observed_2d':
            raise Refusal('v0.1 supports specimen-cluster inference for observed 2D only')


def geometry(points, sx=1., sy=1., angles=180):
    """Exact polygon-hull extrema/mean; quadrature for angular distribution.

    Coordinates may be unordered boundary samples; area/perimeter explicitly
    refer to their convex hull, never the potentially concave original object.
    """
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2 or len(p)<3 or not np.isfinite(p).all():
        raise Refusal('object requires >=3 finite 2D boundary points')
    if not np.isfinite([sx,sy]).all() or min(sx,sy)<=0:
        raise Refusal('calibration must be positive and finite')
    p = (p-p[0])*[sx,sy]
    try:
        h = ConvexHull(p)
    except QhullError as e:
        raise Refusal('degenerate object boundary') from e
    v = p[h.vertices]
    edges = np.roll(v,-1,axis=0)-v
    normals = np.column_stack([-edges[:,1],edges[:,0]])
    normals /= np.linalg.norm(normals,axis=1)[:,None]
    fmin = float(np.ptp(v@normals.T,axis=0).min())
    fmax = float(pdist(v).max())
    theta = np.arange(angles)*np.pi/angles
    widths = np.ptp(v@np.array([np.cos(theta),np.sin(theta)]),axis=0)
    return {'feret_mean': float(h.area/np.pi), 'feret_min':fmin,
            'feret_max':fmax, 'feret_ratio':fmax/fmin,
            'hull_area':float(h.volume), 'hull_perimeter':float(h.area)}, widths


def require_keys(d, required, optional=()):
    if not isinstance(d,dict) or set(d)-set(required)-set(optional) or set(required)-set(d):
        raise Refusal(f'expected keys {required}; optional {optional}')


def identifier(x):
    if not isinstance(x,str) or not x.strip():
        raise Refusal('IDs and units must be nonempty strings')
    return x


def positive(x):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not np.isfinite(x) or x<=0:
        raise Refusal('areas and calibrations must be positive finite numbers')
    return float(x)


def ingest(data, config):
    require_keys(data, ['schema','observation','unit','fields','objects'])
    if data['schema']!='sangha.geometry.dataset.v1' or data['observation'] not in ('section_2d','projection_2d'):
        raise Refusal('unsupported dataset schema or observation type')
    identifier(data['unit'])
    if not isinstance(data['fields'],list) or not data['fields'] or not isinstance(data['objects'],list):
        raise Refusal('fields and objects must be lists; fields must not be empty')
    fields={}; rows=[]; widths=[]; seen=set()
    for f in data['fields']:
        require_keys(f,['specimen_id','field_id','area','pixel_size_x','pixel_size_y','sampling'])
        key=(identifier(f['specimen_id']),identifier(f['field_id']))
        if key in fields: raise Refusal('duplicate field')
        if f['sampling']!='complete_objects_by_reference_point':
            raise Refusal('declare complete objects selected by an unbiased reference-point frame')
        fields[key]={**f, **{k:positive(f[k]) for k in ('area','pixel_size_x','pixel_size_y')}}
    for obj in data['objects']:
        require_keys(obj,['specimen_id','field_id','object_id','points','complete'])
        key=(identifier(obj['specimen_id']),identifier(obj['field_id']))
        oid=key+(identifier(obj['object_id']),)
        if oid in seen or key not in fields: raise Refusal('duplicate object or unknown field')
        if obj['complete'] is not True: raise Refusal('truncated object: acquire full boundary with guard region')
        seen.add(oid); f=fields[key]
        g,w=geometry(obj['points'],f['pixel_size_x'],f['pixel_size_y'],config.angles)
        rows.append({'specimen_id':key[0],'field_id':key[1],'object_id':obj['object_id'],**g})
        widths.append(w)
    if not rows: raise Refusal('no objects: distribution moments undefined')
    # Canonical ordering makes object/field input order irrelevant.
    order=sorted(range(len(rows)),key=lambda i:tuple(rows[i][k] for k in ('specimen_id','field_id','object_id')))
    return dict(sorted(fields.items())), [rows[i] for i in order], np.array([widths[i] for i in order])


def moments(x,w,order):
    """Weighted empirical plug-in moments; not finite-sample unbiased estimators."""
    x=np.asarray(x,float); w=np.asarray(w,float); w=w/w.sum()
    mu=float(np.sum(w*x)); d=x-mu
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        out={f'raw_{k}':float(np.sum(w*(x**k))) for k in range(1,order+1)}
        out.update({f'central_{k}':float(np.sum(w*(d**k))) for k in range(2,order+1)})
        var=out['central_2']; sd=np.sqrt(var)
        out.update(mean=mu,variance=var,sd=float(sd),cv=float(sd/mu) if mu else np.nan,
                   skewness=out['central_3']/sd**3 if sd else np.nan,
                   excess_kurtosis=out['central_4']/var**2-3 if var else np.nan)
    sort=np.argsort(x,kind='stable'); xs=x[sort]; cw=np.cumsum(w[sort]); cw[-1]=1.
    for q in (.05,.25,.5,.75,.95):
        out[f'q{round(100*q):02}']=float(xs[min(np.searchsorted(cw,q),len(xs)-1)])
    out['iqr']=out['q75']-out['q25']
    return out


def summarize(fields,rows,widths,config,draw):
    # draw contains specimen IDs WITH multiplicity. Never collapse duplicates.
    idx=[]; weights=[]; densities=[]; burdens=[]
    total_area=0.; total_count=0; total_burden=np.zeros(config.max_order)
    for s in draw:
        ii=[i for i,r in enumerate(rows) if r['specimen_id']==s]
        area=sum(f['area'] for f in fields.values() if f['specimen_id']==s)
        total_area+=area; total_count+=len(ii)
        if ii:
            idx.extend(ii)
            weights.extend([1/len(ii) if config.weighting=='specimen_equal' else 1.]*len(ii))
        b=np.array([sum(rows[i]['feret_mean']**k for i in ii)/area for k in range(1,config.max_order+1)])
        densities.append(len(ii)/area); burdens.append(b); total_burden+=b*area
    out={}
    if idx:
        for feature in ('feret_mean','feret_min','feret_max','feret_ratio','hull_area','hull_perimeter'):
            for stat,val in moments([rows[i][feature] for i in idx],weights,config.max_order).items():
                out[f'{feature}/{stat}']=val
        # This is the joint object x direction distribution, not moments of mean width.
        ww=np.asarray(weights)
        for stat,val in moments(widths[idx].ravel(),np.repeat(ww/config.angles,config.angles),config.max_order).items():
            out[f'directional_feret/{stat}']=val
        means=widths[idx].mean(axis=1); within=widths[idx].var(axis=1)
        normalized=ww/ww.sum()
        out['directional_feret/within_object_variance']=float(normalized@within)
        out['directional_feret/between_object_variance']=float(normalized@((means-normalized@means)**2))
    out['field/object_density']=float(np.mean(densities) if config.weighting=='specimen_equal' else total_count/total_area)
    burden=np.mean(burdens,axis=0) if config.weighting=='specimen_equal' else total_burden/total_area
    out.update({f'field/feret_mean_power_{k}_per_area':float(burden[k-1]) for k in range(1,config.max_order+1)})
    return out


def analyze(data,config):
    fields,rows,widths=ingest(data,config)
    specimens=sorted({k[0] for k in fields}); n=len(specimens)
    empty=[s for s in specimens if not any(r['specimen_id']==s for r in rows)]
    if empty and config.weighting=='specimen_equal':
        raise Refusal('specimen_equal object distribution undefined for empty specimens; use object_pooled or revise estimand')
    point=summarize(fields,rows,widths,config,specimens)
    reps=[]
    if n>=2:
        rng=np.random.default_rng(config.seed)
        for _ in range(config.bootstrap):
            reps.append(summarize(fields,rows,widths,config,rng.choice(specimens,n).tolist()))
    result=[]; tail=(1-config.confidence)/2
    for key,value in point.items():
        vals=np.array([r.get(key,np.nan) for r in reps]); valid=np.isfinite(vals)
        status='OK'; lo=hi=None
        if not np.isfinite(value): status='UNDEFINED_OR_OVERFLOW'
        elif n<2: status='NO_INDEPENDENT_REPLICATION'
        elif not valid.all(): status='UNDEFINED_BOOTSTRAP_REPLICATES'
        else:
            lo,hi=map(float,np.quantile(vals,[tail,1-tail]))
            if lo==hi: status='DEGENERATE_INTERVAL'
        feature,stat=key.split('/')
        power=2 if feature=='hull_area' else (0 if feature=='feret_ratio' else 1)
        if feature=='field':
            power=-2 if stat=='object_density' else int(stat.split('_')[3])-2
        elif stat in ('cv','skewness','excess_kurtosis'): power=0
        elif stat in ('variance','within_object_variance','between_object_variance'): power*=2
        elif stat.startswith(('raw_','central_')): power*=int(stat.split('_')[1])
        unit='1' if power==0 else f"{data['unit']}^{power}"
        result.append({'feature':feature,'statistic':stat,'unit':unit,
                       'estimate':float(value) if np.isfinite(value) else None,
                       'ci_low':lo,'ci_high':hi,'confidence':config.confidence,
                       'ci_scope':'pointwise','status':status,'valid_replicates':int(valid.sum())})
    warnings=['Intervals are approximate percentile bootstrap intervals, not guaranteed finite-sample coverage.',
              'Sampling uncertainty excludes segmentation and calibration error.',
              'Higher moments may be unstable; no finite-moment or tail adequacy claim is made.',
              'Angular moments use deterministic quadrature; angular error is not a sampling CI.',
              'Frame/complete-boundary sampling declarations are user assertions, not verified from images.']
    if n<10: warnings.append('Fewer than 10 independent specimens: population interval reliability is limited.')
    if (1-config.confidence)*config.bootstrap/2<20: warnings.append('Few bootstrap replicates in each confidence tail.')
    if empty: warnings.append('Empty specimens retained in density; pooled object distribution has no objects from them.')
    return {'schema':'sangha.geometry.result.v1','config':asdict(config),'observation':data['observation'],
            'unit':data['unit'],'n_specimens':n,'n_fields':len(fields),'n_objects':len(rows),
            'warnings':warnings,'statistics':result},rows,fields,widths
