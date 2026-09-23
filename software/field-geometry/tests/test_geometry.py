import copy
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from sangha_geometry.core import Config, Refusal, geometry, moments, analyze, ingest, summarize
from sangha_geometry.cli import run, verify, from_labels


def fixture(n=4):
    fields=[]; objects=[]
    for i in range(n):
        for j in range(2):
            s=str(i); f=str(j)
            fields.append(dict(specimen_id=s,field_id=f,area=100,pixel_size_x=1,pixel_size_y=1,sampling='complete_objects_by_reference_point'))
            a=1+i+j
            objects.append(dict(specimen_id=s,field_id=f,object_id='a',complete=True,points=[[0,0],[a,0],[a,2],[0,2]]))
    return dict(schema='sangha.geometry.dataset.v1',observation='section_2d',unit='um',fields=fields,objects=objects)

class Tests(unittest.TestCase):
    def test_rectangle_analytic(self):
        g,_=geometry([[0,0],[3,0],[3,4],[0,4]])
        self.assertAlmostEqual(g['feret_mean'],14/np.pi)
        self.assertEqual(g['feret_max'],5); self.assertEqual(g['feret_min'],3)
        self.assertEqual(g['hull_area'],12)
    def test_anisotropic_calibration(self):
        g,_=geometry([[0,0],[1,0],[1,1],[0,1]],3,4)
        self.assertEqual(g['feret_max'],5)
    def test_rotation_translation(self):
        p=np.array([[0,0],[3,0],[3,4],[0,4]]); a=.37
        q=p@np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])+200
        g,_=geometry(p); h,_=geometry(q)
        for k in g:self.assertAlmostEqual(g[k],h[k])
    def test_reference_directional_parity(self):
        p=np.random.default_rng(10).normal(size=(50,2))
        _,w=geometry(p,angles=180)
        reference=[]
        for theta in np.arange(180)*np.pi/180:
            projections=[x*np.cos(theta)+y*np.sin(theta) for x,y in p]
            reference.append(max(projections)-min(projections))
        np.testing.assert_allclose(w,reference,rtol=1e-13,atol=1e-13)
    def test_circle_polygon_convergence(self):
        t=np.arange(4096)*2*np.pi/4096
        g,_=geometry(np.column_stack([np.cos(t),np.sin(t)]))
        self.assertAlmostEqual(g['feret_mean'],2,places=6)
    def test_known_moments(self):
        m=moments([1,2,3],[1,1,1],4)
        self.assertEqual(m['mean'],2); self.assertAlmostEqual(m['variance'],2/3)
        self.assertAlmostEqual(m['skewness'],0); self.assertAlmostEqual(m['excess_kurtosis'],-1.5)
    def test_degenerate(self):
        with self.assertRaises(Refusal): geometry([[0,0],[1,1],[2,2]])
        self.assertTrue(np.isnan(moments([1,1],[1,1],4)['skewness']))
    def test_contracts(self):
        for kw in ({'target':'3d'},{'max_order':3},{'angles':True},{'weighting':'area'}):
            with self.assertRaises(Refusal): Config(**kw)
        d=fixture(); d['objects'][0]['complete']=False
        with self.assertRaises(Refusal): ingest(d,Config())
    def test_duplicate(self):
        d=fixture(); d['objects'].append(d['objects'][0])
        with self.assertRaises(Refusal): ingest(d,Config())
    def test_seed_and_input_order(self):
        d=fixture(); c=Config(bootstrap=200,angles=16)
        a=analyze(d,c)[0]; d['objects'].reverse(); d['fields'].reverse()
        self.assertEqual(a,analyze(d,c)[0])
    def test_single_specimen(self):
        r=analyze(fixture(1),Config(bootstrap=200,angles=16))[0]
        self.assertTrue(all(x['ci_low'] is None for x in r['statistics']))
    def test_variance_decomposition(self):
        c=Config(angles=32); f,r,w=ingest(fixture(),c)
        s=summarize(f,r,w,c,['0','1','2','3'])
        self.assertAlmostEqual(s['directional_feret/variance'],s['directional_feret/within_object_variance']+s['directional_feret/between_object_variance'])
    def test_weighting_and_multiplicity(self):
        d=fixture(2); d['objects'].pop()
        c=Config(); f,r,w=ingest(d,c)
        s=summarize(f,r,w,c,['0','1','1'])
        means=[np.mean([x['feret_mean'] for x in r if x['specimen_id']==i]) for i in ['0','1']]
        self.assertAlmostEqual(s['feret_mean/mean'],(means[0]+2*means[1])/3)
    def test_receipt_tamper_and_replay(self):
        with tempfile.TemporaryDirectory() as t:
            t=Path(t); (t/'data.json').write_text(json.dumps(fixture()))
            (t/'config.json').write_text(json.dumps({'bootstrap':200,'angles':16}))
            run(t/'data.json',t/'config.json',t/'a'); run(t/'data.json',t/'config.json',t/'b')
            self.assertEqual((t/'a/receipt.json').read_bytes(),(t/'b/receipt.json').read_bytes())
            self.assertEqual(verify(t/'a')['status'],'PASS')
            (t/'a/objects.csv').write_text('tamper')
            with self.assertRaises(Refusal): verify(t/'a')
    def test_duplicating_objects_does_not_add_independence(self):
        d=fixture(); c=Config(bootstrap=200,angles=16)
        original=analyze(d,c)[0]
        for obj in list(d['objects']):
            clone=copy.deepcopy(obj); clone['object_id']='duplicate'; d['objects'].append(clone)
        repeated=analyze(d,c)[0]
        a=[r for r in original['statistics'] if r['feature']=='feret_mean']
        b=[r for r in repeated['statistics'] if r['feature']=='feret_mean']
        for x,y in zip(a,b):
            for key in ('estimate','ci_low','ci_high'):
                self.assertAlmostEqual(x[key],y[key])
    def test_units_and_empty_field(self):
        d=fixture(); d['fields'].append(dict(specimen_id='0',field_id='empty',area=100,pixel_size_x=1,pixel_size_y=1,sampling='complete_objects_by_reference_point'))
        r=analyze(d,Config(bootstrap=200,angles=16))[0]
        entries={(x['feature'],x['statistic']):x for x in r['statistics']}
        self.assertEqual(entries['hull_area','raw_4']['unit'],'um^8')
        self.assertEqual(entries['field','object_density']['unit'],'um^-2')
        self.assertEqual(entries['field','feret_mean_power_2_per_area']['unit'],'1')
        self.assertAlmostEqual(entries['field','object_density']['estimate'],(2/300+3*2/200)/4)
    def test_pixel_cell_adapter_and_border(self):
        with tempfile.TemporaryDirectory() as t:
            t=Path(t); a=np.zeros((10,10),int); a[4,4]=1; np.save(t/'a.npy',a)
            m={'schema':'sangha.geometry.labels.v1','observation':'section_2d','unit':'um','fields':[dict(specimen_id='s',field_id='f',image='a.npy',pixel_size_x=2,pixel_size_y=3,frame=[2,2,8,8])]}
            (t/'m.json').write_text(json.dumps(m)); from_labels(t/'m.json',t/'d.json')
            _,r,_=ingest(json.loads((t/'d.json').read_text()),Config())
            self.assertEqual(r[0]['hull_area'],6)
            a[:5,4]=1; np.save(t/'a.npy',a)
            m['fields'][0]['frame']=[0,0,10,10]; (t/'m.json').write_text(json.dumps(m))
            with self.assertRaises(Refusal): from_labels(t/'m.json',t/'e.json')

if __name__=='__main__': unittest.main()
