import json
from pathlib import Path
import numpy as np
root=Path(__file__).resolve().parents[1]; rng=np.random.default_rng(1123)
fields=[]; objects=[]
for s in range(12):
    scale=float(rng.lognormal(0,.2))
    for f in range(3):
        fields.append(dict(specimen_id=f's{s:02}',field_id=f'f{f}',area=10000,pixel_size_x=.3,pixel_size_y=.4,sampling='complete_objects_by_reference_point'))
        for i in range(5):
            a,b=scale*rng.uniform(10,25,2); t=np.arange(24)*2*np.pi/24; angle=rng.uniform(0,np.pi)
            p=np.column_stack([a*np.cos(t),b*np.sin(t)])@np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])+[100,100]
            objects.append(dict(specimen_id=f's{s:02}',field_id=f'f{f}',object_id=f'o{i}',points=p.tolist(),complete=True))
data=dict(schema='sangha.geometry.dataset.v1',observation='section_2d',unit='um',fields=fields,objects=objects)
(root/'examples/demo.dataset.json').write_text(json.dumps(data,indent=2)+'\n')
(root/'examples/default.config.json').write_text(json.dumps(dict(schema='sangha.geometry.config.v1',mode='EXACT',max_order=4,angles=180,bootstrap=2000,seed=1729,confidence=.95,weighting='specimen_equal',resampling='specimen_cluster',target='observed_2d'),indent=2)+'\n')
a=np.zeros((32,32),np.uint16); a[10:15,10:18]=1; a[20:23,20:25]=2; np.save(root/'examples/labels.npy',a)
m=dict(schema='sangha.geometry.labels.v1',unit='um',observation='section_2d',fields=[dict(specimen_id='s01',field_id='f01',image='labels.npy',pixel_size_x=.3,pixel_size_y=.4,frame=[4,4,28,28])])
(root/'examples/labels.manifest.json').write_text(json.dumps(m,indent=2)+'\n')
