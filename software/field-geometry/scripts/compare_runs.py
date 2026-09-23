import hashlib
import json
import sys
from pathlib import Path
left,right=map(Path,sys.argv[1:])
files=['results.json','statistics.csv','objects.csv','fields.csv','directional_widths.npy','gates.json','report.html']
checks={f:hashlib.sha256((left/f).read_bytes()).digest()==hashlib.sha256((right/f).read_bytes()).digest() for f in files}
print(json.dumps({'status':'PASS' if all(checks.values()) else 'FAIL','scientific_output_byte_parity':checks},indent=2))
sys.exit(0 if all(checks.values()) else 1)
