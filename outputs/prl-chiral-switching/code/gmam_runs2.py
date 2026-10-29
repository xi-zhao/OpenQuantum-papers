import json, sys
from gmam_runs import val, ghost
kind,x=sys.argv[1],sys.argv[2]
r=val(float(x)) if kind=='val' else ghost(int(x))
print(json.dumps(r),flush=True)
open('../data/gmam_runs.jsonl','a').write(json.dumps(r)+'\n')
