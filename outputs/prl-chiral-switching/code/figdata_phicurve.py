import numpy as np, json, sys
from phiwkb2 import phi
D=float(sys.argv[1]); out=[]
for n in [0.12,0.157,0.20,0.25,0.30,0.35]:
    o=phi(D,2.0,n); out.append([n,o['Phi']]); print(n,o['Phi'],flush=True)
json.dump(out,open('../data/fig_phicurve.json','w'))
