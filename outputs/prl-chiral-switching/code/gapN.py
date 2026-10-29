import numpy as np, json
from ex import eig
V=10.
pts=[(4.0,2.0),(4.3,2.0),(3.8,2.0),(4.6,2.8),(4.4,2.5)]
res={}
for (D,Om) in pts:
    row=[]
    for N in [10,15,20,25,30,35,40,45,50,60]:
        w=eig(N,D,Om,V,k=3,sigma=0.0); l1=-w[1].real; row.append([N,l1]); print(D,Om,N,l1,flush=True)
    res[f"{D},{Om}"]=row
json.dump(res,open('../data/gapN.json','w'))
