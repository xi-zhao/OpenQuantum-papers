import numpy as np, json
from ex import eig
from mf import roots
V=10.; Om=2.0
out={}
for N in [10,20,30,40,50]:
    rows=[]
    for D in np.linspace(3.0,5.2,23):
        w=eig(N,D,Om,V,k=5,sigma=0.0)
        rows.append([D,[[x.real,x.imag] for x in w[1:4]]])
        print(N,round(D,2),np.round(w[1:4],4),len(roots(D,Om,V)),flush=True)
    out[N]=rows
json.dump(out,open('../data/scanN_Om2.json','w'))
