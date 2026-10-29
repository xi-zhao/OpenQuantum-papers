import numpy as np, json, scipy.sparse.linalg as sla
from cep2 import newton
from ex import liou, eig
V=10.; Om=2.0
res=[]
D,t=3.288842+0.0004477j,-0.75913+0.00082j
for N in range(52,72,2):
    D,t,r=newton(N,D,t); res.append([N,D.real,D.imag,t.real,t.imag,r]); print('EP',N,D,t,r,flush=True)
    json.dump(res,open('../data/cep_big.json','w'))
