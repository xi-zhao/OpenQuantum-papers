import numpy as np, json
from pblock import Nmat
from math import factorial
from t6 import perron
out={}
for d in (2,3,4):
    ns=list(range(d,321)); mu=[perron(Nmat(n,d,1.0))[0] for n in ns]
    C=(d+1)*d**(-d/(d+1))*factorial(d)**(-1/(d+1))*4**(-d/(d+1))
    out[d]=dict(n=ns,mu=mu,C=C)
    m=np.array(mu); n=np.array(ns)
    print(d,'C=%.5f'%C,'mu/(C n^a) at n=40,160,320:',[round(m[n==k][0]/(C*k**(2*d/(d+1))),4) for k in (40,160,320)],
          ' gap/(C a n^{a-1}):',[round((m[n==k][0]-m[n==k-1][0])/(C*2*d/(d+1)*k**((d-1)/(d+1))),4) for k in (40,160,320)])
json.dump(out,open('../data/perron_mu.json','w'))
