import numpy as np
from math import comb, lgamma
from pblock import Nmat
from t17 import perronvec
for d in (1,2,3):
    for zeta in (1e-3,1e-1,1,10):
        row=[]
        for n in (20,40,80,160):
            c=perronvec(Nmat(n,d,zeta/n**(d-1))); k=np.arange(n+1)
            lw=np.log(np.maximum(c,1e-300))+np.array([lgamma(n+1)-lgamma(j+1)-lgamma(n-j+1) for j in k])
            w=np.exp(lw-lw.max()); w/=w.sum(); row.append(round((k*w).sum()/n,4))
        print('d=%d zeta=%g filling:'%(d,zeta),row)
