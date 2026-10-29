# quasi-steady filling <k>/n of the Perron state of N_n, and whether it collapses in zeta=z n^{d-1}
import numpy as np
from pblock import Nmat
def perronvec(N,tol=1e-14,it=500000):
    A=N+np.eye(len(N))*(N.max()+1); v=np.ones(len(N))
    for i in range(it):
        w=A@v; w/=w.sum()
        if np.max(np.abs(w-v))<tol: break
        v=w
    return v
if __name__=="__main__":
 for d in (1,2,3):
     for zeta in (1e-3,1e-1,10):
         row=[]
         for n in (20,40,80,160):
             v=perronvec(Nmat(n,d,zeta/n**(d-1))); k=np.arange(n+1); row.append(round((k*v).sum()/n,4))
         print('d=%d zeta=%g <k>/n:'%(d,zeta),row)
 