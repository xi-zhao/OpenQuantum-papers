# N is entrywise nonnegative & irreducible -> mu_max = Perron root, computed stably by power iteration (Collatz-Wielandt bounds)
import numpy as np
from pblock import Nmat
def perron(N,it=200000,tol=1e-12):
    v=np.ones(N.shape[0]); c=N.max()+1.0; A=N+c*np.eye(N.shape[0])
    for i in range(it):
        w=A@v; r=w/v; lo,hi=r.min(),r.max(); v=w/np.linalg.norm(w)
        if hi-lo<tol*hi: break
    return 0.5*(lo+hi)-c, (hi-lo)
if __name__=="__main__":
 for d in (2,3):
     for zeta in (1e-4,1e-2,1,100):
         out=[]
         for n in (20,40,80,160,320):
             m,err=perron(Nmat(n,d,zeta/n**(d-1))); out.append(round(m/n,4))
         print('d=%d zeta=%g mu/n:'%(d,zeta),out)
 