import numpy as np
from pblock import Nmat
from math import factorial
def perron(n,d,tol=1e-13,it=400000):
    # homogeneity: mu(z) = z^{1/(d+1)} mu(1). Choose z=n^{-(d-1)} (zeta=1) and diagonally rescale by WKB p*=const -> well conditioned
    z=float(n)**(-(d-1)); N=Nmat(n,d,z)
    k=np.arange(n+1); s=np.exp(-0.0*k)
    A=N+np.eye(n+1)*(N.max()+1); v=np.ones(n+1)
    for i in range(it):
        w=A@v; r=w/v; lo,hi=r.min(),r.max(); v=w/w.max()
        if hi-lo<tol*hi: break
    c=N.max()+1; return (0.5*(lo+hi)-c)/z**(1/(d+1)), (hi-lo)
if __name__=="__main__":
    import json
    out={}
    for d in (2,3,4):
        C=(d+1)*d**(-d/(d+1))*factorial(d)**(-1/(d+1))*4**(-d/(d+1)); a=2*d/(d+1)
        ns=list(range(max(d,4),321)); mu=[perron(n,d)[0] for n in ns]; m=np.array(mu); n=np.array(ns)
        out[d]=dict(n=ns,mu=mu,C=C)
        print(d,'C=%.5f'%C,'mu/(C n^a):',[round(m[n==k][0]/(C*k**a),5) for k in (20,40,80,160,320)],
          ' gap/(C a n^{a-1}):',[round((m[n==k][0]-m[n==k-1][0])/(C*a*k**(a-1)),5) for k in (20,40,80,160,320)])
    json.dump(out,open('../data/perron_mu.json','w'))
