import numpy as np, sys
from math import comb, log2
rng=np.random.default_rng(int(sys.argv[3]) if len(sys.argv)>3 else 1)
def haar(d):
    v=rng.normal(size=d)+1j*rng.normal(size=d); return v/np.linalg.norm(v)
def depol(D,n,r):
    # apply (1-r) rho + r Tr_k rho (x) I/2 on each qubit
    for k in range(n):
        X=D.reshape(2**k,2,2**(n-k-1),2**k,2,2**(n-k-1))
        tr=X[:,0,:,:,0,:]+X[:,1,:,:,1,:]
        Y=(1-r)*X
        Y[:,0,:,:,0,:]+=r*tr/2; Y[:,1,:,:,1,:]+=r*tr/2
        D=Y.reshape(2**n,2**n)
    return D
def T(n,r):
    d=2**n; a,b=haar(d),haar(d)
    D=np.outer(a,a.conj())-np.outer(b,b.conj())
    D=depol(D,n,r); return 0.5*np.abs(np.linalg.eigvalsh(D)).sum()
def pred(n,r):
    q=3*r/4; s0=-log2(1-q); s1=-log2(r/4) if r>0 else 1e9
    return sum(comb(n,k)*q**k*(1-q)**(n-k) for k in range(n+1) if (n-k)*s0+k*s1<n)
n=int(sys.argv[1]); ns=int(sys.argv[2])
for r in np.linspace(0.10,0.50,9):
    v=[T(n,r) for _ in range(ns)]
    print(n, round(r,3), round(np.mean(v),4), round(np.std(v)/np.sqrt(ns),4), round(pred(n,r),4))
