import numpy as np, sys
sys.argv=["x","2","1"]
from sim import haar, depol
from bound import comb
rng=np.random.default_rng(5)
def lo_hi(n,r):
    q=3*r/4; d=2**n; ks=range(n+1)
    N=[comb(n,k)*3**k for k in ks]; p=[(1-q)**(n-k)*(q/3)**k for k in ks]
    lo=max(sum(N[k]*p[k] for k in ks if k<=K)-sum(N[k] for k in ks if k<=K)/d for K in range(n+1))
    hi=min(sum(N[k]*p[k] for k in ks if p[k]>t)+0.5*np.sqrt(d*sum(N[k]*p[k]**2 for k in ks if p[k]<=t)) for t in sorted(set(p))+[0])
    return max(lo,0),min(hi,1)
n=8; d=2**n
for r in [0.1,0.25,0.4,0.6]:
    v=[]
    for _ in range(30):
        a=haar(d); D=depol(np.outer(a,a.conj()),n,r)-np.eye(d)/d
        v.append(0.5*np.abs(np.linalg.eigvalsh(D)).sum())
    # product state (most robust?) for comparison
    D=depol(np.diag([1]+[0]*(d-1)).astype(complex),n,r)-np.eye(d)/d
    print(r, "haar",round(np.mean(v),4), "bounds",[round(x,4) for x in lo_hi(n,r)], "|0..0>",round(0.5*np.abs(np.linalg.eigvalsh(D)).sum(),4))
