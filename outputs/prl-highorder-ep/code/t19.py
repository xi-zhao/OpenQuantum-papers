import numpy as np, itertools
from t11 import Lfull, G
from pblock import Nmat
from math import comb, factorial
from scipy.linalg import eigvals
from scipy.optimize import brentq
rng=np.random.default_rng(7)
# (1) exact decomposition for more (n,d) with random H
for n,d in [(5,1),(5,3),(6,2),(6,4)]:
    A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); H=(A+A.conj().T)/2; eps=np.linalg.eigvalsh(H); z=0.05
    full=np.linalg.eigvals(Lfull(n,d,z,H))+n
    pred=[]
    for Aset in itertools.chain.from_iterable(itertools.combinations(range(n),k) for k in range(n+1)):
      for Bset in itertools.chain.from_iterable(itertools.combinations(range(n),k) for k in range(n+1)):
        if set(Aset)&set(Bset): continue
        m=n-len(Aset)-len(Bset); dE=eps[list(Aset)].sum()-eps[list(Bset)].sum()
        for j2 in range(m%2,m+1,2):  # 2j = m, m-2, ...
            k=(m-j2)//2; mult=comb(m,k)-(comb(m,k-1) if k else 0)
            pred+=list(-1j*dE+eigvals(Nmat(j2,d,z if j2>=d else 0.0)))*mult
    pred=np.array(pred)
    used=np.zeros(len(pred),bool); worst=0
    for e in full:
        i=np.argmin(np.where(used,np.inf,np.abs(pred-e))); used[i]=True; worst=max(worst,abs(pred[i]-e))
    print(f'n={n} d={d}: dim {len(full)} vs {len(pred)}, worst |diff|={worst:.1e}')
# (2) filling closed form vs numerics
from t17 import perronvec
from math import lgamma
for d in (1,2,3):
    C=(d+1)*d**(-d/(d+1))*factorial(d)**(-1/(d+1))*4**(-d/(d+1))
    for zeta in (1e-2,1,10):
        nu=brentq(lambda x:x+zeta*(1-x)**d/factorial(d)-C*zeta**(1/(d+1)),0,1)
        n=320; c=perronvec(Nmat(n,d,zeta/n**(d-1))); k=np.arange(n+1)
        lw=np.log(np.maximum(c,1e-300))+np.array([lgamma(n+1)-lgamma(j+1)-lgamma(n-j+1) for j in k]); w=np.exp(lw-lw.max()); w/=w.sum()
        print(f'd={d} zeta={zeta}: nu_WKB={nu:.4f}  nu(n=320)={(k*w).sum()/n:.4f}')
