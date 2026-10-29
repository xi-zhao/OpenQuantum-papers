import numpy as np
from full import L as Lfull
from pblock import Nmat
from scipy.linalg import eigvals
n,d,z=5,2,1e-3
ev=np.linalg.eigvals(Lfull(n,d,z).toarray())+n
re=np.round(ev.real,6); vals,cnt=np.unique(re,return_counts=True)
order=np.argsort(-vals)
chains={m:np.round(np.sort(eigvals(Nmat(n-m,d,z)).real)[::-1]-m*0,6) for m in range(n+1)}
print('distinct Re (top 14) with multiplicity:')
for i in order[:14]:
    v=vals[i]; hit=[ (m,np.argmin(abs(chains[m]-v))) for m in chains if np.min(abs(chains[m]-v))<1e-5]
    print('  %.6f  x%3d  chain(n-m) match:'%(v,cnt[i]),hit[:3])
