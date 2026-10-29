import numpy as np
from full import L as Lfull
from pblock import Nmat
from scipy.linalg import eigvals
mu=lambda n,d,z: np.max(eigvals(Nmat(n,d,z)).real)
n,d,z=5,2,1e-3
ev=np.linalg.eigvals(Lfull(n,d,z).toarray())+n
m1=mu(n-1,d,z); sel=ev[np.abs(ev.real-m1)<1e-7]
print('count at Re=mu(n-1):',len(sel)); print(np.round(np.sort_complex(sel),5))
print('2cos q:',np.round(2*np.cos(2*np.pi*np.arange(n)/n),5))
