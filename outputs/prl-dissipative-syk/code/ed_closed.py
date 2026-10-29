import sys, numpy as np
sys.path.insert(0,'.')
from syk_ed import majoranas
import itertools
from math import factorial
N=int(sys.argv[1]); q=4; rng=np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 0)
psi=majoranas(N); d=2**(N//2)
from math import comb
sig=np.sqrt(1/comb(N-1,q-1))
H=np.zeros((d,d),complex)
for idx in itertools.combinations(range(N),q):
    M=psi[idx[0]]
    for i in idx[1:]: M=M@psi[i]
    H+=rng.normal(0,sig)*M.toarray()
H=(H+H.conj().T)/2
E,V=np.linalg.eigh(H)
ts=np.linspace(0,8,161); C=np.zeros_like(ts)
for i in range(4):
    p=V.conj().T@(psi[i].toarray())@V; w=np.abs(p)**2
    dE=E[:,None]-E[None,:]
    for k,t in enumerate(ts): C[k]+= (w*np.cos(dE*t)).sum()/d/4
np.save(f'../data/Cclosed_N{N}.npy',np.vstack([ts,C]))
for k in range(0,161,10): print(ts[k],C[k])
