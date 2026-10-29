import numpy as np, itertools, mpmath as mp
from math import factorial, comb
from t11 import G
from pblock import Nmat
def spinops(n):
    # Boolean lattice basis ordered as in G: subsets by size then lexicographic
    subs=[frozenset(s) for k in range(n+1) for s in itertools.combinations(range(n),k)]
    idx={s:i for i,s in enumerate(subs)}; Jp=np.zeros((len(subs),)*2); Jm=np.zeros_like(Jp)
    for s in subs:
        for q in range(n):
            if q in s: Jm[idx[s-{q}],idx[s]]+=1
            else: Jp[idx[s|{q}],idx[s]]+=1
    return Jp,Jm
for n,d in [(6,3),(8,3),(7,4)]:
    Jp,Jm=spinops(n); z=0.37
    diff=np.max(np.abs(G(n,d,z)-(Jm+z*np.linalg.matrix_power(Jp,d)/factorial(d))))
    print(f'n={n} d={d}: |G - (J- + z J+^d/d!)|max = {diff:.1e}')
mp.mp.dps=50
n,d,z=6,3,mp.mpf('0.37')
ev=mp.eig(mp.matrix(G(n,d,float(z)).tolist()))[0]
pred=[]
for j in range(n//2+1):
    m=n-2*j; mult=comb(n,j)-(comb(n,j-1) if j else 0)
    pred+=list(mp.eig(mp.matrix(Nmat(m,d,float(z) if m>=d else 0.0).tolist()))[0])*mult
used=[False]*len(pred); worst=0
for e in ev:
    k=min((i for i in range(len(pred)) if not used[i]),key=lambda i:abs(pred[i]-e)); used[k]=True; worst=max(worst,abs(pred[k]-e))
print('n=6 d=3, 50-digit: worst |diff| =',mp.nstr(worst,3))
