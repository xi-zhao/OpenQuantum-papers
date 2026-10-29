# exact decomposition check: full spectrum of L' (+n) == union over frozen sets S of { i(E_B-E_A) + spec G_{n-|S|} }
import numpy as np, itertools, scipy.sparse as sp
from full import ops
from math import comb
def Lfull(n,d,z,H):
    c=ops(n); D=2**n; Id=sp.identity(D,format='csr')
    Hm=sum(H[i,j]*c[i].T@c[j] for i in range(n) for j in range(n))
    Lm=-1j*(sp.kron(Hm,Id)-sp.kron(Id,Hm.T.conj().T))
    for j in range(n): Lm=Lm+sp.kron(c[j],c[j].conj())
    Lm=Lm-n*sp.identity(D*D)
    for xi in itertools.combinations(range(n),d):
        o=Id
        for j in xi: o=o@c[j].T
        Lm=Lm+z*sp.kron(o,o.conj())
    return Lm.toarray()
def G(m,d,z):
    subs=[frozenset(s) for k in range(m+1) for s in itertools.combinations(range(m),k)]
    idx={s:i for i,s in enumerate(subs)}; M=np.zeros((len(subs),)*2)
    for s in subs:
        for q in s: M[idx[s-{q}],idx[s]]+=1
        for xi in itertools.combinations(sorted(set(range(m))-s),d): M[idx[s|set(xi)],idx[s]]+=z
    return M
if __name__=="__main__":
    rng=np.random.default_rng(1); n,d,z=5,2,0.05
    A=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); H=(A+A.conj().T)/2   # random hopping, no translation symmetry
    eps=np.linalg.eigvalsh(H)
    full=np.linalg.eigvals(Lfull(n,d,z,H))+n
    pred=[]
    for Aset in itertools.chain.from_iterable(itertools.combinations(range(n),k) for k in range(n+1)):
      for Bset in itertools.chain.from_iterable(itertools.combinations(range(n),k) for k in range(n+1)):
        S=set(Aset)^set(Bset)
        if set(Aset)&set(Bset): continue      # enumerate frozen pair (A\B,B\A) once: A,B disjoint
        m=n-len(S); dE=eps[list(Aset)].sum()-eps[list(Bset)].sum()
        pred+=list(-1j*dE+np.linalg.eigvals(G(m,d,z)))
    pred=np.array(pred)
    print(len(full),len(pred))
    fs=np.sort_complex(np.round(full,6)); ps=np.sort_complex(np.round(pred,6))
    print('max |diff| after sorting:',np.max(np.abs(fs-ps)))
    # also exact z-scaling of the leading real parts
    for zz in (1e-4,1e-2,1.0):
        ev=np.linalg.eigvals(Lfull(n,d,zz,H))+n
        top=np.sort(np.unique(np.round(ev.real,8)))[::-1][:2]
        print(zz, top/zz**(1/3))
    