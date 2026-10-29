import numpy as np
from ryd import *
L=10; Delta=-0.3; V2=0.4
st=np.array([s for s in range(1<<L) if all(not((s>>j)&1 and (s>>((j+1)%L))&1) for j in range(L))])
idx={s:i for i,s in enumerate(st)}; n=len(st); H=np.zeros((n,n)); T=np.zeros((n,n))
for i,s in enumerate(st):
    b=[(s>>j)&1 for j in range(L)]
    d=-Delta*sum(b)
    for j in range(L):
        for r in range(2,L//2+1):
            if r==L//2 and j>=L//2: continue
            d+=V2*(2/r)**6*b[j]*b[(j+r)%L]
    H[i,i]=d
    for j in range(L):
        c=s^(1<<j)
        if c in idx: H[idx[c],i]+=0.5
    T[idx[((s<<1)|(s>>(L-1)))&((1<<L)-1)],i]=1
pc=PeriodicChain(L)
for kpi in (0,1):
    P=sum(np.linalg.matrix_power(T,s)*((-1)**(s*kpi)) for s in range(L))/L
    U,S,_=np.linalg.svd(P); Q=U[:,S>1e-8]
    eb=np.linalg.eigvalsh(Q.T@H@Q)
    e=np.linalg.eigvalsh(pc.H(Delta,V2,kpi).toarray())
    print(kpi,e.size,eb.size,np.max(np.abs(np.sort(e)-np.sort(eb))))
# check O_pi matrix element: <n,pi| sum (-1)^j n_j |0,0>
H0=pc.H(Delta,V2,0).toarray(); Hp=pc.H(Delta,V2,1).toarray()
e0,v0=np.linalg.eigh(H0); ep,vp=np.linalg.eigh(Hp)
coeff,tgt=pc.ops(1); g=np.zeros(Hp.shape[0]); m=tgt>=0
np.add.at(g,tgt[m],coeff[m]*v0[m,0])
amp=vp.T@g
# brute
eall,vall=np.linalg.eigh(H); gs=vall[:,0]
stag=np.array([sum(((s>>j)&1)*(-1)**j for j in range(L)) for s in st])
w=(vall.T@(stag*gs))**2
print('sum weights ksector',np.sum(amp**2),'brute',np.sum(w))
for k in range(4): 
    j=np.argmin(abs(eall-ep[k])); print(ep[k]-e0[0],amp[k]**2, eall[j]-eall[0], w[j])
