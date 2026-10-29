import numpy as np, json
from mf import roots
d=json.load(open('../data/g1_modes.json')); N=40; k=np.arange(N+1)/N
for key,res in d.items():
    D=float(key.split(',')[1]); r=roots(D,2.0,10.)
    for m in res:
        P=np.array(m['P']); P=P/np.max(abs(P))
        # number of sign changes (nodes) in P(k) above 2% of max, and centre of |P|
        s=np.sign(P[abs(P)>0.02]); nodes=int((np.diff(s)!=0).sum())
        lo=(abs(P[k<0.25])).sum()/abs(P).sum()
        print(key,'lam=%.4f nodes=%d weight(n<0.25)=%.3f  MF roots %s'%(m['lam'],nodes,lo,np.round(r,3)))
