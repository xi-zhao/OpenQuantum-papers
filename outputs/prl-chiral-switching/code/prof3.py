import numpy as np, scipy.sparse.linalg as sla, json
from ex import liou, obs_vec
from xel import diagk
V=10.;Om=2.0
av=dict((int(r[0]),r) for r in json.load(open('../data/gN.json')))
# exploration: where along n does the ghost mode's weight sit, and where does the hi-mode's weight decay, AWAY from the crossing
for N in [40,48]:
    for off in (-0.06,+0.04):
        D=av[N][1]+off; L=liou(N,D,Om,V)
        w,R=sla.eigs(L,k=4,sigma=-0.76,tol=1e-13)
        k=np.arange(N+1)/N
        for i in np.argsort(abs(w+0.76))[:2]:
            P,A=diagk(N,R[:,i],'right'); la=np.log(A+1e-300); la-=la.max(); im=np.argmax(A)
            print('N=%d D=%.4f lam=%.4f  peak n=%.3f  ln|r|/N @ n=.10,.157,.20,.30,.40: %s'%(N,D,w[i].real,k[im],np.round(np.interp([0.10,0.157,0.2,0.3,0.4],k,la)/N,4)))
