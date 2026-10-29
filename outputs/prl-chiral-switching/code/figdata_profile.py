import numpy as np, json, scipy.sparse.linalg as sla
from ex import liou, obs_vec, ops
from xel import diagk
V=10.;Om=2.0
g=dict((int(r[0]),r) for r in json.load(open('../data/gN.json')))
out={}
for N in (32,48):
    D=g[N][1]; L=liou(N,D,Om,V); tr,nv=obs_vec(N)
    w0,v0=sla.eigs(L,k=1,sigma=1e-9,tol=1e-14); P,_=diagk(N,v0[:,0]/(tr@v0[:,0]),'right'); P=P.real
    # modes slightly off the crossing (adiabatic identity), D_off = D - 0.06
    Lo=liou(N,D-0.06,Om,V); w,R=sla.eigs(Lo,k=4,sigma=-0.76,tol=1e-13)
    prof=[]
    for i in np.argsort(abs(w+0.76))[:2]:
        rr,A=diagk(N,R[:,i],'right'); prof.append(dict(lam=w[i].real,A=(A/A.max()).tolist()))
    out[N]=dict(D=D,Pss=(P/P.max()).tolist(),modes=prof)
json.dump(out,open('../data/fig_profile.json','w')); print('ok')
