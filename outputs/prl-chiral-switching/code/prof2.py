import numpy as np, scipy.sparse.linalg as sla, json
from ex import liou, obs_vec, ops
from xel import diagk
V=10.;Om=2.0
av=dict((int(r[0]),r) for r in json.load(open('../data/gN.json')))
for N in [32,48]:
    D=av[N][1]; L=liou(N,D,Om,V); tr,nv=obs_vec(N)
    w0,v0=sla.eigs(L,k=1,sigma=1e-9,tol=1e-14); P,_=diagk(N,v0[:,0]/(tr@v0[:,0]),'right'); P=P.real
    w,R=sla.eigs(L,k=4,sigma=-0.76,tol=1e-13); wl,Lv=sla.eigs(L.T.tocsc(),k=4,sigma=-0.76,tol=1e-13)
    k=np.arange(N+1)/N; lp=np.log(np.abs(P)+1e-300)
    print('N=%d D=%.4f  ln P_ss rel. peak at n=0.10,0.157,0.25,0.4:'%(N,D),np.round(np.interp([0.10,0.157,0.25,0.4],k,lp-lp.max())/N,4))
    for i in np.argsort(abs(w+0.76))[:2]:
        rr,A=diagk(N,R[:,i],'right'); j=np.argmin(abs(wl-w[i])); ll,B=diagk(N,Lv[:,j],'left')
        lr=np.log(A+1e-300); lr-=lr.max(); lb=np.log(B+1e-300); lb-=lb.max()
        print('  lam=%.4f  ln|r| /N at n=.10,.157,.25,.4: %s   ln|l| /N: %s'%(w[i].real,np.round(np.interp([0.10,0.157,0.25,0.4],k,lr)/N,4),np.round(np.interp([0.10,0.157,0.25,0.4],k,lb)/N,4)))
