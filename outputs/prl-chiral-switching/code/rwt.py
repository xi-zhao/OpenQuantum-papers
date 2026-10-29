import numpy as np, scipy.sparse.linalg as sla
from ex import liou, obs_vec
from dyn import ems
V=10.;Om=2.0;N=32
# real-parameter EP pair at the right fold (complex-pair window D in [4.9779,4.9791]); spectral weight of N_e
# in relaxation from the MF high-branch-like metastable state (quench from inside the wedge)
D=4.9785
L=liou(N,D,Om,V); tr,nv=obs_vec(N)
w,R=sla.eigs(L,k=8,sigma=-0.894,tol=1e-12); wl,Lv=sla.eigs(L.T.tocsc(),k=8,sigma=-0.894,tol=1e-12)
g,rss,(hi,lo)=ems(N,4.4,Om,V)   # initial: high metastable state prepared inside wedge, then quench to D
tot=[]
for i in np.argsort(-w.real):
    j=np.argmin(abs(wl-w[i])); l=Lv[:,j]; r=R[:,i]; c=(l@hi)/(l@r); A=(nv@r)*c/N
    tot.append((w[i],A)); print('lam=%s  weight in n(t)=%s'%(np.round(w[i],5),np.round(A,5)))
