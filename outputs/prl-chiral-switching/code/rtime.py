import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from ex import liou, obs_vec
from dyn import ems
V=10.;Om=2.0;N=32
L=liou(N,4.9785,Om,V); tr,nv=obs_vec(N); g,rss,(hi,lo)=ems(N,4.4,Om,V)
# exact n(t) via expm_multiply; compare with fit e^{lam t}(a+bt) vs pure exponential on t in [4,14]
ts=np.linspace(0,16,161); vs=sla.expm_multiply(L,hi,start=0,stop=16,num=161,endpoint=True)
n=np.array([(nv@v).real/N for v in vs]); ninf=n[-1]
w=sla.eigs(L,k=3,sigma=-0.894,return_eigenvectors=False); lam=w[np.argmin(abs(w+0.894))].real
m=(ts>5)&(ts<14); y=(n[m]-0.04391)*np.exp(-lam*ts[m])
b,a=np.polyfit(ts[m],y,1); print('lam=%.4f  (n-n_ss)e^{-lam t} = a+bt : a=%.4f b=%.4f  rel. curvature check:'%(lam,a,b),np.round(np.polyfit(ts[m],y,2),5))
print('n(t) samples:',np.round(n[::20],5))
