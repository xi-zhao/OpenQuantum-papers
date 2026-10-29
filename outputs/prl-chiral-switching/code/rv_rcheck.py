import numpy as np, scipy.sparse.linalg as sla
from ex import liou
V=10.;Om=2.0;N=32; D=4.97911683469913
w,v=sla.eigs(liou(N,D,Om,V),k=6,sigma=-0.8945,tol=1e-13); o=np.argsort(abs(w+0.8945))[:2]
a,b=v[:,o[0]],v[:,o[1]]; print('eigs',w[o],'overlap',abs(np.vdot(a,b))/np.linalg.norm(a)/np.linalg.norm(b))
# fine real scan of the gap
for d in np.linspace(D-0.004,D+0.004,9):
    w=sla.eigs(liou(N,d,Om,V),k=4,sigma=-0.8945,return_eigenvectors=False,tol=1e-13); w=np.sort(w[np.argsort(abs(w+0.8945))][:2].real)
    print(round(d,5),np.round(w,6),'gap %.2e'%(w[1]-w[0]))
