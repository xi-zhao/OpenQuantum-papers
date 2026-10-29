import numpy as np, sys
from ex import liou
import scipy.sparse.linalg as sla
V=10.; Om=float(sys.argv[2]); N=int(sys.argv[1])
D0,D1,n=float(sys.argv[3]),float(sys.argv[4]),int(sys.argv[5]); sig=float(sys.argv[6])
best=(9,0)
for D in np.linspace(D0,D1,n):
    w=sla.eigs(liou(N,D,Om,V),k=6,sigma=sig,return_eigenvectors=False,tol=1e-12)
    w=w[np.argsort(abs(w-sig))][:3]; w=w[np.argsort(-w.real)]
    gap=abs(w[0]-w[1]); best=min(best,(gap,D))
    print(round(D,5),' '.join(f"{x.real:.6f}{x.imag:+.6f}j" for x in w),flush=True)
print('min gap',best)
