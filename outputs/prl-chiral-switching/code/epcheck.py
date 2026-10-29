import numpy as np, sys
from ex import liou
import scipy.sparse.linalg as sla
V=10.; Om=2.0; N=int(sys.argv[1])
for D in np.linspace(float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4])):
    L=liou(N,D,Om,V)
    w=sla.eigs(L,k=4,sigma=float(sys.argv[5]),which='LM',return_eigenvectors=False,tol=1e-12)
    w=w[np.argsort(-w.real)]
    print(round(D,4),' '.join(f"{x.real:.6f}{x.imag:+.6f}j" for x in w),flush=True)
