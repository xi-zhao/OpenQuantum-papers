import numpy as np, scipy.sparse.linalg as sla
from ex import liou
# encircle the complex EP at N=30 in the complex-Delta plane: check eigenvalue swap (square-root branch point)
V=10.; Om=2.0
D0=3.1914927823294166+0.010227874812525385j; lam0=-0.7542051742339325+0.014038503952625566j
def two(D,guess):
    w=sla.eigs(liou(30,D,Om,V),k=4,sigma=guess,return_eigenvectors=False,tol=1e-12); return w
r=0.004; phis=np.linspace(0,2*np.pi,81)
a=None
track=[]
for ph in phis:
    D=D0+r*np.exp(1j*ph); w=two(D,lam0); w=w[np.argsort(abs(w-lam0))][:2]
    if a is None: a=w[np.argmax(w.real)]
    else: a=w[np.argmin(abs(w-a))]
    track.append(a)
print('start',np.round(track[0],6),'end',np.round(track[-1],6),'other at start',np.round(two(D0+r,lam0)[np.argsort(abs(two(D0+r,lam0)-lam0))][:2],6))
