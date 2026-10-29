import numpy as np, scipy.sparse.linalg as sla
from ex import liou
from mfdyn import rhs, steady
from mf import roots
V=10.; Om=2.0
# (1) reviewer counterexample, N=30, Delta=3.19149277256 : slowest modes by real part
L=liou(30,3.19149277256,Om,V)
w=[]
for sig in [-0.5,-0.7+2.9j,-0.7-2.9j,-0.7+1.5j,-0.7-1.5j,-0.7+4.5j,-0.7-4.5j]:
    w+=list(sla.eigs(L,k=8,sigma=sig,return_eigenvectors=False,tol=1e-12))
w=np.unique(np.round(np.array(w),9)); w=w[np.argsort(-w.real)]
print('N=30 slowest by Re:',np.round(w[:7],6))
# (2) MF Jacobian at fold points: compare with lambda[(lambda+1)^2+Om^2/2]
def J(D,Om,n):
    Dt=D-V*n; c=(-1j*Om/2*(1-2*n))/(0.5-1j*Dt); y=np.array([n,c.real,c.imag])
    f=lambda y: np.array(rhs(0,y,lambda t:(D,Om,V))); h=1e-7
    return np.array([(f(y+h*e)-f(y-h*e))/(2*h) for e in np.eye(3)]).T
from folds import fold_curve
fc=fold_curve(); 
for k in [500,3000,9000,15000]:
    n,D,O=fc[k]; ev=np.sort_complex(np.linalg.eigvals(J(D,O,n)))
    print('fold pt (D,Om)=(%.3f,%.3f) eig='%(D,O),np.round(ev,5),' predicted -1±i Om/sqrt2 =',np.round(-1+1j*O/np.sqrt(2),5))
