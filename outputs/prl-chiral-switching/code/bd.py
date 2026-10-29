import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.linalg import eig_banded, eigvals
# incoherent (rate-equation) reduction of the same model: k = number of excited atoms
#   r(n) = (Om^2/4)/((Delta - V n)^2 + 1/4);   k->k+1 : (N-k) r ;  k->k-1 : k (1 + r)
# MF: dn/dt = r(1-2n) - n  -> same fixed points as the Bloch model.
V=10.; Om=2.0
def r(n,D): return (Om**2/4)/((D-V*n)**2+0.25)
def gen(N,D):
    k=np.arange(N+1); n=k/N
    up=(N-k)*r(n,D); dn=k*(1+r(n,D))
    M=sp.diags([-(up+dn), up[:-1], dn[1:]],[0,-1,1],shape=(N+1,N+1),format='csc')  # column-stochastic: dp/dt = M p
    return M
def slow(N,D,k=4,sigma=-0.3):
    w=sla.eigs(gen(N,D),k=k,sigma=sigma,return_eigenvectors=False,tol=1e-12)
    return w[np.argsort(-w.real)]
if __name__=="__main__":
    for N in [50,200,800]:
        print(N,[np.round(slow(N,D,6,-0.3)[:4],4) for D in (3.0,3.4,3.55,3.8)])
