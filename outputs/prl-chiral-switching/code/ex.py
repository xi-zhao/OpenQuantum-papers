import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from qutip import piqs
# Exact permutation-invariant Liouvillian of N two-level Rydberg atoms, all-to-all V/N:
# H = -Delta Ne + Om Jx + (V/2N) Ne(Ne-1), local decay gamma=1 (PIQS, Shammah et al. PRA 98, 063815)
_cache={}
def ops(N):
    if N not in _cache:
        jx,jy,jz=piqs.jspin(N,basis='dicke')
        Ne=(jz+N/2).full().real
        Jx=jx.full()
        B=np.asarray(piqs.block_matrix(N).toarray()).astype(bool)
        idx=np.flatnonzero(B.flatten(order='F'))   # qutip column-stacking
        s=piqs.Dicke(N=N,emission=1.0)
        D=sp.csr_matrix(s.lindbladian().data.as_scipy())[idx][:,idx]
        _cache[N]=(Ne,Jx,B,idx,D)
    return _cache[N]
def liou(N,Delta,Om,V):
    Ne,Jx,B,idx,D=ops(N)
    H=-Delta*Ne+Om*Jx+(V/(2*N))*(Ne@Ne-Ne)
    H=sp.csr_matrix(H*B)   # H is block diagonal anyway
    I=sp.identity(H.shape[0],format='csr')
    Lh=-1j*(sp.kron(I,H)-sp.kron(H.T,I))
    Lh=sp.csr_matrix(Lh)[idx][:,idx]
    return (Lh+D).tocsc()
def obs_vec(N):
    Ne,Jx,B,idx,D=ops(N)
    tr=np.eye(Ne.shape[0]).flatten(order='F')[idx]
    nv=Ne.T.flatten(order='F')[idx]      # Tr(Ne rho) = sum_ij Ne_ji rho_ij
    return tr,nv
def eig(N,Delta,Om,V,k=6,sigma=0.0):
    L=liou(N,Delta,Om,V)
    w=sla.eigs(L,k=k,sigma=sigma,which='LM',return_eigenvectors=False)
    return w[np.argsort(-w.real)]
