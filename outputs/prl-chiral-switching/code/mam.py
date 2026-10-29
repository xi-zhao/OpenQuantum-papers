"""Minimum-action (heteroclinic) solver for the reduced 3-DOF WKB Hamiltonian, via BVP + continuation.
   Instanton A -> S: leaves A in its 3D unstable subspace, enters S in its 3D stable subspace; h=0 automatically.
   Action S = int p.dq  (switching exponent: rate ~ exp(-N S))."""
import numpy as np, json, sys
from scipy.integrate import solve_bvp
from scipy.linalg import null_space
from wkbr import F, JAC, H, fixed
from mf import roots
Om=2.0; V=10.
def subspaces(y,D,stable):
    w,R=np.linalg.eig(JAC(y,D,Om,V)); sel=(w.real<0) if stable else (w.real>0)
    B=[]
    for i in np.flatnonzero(sel):
        B+=[R[:,i].real,R[:,i].imag] if abs(w[i].imag)>1e-10 else [R[:,i].real]
    U=np.linalg.qr(np.array(B).T)[0][:,:3]; return U, null_space(U.T)   # basis, complement (for BCs)
def guess(D,src,dst,m=400):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[src]); S=fixed(D,Om,V,r[dst]); t=np.linspace(0,1,m)
    s=0.5*(1+np.tanh(6*(t-0.5)))/np.tanh(3)*np.tanh(3); s=(s-s[0])/(s[-1]-s[0])
    n=A[0]+(S[0]-A[0])*s; Y=np.zeros((6,m))
    for i,nn in enumerate(n): Y[:,i]=fixed(D,Om,V,nn)
    Y[3]=2.0*np.sin(np.pi*t)*np.sign(S[0]-A[0])*np.sign(r[1]-r[src]) # pn uphill
    return t,Y,A,S
def solve(D,src,dst,init=None,eps=1e-4,T0=30.,tol=1e-7):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[src]); S=fixed(D,Om,V,r[dst])
    UA,WA=subspaces(A,D,False); US,WS=subspaces(S,D,True)
    if init is None: t,Y,_,_=guess(D,src,dst); T=T0
    else: t,Y,T=init
    f=lambda t,y,p: p[0]*F(y,D,Om,V)
    fj=lambda t,y,p: (p[0]*JAC(y,D,Om,V), F(y,D,Om,V)[:,None,:])
    def bc(ya,yb,p): return np.concatenate([WA.T@(ya-A),[np.sum((ya-A)**2)-eps**2],WS.T@(yb-S)])
    sol=solve_bvp(f,bc,t,Y,p=[T],fun_jac=fj,tol=tol,max_nodes=200000,verbose=0)
    y=sol.sol(np.linspace(0,1,20001)); dq=np.diff(y[:3],axis=1); pm=(y[3:,1:]+y[3:,:-1])/2
    Sact=float(np.sum(pm*dq)); hmax=float(np.abs(H(y,D,Om,V)).max())
    return dict(D=D,src=src,dst=dst,S=Sact,T=float(sol.p[0]),ok=bool(sol.success),msg=sol.message,hmax=hmax,
                endgap=float(np.linalg.norm(y[:,-1]-S)),sol=sol)
if __name__=="__main__":
    D=float(sys.argv[1]) if len(sys.argv)>1 else 3.9
    for src,lab in [(2,'hi'),(0,'lo')]:
        R=solve(D,src,1); print('D=%.3f %s->saddle S=%.5f T=%.1f ok=%s hmax=%.1e endgap=%.1e %s'%(D,lab,R['S'],R['T'],R['ok'],R['hmax'],R['endgap'],R['msg']),flush=True)
