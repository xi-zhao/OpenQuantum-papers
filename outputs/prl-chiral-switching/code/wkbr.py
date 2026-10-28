"""Reduced (gauge-fixed, K=Tr(X rho)=1) 3-DOF WKB Hamiltonian.
   q=(n,u,v): rho_hat=[[1-n,u-iv],[u+iv,n]];  p=(pn,pu,pv): X_hat=[[pgg,(pu-ipv)/2],[(pu+ipv)/2,pee]],
   pgg = 1 - n pn - u pu - v pv,  pee = pgg + pn   (so Tr(X_hat rho_hat)=1, Tr(X_hat d rho_hat)=p.dq)."""
import numpy as np, sympy as sp
from wkbq import h as h8, q as Q8, p as P8, D_, O_, V_
n,u,v,pn,pu,pv=sp.symbols('n u v pn pu pv',real=True)
pgg=1-n*pn-u*pu-v*pv; pee=pgg+pn
sub={Q8[0]:1-n,Q8[1]:n,Q8[2]:u,Q8[3]:v,P8[0]:pgg,P8[1]:pee,P8[2]:pu,P8[3]:pv}
hr=sp.expand(h8.subs(sub))
qs=[n,u,v]; ps=[pn,pu,pv]; z=qs+ps
rhs_expr=[sp.diff(hr,s) for s in ps]+[-sp.diff(hr,s) for s in qs]
jac_expr=sp.Matrix(rhs_expr).jacobian(z)
_f=sp.lambdify([z,D_,O_,V_],rhs_expr,'numpy'); _Je=[[sp.lambdify([z,D_,O_,V_],jac_expr[i,j],'numpy') for j in range(6)] for i in range(6)]; _h=sp.lambdify([z,D_,O_,V_],hr,'numpy')
def F(y,D,Om,V):
    out=_f(list(y),D,Om,V); return np.array(np.broadcast_arrays(*out),float)
def JAC(y,D,Om,V):
    m=y.shape[1] if y.ndim==2 else None; yl=list(y)
    A=np.empty((6,6)+((m,) if m else ()))
    for i in range(6):
        for j in range(6):
            e=_Je[i][j](yl,D,Om,V); A[i,j]=e if m is None else np.broadcast_to(e,(m,))
    return A
def H(y,D,Om,V): return np.broadcast_to(_h(list(y),D,Om,V),np.shape(y[0])).astype(float)
def fixed(D,Om,V,n0):
    Dt=D-V*n0; c=(-1j*Om/2*(1-2*n0))/(0.5-1j*Dt); return np.array([n0,c.real,c.imag,0,0,0.])
if __name__=="__main__":
    from mfdyn import rhs as mfrhs
    from wkbf import flow as flow8
    from mf import roots
    D,Om,V=3.9,2.0,10.
    rng=np.random.default_rng(0)
    # (a) p=0 reproduces MF
    y=np.array([0.3,0.1,-0.2,0,0,0.]); print('MF check:',np.abs(F(y,D,Om,V)[:3]-np.array(mfrhs(0,y[:3],lambda t:(D,Om,V)))).max())
    # (b) reduced flow = projection of 8D flow (q_hat = rho/Tr rho, p from X_hat = X Tr rho)
    y=np.concatenate([rng.normal(size=3)*0.3+[0.3,0,0],rng.normal(size=3)*0.5]); nn,uu,vv,a,b,c=y
    g=1-nn*a-uu*b-vv*c; y8=np.array([1-nn,nn,uu,vv,g,g+a,b,c]); d8=flow8(0,y8,D,Om,V)
    tr=d8[0]+d8[1]   # d(Tr rho)/dt ; with Tr rho=1 initially
    dq=np.array([d8[1]-nn*tr,d8[2]-uu*tr,d8[3]-vv*tr])
    dX=np.array([d8[4]+g*tr,d8[5]+(g+a)*tr,d8[6]+b*tr,d8[7]+c*tr]); dp=np.array([dX[1]-dX[0],dX[2],dX[3]])
    print('projection check q:',np.abs(dq-F(y,D,Om,V)[:3]).max(),' p:',np.abs(dp-F(y,D,Om,V)[3:]).max())
    for n0 in roots(D,Om,V):
        w=np.linalg.eigvals(JAC(fixed(D,Om,V,n0),D,Om,V)); print('n=%.4f eig'%n0,np.round(np.sort_complex(w),4))
