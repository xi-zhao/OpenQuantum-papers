"""Exact WKB action (reduced quadratic-in-p Hamiltonian, indefinite) by BVP, seeded with the gMAM path."""
import numpy as np, json, sys, warnings, sympy as sp
from scipy.integrate import solve_bvp
from scipy.linalg import null_space
from wkbr import F, JAC, H, fixed, hr, ps
from gmamv import gmam, drift
from mf import roots
warnings.filterwarnings('ignore')
Om=2.0; V=10.
n_,u_,v_,D_,O_,V_=sp.symbols('n u v Delta Omega V',real=True)
Gf=sp.lambdify([n_,u_,v_,D_,O_,V_],sp.Matrix([[sp.diff(hr,a,b) for b in ps] for a in ps]),'numpy')
def sub(y,D,stable):
    w,R=np.linalg.eig(JAC(y,D,Om,V)); sel=(w.real<0) if stable else (w.real>0); B=[]
    for i in np.flatnonzero(sel): B+=[R[:,i].real,R[:,i].imag] if abs(w[i].imag)>1e-10 else [R[:,i].real]
    U=np.linalg.qr(np.array(B).T)[0][:,:3]; return U,null_space(U.T)
def seed(D,P,T,m=600,cut=0.02):
    """time-parameterise the gMAM path by the reversed relaxation speed; p from least squares G p = q' - f."""
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(P,axis=0),axis=1))]
    sp_=np.linspace(s[0],s[-1],m); Q=np.array([np.interp(sp_,s,P[:,k]) for k in range(3)]).T
    f=drift(Q,D); sp_speed=np.maximum(np.linalg.norm(f,axis=1),cut)
    t=np.r_[0,np.cumsum(np.diff(sp_)/((sp_speed[1:]+sp_speed[:-1])/2))]; t=t/t[-1]
    dq=np.gradient(Q,t,axis=0)/T
    Pm=np.array([np.linalg.lstsq(np.array(Gf(*q,D,Om,V),float),d-fi,rcond=None)[0] for q,d,fi in zip(Q,dq,f)])
    return t,np.vstack([Q.T,Pm.T])
def solve(D,src,T,P=None,eps=1e-3,tol=1e-6):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[src]); S=fixed(D,Om,V,r[1])
    UA,WA=sub(A,D,False); US,WS=sub(S,D,True)
    if P is None:
        P=np.array([fixed(D,Om,V,x)[:3] for x in np.linspace(r[src],r[1],40)])
        for m in (40,80): _,P=gmam(D,A[:3],S[:3],P,m)
    t,Y=seed(D,P,T)
    fun=lambda t,y: T*F(y,D,Om,V); jac=lambda t,y: T*JAC(y,D,Om,V)
    bc=lambda ya,yb: np.concatenate([WA.T@(ya-A),WS.T@(yb-S)])
    sol=solve_bvp(fun,bc,t,Y,fun_jac=jac,tol=tol,max_nodes=100000)
    y=sol.sol(np.linspace(0,1,40001)); act=float(np.sum((y[3:,1:]+y[3:,:-1])/2*np.diff(y[:3],axis=1)))
    return dict(S=act,ok=bool(sol.success),hmax=float(np.abs(H(y,D,Om,V)).max()),dA=float(np.linalg.norm(y[:,0]-A)),dS=float(np.linalg.norm(y[:,-1]-S)),P=P,msg=sol.message)
if __name__=="__main__":
    D=float(sys.argv[1]); P=None; out=[]
    for T in (15.,25.,40.,60.):
        R=solve(D,2,T,P); P=R['P']; out.append([T,R['S'],R['ok'],R['hmax'],R['dA'],R['dS']])
        print('D=%.3f T=%4.0f  S=%.5f ok=%s hmax=%.1e dA=%.1e dS=%.1e %s'%(D,T,R['S'],R['ok'],R['hmax'],R['dA'],R['dS'],R['msg']),flush=True)
    json.dump(out,open(f'../data/bvp_exact_{D}.json','w'))
