import numpy as np, json, sys, warnings
from scipy.integrate import solve_bvp
from scipy.linalg import null_space
from scipy.optimize import brentq
from wkbr import F, JAC, H, fixed
from mf import roots
from lna import D1
warnings.filterwarnings('ignore')
Om=2.0; V=10.
def sub(y,D,stable):
    w,R=np.linalg.eig(JAC(y,D,Om,V)); sel=(w.real<0) if stable else (w.real>0); B=[]
    for i in np.flatnonzero(sel): B+=[R[:,i].real,R[:,i].imag] if abs(w[i].imag)>1e-10 else [R[:,i].real]
    U=np.linalg.qr(np.array(B).T)[0][:,:3]; return null_space(U.T)
def kramers_guess(D,src,dst,m=300):
    """q: straight MF-fixed-point-like interpolation; p: 1D Kramers p=-2 f(x) l0/Q along slow left eigvec at saddle."""
    r=roots(D,Om,V); A=fixed(D,Om,V,r[src]); S=fixed(D,Om,V,r[dst])
    Js=JAC(S,D,Om,V)[:3,:3]; w,Rv=np.linalg.eig(Js.T); l0=Rv[:,np.argmax(w.real)].real   # left eigvec of unstable MF dir
    Bm,_=D1(D,Om,S[0],S[1]+1j*S[2]); Q=l0@Bm@l0
    t=np.linspace(0,1,m); s=t; Y=np.zeros((6,m))
    for i,si in enumerate(s):
        q=A[:3]+(S[:3]-A[:3])*si; f=F(np.r_[q,0,0,0],D,Om,V)[:3]
        Y[:3,i]=q; Y[3:,i]=-2*(l0@f)/Q*l0/(l0@l0)
    return t,Y
def solve(D,src,dst,T,init=None,tol=1e-6,m=300):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[src]); S=fixed(D,Om,V,r[dst]); WA=sub(A,D,False); WS=sub(S,D,True)
    t,Y=init if init is not None else kramers_guess(D,src,dst,m)
    f=lambda t,y: T*F(y,D,Om,V); fj=lambda t,y: T*JAC(y,D,Om,V)
    bc=lambda ya,yb: np.concatenate([WA.T@(ya-A),WS.T@(yb-S)])
    sol=solve_bvp(f,bc,t,Y,fun_jac=fj,tol=tol,max_nodes=60000)
    y=sol.sol(np.linspace(0,1,40001)); Sact=float(np.sum((y[3:,1:]+y[3:,:-1])/2*np.diff(y[:3],axis=1)))
    return dict(S=Sact,ok=bool(sol.success),hmax=float(np.abs(H(y,D,Om,V)).max()),gapA=float(np.linalg.norm(y[:,0]-A)),
                gapS=float(np.linalg.norm(y[:,-1]-S)),sol=sol,msg=sol.message)
def series(src,Ds,Ts=(20.,40.,60.)):
    out=[]; init=None
    for D in Ds:
        res={}
        for T in Ts:
            R=solve(D,src,1,T,init=init if (init is not None and T==Ts[0]) else (None if T==Ts[0] else (R['sol'].x,R['sol'].y)))
            res[T]=R
        R=res[Ts[-1]]; init=(res[Ts[0]]['sol'].x,res[Ts[0]]['sol'].y)
        out.append([D]+[res[T]['S'] for T in Ts]+[R['ok'],R['hmax'],R['gapA'],R['gapS']])
        print('D=%.4f  S(T=%s)=%s ok=%s hmax=%.1e gapA=%.1e gapS=%.1e'%(D,Ts,np.round([res[T]['S'] for T in Ts],5),R['ok'],R['hmax'],R['gapA'],R['gapS']),flush=True)
    return out
if __name__=="__main__":
    which=sys.argv[1]
    if which=='hi':   # hi -> saddle, start near right fold (4.5750) and continue downward
        out=series(2,[4.56,4.54,4.5,4.45,4.4,4.3,4.2,4.1,4.0,3.9,3.8,3.7,3.65])
    else:             # lo -> saddle, start near left fold (3.5996) and continue upward
        out=series(0,[3.61,3.63,3.66,3.7,3.75,3.8,3.9,4.0,4.1,4.2,4.3])
    json.dump(out,open(f'../data/mam2_{which}.json','w'))
