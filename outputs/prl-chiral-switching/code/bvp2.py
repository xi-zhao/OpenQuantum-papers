"""Heteroclinic BVP with projection boundary conditions (Beyn): y(0)=A+eps_A*UA c (3 unknown via BC),
   y(1) in stable subspace of S at distance eps_S; T free parameter; phase condition |y(0)-A|=eps_A.
   Action = int p.dq + tail corrections from linearised end pieces (quadratic, O(eps^2), reported)."""
import numpy as np, json, sys, warnings
from scipy.integrate import solve_bvp
from scipy.linalg import null_space
from wkbr import F, JAC, H, fixed
from bvp_exact import seed
from gmamv import gmam
from mf import roots
warnings.filterwarnings('ignore')
Om=2.0; V=10.
def sub(y,D,stable):
    w,R=np.linalg.eig(JAC(y,D,Om,V)); sel=(w.real<0) if stable else (w.real>0); B=[]
    for i in np.flatnonzero(sel): B+=[R[:,i].real,R[:,i].imag] if abs(w[i].imag)>1e-10 else [R[:,i].real]
    U=np.linalg.qr(np.array(B).T)[0][:,:3]; return U,null_space(U.T)
def tailS(eps,y,D):   # action of the linear piece near a fixed point: S_lin = 0.5 * dy^T Omega... approximate by 0 (O(eps^2))
    return 0.0
def hetero(D,src,dst,P,eA=0.02,eS=0.02,T0=20.,tol=1e-6,m=400):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[src]); S=fixed(D,Om,V,r[dst])
    UA,WA=sub(A,D,False); US,WS=sub(S,D,True)
    t,Y=seed(D,P,T0,m)
    # trim seed to the eps-balls
    dA=np.linalg.norm(Y-A[:,None],axis=0); dS=np.linalg.norm(Y-S[:,None],axis=0)
    i0=np.argmax(dA>eA); i1=len(dS)-1-np.argmax(dS[::-1]>eS); Y=Y[:,i0:i1+1]; t=np.linspace(0,1,Y.shape[1])
    fun=lambda t,y,p: p[0]*F(y,D,Om,V)
    fj=lambda t,y,p: (p[0]*JAC(y,D,Om,V),F(y,D,Om,V)[:,None,:])
    def bc(ya,yb,p): return np.concatenate([WA.T@(ya-A),[np.sum((ya-A)**2)-eA**2],WS.T@(yb-S)])
    sol=solve_bvp(fun,bc,t,Y,p=[T0],fun_jac=fj,tol=tol,max_nodes=200000)
    y=sol.sol(np.linspace(0,1,40001)); act=float(np.sum((y[3:,1:]+y[3:,:-1])/2*np.diff(y[:3],axis=1)))
    return dict(S=act,T=float(sol.p[0]),ok=bool(sol.success),hmax=float(np.abs(H(y,D,Om,V)).max()),
                dS=float(np.linalg.norm(y[:,-1]-S)),msg=sol.message,sol=sol)
if __name__=="__main__":
    D=float(sys.argv[1]); src=int(sys.argv[2]) if len(sys.argv)>2 else 2
    r=roots(D,Om,V); A=fixed(D,Om,V,r[src])[:3]; Z=fixed(D,Om,V,r[1])[:3]
    P=np.array([fixed(D,Om,V,x)[:3] for x in np.linspace(r[src],r[1],40)])
    for mm in (40,80): _,P=gmam(D,A,Z,P,mm)
    out=[]
    for e in (0.04,0.02,0.01):
        R=hetero(D,src,1,P,eA=e,eS=e); out.append([e,R['S'],R['T'],R['ok'],R['hmax']])
        print('D=%.3f src=%d eps=%.3f  S=%.5f T=%.1f ok=%s hmax=%.1e dS=%.1e %s'%(D,src,e,R['S'],R['T'],R['ok'],R['hmax'],R['dS'],R['msg']),flush=True)
    json.dump(out,open(f'../data/bvp2_{D}_{src}.json','w'))
