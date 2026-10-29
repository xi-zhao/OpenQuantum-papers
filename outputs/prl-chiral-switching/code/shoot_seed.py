"""Reliable exact action: initial-value shooting from the hi fixed point inside its 3D unstable manifold,
   with the launch direction seeded by the BVP/gMAM solution and refined to hit the saddle's stable manifold.
   Action = int p.dq up to closest approach; reported with miss distance (convergence criterion)."""
import numpy as np, json, sys, warnings
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from wkbr import F, JAC, H, fixed
from bvp_exact import solve as bvp_solve
from mf import roots
warnings.filterwarnings('ignore')
Om=2.0; V=10.
def basis(y,D,stable):
    w,R=np.linalg.eig(JAC(y,D,Om,V)); sel=(w.real<0) if stable else (w.real>0); B=[]
    for i in np.flatnonzero(sel): B+=[R[:,i].real,R[:,i].imag] if abs(w[i].imag)>1e-10 else [R[:,i].real]
    return np.linalg.qr(np.array(B).T)[0][:,:3]
def run(c,A,UA,D,S,WS,T=80,eps=1e-6):
    c=np.asarray(c); y0=A+eps*UA@(c/np.linalg.norm(c))
    def ev(t,y,*a): return np.abs(y).max()-20
    ev.terminal=True
    sol=solve_ivp(lambda t,y: F(y,D,Om,V),[0,T],y0,rtol=1e-11,atol=1e-13,events=ev,dense_output=True)
    Y=sol.y; d=np.linalg.norm(Y-S[:,None],axis=0); k=int(np.argmin(d))
    act=float(np.sum((Y[3:,1:k+1]+Y[3:,:k])/2*np.diff(Y[:3,:k+1],axis=1)))
    return d[k],act,k,sol
if __name__=="__main__":
    D=float(sys.argv[1]); r=roots(D,Om,V); A=fixed(D,Om,V,r[2]); S=fixed(D,Om,V,r[1])
    UA=basis(A,D,False)
    R=bvp_solve(D,2,40.); # seed: take the BVP trajectory point at small distance from A, project on UA
    from bvp_exact import seed
    P=R['P']; t,Y=seed(D,P,40.)
    best=None
    for frac in (0.02,0.05,0.1):
        k=np.argmin(abs(np.linalg.norm(Y-A[:,None],axis=0)-frac)); c0=UA.T@(Y[:,k]-A)
        res=least_squares(lambda c: [run(c,A,UA,D,S,None)[0]],c0,diff_step=1e-6,max_nfev=300)
        d,act,kk,sol=run(res.x,A,UA,D,S,None)
        if best is None or d<best[0]: best=(d,act)
        print('D=%.3f seed frac %.2f -> miss=%.2e  S=%.5f'%(D,frac,d,act),flush=True)
    json.dump(dict(D=D,miss=best[0],S=best[1]),open(f'../data/shoot_{D}.json','w'))
