"""Gate B, exactly as pre-registered (data/PREREGISTRATION_v2.md):
   at D = Re Delta_EP(N): minimum action from the hi fixed point to the hyperplane n = n_b(D),
   n_b = argmin |g(n)| of the MF slow-manifold speed g(n)=r(n)(1-2n)-n on n<n_hi, transversality p_u=p_v=0.
   Solver: forward rays in W^u(hi) (2 angles) to the event n=n_b; unknowns 2 angles, residual (p_u,p_v)=0."""
import numpy as np, json, sys, warnings
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, minimize_scalar
import lin
from wkbr import F, H
from mf import roots
warnings.filterwarnings('ignore')
Om=2.0; V=10.
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
def nb(D):
    r=roots(D,Om,V); g=lambda n: abs((Om**2/4)/((D-V*n)**2+0.25)*(1-2*n)-n)
    return minimize_scalar(g,bounds=(0.02,r[-1]-0.05),method='bounded').x
def run(N):
    D=EP[N][1]; r=roots(D,Om,V); A=lin.fixed(D,Om,V,r[-1]); UA,_=lin.basis(A,D,False); n_b=nb(D)
    grid=[(a,b) for a in np.linspace(0.05,np.pi-0.05,16) for b in np.linspace(0,2*np.pi,32,endpoint=False)]
    hits=[(g,lin.ray(g,A,UA,D,n_b,+1)) for g in grid]; hits=[(g,y) for g,y in hits if y is not None]
    if not hits: return dict(N=N,D=D,nb=n_b,ok=False)
    hits.sort(key=lambda gy: np.hypot(gy[1][4],gy[1][5])); best=None
    for g,y in hits[:6]:
        res=least_squares(lambda a: (lambda yy: np.array([10.,10.]) if yy is None else yy[4:6])(lin.ray(a,A,UA,D,n_b,+1)),
                          g,xtol=1e-14,ftol=1e-14,gtol=1e-14,diff_step=1e-7,max_nfev=300)
        y=lin.ray(res.x,A,UA,D,n_b,+1)
        if y is not None and (best is None or np.hypot(y[4],y[5])<best[1]): best=(y[6],float(np.hypot(y[4],y[5])),y)
    return dict(N=N,D=D,nb=n_b,ok=True,S=float(best[0]),resid=best[1],end=best[2][:6].tolist(),nhits=len(hits))
if __name__=="__main__":
    o=run(int(sys.argv[1])); print(json.dumps(o),flush=True); open('../data/gateB.jsonl','a').write(json.dumps(o)+'\n')
