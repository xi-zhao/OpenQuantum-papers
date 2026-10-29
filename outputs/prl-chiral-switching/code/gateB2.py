import numpy as np, json, sys, warnings
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares, brentq
from scipy.spatial import cKDTree
import lin
from wkbr import F, H
from gateB import nb, EP
from mf import roots
warnings.filterwarnings('ignore')
Om=2.0; V=10.
def endpoint(up,n_b,D):
    """transversality p=(pn,0,0) and h=-pn(Om v+n)=0 with pn!=0  ->  v=-n_b/Om; free (u,pn)."""
    u,pn=up; return np.array([n_b,u,-n_b/Om,pn,0,0.])
def back(uv,n_b,D,nsec,T=80.):
    y0=endpoint(uv,n_b,D)
    if y0 is None: return None
    def rhs(t,z): f=F(z[:6],D,Om,V); return -np.r_[f,z[3:6]@f[:3]]
    def sec(t,z): return z[0]-nsec
    sec.terminal=True
    def blow(t,z): return np.abs(z[:6]).max()-15
    blow.terminal=True
    s=solve_ivp(rhs,[0,T],np.r_[y0,0.],rtol=1e-10,atol=1e-12,events=[sec,blow])
    return s.y_events[0][0] if (s.status==1 and len(s.t_events[0])) else None
def run(N):
    D=EP[N][1]; r=roots(D,Om,V); A=lin.fixed(D,Om,V,r[-1]); UA,_=lin.basis(A,D,False); n_b=nb(D); nsec=(n_b+r[-1])/2
    ga=[(a,b) for a in np.linspace(0.05,np.pi-0.05,16) for b in np.linspace(0,2*np.pi,32,endpoint=False)]
    Fw=[(g,lin.ray(g,A,UA,D,nsec,+1)) for g in ga]; Fw=[(g,y) for g,y in Fw if y is not None]
    gb=[(u,pn) for u in np.linspace(-0.3,0.3,25) for pn in np.linspace(-3,3,25) if abs(pn)>1e-3]
    Bw=[(g,back(g,n_b,D,nsec)) for g in gb]; Bw=[(g,y) for g,y in Bw if y is not None]
    if not Fw or not Bw: return dict(N=N,D=D,nb=n_b,ok=False,nf=len(Fw),nbk=len(Bw))
    tr=cKDTree(np.array([y[1:6] for g,y in Bw])); d,idx=tr.query(np.array([y[1:6] for g,y in Fw]))
    def res(x):
        yf=lin.ray(x[:2],A,UA,D,nsec,+1); yb=back(x[2:],n_b,D,nsec)
        return np.full(5,10.) if (yf is None or yb is None) else yf[1:6]-yb[1:6]
    best=None
    for i in np.argsort(d)[:8]:
        R=least_squares(res,np.r_[Fw[i][0],Bw[idx[i]][0]],xtol=1e-14,ftol=1e-14,gtol=1e-14,diff_step=1e-7,max_nfev=400)
        if best is None or R.cost<best.cost: best=R
    yf=lin.ray(best.x[:2],A,UA,D,nsec,+1); yb=back(best.x[2:],n_b,D,nsec)
    return dict(N=N,D=D,nb=n_b,nsec=nsec,ok=True,S=float(yf[6]-yb[6]),miss=float(np.linalg.norm(yf[1:6]-yb[1:6])),
                nf=len(Fw),nbk=len(Bw),seed=float(d.min()),x=best.x.tolist())
if __name__=="__main__":
    o=run(int(sys.argv[1])); print(json.dumps(o),flush=True); open('../data/gateB2.jsonl','a').write(json.dumps(o)+'\n')
