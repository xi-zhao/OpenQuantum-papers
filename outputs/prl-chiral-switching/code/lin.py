"""Lin's method for the hi -> saddle heteroclinic of the reduced (gauge-fixed) WKB Hamiltonian h(q,p).
   forward rays: y(0)=A+eps*U_A c(a), c on S^2 in A's 3D unstable subspace; backward rays from saddle's 3D stable subspace.
   Both stop at the section n=n_sec. Unknowns (a1,a2,b1,b2); residual = mismatch of (u,v,pn,pu,pv) at the section.
   Action S = int p.dq integrated as an extra ODE component (mesh-independent)."""
import numpy as np, json, sys, time, warnings
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from wkbr import F, JAC, H, fixed
from mf import roots
warnings.filterwarnings('ignore')
Om=2.0; V=10.
def basis(y,D,stable):
    w,R=np.linalg.eig(JAC(y,D,Om,V)); sel=(w.real<0) if stable else (w.real>0); B=[]
    for i in np.flatnonzero(sel):
        if abs(w[i].imag)>1e-10:
            if w[i].imag>0: B+=[R[:,i].real,R[:,i].imag]
        else: B.append(R[:,i].real)
    U=np.linalg.qr(np.array(B).T)[0]; assert U.shape[1]==3
    return U, np.sort(w[sel].real)
def sph(a): return np.array([np.cos(a[0]),np.sin(a[0])*np.cos(a[1]),np.sin(a[0])*np.sin(a[1])])
def ray(a,Y0,U,D,nsec,sgn,eps=1e-6,T=80.):
    y0=np.r_[Y0+eps*U@sph(a),0.0]
    def rhs(t,z):
        f=F(z[:6],D,Om,V); return sgn*np.r_[f,z[3:6]@f[:3]]
    def sec(t,z): return z[0]-nsec
    sec.terminal=True
    def blow(t,z): return np.abs(z[:6]).max()-15
    blow.terminal=True
    s=solve_ivp(rhs,[0,T],y0,rtol=1e-10,atol=1e-12,events=[sec,blow])
    if s.status==1 and len(s.t_events[0]): return s.y_events[0][0]
    return None
def mism(x,A,UA,S,US,D,nsec):
    yf=ray(x[:2],A,UA,D,nsec,+1); yb=ray(x[2:],S,US,D,nsec,-1)
    if yf is None or yb is None: return np.full(5,10.)
    return (yf[1:6]-yb[1:6])
def hetero(D,ngrid=16,nsec_frac=0.5,verbose=True):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[2]); S=fixed(D,Om,V,r[1])
    UA,wA=basis(A,D,False); US,wS=basis(S,D,True); nsec=r[1]+nsec_frac*(r[2]-r[1])
    grid=[(a,b) for a in np.linspace(0.05,np.pi-0.05,ngrid) for b in np.linspace(0,2*np.pi,2*ngrid,endpoint=False)]
    Fw=[(g,ray(g,A,UA,D,nsec,+1)) for g in grid]; Fw=[(g,y) for g,y in Fw if y is not None]
    Bw=[(g,ray(g,S,US,D,nsec,-1)) for g in grid]; Bw=[(g,y) for g,y in Bw if y is not None]
    if not Fw or not Bw: return dict(D=D,ok=False,reason='no section hits',nf=len(Fw),nb=len(Bw))
    tree=cKDTree(np.array([y[1:6] for g,y in Bw])); d,idx=tree.query(np.array([y[1:6] for g,y in Fw]))
    order=np.argsort(d)[:6]; best=None
    for i in order:
        x0=np.r_[Fw[i][0],Bw[idx[i]][0]]
        res=least_squares(mism,x0,args=(A,UA,S,US,D,nsec),xtol=1e-14,ftol=1e-14,gtol=1e-14,diff_step=1e-7,max_nfev=400)
        if best is None or res.cost<best.cost: best=res
    yf=ray(best.x[:2],A,UA,D,nsec,+1); yb=ray(best.x[2:],S,US,D,nsec,-1)
    Sact=yf[6]-yb[6]      # backward ray accumulates int p.dq along reversed time: forward action = S_fwd - S_bwd
    out=dict(D=D,ok=True,S=float(Sact),miss=float(np.linalg.norm(yf[1:6]-yb[1:6])),h_sec=float(H(yf[:6],D,Om,V)),
             nf=len(Fw),nb=len(Bw),seedgap=float(d[order[0]]),x=best.x.tolist(),nsec=nsec)
    if verbose: print(json.dumps(out),flush=True)
    return out
if __name__=="__main__":
    D=float(sys.argv[1]); frac=float(sys.argv[2]) if len(sys.argv)>2 else 0.5
    t=time.time(); o=hetero(D,nsec_frac=frac); o['time']=time.time()-t
    open('../data/lin.jsonl','a').write(json.dumps(o)+'\n')
