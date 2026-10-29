"""A-priori WKB quasipotential Phi(n) for the excitation-number distribution (pre-registration v4).
   Branch A endpoint (pu=pv=0): occupation constraint x_ee = n_rho(1+(1-n_rho)pn) = n, energy h=-pn(Om v + n_rho)=0 -> v=-n_rho/Om.
   Free endpoint parameters (n_rho,u) ; pn fixed by occupation.  Phi = S_Doi + (1-n)ln X_gg + n ln X_ee.
   Branch B endpoint (u=v=0): occupation x_ee=n_rho*(X_ee) with X_ee=1-n_rho pn+pn ... solved numerically with h=0.
   Lin matching: forward rays in W^u(hi) to section n_sec, backward rays from the endpoint family; 4 unknowns, 5 residuals."""
import numpy as np, json, sys, warnings
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
import lin
from wkbr import F, H
from mf import roots
warnings.filterwarnings('ignore')
V=10.
def endA(z,n,Om):
    nr,u=z
    if not (1e-4<nr<0.5): return None
    pn=(n/nr-1)/(1-nr); return np.array([nr,u,-nr/Om,pn,0.,0.])
def endB(z,n,Om,D):
    nr,pv=z           # u=v=0 ; unknown pn,pu from occupation + energy: x_ee=(1-nr*pn+pn)*nr=n -> pn ; h=0 -> pu
    if not (1e-4<nr<0.5): return None
    pn=(n/nr-1)/(1-nr)
    f=lambda pu: H(np.array([nr,0.,0.,pn,pu,pv]),D,Om,V)
    a=f(0.); b=(f(1.)-f(-1.))/2; c=(f(1.)+f(-1.))/2-a
    if abs(c)<1e-14:
        if abs(b)<1e-14: return None
        pu=-a/b
    else:
        d=b*b-4*c*a
        if d<0: return None
        pu=(-b+np.sqrt(d))/(2*c)
    return np.array([nr,0.,0.,pn,pu,pv])
def back(y0,D,Om,nsec,T=80.):
    if y0 is None: return None
    def rhs(t,z): f=F(z[:6],D,Om,V); return -np.r_[f,z[3:6]@f[:3]]
    def sec(t,z): return z[0]-nsec
    sec.terminal=True
    def blow(t,z): return np.abs(z[:6]).max()-25
    blow.terminal=True
    def phys(t,z): return min(z[0]+0.02,1.02-z[0])
    phys.terminal=True
    s=solve_ivp(rhs,[0,T],np.r_[y0,0.],rtol=1e-10,atol=1e-12,events=[sec,blow,phys],max_step=0.5)
    return s.y_events[0][0] if (s.status==1 and len(s.t_events[0])) else None
def fwd(a,A,UA,D,Om,nsec,eps=1e-6,T=80.):
    y0=np.r_[A+eps*UA@lin.sph(a),0.0]
    def rhs(t,z): f=F(z[:6],D,Om,V); return np.r_[f,z[3:6]@f[:3]]
    def sec(t,z): return z[0]-nsec
    sec.terminal=True
    def blow(t,z): return np.abs(z[:6]).max()-25
    blow.terminal=True
    def phys(t,z): return min(z[0]+0.02,1.02-z[0])
    phys.terminal=True
    s=solve_ivp(rhs,[0,T],y0,rtol=1e-10,atol=1e-12,events=[sec,blow,phys],max_step=0.5)
    return s.y_events[0][0] if (s.status==1 and len(s.t_events[0])) else None
def endpoint_term(y,n):
    nr,u,v,pn,pu,pv=y; Xgg=1-nr*pn-u*pu-v*pv; Xee=Xgg+pn
    if Xgg<=0 or Xee<=0: return None
    return (1-n)*np.log(Xgg)+n*np.log(Xee)
def phi(D,Om,n,ngrid=14,verbose=False):
    r=roots(D,Om,V); A=lin.fixed(D,Om,V,r[-1]); UA,_=lin.basis(A,D,False); nhi=r[-1]
    if n>=nhi: return dict(ok=False,reason='n>=nhi')
    sols=[]
    for frac in (0.5,0.3):
        nsec=n+frac*(nhi-n)
        ga=[(a,b) for a in np.linspace(0.05,np.pi-0.05,ngrid) for b in np.linspace(0,2*np.pi,2*ngrid,endpoint=False)]
        Fw=[(g,fwd(g,A,UA,D,Om,nsec)) for g in ga]; Fw=[(g,y) for g,y in Fw if y is not None]
        if not Fw: continue
        for br in ('A','B'):
            if br=='A': gb=[(nr,u) for nr in np.linspace(0.01,0.49,25) for u in np.linspace(-0.45,0.45,25)]; mk=lambda z: endA(z,n,Om)
            else:      gb=[(nr,pv) for nr in np.linspace(0.01,0.49,25) for pv in np.linspace(-4,4,25)]; mk=lambda z: endB(z,n,Om,D)
            Bw=[(g,back(mk(g),D,Om,nsec)) for g in gb]; Bw=[(g,y) for g,y in Bw if y is not None]
            if not Bw: continue
            tr=cKDTree(np.array([y[1:6] for g,y in Bw])); d,idx=tr.query(np.array([y[1:6] for g,y in Fw]))
            def res(x):
                yf=fwd(x[:2],A,UA,D,Om,nsec); yb=back(mk(x[2:]),D,Om,nsec)
                return np.full(5,10.) if (yf is None or yb is None) else yf[1:6]-yb[1:6]
            print(' frac=%.1f br=%s nF=%d nB=%d seedgap=%.3f'%(frac,br,len(Fw),len(Bw),d.min()),flush=True)
            for i in np.argsort(d)[:3]:
                R=least_squares(res,np.r_[Fw[i][0],Bw[idx[i]][0]],xtol=1e-14,ftol=1e-14,gtol=1e-14,diff_step=1e-7,max_nfev=150)
                yf=fwd(R.x[:2],A,UA,D,Om,nsec); e=mk(R.x[2:]); yb=back(e,D,Om,nsec)
                if yf is None or yb is None: continue
                miss=float(np.linalg.norm(yf[1:6]-yb[1:6]))
                if miss>1e-6: continue
                et=endpoint_term(e,n)
                if et is None: continue
                Sd=float(yf[6]-yb[6]); sols.append(dict(br=br,frac=frac,S_doi=Sd,endterm=float(et),Phi=Sd+float(et),miss=miss,end=e.tolist()))
                print('  sol br=%s frac=%.1f S_doi=%.5f end=%.5f Phi=%.5f miss=%.1e'%(br,frac,Sd,et,Sd+et,miss),flush=True)
    pos=[s for s in sols if s['Phi']>0]
    if not pos: return dict(ok=False,reason='no converged positive solution',nsol=len(sols),sols=sols)
    best=min(pos,key=lambda s:s['Phi'])
    return dict(ok=True,Phi=best['Phi'],best=best,nsol=len(sols),allPhi=sorted(round(s['Phi'],5) for s in sols))
if __name__=="__main__":
    D,Om,n=float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]); tag=sys.argv[4] if len(sys.argv)>4 else 'x'
    o=phi(D,Om,n); o.update(D=D,Om=Om,n=n,tag=tag); o.pop('sols',None)
    print(json.dumps(o),flush=True); open('../data/phiwkb.jsonl','a').write(json.dumps(o)+'\n')
