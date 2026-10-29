"""v4 solver (backward-to-A): endpoint family on the occupation constraint x=(1-n,n,0,0), integrate backward in time,
   require the ray to land on W^u(A) (codim 2 in h=0) -> minimise closest-approach distance to A over 2 endpoint params.
   Forward Doi action from A to the endpoint = -(backward-accumulated int p.dq) at closest approach (+O(d^2))."""
import numpy as np, json, sys, warnings
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import lin
from wkbr import F, H
from mf import roots
from phiwkb import endA, endB, endpoint_term, V
warnings.filterwarnings('ignore')
def back_to_A(y0,A,D,Om,T=60.):
    if y0 is None: return None
    def rhs(t,z): f=F(z[:6],D,Om,V); return -np.r_[f,z[3:6]@f[:3]]
    def blow(t,z): return np.abs(z[:6]).max()-25
    blow.terminal=True
    def phys(t,z): return min(z[0]+0.02,1.02-z[0])
    phys.terminal=True
    s=solve_ivp(rhs,[0,T],np.r_[y0,0.],rtol=1e-10,atol=1e-12,events=[blow,phys],max_step=0.5,dense_output=False)
    d=np.linalg.norm(s.y[:6]-A[:,None],axis=0); k=int(np.argmin(d))
    return d[k], -s.y[6,k], k
def phi(D,Om,n):
    r=roots(D,Om,V); A=lin.fixed(D,Om,V,r[-1])
    sols=[]
    for br in ('A','B'):
        mk=(lambda z: endA(z,n,Om)) if br=='A' else (lambda z: endB(z,n,Om,D))
        grid=[(nr,u) for nr in np.linspace(0.01,0.49,30) for u in np.linspace(-0.45,0.45,30)] if br=='A' else \
             [(nr,pv) for nr in np.linspace(0.01,0.49,30) for pv in np.linspace(-4,4,30)]
        vals=[]
        for g in grid:
            o=back_to_A(mk(g),A,D,Om)
            if o is not None: vals.append((o[0],g))
        vals.sort(key=lambda x:x[0]); print(' br=%s  valid=%d  best grid d=%.3e'%(br,len(vals),vals[0][0] if vals else -1),flush=True)
        for d0,g in vals[:4]:
            obj=lambda z: (lambda o: 10. if o is None else np.log10(o[0]+1e-16))(back_to_A(mk(z),A,D,Om))
            R=minimize(obj,np.array(g),method='Nelder-Mead',options=dict(xatol=1e-12,fatol=1e-6,maxiter=600))
            e=mk(R.x); o=back_to_A(e,A,D,Om)
            if o is None: continue
            et=endpoint_term(e,n)
            print('   seed d=%.2e -> d=%.2e  S_doi=%.5f  end=%s  Phi=%s'%(d0,o[0],o[1],None if et is None else '%.5f'%et,None if et is None else '%.5f'%(o[1]+et)),flush=True)
            if o[0]<1e-5 and et is not None: sols.append(dict(br=br,d=o[0],S_doi=o[1],endterm=et,Phi=o[1]+et,z=R.x.tolist()))
    pos=[s for s in sols if s['Phi']>0]
    return dict(ok=bool(pos),Phi=min(s['Phi'] for s in pos) if pos else None,sols=sols)
if __name__=="__main__":
    D,Om,n=float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]); tag=sys.argv[4] if len(sys.argv)>4 else 'x'
    o=phi(D,Om,n); o.update(D=D,Om=Om,n=n,tag=tag)
    print(json.dumps({k:o[k] for k in ('D','Om','n','ok','Phi','tag')}),flush=True); open('../data/phiwkb2.jsonl','a').write(json.dumps(o)+'\n')
