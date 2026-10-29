"""Pre-registered v3 job. Writes g_N to a separate file WITHOUT printing it; prints only Delta_min and Phi."""
import numpy as np, scipy.sparse.linalg as sla, json, sys, warnings
from scipy.optimize import brentq, minimize_scalar
from ex import liou, obs_vec, ops
from mf import roots
from mfdyn import rhs, steady
warnings.filterwarnings('ignore')
V=10.
Om=float(sys.argv[1]); N=int(sys.argv[2])
Df=brentq(lambda D: (len(roots(D,Om,V))==3)-0.5,3.0,(4.2 if Om<2.5 else 4.6),xtol=1e-12)
rr=roots(Df+1e-7,Om,V); n_f=(rr[0]+rr[1])/2
y=np.array(steady(Df-1e-7,Om,V,0)); f=lambda y: np.array(rhs(0,y,lambda t:(Df-1e-7,Om,V)))
Jm=np.array([(f(y+1e-7*e)-f(y-1e-7*e))/2e-7 for e in np.eye(3)]).T; ev=np.linalg.eigvals(Jm); kappa=-ev[np.argmin(abs(ev.imag))].real
def pair(D):
    w=sla.eigs(liou(N,D,Om,V),k=8,sigma=-kappa,return_eigenvectors=False,tol=1e-12)
    w=w[np.abs(w.imag)<1e-6].real; w=w[np.argsort(abs(w+kappa))][:2]
    return np.sort(w) if len(w)==2 else None
def gap(D):
    p=pair(D); return 9.0 if p is None else p[1]-p[0]
Ds=np.linspace(Df-1.0,Df-0.02,50); G=np.array([gap(D) for D in Ds]); i=int(np.argmin(G))
lo,hi=Ds[max(i-1,0)],Ds[min(i+1,len(Ds)-1)]
res=minimize_scalar(gap,bounds=(lo,hi),method='bounded',options=dict(xatol=1e-10))
Dmin=float(res.x); g=float(res.fun/2)
def Pk(D):
    L=liou(N,D,Om,V); tr,nv=obs_vec(N); Ne,Jx,B,idx,Dd=ops(N)
    w,v=sla.eigs(L,k=1,sigma=1e-9,tol=1e-15); r=v[:,0]/(tr@v[:,0])
    d=Ne.shape[0]; full=np.zeros(d*d,complex); full[idx]=r; R=full.reshape(d,d,order='F')
    ne=np.round(np.diag(Ne).real).astype(int); P=np.zeros(N+1)
    for j in range(d): P[ne[j]]+=R[j,j].real
    return P
P=Pk(Dmin); lp=np.log(np.clip(P,1e-300,None)); k=np.arange(N+1)/N
Phi=float(-(np.interp(n_f,k,lp)-lp.max())/N)
rec=dict(Om=Om,N=N,Df=Df,n_f=n_f,kappa=kappa,Dmin=Dmin,Phi=Phi,edge=bool(i in (0,len(Ds)-1)),scan_min=float(G.min()))
open("../data/v4_dmin.jsonl","a").write(json.dumps(dict(Om=Om,N=N,Dmin=Dmin,n_f=n_f,Df=Df,edge=rec["edge"]))+"\n")
open("../data/v4_SEALED.jsonl","a").write(json.dumps(dict(Om=Om,N=N,g=g,PhiExact=Phi))+"\n")
print('Om=%.2f N=%d Dmin=%.5f edge=%s'%(Om,N,Dmin,rec['edge']),flush=True)
