"""Geometric minimum action (gMAM, Heymann–Vanden-Eijnden) for dq = f(q)dt + sqrt(B(q)/N) dW in q=(n,u,v):
   S[phi] = int ( sqrt(<phi',B^-1 phi'><f,B^-1 f>) - <phi',B^-1 f> ) ds ,  rate ~ exp(-N S).
   f: MF Bloch drift; B: single-atom symmetric-ordered (LNA) noise at the current single-atom state."""
import numpy as np, json, sys, warnings
from scipy.optimize import minimize
from mfdyn import rhs
from lna import D1
from wkbr import fixed
from mf import roots
warnings.filterwarnings('ignore')
Om=2.0; V=10.
def drift(q,D): return np.array(rhs(0,q,lambda t:(D,Om,V)))
def noise(q,D):
    B,_=D1(D,Om,q[0],q[1]+1j*q[2]); B=(B+B.T)/2
    w,R=np.linalg.eigh(B); w=np.clip(w,1e-6,None); return R@np.diag(w)@R.T
def action(X,D,A,S):
    P=np.vstack([A,X.reshape(-1,3),S]); dP=np.diff(P,axis=0); Pm=(P[1:]+P[:-1])/2; tot=0.0
    for d,q in zip(dP,Pm):
        Bi=np.linalg.inv(noise(q,D)); f=drift(q,D)
        tot+=np.sqrt(max(d@Bi@d,0)*max(f@Bi@f,0))-d@Bi@f
    return tot
def reparam(P,M):
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(P,axis=0),axis=1))]; s/=s[-1]; t=np.linspace(0,1,M)
    return np.array([np.interp(t,s,P[:,k]) for k in range(3)]).T
def gmam(D,src,dst,M=40,init=None,rounds=6):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[src])[:3]; S=fixed(D,Om,V,r[dst])[:3]
    if init is None:   # start on the MF slaved manifold (fixed-point curve) between A and S
        P=np.array([fixed(D,Om,V,x)[:3] for x in np.linspace(r[src],r[dst],M)])
    else: P=reparam(init,M)
    for k in range(rounds):
        res=minimize(action,P[1:-1].ravel(),args=(D,A,S),method='L-BFGS-B',options=dict(maxiter=400))
        P=reparam(np.vstack([A,res.x.reshape(-1,3),S]),M)
    Sf=action(P[1:-1].ravel(),D,A,S)
    Pf=reparam(P,2*M); S2=action(minimize(action,Pf[1:-1].ravel(),args=(D,A,S),method='L-BFGS-B',options=dict(maxiter=600)).x,D,A,S)
    return Sf,S2,P
if __name__=="__main__":
    out=[]
    for D in [float(x) for x in sys.argv[1:]]:
        S1,S2,P=gmam(D,2,1); out.append([D,S1,S2]); print('D=%.3f  S_FW(hi->saddle) M=40: %.5f  M=80: %.5f'%(D,S1,S2),flush=True)
    json.dump(out,open('../data/gmam_%s.json'%('_'.join(sys.argv[1:])),'w'))
