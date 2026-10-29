"""Vectorised gMAM (Heymann–Vanden-Eijnden) for the coherent LNA Freidlin–Wentzell action.
   dq = f(q) dt + sqrt(B(q)/N) dW,  q=(n,u,v)=(rho_ee, Re rho_eg, Im rho_eg)
   f: mean-field Bloch drift;  B = [[n,u/2,v/2],[u/2,1/4,0],[v/2,0,1/4]] (single-atom symmetric noise, checked vs lna.D1)
   S[phi] = int ( |phi'|_{B^-1} |f|_{B^-1} - <phi',f>_{B^-1} ) ds ;  rate ~ exp(-N S)."""
import numpy as np, json, sys, warnings
from scipy.optimize import minimize
from wkbr import fixed
from mf import roots
warnings.filterwarnings('ignore')
Om=2.0; V=10.
def drift(Q,D):
    n,u,v=Q.T; Dt=D-V*n
    return np.stack([-Om*v-n, -Dt*v-u/2, Dt*u-Om/2*(1-2*n)-v/2],1)
def Binv(Q):
    n,u,v=Q.T; m=len(n); B=np.zeros((m,3,3)); B[:,0,0]=n; B[:,0,1]=B[:,1,0]=u/2; B[:,0,2]=B[:,2,0]=v/2; B[:,1,1]=B[:,2,2]=0.25
    return np.linalg.inv(B)
def S_path(x,D,A,Z,grad=True,h=1e-7):
    P=np.vstack([A,x.reshape(-1,3),Z]); d=np.diff(P,axis=0); q=(P[1:]+P[:-1])/2
    def seg(d,q):
        Bi=Binv(q); f=drift(q,D)
        a=np.einsum('mi,mij,mj->m',d,Bi,d); b=np.einsum('mi,mij,mj->m',f,Bi,f); c=np.einsum('mi,mij,mj->m',d,Bi,f)
        return np.sqrt(np.clip(a*b,1e-300,None))-c
    L=seg(d,q); val=L.sum()
    if not grad: return val
    # analytic-structure gradient by per-segment finite differences (each interior point touches 2 segments)
    m=len(x)//3; g=np.zeros((m,3))
    for k in range(3):
        e=np.zeros(3); e[k]=h
        Lp_d=seg(d+e,q+e/2); Lm_d=seg(d-e,q-e/2)      # moving the segment's end point
        Lp_s=seg(d-e,q+e/2); Lm_s=seg(d+e,q-e/2)      # moving the segment's start point
        dend=(Lp_d-Lm_d)/(2*h); dstart=(Lp_s-Lm_s)/(2*h)
        g[:,k]=dend[:-1]+dstart[1:]
    return val,g.ravel()
def reparam(P,M):
    s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(P,axis=0),axis=1))]; s/=s[-1]; t=np.linspace(0,1,M)
    return np.array([np.interp(t,s,P[:,k]) for k in range(3)]).T
def gmam(D,A,Z,init,M):
    P=reparam(init,M)
    for k in range(8):
        r=minimize(S_path,P[1:-1].ravel(),args=(D,A,Z),jac=True,method='L-BFGS-B',options=dict(maxiter=3000))
        P=reparam(np.vstack([A,r.x.reshape(-1,3),Z]),M)
    return S_path(P[1:-1].ravel(),D,A,Z,grad=False),P
def barrier(D,src,M=(40,80,160)):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[src])[:3]; Z=fixed(D,Om,V,r[1])[:3]
    P=np.array([fixed(D,Om,V,x)[:3] for x in np.linspace(r[src],r[1],M[0])]); out=[]
    for m in M:
        s,P=gmam(D,A,Z,P,m); out.append(s)
    return out,P
if __name__=="__main__":
    res=[]
    for D in [float(x) for x in sys.argv[1:]]:
        hl,P1=barrier(D,2); lh,P2=barrier(D,0)
        res.append(dict(D=D,hl=hl,lh=lh,Phl=P1.tolist(),Plh=P2.tolist()))
        print('D=%.3f  S_hl(M=40,80,160)=%s   S_lh=%s'%(D,np.round(hl,5),np.round(lh,5)),flush=True)
    json.dump(res,open('../data/gmamv_%s.json'%sys.argv[1],'w'))
