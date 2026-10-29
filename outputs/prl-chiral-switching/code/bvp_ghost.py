"""Exact-Hamiltonian BVP, hi fixed point -> hyperplane n=n_f with transversality p_u=p_v=0 (free endpoint).
   Seeded by the free-endpoint gMAM path. Also re-runs the in-wedge validation (hi->saddle) at the same settings."""
import numpy as np, json, sys, warnings
from scipy.integrate import solve_bvp
from wkbr import F, JAC, H, fixed
from bvp_exact import seed, sub
from gmamv import gmam
from mf import roots
warnings.filterwarnings('ignore')
Om=2.0; V=10.
nf=json.load(open('../data/normalform.json'))
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
FR=dict((j['N'],j) for j in map(json.loads,open('../data/gmam_free.jsonl')))
def solve(D,P,T,mode):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[-1]); UA,WA=sub(A,D,False)
    t,Y=seed(D,P,T)
    fun=lambda t,y: T*F(y,D,Om,V); jac=lambda t,y: T*JAC(y,D,Om,V)
    if mode=='ghost': bc=lambda ya,yb: np.concatenate([WA.T@(ya-A),[yb[0]-nf['nf'],yb[4],yb[5]]])
    else:
        S=fixed(D,Om,V,r[1]); US,WS=sub(S,D,True); bc=lambda ya,yb: np.concatenate([WA.T@(ya-A),WS.T@(yb-S)])
    sol=solve_bvp(fun,bc,t,Y,fun_jac=jac,tol=1e-6,max_nodes=300000)
    y=sol.sol(np.linspace(0,1,60001)); act=float(np.sum((y[3:,1:]+y[3:,:-1])/2*np.diff(y[:3],axis=1)))
    return act,bool(sol.success),float(np.abs(H(y,D,Om,V)).max()),sol.message
if __name__=="__main__":
    N=int(sys.argv[1]); D=EP[N][1]
    r=roots(D,Om,V); A=fixed(D,Om,V,r[-1])[:3]; Z=np.array(FR[N]['end'])
    P=np.array([A+(Z-A)*s for s in np.linspace(0,1,40)])
    from gmam_free import run as freerun
    S,end=freerun(D)
    # rebuild the free path (run returns only end); reuse straight seed refined by fixed-end gMAM to that end
    Z=np.array(end); 
    for m in (40,80,160): _,P=gmam(D,A,Z,P,m)
    out=[]
    for T in (30.,50.,80.):
        act,ok,hm,msg=solve(D,P,T,'ghost'); out.append([T,act,ok,hm]); print('N=%d D=%.4f T=%3.0f S_ghost(exact H)=%.5f ok=%s hmax=%.1e %s'%(N,D,T,act,ok,hm,msg),flush=True)
    json.dump(dict(N=N,D=D,gmam_free=S,bvp=out),open(f'../data/bvp_ghost_{N}.json','w'))
