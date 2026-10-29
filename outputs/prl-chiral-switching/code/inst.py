import numpy as np, json, sys
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from wkbf import flow, lin
from wkbq_pt import mf_point
from mf import roots
def unstable_basis(y,D,Om,V):
    w,R=np.linalg.eig(lin(y,D,Om,V)); B=[]
    for i in np.argsort(-w.real):
        if w[i].real>1e-6:
            v=R[:,i]
            if abs(w[i].imag)>1e-9:
                if w[i].imag>0: B+= [v.real/np.linalg.norm(v.real), v.imag/np.linalg.norm(v.imag)]
            else: B.append(v.real/np.linalg.norm(v.real))
    return np.array(B), w
def dist_all(Y,ys):
    # gauge-invariant distance (rho->lam rho, X->X/lam): normalise by Tr rho
    tr=Y[0]+Y[1]; q=Y[:4]/tr; p=Y[4:]*tr
    return np.linalg.norm(q-ys[:4,None],axis=0)+np.linalg.norm(p-ys[4:,None],axis=0)
def shoot(ang,y0,B,ys,D,Om,V,eps=1e-5,T=45):
    th,ph=ang; c=np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
    y=y0+eps*(c@B)
    best=[1e9,None,None]
    def ev(t,y,*a): return np.abs(y).max()-50
    ev.terminal=True
    S=solve_ivp(flow,[0,T],y,args=(D,Om,V),rtol=1e-9,atol=1e-11,events=ev,max_step=0.1)
    d=dist_all(S.y,ys); k=int(np.argmin(d)); return [float(d[k]),k,S]
def action(S,k):
    q=S.y[:4,:k+1]; p=S.y[4:,:k+1]
    dq=np.diff(q,axis=1); pm=(p[:,1:]+p[:,:-1])/2
    return float(np.sum(pm*dq))   # int Tr(X d rho)
def instanton(D,Om,V,src,dst,ngrid=7):
    r=roots(D,Om,V); y0=mf_point(D,Om,V,r[src]); ys=mf_point(D,Om,V,r[dst])
    B,w=unstable_basis(y0,D,Om,V)
    grid=[(t,p) for t in np.linspace(0.05,np.pi-0.05,ngrid) for p in np.linspace(0,2*np.pi,2*ngrid,endpoint=False)]
    vals=[(shoot(a,y0,B,ys,D,Om,V,T=40)[0],a) for a in grid]
    vals.sort(key=lambda x:x[0]); a0=vals[0][1]
    res=minimize(lambda a: shoot(a,y0,B,ys,D,Om,V)[0],a0,method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-12,maxiter=300))
    d,k,S=shoot(res.x,y0,B,ys,D,Om,V)
    return action(S,k), d, res.x
def job(a):
    D,src=a; V=10.;Om=2.0
    S,d,ang=instanton(D,Om,V,src,1); return [D,src,S,d,list(ang)]
if __name__=="__main__":
    from multiprocessing import Pool
    Ds=[float(x) for x in sys.argv[1:]]
    with Pool(8) as P:
        for r in P.imap_unordered(job,[(D,s) for D in Ds for s in (2,0)]):
            print('D=%.3f %s->saddle  S=%.5f  miss=%.1e'%(r[0],'hi' if r[1]==2 else 'lo',r[2],r[3]),flush=True)
            with open('../data/inst.jsonl','a') as f: f.write(json.dumps(r)+'\n')
