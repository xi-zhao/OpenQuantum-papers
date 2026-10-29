import numpy as np, json, sys
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from wkbf import flow, lin
from wkbq_pt import mf_point
from mf import roots
V=10.;Om=2.0
def ubasis(y,D):
    w,R=np.linalg.eig(lin(y,D,Om,V)); B=[]
    for i in np.argsort(-w.real):
        if w[i].real>1e-6:
            v=R[:,i]
            if abs(w[i].imag)>1e-9:
                if w[i].imag>0: B+=[v.real/np.linalg.norm(v.real),v.imag/np.linalg.norm(v.imag)]
            else: B.append(v.real/np.linalg.norm(v.real))
    return np.array(B)
def ev(t,y,*a): return np.abs(y).max()-30
ev.terminal=True
def ray(ang,y0,B,D,T=45,eps=1e-5):
    th,ph=ang; c=np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
    return solve_ivp(flow,[0,T],y0+eps*(c@B),args=(D,Om,V),rtol=1e-9,atol=1e-11,events=ev,max_step=0.1)
def qn(Y): tr=Y[0]+Y[1]; return Y[:4]/tr
def miss(S,qt):
    d=np.linalg.norm(qn(S.y)-qt[:,None],axis=0); k=int(np.argmin(d)); return d[k],k
def act(S,k):
    q=S.y[:4,:k+1]; p=S.y[4:,:k+1]; return float(np.sum((p[:,1:]+p[:,:-1])/2*np.diff(q,axis=1)))
def to_point(D,src_n,qt):
    y0=mf_point(D,Om,V,src_n); B=ubasis(y0,D)
    grid=[(t,p) for t in np.linspace(0.15,np.pi-0.15,6) for p in np.linspace(0,2*np.pi,12,endpoint=False)]
    vals=sorted([(miss(ray(a,y0,B,D),qt)[0],a) for a in grid],key=lambda x:x[0])
    best=None
    for d0,a0 in vals[:3]:
        r=minimize(lambda a: miss(ray(a,y0,B,D),qt)[0],a0,method='Nelder-Mead',options=dict(xatol=1e-8,fatol=1e-11,maxiter=250))
        if best is None or r.fun<best.fun: best=r
    S=ray(best.x,y0,B,D); d,k=miss(S,qt); return act(S,k),d
def job(a):
    kind,D=a; r=roots(D,Om,V); hi=r[-1]
    if kind=='saddle': qt=mf_point(D,Om,V,r[1])[:4]
    else:  # frozen fold point rho_f (MF state at the left fold)
        nf=json.load(open('../data/normalform.json')); qt=mf_point(nf['Df'],Om,V,nf['nf'])[:4]
    S,d=to_point(D,hi,qt); return [kind,D,S,d]
if __name__=="__main__":
    from multiprocessing import Pool
    jobs=[('saddle',D) for D in (3.62,3.7,3.9)]+[('fold',D) for D in (3.1,3.2,3.3,3.45,3.6)]
    with Pool(8) as P:
        for r in P.imap_unordered(job,jobs):
            print('%s D=%.3f  S=%.5f  miss=%.1e'%tuple(r),flush=True)
            with open('../data/inst2.jsonl','a') as f: f.write(json.dumps(r)+'\n')
