"""POST-HOC check (not pre-registered): gMAM from hi to the hyperplane n = n_f (u,v free)."""
import numpy as np, json, sys, warnings
from scipy.optimize import minimize
from gmamv import S_path, reparam, fixed, roots, Om, V
warnings.filterwarnings('ignore')
nf=json.load(open('../data/normalform.json'))
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
def Sfree(z,D,A,m):
    x=z[:-2]; Z=np.array([nf['nf'],z[-2],z[-1]]); v,g=S_path(x,D,A,Z)
    h=1e-7; gz=[]
    for k in (1,2):
        e=np.zeros(3); e[k]=h; gz.append((S_path(x,D,A,Z+e,grad=False)-S_path(x,D,A,Z-e,grad=False))/(2*h))
    return v,np.r_[g,gz]
def run(D):
    r=roots(D,Om,V); A=fixed(D,Om,V,r[-1])[:3]; Z0=fixed(nf['Df'],Om,V,nf['nf'])[:3]
    P=np.array([A+(Z0-A)*s for s in np.linspace(0,1,40)]); out=[]
    for m in (40,80,160):
        P=reparam(P,m)
        for it in range(6):
            z=np.r_[P[1:-1].ravel(),P[-1,1:]]
            res=minimize(Sfree,z,args=(D,A,m),jac=True,method='L-BFGS-B',options=dict(maxiter=3000))
            Z=np.array([nf['nf'],res.x[-2],res.x[-1]]); P=reparam(np.vstack([A,res.x[:-2].reshape(-1,3),Z]),m)
        out.append(float(res.fun))
    return out,P[-1].tolist()
if __name__=="__main__":
    N=int(sys.argv[1]); D=EP[N][1]; S,end=run(D)
    r=dict(N=N,D=D,S=S,end=end); print(json.dumps(r),flush=True); open('../data/gmam_free.jsonl','a').write(json.dumps(r)+'\n')
