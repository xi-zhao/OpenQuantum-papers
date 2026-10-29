import numpy as np, scipy.sparse.linalg as sla, json, sys, warnings
from ex import liou, obs_vec, ops
warnings.filterwarnings('ignore')
V=10.
P={(o['Om'],o['N']):o for o in map(json.loads,open('../data/v3_phi.jsonl'))}
G={(o['Om'],o['N']):o['g'] for o in map(json.loads,open('../data/v3_g_SEALED.jsonl'))}
def Pk(N,D,Om):
    L=liou(N,D,Om,V); tr,nv=obs_vec(N); Ne,Jx,B,idx,Dd=ops(N)
    w,v=sla.eigs(L,k=1,sigma=1e-9,tol=1e-15); r=v[:,0]/(tr@v[:,0])
    d=Ne.shape[0]; full=np.zeros(d*d,complex); full[idx]=r; R=full.reshape(d,d,order='F')
    ne=np.round(np.diag(Ne).real).astype(int); out=np.zeros(N+1)
    for j in range(d): out[ne[j]]+=R[j,j].real
    return out
Om=float(sys.argv[1]); Ns=[24,32,40,48]
dist={N:Pk(N,P[(Om,N)]['Dmin'],Om) for N in Ns}
nf=P[(Om,24)]['n_f']; Y=np.log([G[(Om,N)] for N in Ns])
print('Om=%.1f n_f=%.4f'%(Om,nf))
for n in np.round(np.linspace(0.06,0.34,15),3):
    X=[]
    for N in Ns:
        p=dist[N]; lp=np.log(np.clip(p,1e-300,None)); k=np.arange(N+1)/N; X.append(-(np.interp(n,k,lp)-lp.max()))
    X=np.array(X); sl=-np.polyfit(X,Y,1)[0]; loc=(X[2]-X[1])/8; obs=-(Y[2]-Y[1])/8
    print('  n=%.3f  slope=%.3f  local pred/obs=%.3f %s'%(n,sl,loc/obs,'<- n_f' if abs(n-nf)<0.011 else ''))
