import numpy as np, json, sys, scipy.sparse.linalg as sla
from ex import liou
V=10.
def pair(N,D,Om,t):
    w=sla.eigs(liou(N,D,Om,V),k=4,sigma=t,return_eigenvectors=False,tol=1e-13)
    return w[np.argsort(abs(w-t))][:2]
def newton(N,D,Om,t):
    h=1e-4*(1+1j)
    for it in range(40):
        a,b=pair(N,D,Om,t); f0=(a-b)**2; t=(a+b)/2
        a,b=pair(N,D+h,Om,t); f1=(a-b)**2
        dD=-f0*h/(f1-f0)
        if abs(dD)>0.05: dD*=0.05/abs(dD)
        D+=dD
        if abs(dD)<1e-10: break
    a,b=pair(N,D,Om,t); return D,(a+b)/2,abs(a-b)**2
tag=sys.argv[1]; Om=float(sys.argv[2]); D=complex(sys.argv[3]); t=complex(sys.argv[4])
res=[]
for N in [20,24,28,32,36,40]:
    D,t,r=newton(N,D,Om,t); res.append([N,D.real,D.imag,t.real,t.imag,r]); print(tag,N,D,t,r,flush=True)
    json.dump(res,open(f'../data/cep_{tag}.json','w'))
