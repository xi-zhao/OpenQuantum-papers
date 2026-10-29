import numpy as np, scipy.sparse.linalg as sla, json, sys
from ex import liou
V=10.; Om=2.0
def pair(N,D,target):
    w=sla.eigs(liou(N,D,Om,V),k=4,sigma=target,return_eigenvectors=False,tol=1e-13)
    w=w[np.argsort(abs(w-target))][:2]; return w
def F(N,D,t):
    a,b=pair(N,D,t); return (a-b)**2,(a+b)/2
def newton(N,D,t):
    # secant on complex D for disc(D)=0, keep target updated
    h=1e-4*(1+1j)
    for it in range(30):
        f0,c=F(N,D,t); t=c
        f1,_=F(N,D+h,t)
        dD=-f0*h/(f1-f0)
        if abs(dD)>0.05: dD*=0.05/abs(dD)
        D=D+dD
        if abs(dD)<1e-10: break
    f0,c=F(N,D,t)
    return D,c,abs(f0)
if __name__=="__main__":
    D,t=3.2591257552971+0.001555639344660722j,-0.7581141031979668+0.002579229367763204j
    res=[]
    for N in list(range(42,52,2)):
        D,t,r=newton(N,D,t); res.append([N,D.real,D.imag,t.real,t.imag,r]); print(N,D,t,r,flush=True)
    D,t=3.0843622141439826+0.04969689225682156j,-0.749282551603581+0.0527139843008112j
    for N in range(18,7,-2):
        D,t,r=newton(N,D,t); res.append([N,D.real,D.imag,t.real,t.imag,r]); print(N,D,t,r,flush=True)
    json.dump(res,open('../data/cep_left.json','w'))
