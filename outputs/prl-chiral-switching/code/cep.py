import numpy as np, scipy.sparse.linalg as sla, sys
from ex import liou
# locate genuine EP of the exact (linear) Liouvillian at complex Delta, Om fixed
V=10.; Om=2.0
def pair(N,D,target):
    w=sla.eigs(liou(N,D,Om,V),k=4,sigma=target,return_eigenvectors=False,tol=1e-13)
    w=w[np.argsort(abs(w-target))][:2]; return w
def f(N,D,target):
    a,b=pair(N,D,target); return (a-b)**2,(a+b)/2
def find(N,D0,target):
    D1=D0+1e-3j; f0,c=f(N,D0,target); f1,c=f(N,D1,target)
    for it in range(40):
        D2=D1-f1*(D1-D0)/(f1-f0)
        D0,f0=D1,f1; D1=D2; f1,c=f(N,D1,target); target=c
        if abs(D1-D0)<1e-11: break
    return D1,c,abs(f1)
if __name__=="__main__":
    for N in [int(x) for x in sys.argv[1:]]:
        # start at real avoided crossing: scan for min gap near left fold
        Ds=np.linspace(3.0,3.5,51); g=[]
        for D in Ds:
            w=np.sort(sla.eigs(liou(N,D,Om,V),k=6,sigma=-0.75,return_eigenvectors=False).real)[::-1]
            w=w[w<-1e-6]; g.append(min(np.diff(np.sort(w)))) 
        Dm=Ds[int(np.argmin(g))]
        D,c,res=find(N,Dm+0.01j,-0.76)
        print(N,'avoided at',round(Dm,3),'mingap',round(min(g),5),'EP D=',D,'lambda_EP=',c,'res',res,flush=True)
