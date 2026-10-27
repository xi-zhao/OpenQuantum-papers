import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.optimize import brentq
LN2=np.log(2)
def classes(n,r):
    q=3*r/4; k=np.arange(n+1)
    lna=gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+k*np.log(3)-n*LN2
    lnp=(n-k)*np.log1p(-q)+k*np.log(q/3)
    return lna,lnp
def lnT(n,r,N=20001):
    lna,lnp=classes(n,r)
    L=lambda x: np.logaddexp(0,2*(lnp+x))
    lnS=lambda x: logsumexp(lna+2*lnp-L(x))
    x0=brentq(lambda x: LN2+2*x+lnS(x),-5*n-50,5*n+50)
    xs=np.linspace(x0-3*n-80,x0,N)
    out=[]
    for x in xs:
        a=LN2+2*x+lnS(x)                      # ln 2mu^2 S
        # -ds/dmu = mu^-2 + 2 sum a p^2 (1-p^2mu^2)/(1+p^2mu^2)^2
        y=2*(lnp+x)
        lab=np.where(y>0, y+np.log(-np.expm1(-np.abs(y))+1e-300), np.log(-np.expm1(-np.abs(y))+1e-300))
        terms=np.concatenate([[-2*x], LN2+lna+2*lnp-2*L(x)+lab])
        sg=np.concatenate([[1.0],np.where(y>0,-1.0,1.0)])
        b,s=logsumexp(terms,b=sg,return_sign=True)
        out.append(a+b+x)
    out=np.array(out); dx=xs[1]-xs[0]
    return n*LN2-np.log(np.pi)+logsumexp(out)+np.log(dx)
if __name__=="__main__":
    for n in [8,12]:
        print(n,[round(np.exp(lnT(n,r)),4) for r in [0.1,0.25,0.4,0.5]])
    rc=None
    for r in [0.2524,0.30,0.35,0.40,0.45,0.50,0.6]:
        v=[lnT(n,r)/LN2/n for n in [100,200,400,800]]
        print(r,[round(-x,5) for x in v])
