import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.optimize import brentq
from scipy.integrate import quad
# Deterministic-equivalent (free-probability) prediction for T = (1/2)||N(psi1)-N(psi2)||_1
# under local depolarizing N(rho)=(1-r)rho + r I/2 per qubit, psi_i Haar.
def classes(n,r):
    q=3*r/4; k=np.arange(n+1)
    lnN=gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+k*np.log(3)   # multiplicity
    lnp=(n-k)*np.log1p(-q)+k*np.log(q/3)                        # weight
    return lnN-n*np.log(2), lnp                                 # ln a_k, ln p_k
def T(n,r):
    lna,lnp=classes(n,r)
    # work in x=ln mu; S(mu)=sum a p^2/(1+p^2 mu^2)
    def lnS(x): return logsumexp(lna+2*lnp-np.logaddexp(0,2*(lnp+x)))
    def g(x):  # 2 mu^2 S - 1
        return np.log(2)+2*x+lnS(x)
    x0=brentq(g,-200*n-50,200*n+50)   # 2mu0^2 S(mu0)=1
    # T = (d/pi) int_0^{mu0} 2 mu^2 S (1/mu^2 + 2S + 2 mu S') dmu ; S'=-sum a p^4 2mu/(1+p^2mu^2)^2
    # d/dmu[1/mu - 2 mu S] = -(1/mu^2 + 2S + 2 mu S'), so T=(d/pi) int 2mu^2 S |dt/dmu| dmu
    def integrand(x):
        mu=np.exp(x)
        S=np.exp(lnS(x))
        Sp=-2*mu*np.exp(logsumexp(lna+4*lnp-2*np.logaddexp(0,2*(lnp+x))))
        val=2*mu**2*S*(1/mu**2+2*S+2*mu*Sp)*mu   # dmu = mu dx
        return val
    xs=np.linspace(x0-60-2*n,x0,4001)
    v=np.array([integrand(x) for x in xs])
    return 2**n/np.pi*np.trapezoid(v,xs)
if __name__=="__main__":
    import sys
    for n in [4,6,8,10,12]:
        print(n,[round(T(n,r),4) for r in np.linspace(0.1,0.5,9)])
