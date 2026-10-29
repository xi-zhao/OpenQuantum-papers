import numpy as np
from math import comb
def Nmat(n,d,z):
    # coefficient-space generator on span{P_0..P_n}, shifted by +n gamma (gamma=1)
    N=np.zeros((n+1,n+1))
    for k in range(1,n+1): N[k-1,k]=n-k+1
    for k in range(0,n+1-d): N[k+d,k]=z*comb(k+d,d)
    return N
def semiclassical(n,d,z):
    from math import factorial
    a=d/(d+1)
    C=(d+1)*d**(-a)*factorial(d)**(-1/(d+1))*z**(1/(d+1))
    Hm=C*(n*n/4)**a
    return Hm, Hm*np.sqrt(8*d*d/(d+1))/n
