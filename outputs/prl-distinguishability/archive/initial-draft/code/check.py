import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.optimize import brentq
from theory2 import lnT, LN2
h=lambda f: -f*np.log2(f)-(1-f)*np.log2(1-f)
qc=brentq(lambda f: h(f)+f*np.log2(3)-1,0.01,0.5)
D=lambda a,b: a*np.log2(a/b)+(1-a)*np.log2((1-a)/(1-b))
print("qc",qc,"rc",4*qc/3)
for r in [0.3,0.4,0.5,0.6,0.8,0.9,1.0]:
    q=3*r/4
    H2=-np.log2((1-q)**2+q*q/3)
    v=[-lnT(n,r)/LN2/n for n in [200,400,800,1600]]
    # Richardson in 1/n assuming a + b log n / n + c/n
    print(r, "sanov",round(D(qc,q),4),"coll",round((H2-1)/2,4), [round(x,4) for x in v])
