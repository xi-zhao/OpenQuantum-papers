# Semiclassics for the P-chain. Right eigvec ansatz c_k ~ exp(n S(x)), x=k/n.
# Generator acting on coefficient vector (N v)_k = (n-k) v_{k+1} + z C(k,d) v_{k-d}
# With v_{k+1}/v_k = e^{p}:  mu = n[(1-x) e^{p} + zeta x^d/d! e^{-dp}] , zeta = z n^{d-1}
# Max over p at fixed x then max over x gives mu/n = max_x min_p h  (saddle). Compare to Perron roots.
import numpy as np
from scipy.optimize import minimize_scalar
from math import factorial
def wkb(zeta,d):
    def hmin(x):
        a=(1-x); b=zeta*x**d/factorial(d)
        # min over p of a e^p + b e^{-dp}: e^{(d+1)p}=d b/a
        return (d+1)*a**(d/(d+1))*(b/d)**(1/(d+1))
    r=minimize_scalar(lambda x:-hmin(x),bounds=(1e-12,1-1e-12),method='bounded')
    return -r.fun, r.x
for d in (2,3):
    for zeta in (1e-4,1e-2,1,100):
        print(d,zeta,'WKB mu/n=%.4f  x*=%.3f'%wkb(zeta,d))
# closed form: maximize (1-x)^{d/(d+1)} x^{d/(d+1)} -> x*=1/2
for d in (1,2,3,4):
    C=(d+1)*(0.25)**(d/(d+1))*(1/(d*factorial(d)))**(1/(d+1))
    print('d=%d  mu/n -> C_d zeta^{1/(d+1)}, C_d=%.5f'%(d,C))
