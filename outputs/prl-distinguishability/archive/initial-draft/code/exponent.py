import numpy as np
from scipy.optimize import brentq
L=np.log2
h=lambda f: -f*L(f)-(1-f)*L(1-f)
H=lambda f: h(f)+f*L(3)
qc=brentq(lambda f:H(f)-1,1e-3,.5)
D=lambda a,b: a*L(a/b)+(1-a)*L((1-a)/(1-b))
def E(r):
    q=3*r/4
    if q<=qc: return 0.0
    fs=(q*q/3)/((1-q)**2+q*q/3)
    return D(qc,q) if fs<qc else 0.5*(-L((1-q)**2+q*q/3)-1)
rs=brentq(lambda r:(lambda q:(q*q/3)/((1-q)**2+q*q/3)-qc)(3*r/4),.3,.9)
print("r_c",4*qc/3,"r_*",rs)
# quadratic coefficient near threshold
for dr in [1e-2,1e-3]: print(dr,E(4*qc/3+dr)/dr**2)
np.savetxt("exponent.dat",[[r,E(r)] for r in np.linspace(0,1,201)])
