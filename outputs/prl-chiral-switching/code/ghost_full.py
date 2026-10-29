import numpy as np, json
from scipy.optimize import minimize_scalar, brentq
from mfdyn import rhs
from lna import D1, J
V=10.;Om=2.0
nf=json.load(open('../data/normalform.json')); Df=nf['Df']
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
# beyond the normal form: bottleneck speed of the full 3D MF flow along its slow manifold, outside the wedge
# slow manifold param. by n: coherence slaved (fast, rate 1/2) -> dn/dt = g(n) = r(n)(1-2n)-n,  r=(Om^2/4)/((D-Vn)^2+1/4)
def g(n,D): r=(Om**2/4)/((D-V*n)**2+0.25); return r*(1-2*n)-n
def bottleneck(D):
    res=minimize_scalar(lambda n: g(n,D),bounds=(0.05,0.3),method='bounded'); return res.x,res.fun   # g<0 everywhere outside wedge? sign
for D in [3.0,3.2,3.3,3.5,3.59]:
    print(D,'g(n) max/min on (0.05,0.3):',round(max(g(n,D) for n in np.linspace(.05,.3,500)),5),round(min(g(n,D) for n in np.linspace(.05,.3,500)),5))
