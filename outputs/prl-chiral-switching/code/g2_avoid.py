import numpy as np, scipy.sparse.linalg as sla, json
from scipy.optimize import minimize_scalar
from ex import liou
V=10.;Om=2.0
# Test A: is the complex EP the analytic continuation of a Hermitian-type avoided crossing?
#  2x2 model  [[l_a(D), g],[g, l_b(D)]]  ->  EP at D* = D0 +- 2i g/|l_a'-l_b'|
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
def pair(N,D,t):
    w=sla.eigs(liou(N,D,Om,V),k=4,sigma=t,return_eigenvectors=False,tol=1e-13)
    w=np.sort(w[np.argsort(abs(w-t))][:2].real); return w
out=[]
for N in [20,30,40]:
    D0=EP[N][1]; t=EP[N][3]
    res=minimize_scalar(lambda D: np.diff(pair(N,D,t))[0],bracket=(D0-0.01,D0,D0+0.01),tol=1e-10)
    Dm=res.x; gap=res.fun
    h=0.03; a1=pair(N,Dm-h,t); a2=pair(N,Dm+h,t)
    # diabatic slopes: far from the crossing the branches exchange; slope difference from asymptotes
    s=abs(((a2[1]-a1[0])-(a2[0]-a1[1]))/(2*h))
    pred=gap/s   # Im D_EP = 2g/|ds| with gap=2g
    out.append([N,Dm,gap,s,pred,EP[N][2]])
    print(N,'Dmin=%.5f (ReD_EP=%.5f) gap=%.3e |dslope|=%.3f  pred Im D_EP=%.3e  actual=%.3e  ratio=%.3f'%(Dm,D0,gap,s,pred,EP[N][2],pred/EP[N][2]),flush=True)
json.dump(out,open('../data/g2_avoid.json','w'))
