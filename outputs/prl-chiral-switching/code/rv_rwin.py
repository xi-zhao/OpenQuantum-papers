import numpy as np, scipy.sparse.linalg as sla, json, sys
from ex import liou
V=10.;Om=2.0;N=int(sys.argv[1]); c=float(sys.argv[2]); hw=float(sys.argv[3]); out=[]
for D in np.linspace(c-hw,c+hw,41):
    w=sla.eigs(liou(N,D,Om,V),k=4,sigma=-0.894,return_eigenvectors=False,tol=1e-13); w=w[np.argsort(abs(w+0.894))][:2]
    out.append([D,abs(w[0].imag),abs(w[0]-w[1])])
out=np.array(out); cx=out[out[:,1]>1e-9]
print(N,'complex-pair window:',(cx[:,0].min(),cx[:,0].max()) if len(cx) else None,'min |gap| %.2e'%out[:,2].min(),flush=True)
json.dump(out.tolist(),open(f'../data/rwin_N{N}.json','w'))
