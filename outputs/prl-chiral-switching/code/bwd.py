import numpy as np, warnings, json, sys
from scipy.integrate import solve_ivp
import lin
from wkbr import F
from mf import roots
warnings.filterwarnings('ignore')
D=float(sys.argv[1]); r=roots(D,2.,10.); S=lin.fixed(D,2.,10.,r[1]); US,w=lin.basis(S,D,True); A=lin.fixed(D,2.,10.,r[2])
best=(9,None,None,None); nmax=0
for a1 in np.linspace(0.05,np.pi-0.05,12):
    for a2 in np.linspace(0,2*np.pi,12,endpoint=False):
        y0=S+1e-6*US@lin.sph((a1,a2))
        def blow(t,z,*k): return np.abs(z).max()-15
        blow.terminal=True
        sol=solve_ivp(lambda t,y: -F(y,D,2.,10.),[0,30],y0,rtol=1e-8,atol=1e-10,events=blow)
        nmax=max(nmax,sol.y[0].max())
        d=np.linalg.norm(sol.y-A[:,None],axis=0); k=int(np.argmin(d))
        if d[k]<best[0]: best=(d[k],(a1,a2),sol.y[:,k],sol.t[k])
print('D=%.2f closest backward approach to hi: dist=%.4f t=%.1f  max n reached %.4f (n_hi %.4f)'%(D,best[0],best[3],nmax,r[2]),flush=True)
