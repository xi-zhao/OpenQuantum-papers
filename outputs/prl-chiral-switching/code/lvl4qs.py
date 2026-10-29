import numpy as np, sys
from scipy.optimize import fsolve
from lvl4 import mkpath
from lvl4ss import G, all_roots
# quasi-static (T->inf, MF) outcome: follow branch by continuation; at loss jump to the other stable branch
def qs(name,sgn,sign,x0,M=600):
    p=mkpath(name,1.0,sgn,sign); x=np.array(x0)
    for k in range(1,M+1):
        D,O=p(k/M)
        y=fsolve(G,x,args=(D,O),xtol=1e-12)
        if np.max(abs(G(y,D,O)))>1e-9 or np.hypot(*(y-x))>0.03:
            rs=all_roots(D,O,6)
            y=rs[0] if (x[0]+x[1])>0.2 else rs[-1]   # fold: high branch lost -> low, and vice versa
        x=y
    return x[0]+x[1]
name=sys.argv[1]; sign=float(sys.argv[2])
p=mkpath(name,1.0,1,sign); D0,O0=p(0); rs=all_roots(D0,O0,8)
for sgn,lab in [(1,'CW'),(-1,'CCW')]:
    print(name,'sign',sign,lab,[ '%.4f->%.4f'%(r[0]+r[1],qs(name,sgn,sign,r)) for r in (rs[0],rs[-1])],flush=True)
