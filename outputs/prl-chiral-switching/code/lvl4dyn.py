import numpy as np, sys, json, time
from scipy.integrate import solve_ivp
from lvl4 import rhs, mkpath, LOOPS, n
from lvl4ss import all_roots, Hlin
from lvl import ss_lin
from lvl4 import LS
def vec(r): v=r.flatten(order='F'); return np.concatenate([v.real,v.imag])
name=sys.argv[1]; T=float(sys.argv[2]); sign=float(sys.argv[3]) if len(sys.argv)>3 else -1.0
O0,Ot,D0,dl,cut=LOOPS[name]; Ds,Os=D0,O0-Ot
rs=all_roots(Ds,Os,8); print(name,'start',(Ds,Os),'roots n=',[round(s[0]+s[1],4) for s in rs],flush=True)
res={}
for k in (0,-1):
    y0=vec(ss_lin(Hlin(Ds,Os,*rs[k]),LS))
    for sgn in (1,-1):
        t=time.time()
        S=solve_ivp(rhs,[0,T],y0,args=(mkpath(name,T,sgn,sign),),method='LSODA',rtol=1e-7,atol=1e-9)
        lab='CW' if sgn>0 else 'CCW'; res[f'{lab},{k}']=n(S.y[:,-1])
        print('  ',lab,'init n=%.4f -> %.4f  (%.0fs)'%(n(y0),n(S.y[:,-1]),time.time()-t),flush=True)
json.dump(res,open(f'../data/lvl4_{name}_T{int(T)}_s{int(sign)}.json','w'))
