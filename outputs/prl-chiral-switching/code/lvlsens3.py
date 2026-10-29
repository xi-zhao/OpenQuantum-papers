import numpy as np
from multiprocessing import Pool
from lvl import roots3, P3
xs=np.concatenate([np.linspace(0,0.02,150),np.linspace(0.02,0.3,200)[1:]])
def job(a):
    k,f=a; p=dict(P3); p[k]*=f
    hits=[(O,D) for O in np.arange(2,30.1,2.0) for D in np.arange(-130,0.1,2.5) if len(roots3(D,O,p,xs=xs))>=3]
    return [k,f,(min(h[0] for h in hits),max(h[0] for h in hits)) if hits else None, len(hits)]
if __name__=="__main__":
    with Pool(5) as P: res=P.map(job,[('g21',0.8),('g21',0.9),('g3',1.2),('alpha',1.1),('alpha',1.2)])
    for r in res: print(r,flush=True)
