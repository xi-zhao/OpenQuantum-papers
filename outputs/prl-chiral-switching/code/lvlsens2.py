import numpy as np, json
from multiprocessing import Pool
from lvl import roots3, P3
xs=np.concatenate([np.linspace(0,0.02,150),np.linspace(0.02,0.3,200)[1:]])
def has(O,p): return any(len(roots3(D,O,p,xs=xs))>=3 for D in np.linspace(-60,-10,151))
def cusp_Om(p,lo=2.0,hi=16.0):
    if not has(hi,p): return ('>',hi)
    if has(lo,p): return ('<',lo)
    for _ in range(9):
        m=(lo+hi)/2; (hi:=m) if has(m,p) else (lo:=m)
    return hi
def job(a):
    k,f=a; p=dict(P3); p[k]*=f; return [k,f,cusp_Om(p)]
if __name__=="__main__":
    args=[('g21',0.8),('g3',1.2),('alpha',1.2),('g21',0.9),('g3',1.1),('alpha',1.1),('alpha',0.9),('g3',0.9)]
    with Pool(8) as P: res=P.map(job,args)
    for r in res: print(r,flush=True)
    json.dump(res,open('../data/lvlsens2.json','w'))
