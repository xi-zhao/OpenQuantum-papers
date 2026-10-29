import numpy as np, json
from multiprocessing import Pool
from lvl import roots3, P3
xs=np.concatenate([np.linspace(0,0.02,150),np.linspace(0.02,0.3,200)[1:]])
def has(O,p):
    return any(len(roots3(D,O,p,xs=xs))>=3 for D in np.linspace(-45,-15,121))
def cusp_Om(p,lo=4.0,hi=11.0):
    if not has(hi,p): return None
    for _ in range(9):
        m=(lo+hi)/2; (hi:=m) if has(m,p) else (lo:=m)
    return hi
def job(a):
    k,f=a; p=dict(P3); p[k]*=f; return [k,f,cusp_Om(p)]
if __name__=="__main__":
    args=[(k,f) for k in ['g21','g32','g3','alpha','Oc'] for f in (0.8,1.2)]
    with Pool(8) as P: res=P.map(job,args)
    for r in res: print(r,flush=True)
    json.dump(res,open('../data/lvlsens.json','w'))
