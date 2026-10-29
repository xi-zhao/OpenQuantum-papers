import numpy as np, json
from scipy.optimize import fsolve
from lvl4ss import G
from multiprocessing import Pool
def sweep(Omw,Ds):
    x=np.array([0.0,0.0]); out=[]
    for D in Ds:
        for _ in range(3):
            y=fsolve(G,x,args=(D,Omw),xtol=1e-12)
            if np.max(abs(G(y,D,Omw)))<1e-9: break
        # follow branch; if residual bad do small relaxation by restarting from previous
        x=y; out.append(x[0]+x[1])
    return np.array(out)
def job(O):
    Ds=np.linspace(-90,-30,241)
    f=sweep(O,Ds); b=sweep(O,Ds[::-1])[::-1]
    bi=np.abs(f-b)>1e-3
    return [O,(Ds[bi].min(),Ds[bi].max()) if bi.any() else None]
if __name__=="__main__":
    Os=np.arange(0,70.1,2.5)
    with Pool(8) as P: res=P.map(job,Os)
    for r in res: print(r)
    json.dump([[r[0],list(r[1]) if r[1] else None] for r in res],open('../data/lvl4map.json','w'))
