import numpy as np, json
from multiprocessing import Pool
from lvl4ss import all_roots
def job(O):
    Ds=np.arange(-100,-29.9,1.0); b=[]
    for D in Ds:
        r=all_roots(D,O,6)
        if len(r)>=3: b.append(D)
    return [float(O),(min(b),max(b)) if b else None]
if __name__=="__main__":
    Os=[5,10,15.67,20,25,30,35,40,45,50,60,70,80,90,100]
    with Pool(8) as P: res=P.map(job,Os)
    for r in res: print(r,flush=True)
    json.dump(res,open('../data/lvl4wedge.json','w'))
