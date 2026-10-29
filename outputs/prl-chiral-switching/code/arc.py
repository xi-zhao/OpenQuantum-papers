import numpy as np, json
from multiprocessing import Pool
from dyn import ems
from ex import obs_vec
from loops import LOOPS, start_angle, P0, V
def job(a):
    N,loop,th=a; L=LOOPS[loop]; C,R=L['C'],L['R']
    D,Om=C[0]+R[0]*np.cos(th),C[1]+R[1]*np.sin(th)
    try:
        g,rss,(hi,lo)=ems(N,D,Om,V); tr,nv=obs_vec(N)
        return [N,loop,th,-g,(nv@rss).real/N,(nv@hi).real/N,(nv@lo).real/N]
    except Exception as e: return [N,loop,th,None,None,None,None]
if __name__=="__main__":
    args=[]
    for loop in 'AB':
        t0=start_angle(LOOPS[loop])
        for N in [10,20,30,40]:
            for th in t0+np.linspace(-np.pi,np.pi,121): args.append((N,loop,th))
    with Pool(9) as P: res=P.map(job,args,chunksize=4)
    json.dump(res,open('../data/arc.json','w')); print('done',len(res))
