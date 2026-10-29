import numpy as np, json
from multiprocessing import Pool
from dyn import ems
from ex import obs_vec
from loops import LOOPS, start_angle, P0, V
from mf import roots
def seg(loop):
    # bistable angular window containing th0 (MF), returned as (th_lo, th_hi)
    L=LOOPS[loop]; C,R=L['C'],L['R']; t0=start_angle(L)
    f=lambda th: len(roots(C[0]+R[0]*np.cos(th),C[1]+R[1]*np.sin(th),V))==3
    a=t0
    while f(a-1e-3): a-=1e-3
    b=t0
    while f(b+1e-3): b+=1e-3
    return a,b,t0
def job(x):
    N,loop,th=x; L=LOOPS[loop]; C,R=L['C'],L['R']
    D,Om=C[0]+R[0]*np.cos(th),C[1]+R[1]*np.sin(th)
    g,rss,(hi,lo)=ems(N,D,Om,V); tr,nv=obs_vec(N)
    return [N,loop,th,-g,(nv@rss).real/N,(nv@hi).real/N,(nv@lo).real/N]
if __name__=="__main__":
    args=[]; segs={}
    for loop in 'AB':
        a,b,t0=seg(loop); segs[loop]=(a,b,t0); print(loop,a,b,t0)
        for N in [10,20,30,40]:
            for th in np.linspace(a-0.02,b+0.02,81): args.append((N,loop,th))
    with Pool(5) as P: res=P.map(job,args,chunksize=4)
    json.dump(dict(segs=segs,res=res),open('../data/twostate.json','w')); print('done')
