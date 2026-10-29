import numpy as np, json, sys
from multiprocessing import Pool
Dc,Vc=-1.5011,-1.6260
L=int(sys.argv[1])
pc=None
def init():
    global pc
    from ryd import PeriodicChain
    pc=PeriodicChain(L)
def job(p):
    from masses import levels
    th,t=p; a=np.deg2rad(th)
    e0,w0,ep,wp=levels(pc,Dc+t*np.cos(a),Vc+t*np.sin(a),n0=8,npi=6)
    return dict(L=L,th=th,t=t,e0=e0.tolist(),w0=w0.tolist(),ep=ep.tolist(),wp=wp.tolist())
if __name__=="__main__":
    pts=[(th,t) for t in (0.1,0.15,0.2,0.3) for th in range(60,215,10)]
    with Pool(int(sys.argv[2]),initializer=init) as P: res=P.map(job,pts)
    json.dump(res,open(f'../data/thscan_L{L}.json','w'))
