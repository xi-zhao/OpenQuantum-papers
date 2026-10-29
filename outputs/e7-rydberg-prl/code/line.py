import numpy as np, json, sys
from multiprocessing import Pool
Dc,Vc=-1.5011,-1.6260
L=int(sys.argv[1]); pc=None
def init():
    global pc
    from ryd import PeriodicChain
    pc=PeriodicChain(L)
def job(p):
    from masses import levels
    th,t=p; a=np.deg2rad(th)
    e0,w0,ep,wp=levels(pc,Dc+t*np.cos(a),Vc+t*np.sin(a),n0=6,npi=4)
    return dict(L=L,th=th,t=t,e0=e0.tolist(),w0=w0.tolist(),ep=ep.tolist(),wp=wp.tolist())
if __name__=="__main__":
    pts=[(th,t) for th in (90,110,130,180) for t in (0.05,0.075,0.1,0.15,0.2,0.3,0.4)]
    with Pool(int(sys.argv[2]),initializer=init) as P: res=P.map(job,pts)
    json.dump(res,open(f'../data/line_L{L}.json','w'))
