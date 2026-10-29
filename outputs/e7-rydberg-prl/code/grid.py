import numpy as np, json, sys
from multiprocessing import Pool
L=int(sys.argv[1]); pc=None
def init():
    global pc
    from ryd import PeriodicChain
    pc=PeriodicChain(L)
def job(p):
    from masses import levels
    D,V=p
    try: e0,w0,ep,wp=levels(pc,D,V,n0=6,npi=4)
    except Exception as e: return None
    return dict(D=D,V=V,e0=e0.tolist(),w0=w0.tolist(),ep=ep.tolist(),wp=wp.tolist())
if __name__=="__main__":
    Ds=np.linspace(-2.0,-1.2,21); Vs=np.linspace(-1.75,-1.05,21)
    pts=[(float(D),float(V)) for D in Ds for V in Vs]
    with Pool(int(sys.argv[2]),initializer=init) as P: res=P.map(job,pts,chunksize=4)
    json.dump(res,open(f'../data/grid_L{L}.json','w'))
