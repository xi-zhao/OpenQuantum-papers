import numpy as np
from multiprocessing import Pool
from lvl4ss import all_roots
from matplotlib.path import Path
from lvl4 import mkpath
def width(O):
    b=[D for D in np.arange(-80,-64.9,0.25) if len(all_roots(D,O,6))>=3]
    return (O,(min(b),max(b)) if b else None)
if __name__=="__main__":
    with Pool(8) as P: res=P.map(width,np.arange(66,78.1,1.0))
    for r in res: print(r)
    last=[r for r in res if r[1]][-1]; cusp=(np.mean(last[1]),last[0]); print('cusp approx (Delta,Omw)=',cusp)
    for name in ['black','cyan','magenta']:
        p=mkpath(name,1.0,1,1.0); xy=np.array([p(s) for s in np.linspace(0,1,2001)])
        print(name,'encloses cusp:',Path(xy).contains_point(cusp),' Omw range',xy[:,1].min().round(2),xy[:,1].max().round(2),' Delta range',xy[:,0].min().round(2),xy[:,0].max().round(2))
