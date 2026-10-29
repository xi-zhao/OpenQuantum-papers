import numpy as np, json, sys
from ex import eig
from mf import roots
from multiprocessing import Pool
V=10.
# points along the wedge axis toward the upper cusp (D_c,Om_c)=(5.4848,4.4221) and across folds at Om=2
Dc,Oc=5.484761,4.422111
def job(p):
    D,Om=p; out=[]
    for N in [20,30,40,50,60]:
        w=eig(N,D,Om,V,k=3,sigma=0.0); out.append(-w[1].real)
    return [D,Om,out,len(roots(D,Om,V))]
if __name__=="__main__":
    pts=[(D,2.0) for D in np.linspace(3.55,4.6,22)]
    # axis of wedge: midpoints of bistable interval at several Om
    for Om in [2.4,2.8,3.2,3.6,4.0]:
        Ds=np.linspace(0,7,14001); bi=[D for D in Ds if len(roots(D,Om,V))==3]
        # point of max barrier ~ near Maxwell; scan across interval
        for D in np.linspace(min(bi),max(bi),9): pts.append((D,Om))
    with Pool(9) as P: res=P.map(job,pts)
    json.dump(res,open('../data/sfit.json','w'))
    for r in res: print(r[0],r[1],np.round(r[2],6),r[3])
