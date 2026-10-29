import numpy as np, json, scipy.sparse.linalg as sla
from multiprocessing import Pool
from dyn import ems
from ex import obs_vec
V=10.;Om=2.0
def job(a):
    N,D=a; g,rss,(hi,lo)=ems(N,D,Om,V); tr,nv=obs_vec(N)
    nh,nl,ns=[(nv@x).real/N for x in (hi,lo,rss)]; q=(ns-nl)/(nh-nl); lam=-g
    return [N,D,lam,q,lam*(1-q),lam*q]
if __name__=="__main__":
    args=[(N,D) for D in [3.62,3.65,3.7,3.75,3.8,3.9,4.0] for N in [20,30,40,50,60]]
    with Pool(8) as P: res=P.map(job,args)
    json.dump(res,open('../data/g2_rates.json','w'))
    for D in sorted(set(r[1] for r in res)):
        rr=sorted([r for r in res if r[1]==D]); Ns=np.array([r[0] for r in rr]); khl=np.array([r[4] for r in rr]); klh=np.array([r[5] for r in rr])
        print('D=%.2f  s_hl(40-60)=%.4f  s_hl(20-40)=%.4f  s_lh(40-60)=%.4f  q60=%.3e'%(D,-np.polyfit(Ns[2:],np.log(khl[2:]),1)[0],-np.polyfit(Ns[:3],np.log(khl[:3]),1)[0],-np.polyfit(Ns[2:],np.log(klh[2:]),1)[0],rr[-1][3]),flush=True)
