import numpy as np, json, sys
from multiprocessing import Pool
from dyn import ems, evolve
from loops import mkpath, P0, V
def job(a):
    N,loop,T,sgn,br,steps=a
    g,rss,(hi,lo)=ems(N,P0[0],P0[1],V)
    rho0={'hi':hi,'lo':lo,'ss':rss}[br]
    v,tr=evolve(N,rho0,mkpath(loop,T,sgn),T,steps)
    return dict(N=N,loop=loop,T=T,sgn=sgn,br=br,steps=steps,nf=float(tr[-1]),traj=tr[::max(1,steps//400)].tolist())
if __name__=="__main__":
    N=int(sys.argv[1]); loop=sys.argv[2]; steps=int(sys.argv[3])
    Ts=[float(x) for x in sys.argv[4].split(',')]
    brs=sys.argv[5].split(',') if len(sys.argv)>5 else ['ss']
    args=[(N,loop,T,s,b,steps) for T in Ts for s in (1,-1) for b in brs]
    with Pool(8) as P: res=P.map(job,args)
    json.dump(res,open(f'../data/chir_N{N}_{loop}_s{steps}.json','w'))
    for r in res: print(r['T'],r['sgn'],r['br'],round(r['nf'],4))
