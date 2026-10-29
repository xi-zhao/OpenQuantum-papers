import numpy as np, json
from predict import tables, segs
# path functions r_N(u)=lambda_N/lambda_min and q_N(u) along the switching arc; show convergence of the
# Maxwell-step structure (reviewer point 3.1)
out={}
for loop in 'AB':
    sw=1 if loop=='A' else -1
    for N in [10,20,30,40]:
        th,lam,pss,nh,nl=tables(N,loop); t0=segs[loop][2]
        x=np.linspace(th[0],t0,300) if sw==1 else np.linspace(t0,th[-1],300)
        out[f'{loop}{N}']=dict(x=x.tolist(),lam=lam(x).tolist(),q=pss(x).tolist(),t0=t0)
json.dump(out,open('../data/pathfun.json','w')); print('ok')
