import numpy as np, json
from predict import predict, tables, segs
from scipy.integrate import quad
R=[json.loads(l) for l in open('../data/batch.jsonl')]
def Lam(N,loop):
    th,lam,pss,nh,nl=tables(N,loop); a,b,t0=segs[loop]
    # arc traversed inside bistable window before returning to th0: CCW enters at a for A? use both
    return quad(lambda x: lam(x),th[0],t0)[0]/(2*np.pi), quad(lambda x: lam(x),t0,th[-1])[0]/(2*np.pi)
for loop in 'AB':
  for N in [10,20,30]:
    rs={(round(r['T'],3),r['sgn']):r['nf'] for r in R if r['N']==N and r['loop']==loop}
    Ts=sorted(set(k[0] for k in rs))
    th,lam,pss,nh,nl=tables(N,loop); t0=segs[loop][2]; dn=float(nh(t0)-nl(t0))
    La=Lam(N,loop)
    for T in Ts:
        if (T,1) in rs and (T,-1) in rs:
            Cex=(rs[(T,1)]-rs[(T,-1)])/dn
            Cp=(predict(N,loop,T,1)[0]-predict(N,loop,T,-1)[0])/dn
            print(loop,N,T,'C_exact=%.4f C_2state=%.4f'%(Cex,Cp))
    print(loop,N,'Lambda arcs',La)
