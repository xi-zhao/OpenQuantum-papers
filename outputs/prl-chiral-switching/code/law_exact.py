import numpy as np, json
from scipy.integrate import quad
from scipy.optimize import brentq
from predict import tables, segs
R=[json.loads(l) for l in open('../data/batch.jsonl')]
out=[]
for loop in 'AB':
    sw=1 if loop=='A' else -1
    for N in [10,20,30,40]:
        th,lam,pss,nh,nl=tables(N,loop); t0=segs[loop][2]; dn=float(nh(t0)-nl(t0))
        lo,hi=(th[0],t0) if sw==1 else (t0,th[-1])
        xm=brentq(lambda x: float(pss(x))-0.5,lo+1e-6,hi-1e-6); a,b=sorted([xm,t0])
        Lam=quad(lambda z: float(lam(z)),a,b,limit=200)[0]/(2*np.pi)
        rs={(round(r['T'],3),r['sgn']):r['nf'] for r in R if r['N']==N and r['loop']==loop}
        ex=[[T,abs(rs[(T,1)]-rs[(T,-1)])/dn] for T in sorted(set(k[0] for k in rs)) if (T,1) in rs and (T,-1) in rs]
        out.append(dict(loop=loop,N=N,Lam=Lam,Thalf=np.log(2)/Lam,ex=ex))
        print(loop,N,'Lambda=%.3e T_half=ln2/Lambda=%.4g'%(Lam,np.log(2)/Lam),' '.join('x=%.2f:%.3f/%.3f'%(T*Lam,c,np.exp(-T*Lam)) for T,c in ex if 0.05<T*Lam<5))
json.dump(out,open('../data/law_exact.json','w'))
