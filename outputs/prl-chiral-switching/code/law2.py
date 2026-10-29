import numpy as np, json
from scipy.integrate import quad
from predict import predict, tables, segs
def theta(N,loop):
    """Theta_N = (1/T) * integral of k_exit dt over the switching arc, k_exit = lambda_sw (1 - q) for the hi state
       (state enters the bistable arc on hi through the 'switching' fold and must survive to th0)."""
    th,lam,pss,nh,nl=tables(N,loop); a,b,t0=segs[loop]
    sw=1 if loop=='A' else -1                  # direction that enters on the branch that is metastable at P0
    lo,hi=(th[0],t0) if sw==1 else (t0,th[-1])
    f=lambda x: float(lam(x)*(1-pss(x)))
    return quad(f,lo,hi,limit=200)[0]/(2*np.pi), sw
if __name__=="__main__":
    R=[json.loads(l) for l in open('../data/batch.jsonl')]
    rows=[]
    for loop in 'AB':
        for N in [10,20,30,40]:
            Th,sw=theta(N,loop); th,lam,pss,nh,nl=tables(N,loop); t0=segs[loop][2]; dn=float(nh(t0)-nl(t0))
            x=np.linspace(th[0],t0,400) if sw==1 else np.linspace(t0,th[-1],400); lm=float(lam(x).min())
            print(loop,N,'Theta=%.3e lam_min=%.3e ratio=%.3f'%(Th,lm,Th/lm))
            rs={(round(r['T'],3),r['sgn']):r['nf'] for r in R if r['N']==N and r['loop']==loop}
            for T in sorted(set(k[0] for k in rs)):
                if (T,1) in rs and (T,-1) in rs:
                    Cex=abs(rs[(T,1)]-rs[(T,-1)])/dn
                    C2=abs(predict(N,loop,T,1)[0]-predict(N,loop,T,-1)[0])/dn
                    rows.append([loop,N,T,Th,lm,Cex,C2]); print('   T=%-9.1f x=T*Theta=%.3f  C_exact=%.4f  C_2state=%.4f  exp(-x)=%.4f'%(T,T*Th,Cex,C2,np.exp(-T*Th)))
    json.dump(rows,open('../data/law2.json','w'))
