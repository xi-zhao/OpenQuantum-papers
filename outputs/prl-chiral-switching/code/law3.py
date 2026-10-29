import numpy as np, json
from scipy.integrate import quad
from scipy.optimize import brentq
from predict import predict, tables, segs
# exact two-state solution, entering at a fold with p=q:
#   p(1)-q(1) = -int_0^1 q'(v) exp[-T int_v^1 lambda(w) dw] dv
# large-N: q_N(u) -> step at the Maxwell point u_M (width ~1/N)  =>  C -> exp(-T Lambda_N),
#   Lambda_N = int_{u_M}^{1} lambda_N(u) du  (arc from Maxwell point to the end point)
out=[]
for loop in 'AB':
    sw=1 if loop=='A' else -1
    for N in [10,20,30,40]:
        th,lam,pss,nh,nl=tables(N,loop); t0=segs[loop][2]
        lo,hi=(th[0],t0) if sw==1 else (t0,th[-1])
        g=lambda x: float(pss(x))-0.5; xm=brentq(g,lo+1e-6,hi-1e-6)
        x=np.linspace(lo,hi,4001); q=pss(x)
        q75=x[np.argmin(abs(q-0.75))]; q25=x[np.argmin(abs(q-0.25))]; width=abs(q75-q25)
        a,b=sorted([xm,t0]); Lam=quad(lambda z: float(lam(z)),a,b,limit=200)[0]/(2*np.pi)
        Ts=np.logspace(1,6,41); C2=np.array([abs(predict(N,loop,T,1)[1]-predict(N,loop,T,-1)[1]) for T in Ts])
        dev=np.abs(C2-np.exp(-Ts*Lam)); m=(Ts*Lam<3)
        print(loop,N,'theta_M=%.4f width(q 0.75->0.25)=%.4f  Lambda=%.3e  max|C-e^{-TLam}| (x<3)=%.3f'%(xm,width,Lam,dev[m].max()))
        out.append(dict(loop=loop,N=N,thM=xm,width=width,Lam=Lam,Ts=Ts.tolist(),C2=C2.tolist()))
json.dump(out,open('../data/law3.json','w'))
