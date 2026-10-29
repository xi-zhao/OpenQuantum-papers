import numpy as np, json
from scipy.integrate import quad
from predict import predict, tables, segs
d=json.load(open('../data/twostate.json'))
def info(N,loop):
    th,lam,pss,nh,nl=tables(N,loop); a,b,t0=segs[loop]
    out={}
    for sgn,lo,hi in [(1,th[0],t0),(-1,t0,th[-1])]:
        x=np.linspace(lo,hi,2001); q=pss(x); l=lam(x)
        ent=x[0] if sgn>0 else x[-1]; pe=float(pss(ent)); q0=float(pss(t0))
        # Maxwell crossing: where q crosses 1/2 on this arc (if any)
        cr=np.flatnonzero(np.diff(np.sign(q-0.5)))
        if len(cr):
            xm=x[cr[-1]] if sgn>0 else x[cr[0]]
            Lam=quad(lambda z: float(lam(z)),min(xm,t0),max(xm,t0))[0]/(2*np.pi)
        else: xm=None; Lam=0.0
        dl=(np.log(lam(t0+1e-4))-np.log(lam(t0-1e-4)))/2e-4
        LamL=float(lam(t0))/(2*np.pi*abs(dl))
        out[sgn]=dict(pe=pe,q0=q0,Lam=Lam,LamL=LamL,xm=xm)
    return out
def Cpred(N,loop,T,key='Lam'):
    o=info(N,loop); th,lam,pss,nh,nl=tables(N,loop)
    p={s: o[s]['q0']+(o[s]['pe']-o[s]['q0'])*np.exp(-T*o[s][key]) for s in (1,-1)}
    return p[1]-p[-1], o
if __name__=="__main__":
    for loop in 'AB':
        for N in [10,20,30,40]:
            o=info(N,loop)
            print(loop,N,' '.join(f"[{'CCW' if s>0 else 'CW'} pe={o[s]['pe']:.3f} q0={o[s]['q0']:.4f} Lam={o[s]['Lam']:.3e} LamLaplace={o[s]['LamL']:.3e}]" for s in (1,-1)))
            th,lam,pss,nh,nl=tables(N,loop); t0=segs[loop][2]
            for T in [100,1000,1e4,1e5]:
                ex=predict(N,loop,T,1)[1]-predict(N,loop,T,-1)[1]
                c1,_=Cpred(N,loop,T); c2,_=Cpred(N,loop,T,'LamL')
                print('    T=%g  2state=%.4f  law(Lam)=%.4f  law(Laplace)=%.4f'%(T,ex,c1,c2))
