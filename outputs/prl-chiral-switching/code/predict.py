import numpy as np, json
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
d=json.load(open('../data/twostate.json')); segs=d['segs']
from functools import lru_cache
@lru_cache(None)
def tables(N,loop):
    rows=sorted([r for r in d['res'] if r[0]==N and r[1]==loop],key=lambda r:r[2])
    th,lam,nss,nh,nl=[np.array([r[k] for r in rows]) for k in range(2,7)]
    p=(nss-nl)/(nh-nl)
    f=lambda y: interp1d(th,y,kind='cubic')
    return th,f(lam),f(p),f(nh),f(nl)
def predict(N,loop,T,sgn):
    """two-state reduction: dp/dt=-lam1(th)(p-pss(th)); re-entry at fold, then freeze/relax to th0"""
    th,lam,pss,nh,nl=tables(N,loop); a,b,t0=segs[loop]
    w=2*np.pi/T
    if sgn>0: s0,s1=th[0],t0      # enter through th=a side, move up to t0
    else:     s0,s1=th[-1],t0     # enter through th=b side, move down to t0
    # dp/dth = -(lam/(w)) (p-pss) * sign
    rhs=lambda x,p: -lam(x)/w*(p-pss(x))*np.sign(s1-s0)
    sol=solve_ivp(rhs,[s0,s1],[float(pss(s0))],rtol=1e-9,atol=1e-12)
    p=sol.y[0,-1]; return float(nl(t0)+p*(nh(t0)-nl(t0))), float(p)
if __name__=="__main__":
    for loop in 'AB':
        for N in [10,20,30,40]:
            print(loop,N,[ (round(T),round(predict(N,loop,T,1)[0],4),round(predict(N,loop,T,-1)[0],4)) for T in np.logspace(1,5,9)])
