import numpy as np
from lvl import ops, liou
from scipy.integrate import solve_ivp
P4=dict(Oc=12.21,Op=13.07,g21=7.34,g32=0.70,g31=0.13,g41=0.13,g3=11.78,g4=11.78,alpha=-174.32,beta=-954.93,gam=6843.66,eps=1114.08)
def H4(Dc,Omw,r,p=P4):
    n3,n4=r[2,2].real,r[3,3].real
    D1=Dc-p['alpha']*n3-p['beta']*n4; D2=-p['gam']*n3-p['eps']*n4
    return np.array([[0,p['Op']/2,0,0],[p['Op']/2,0,p['Oc']/2,0],[0,p['Oc']/2,-D1,Omw/2],[0,0,Omw/2,-D2]],complex)
s=ops(4); P=P4
# paper writes L_i = gamma |i><j| ; we take rates as gamma (amplitude sqrt(gamma)), as in their 3-level Bloch matrix
LS=[np.sqrt(P['g21'])*s(0,1),np.sqrt(P['g32'])*s(1,2),np.sqrt(P['g31'])*s(0,2),np.sqrt(P['g3'])*s(2,2),np.sqrt(P['g41'])*s(0,3),np.sqrt(P['g4'])*s(3,3)]
def rhs(t,y,path):
    r=(y[:16]+1j*y[16:]).reshape(4,4,order='F'); Dc,Omw=path(t)
    dr=liou(H4(Dc,Omw,r),LS)@r.flatten(order='F'); return np.concatenate([dr.real,dr.imag])
def vec(r): v=r.flatten(order='F'); return np.concatenate([v.real,v.imag])
def n(y): return y[10]+y[15]   # rho33+rho44 (F-order index 2*4+2=10, 3*4+3=15)
LOOPS={'black':(34.48,18.81,-57.80,7.54,None),'cyan':(34.48,18.81,-57.80,7.54,'cyan'),'magenta':(53.27,37.60,-57.80,29.86,'mag')}
SIGN=-1.0  # Delta_c = D0 + SIGN*dl*f*sin(th); paper table formula has SIGN=-1
def mkpath(name,T,sgn,sign=None):
    O0,Ot,D0,dl,cut=LOOPS[name]; sg=SIGN if sign is None else sign
    def path(t):
        th=(2*np.pi*t/T) if sgn>0 else 2*np.pi*(1-t/T)   # paper: theta=2pi n/N is CW
        th=th%(2*np.pi); f=1.0
        if cut=='cyan': f=1.0 if th<np.pi else 0.0
        if cut=='mag': f=1.0 if th>np.pi else 0.0
        return D0+sg*dl*f*np.sin(th), O0-Ot*np.cos(th)
    return path
def relax(Dc,Omw,y0,T=200):
    p=lambda t:(Dc,Omw); S=solve_ivp(rhs,[0,T],y0,args=(p,),rtol=1e-8,atol=1e-10,method='LSODA'); return S.y[:,-1]
