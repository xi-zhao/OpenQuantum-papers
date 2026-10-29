import numpy as np, json
from mf import roots
from mfdyn import rhs
from lna import D1, J
V=10.;Om=2.0
def ystate(D,n):
    Dt=D-V*n; c=(-1j*Om/2*(1-2*n))/(0.5-1j*Dt); return np.array([n,c.real,c.imag])
# locate left fold precisely: double root of the cubic
from scipy.optimize import brentq
g=lambda D: np.diff(roots(D,Om,V)[:2])[0] if len(roots(D,Om,V))==3 else -1.0
Df=brentq(lambda D: (len(roots(D,Om,V))==3)-0.5, 3.5, 3.7, xtol=1e-13)
r=roots(Df+1e-9,Om,V); nf=(r[0]+r[1])/2
y=ystate(Df,nf); Jm=J(Df,Om,y); w,R=np.linalg.eig(Jm)
i0=np.argmin(abs(w)); r0=R[:,i0].real; l0=np.linalg.inv(R)[i0].real; l0/=l0@r0
f=lambda y,D: np.array(rhs(0,y,lambda t:(D,Om,V))); h=1e-4
a=l0@((f(y,Df+h)-f(y,Df-h))/(2*h))
b=l0@((f(y+h*r0,Df)+f(y-h*r0,Df)-2*f(y,Df))/(2*h*h))
B,S=D1(Df,Om,y[0],y[1]+1j*y[2]); Q=l0@B@l0
print('fold Df=%.8f n_f=%.5f  J eig at fold:'%(Df,nf),np.round(w,5))
print('normal form  dx/dt = a mu + b x^2 + sqrt(Q/N) xi :  a=%.5f  b=%.5f  Q=%.5f'%(a,b,Q))
json.dump(dict(Df=Df,nf=nf,a=a,b=b,Q=Q),open('../data/normalform.json','w'))
