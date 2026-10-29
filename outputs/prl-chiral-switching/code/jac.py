import numpy as np
from mfdyn import rhs, steady
from mf import roots
V=10.; Om=2.0
def J(D,br):
    y=np.array(steady(D,Om,V,br)); f=lambda y: np.array(rhs(0,y,lambda t:(D,Om,V)))
    h=1e-7; return np.array([(f(y+h*e)-f(y-h*e))/(2*h) for e in np.eye(3)]).T
for D in np.linspace(3.1,3.7,13):
    r=roots(D,Om,V); out=[]
    for br in range(len(r)):
        out.append((round(r[br],3),np.round(np.sort_complex(np.linalg.eigvals(J(D,br))),4)))
    print(round(D,3),out)
