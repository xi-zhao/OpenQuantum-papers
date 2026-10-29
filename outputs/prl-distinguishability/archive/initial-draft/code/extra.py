import numpy as np, sys
sys.argv=["x","2","1"]
from sim import T as Tsim
from theory2 import lnT
from exponent import E, qc
L=np.log2
for n in [8,10]:
  for r in [0.7,0.9]:
    ns=30 if n==8 else 8
    v=[Tsim(n,r) for _ in range(ns)]
    print("frozen",n,r,round(np.mean(v),5),round(np.std(v)/np.sqrt(ns),5),round(np.exp(lnT(n,r,20001)),5))
rc=4*qc/3; q=qc
s0,s1=-L(1-q),-L(q/3); m=(1-q)*s0+q*s1; s2=(1-q)*s0**2+q*s1**2-m*m
dH=0.75*(-L(1-q)*-1 + 0)  # placeholder
from scipy.misc import derivative
