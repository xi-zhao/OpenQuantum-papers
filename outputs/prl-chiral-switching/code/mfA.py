import numpy as np
from scipy.integrate import solve_ivp
from mfdyn import rhs, steady
from loops import mkpath, P0, V
for loop in 'AB':
  for T in [100,1000]:
    for sgn in (1,-1):
      for br in (0,2):
        s=solve_ivp(rhs,[0,T],steady(*P0,V,br),args=(mkpath(loop,T,sgn),),rtol=1e-8,atol=1e-10,method='LSODA')
        print(loop,T,sgn,br,round(s.y[0,-1],4))
