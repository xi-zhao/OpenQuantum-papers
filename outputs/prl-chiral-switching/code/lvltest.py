import numpy as np, json, sys
from lvl import *
Om0,pp,D0,dl=10.58,0.31,-52.80,14.92
def exp_path(T,s):
    return lambda t:(D0+dl*np.sin(s*2*np.pi*t/T), Om0*np.sqrt(1+pp*np.cos(s*2*np.pi*t/T)))
r=roots3(D0,Om0*np.sqrt(1+pp)); print('start roots',r)
out={}
for T in [5.0,50.0]:   # in units 1/(2pi MHz)^-1 -> microseconds/(2pi)
  for s in (1,-1):
    for k in (0,-1):
        S=run3(exp_path(T,s),T,ss_branch(D0,Om0*np.sqrt(1+pp),r[k]))
        nf=S.y[8,-1]; out[f'{T},{s},{k}']=nf; print('T',T,'dir',s,'init',round(r[k],4),'->',round(nf,4),flush=True)
json.dump(out,open('../data/lvl_exploop.json','w'))
