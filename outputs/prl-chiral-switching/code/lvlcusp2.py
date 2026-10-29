import numpy as np
from lvl import roots3
from scipy.optimize import brentq
xs=np.concatenate([np.linspace(0,0.02,200),np.linspace(0.02,0.3,300)[1:]])
nb=lambda D,O: len(roots3(D,O,xs=xs))
def interval(O):
    Ds=np.linspace(-45,-25,201); b=[D for D in Ds if nb(D,O)>=3]
    return (min(b),max(b)) if b else None
for O in np.arange(7.6,8.61,0.1): print(round(O,2),interval(O),flush=True)
