import numpy as np
from lvl import roots3
res=[]
for Op in np.arange(8,30.1,1.0):
    Ds=np.linspace(-110,-25,171); b=[D for D in Ds if len(roots3(D,Op))>=3]
    res.append((Op,(min(b),max(b)) if b else None)); print(Op,res[-1][1],flush=True)
    if not b and Op>12: break
