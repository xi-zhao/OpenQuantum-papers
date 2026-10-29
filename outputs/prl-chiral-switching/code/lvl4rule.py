import numpy as np
from lvl4 import mkpath
from lvl4ss import all_roots
# fold rule: final branch = the unique root just outside the entry fold (continued inward)
cases={'black':[(0.925,1),(0.042,-1)],'cyan':[(0.783,1),(0.042,-1)],'magenta':[(0.975,1),(0.142,-1)]}
for name,cs in cases.items():
    p=mkpath(name,1.0,1)
    for s,d in cs:
        so=s-0.012 if d>0 else s+0.012   # just outside, on the approach side
        D,O=p(so); r=all_roots(D,O,6)
        print(name,'CW' if d>0 else 'CCW','entry at',np.round(p(s),2),'outside root n=',[round(x[0]+x[1],4) for x in r])
