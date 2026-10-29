import numpy as np, sys
from lvl4 import mkpath
from lvl4ss import all_roots
name=sys.argv[1]; M=int(sys.argv[2])
sign=float(sys.argv[3]) if len(sys.argv)>3 else -1.0
p=mkpath(name,1.0,1,sign)  # sgn=+1 : CW per paper convention
s=np.linspace(0,1,M+1); nb=[]
for x in s:
    D,O=p(x); nb.append(len(all_roots(D,O,6)))
nb=np.array(nb); print(name,''.join('3' if k>=3 else '.' for k in nb))
for i in np.flatnonzero(np.diff(nb>=3)):
    D,O=p(s[i]); print('  crossing at s=%.3f (D=%.2f, Omw=%.2f) %s'%(s[i],D,O,'enter' if nb[i+1]>=3 else 'leave'),flush=True)
