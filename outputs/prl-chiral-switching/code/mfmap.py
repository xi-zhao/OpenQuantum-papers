import numpy as np
from mf import roots
V=10
for Om in [0.3,1,2,3,4,4.4,5]:
    Ds=np.linspace(-2,15,3401); bi=[D for D in Ds if len(roots(D,Om,V))==3]
    print(Om, (round(min(bi),3),round(max(bi),3)) if bi else None)
