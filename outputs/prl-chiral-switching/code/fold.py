import numpy as np
from mf import roots
V=10.
for Om in np.arange(0.4,4.6,0.3):
    Ds=np.linspace(0,7,7001); bi=[D for D in Ds if len(roots(D,Om,V))==3]
    print(round(Om,2), (round(min(bi),3),round(max(bi),3)) if bi else None)
