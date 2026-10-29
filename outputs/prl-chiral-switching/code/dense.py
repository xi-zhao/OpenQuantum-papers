import numpy as np, sys
from ex import liou
V=10.; Om=2.0; N=int(sys.argv[1])
for D in np.concatenate([np.linspace(3.1,3.6,11),np.linspace(4.5,5.1,13)]):
    w=np.linalg.eigvals(liou(N,D,Om,V).toarray()); w=w[np.argsort(-w.real)]
    print(round(D,3),np.round(w[1:6],5))
