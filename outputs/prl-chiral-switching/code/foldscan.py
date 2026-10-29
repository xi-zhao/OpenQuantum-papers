import numpy as np, sys
from ex import eig
V=10.; Om=float(sys.argv[2]); N=int(sys.argv[1]); D0,D1=float(sys.argv[3]),float(sys.argv[4])
for D in np.linspace(D0,D1,int(sys.argv[5])):
    w=eig(N,D,Om,V,k=6,sigma=0.0)[1:]
    w=w[np.abs(w.imag)<1.0]
    print(round(D,4),' '.join(f"{x.real:.5f}{x.imag:+.5f}j" for x in w[:4]),flush=True)
