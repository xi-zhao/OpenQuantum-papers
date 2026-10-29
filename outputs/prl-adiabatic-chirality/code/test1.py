import numpy as np,sys
from model import *
G=1.0
for gphi in [0.05,0.2,1.0]:
    A,w=amp(G,gphi); print('gphi',gphi,'pred C*T=',A,'gaps',np.sort(w.real)[::-1][:2])
    for T in [10,30,100,300,1000,3000]:
        c=tdist(evolve(T,1,G,gphi,nst=max(2000,int(T*20))),evolve(T,-1,G,gphi,nst=max(2000,int(T*20))))
        print('  T=%6g C=%.4e C*T=%.4f'%(T,c,c*T))
