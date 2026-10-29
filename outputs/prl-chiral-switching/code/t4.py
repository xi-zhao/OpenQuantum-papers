import numpy as np, time
from ex import liou, eig, obs_vec
L=liou(4,0,0,0).toarray(); w=np.linalg.eigvals(L); print(L.shape,np.round(np.sort(w.real)[-4:],3))
for N in [10,20,40]:
    t=time.time(); print(N, liou(N,4,2,10).shape, np.round(eig(N,4.1,2,10,sigma=0.01),5), time.time()-t)
