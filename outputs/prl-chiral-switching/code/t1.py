import numpy as np
from ex import liou
L,Ne=liou(6,4,2,10)
M=L.full(); w=np.linalg.eigvals(M); w=w[np.argsort(-w.real)]
print(M.shape, np.round(w[:8],4))
