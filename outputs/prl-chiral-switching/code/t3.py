import numpy as np
from qutip import piqs
N=4
s=piqs.Dicke(N=N,emission=1.0)
D=s.lindbladian().full(); w=np.linalg.eigvals(D)
print(np.round(np.sort(w.real),3)); print(np.abs(D).sum())
