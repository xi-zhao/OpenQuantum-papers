import numpy as np
from qutip import piqs
N=6
s=piqs.Dicke(N=N,emission=1.0)
D=s.lindbladian(); print(type(D), D.shape)
w=np.linalg.eigvals(D.full()); print(np.round(sorted(w.real)[-5:],4))
jx,jy,jz=piqs.jspin(N,basis='dicke'); print(jz.shape)
s2=piqs.Dicke(N=N,hamiltonian=jx,emission=1.0); L=s2.liouvillian(); print(np.round(sorted(np.linalg.eigvals(L.full()).real)[-5:],4))
