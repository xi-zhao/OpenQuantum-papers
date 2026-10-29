import numpy as np
from pblock import Nmat
from scipy.linalg import eigvals
def mu(n,d,z): return np.max(eigvals(Nmat(n,d,z)).real)
# scaling collapse mu/n = F(z n^{d-1})
for d in (2,3):
    print('d=',d)
    for zeta in (1e-4,1e-2,1,100):
        row=[]
        for n in (20,40,80,160):
            z=zeta/n**(d-1); row.append(mu(n,d,z)/n)
        print('  zeta=%g: mu/n ='%zeta, np.round(row,4))
