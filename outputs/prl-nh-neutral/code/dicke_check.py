# exact thermal Dicke ladder (r->52P3/2 channel), start fully inverted (no r-P coherence, as after a gate pulse)
import numpy as np
from scipy.integrate import solve_ivp
G0=7.39; nth=227; tau=1.02e-7
for M in (100,1000,10000,100000):
    k=np.arange(M+1)
    down=(M-k)*(k+1)*G0*(nth+1); up=k*(M-k+1)*G0*nth
    def f(t,P):
        d=-(down+up)*P; d[1:]+=down[:-1]*P[:-1]; d[:-1]+=up[1:]*P[1:]; return d
    P0=np.zeros(M+1); P0[0]=1
    P=solve_ivp(f,(0,tau),P0,method='LSODA',rtol=1e-10,atol=1e-14).y[:,-1]
    m=(k*P).sum(); ind=M*G0*(nth+1)*tau
    print('M=%6d  <k>=%.4e  independent=%.4e  rel.excess=%.2e  (pred M G0 tau/2=%.2e)'%(M,m,ind,m/ind-1,M*G0*tau/2))
