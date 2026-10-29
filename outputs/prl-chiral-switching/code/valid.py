import numpy as np, scipy.sparse.linalg as sla
from ex import liou, obs_vec
# N=1 check vs single-atom formula; N=2 check vs brute-force qutip
def ss(N,D,Om,V):
    L=liou(N,D,Om,V); tr,nv=obs_vec(N)
    w,v=sla.eigs(L,k=1,sigma=1e-6); r=v[:,0]; r=r/(tr@r); return (nv@r).real/N
D,Om=1.3,0.8
print(ss(1,D,Om,10), (Om**2/4)/(D**2+0.25+Om**2/2))
import qutip as q
sm=q.destroy(2); 
def op(o,i,N): return q.tensor([o if k==i else q.qeye(2) for k in range(N)])
N=3;V=4.
ne=[op(q.num(2),i,N) for i in range(N)]; Ne=sum(ne)
H=-D*Ne+Om*sum(op(q.sigmax()/2,i,N) for i in range(N))+(V/(2*N))*(Ne*Ne-Ne)
c=[op(q.destroy(2),i,N) for i in range(N)]
rho=q.steadystate(H,c); print(q.expect(Ne,rho)/N, ss(N,D,Om,V))
