import sys, numpy as np
sys.path.insert(0,'.')
from syk_ed import EvenLindbladian
N=int(sys.argv[1]); L=EvenLindbladian(N, seed=7)
n=2*L.D**2
h=L.vec_of(L.H); one=L.vec_of(np.eye(2*L.D))
H2=L.H@L.H; h2=L.vec_of(H2)
def orth(v,B):
    for b in B: v=v-np.vdot(b,v)/np.vdot(b,b)*b
    return v
h2=orth(h2,[one,h])
# columns of L once (mu-independent pieces)
A=np.empty((n,n),complex); S=np.empty((n,n),complex); e=np.zeros(n,complex)
for j in range(n):
    e[j]=1; A[:,j]=L.apply(e,0.0); S[:,j]=-L.apply(e,1.0)+A[:,j]; e[j]=0
# S = -size superoperator
res=[]
for mu in [0.01,0.03,0.05,0.1,0.2,0.4]:
    w,V=np.linalg.eig(A-mu*S)
    o=np.argsort(-w.real); w=w[o]; V=V[:,o]
    ovH=np.abs(V.conj().T@h)/np.linalg.norm(h); ovH2=np.abs(V.conj().T@h2)/np.linalg.norm(h2)
    iH=np.argmax(ovH); iH2=np.argmax(ovH2)
    perp=[k for k in range(1,n) if k!=iH and k!=iH2]
    res.append([mu, -w[iH].real/mu, ovH[iH], -w[iH2].real/mu, ovH2[iH2], -w[perp[0]].real, -w[1].real])
    print(res[-1], flush=True)
np.save(f'../data/dense_N{N}.npy', np.array(res))
