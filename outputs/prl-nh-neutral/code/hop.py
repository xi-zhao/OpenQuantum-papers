# Levine-Pichler two-pulse CZ under perfect blockade: probability that a leaked |r> on atom A hops to atom B
import numpy as np
from scipy.linalg import expm
# single atom basis 0,1,r ; two-atom 9 states, |rr> removed by V->inf (large V)
def H(Om,De,phi,V=1e4):
    h=np.zeros((3,3),complex); h[1,2]=Om/2*np.exp(1j*phi); h[2,1]=np.conj(h[1,2]); h[2,2]=-De
    I=np.eye(3); HH=np.kron(h,I)+np.kron(I,h); P=np.zeros((3,3)); P[2,2]=1
    return HH+V*np.kron(P,P)
Om=1.0; De=0.377371*Om; tau=4.29268/Om; xi=3.90242
U=expm(-1j*H(Om,De,xi)*tau)@expm(-1j*H(Om,De,0)*tau)
idx=lambda a,b:3*a+b
for a,lab in [(1,'|r1>'),(0,'|r0>')]:
    v=np.zeros(9,complex); v[idx(2,a)]=1; w=U@v; p=np.abs(w)**2
    print(lab,'-> P(rA)=%.3f P(rB)=%.3f P(no r)=%.3f'%(p[idx(2,0)]+p[idx(2,1)],p[idx(0,2)]+p[idx(1,2)],1-p[[idx(2,0),idx(2,1),idx(0,2),idx(1,2)]].sum()))
# check gate: phases
for s in [(0,1),(1,0),(1,1)]:
    v=np.zeros(9,complex); v[idx(*s)]=1; w=U@v; print(s,abs(w[idx(*s)])**2, np.angle(w[idx(*s)]))
