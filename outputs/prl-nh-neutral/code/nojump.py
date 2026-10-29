# Levine-Pichler CZ with Rydberg decay; no-jump Kraus K = T exp(-i H_eff t)
import numpy as np
from scipy.linalg import expm
def H(Om,De,phi,G,V=1e4):
    h=np.zeros((3,3),complex); h[1,2]=Om/2*np.exp(1j*phi); h[2,1]=np.conj(h[1,2]); h[2,2]=-De-0.5j*G
    I=np.eye(3); P=np.zeros((3,3)); P[2,2]=1
    return np.kron(h,I)+np.kron(I,h)+V*np.kron(P,P)
Om=1.0; De=0.377371; tau=4.29268; xi=3.90242
for G in [1e-3,1e-2]:
    K=expm(-1j*H(Om,De,xi,G)*tau)@expm(-1j*H(Om,De,0,G)*tau)
    comp=[0,1,3,4]; Kc=K[np.ix_(comp,comp)]
    print('Gamma/Om=',G,' |K| diag:',np.round(np.abs(np.diag(Kc))**2,6), ' 1-p_noerase per state:',np.round(1-np.abs(np.diag(Kc))**2,6))
    print('   time in r (from 1-|K|^2)/G:',np.round((1-np.abs(np.diag(Kc))**2)/G,3))
