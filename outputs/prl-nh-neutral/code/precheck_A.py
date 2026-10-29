# Candidate A pre-check: collective (Dicke) vs independent Rydberg leakage in a parallel CZ layer, Rb 53S1/2, T=300K
import numpy as np, json
from scipy.special import gammaln
d=json.load(open('../data_nu.json'))
G_sp,G_bbr,nu=d['G_sp'],d['G_bbr'],d['nu_eff']
Om=2*np.pi*4.6e6; tauR=2.96/Om            # mean time in |r> per time-optimal CZ (Jandura-Pupillo), Mok SI params
G0=7.39                                     # bare 53S->52P3/2 rate (ARC, T=0)
p=(G_sp+G_bbr)*tauR; pB=G_bbr*tauR
print('tau_R=%.3g s  p_leak=%.2e  p_BBR=%.2e  nu_eff=%.1f'%(tauR,p,pB,nu))
# (1) superradiant excess: gate states carry no S-P coherence -> Dicke ladder starts at full inversion;
#     excess emission fraction on dominant channel ~ (M-1) G0 tau/2 ; Mok worst case per-atom NG0/4 vs G_ind
for M in (256,1e3,1e4,1.6e4,1e5):
    exc=(M-1)*G0*tauR/2*G0/(G_sp+G_bbr)
    mok=(M-2)*G0/4/(G_sp+G_bbr)
    # (2) common thermal field: Fano factor of leak count, pair excess correlation
    F=1+M*pB/nu
    print('M=%7d  Dicke excess(gate state)=%.1e   Mok worst-case ratio=%.2f   Fano=%.2f'%(M,exc,mok,F))
print('pair g2-1 = 1/nu = %.3f ; M* (Fano=2) = nu/p_BBR = %.0f atoms'%(1/nu,nu/pB))
# (3) k-fold joint leakage enhancement within a layer: <X^k>, X~Gamma(nu)/nu
for k in (2,3,5,8,12):
    print(' k=%2d  enhancement=%.2f'%(k,np.exp(gammaln(nu+k)-gammaln(nu)-k*np.log(nu))))
