import numpy as np
from scipy.integrate import solve_ivp
from mf import roots
# MF Bloch eqs for H = -Dt n + Om sx (sx=Pauli/2*2? use Om Jx -> Om/2 sigma_x per atom), Dt = Delta - V n, decay 1
# rho_ee=n, coherence c=rho_ge ; dn/dt = -n - Om Im? derive via Lindblad numerically
def rhs(t,y,path):
    D,Om,V=path(t)
    n=y[0]; c=y[1]+1j*y[2]           # c = rho_eg
    Dt=D-V*n
    # H = -Dt|e><e| + Om/2 (|e><g|+|g><e|)
    # d rho_eg/dt = -i[H,rho]_eg - c/2 ; [H,rho]_eg = -Dt c + Om/2 (rho_gg - rho_ee)
    dc=-1j*(-Dt*c+Om/2*(1-2*n))-c/2
    # d rho_ee/dt = -i[H,rho]_ee - n ; [H,rho]_ee = Om/2 (rho_ge - rho_eg) = Om/2 (c* - c)
    dn=(-1j*Om/2*(np.conj(c)-c)).real-n
    return [dn,dc.real,dc.imag]
def steady(D,Om,V,branch):
    r=roots(D,Om,V); n=r[branch]
    Dt=D-V*n; c=(-1j*Om/2*(1-2*n))/(0.5-1j*Dt)
    return [n,c.real,c.imag]
