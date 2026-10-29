import numpy as np
def mf_point(D,Om,V,n0):
    Dt=D-V*n0; c=(-1j*Om/2*(1-2*n0))/(0.5-1j*Dt)
    return np.array([1-n0,n0,c.real,c.imag,1,1,0,0],float)
