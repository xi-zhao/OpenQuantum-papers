import numpy as np
from scipy.integrate import solve_ivp
from model import liou, loop, tdist
TR=np.array([1,0,0,1.0])
def evolve_ode(T,s,G,gphi,rho0=None,rtol=1e-10,**kw):
    if rho0 is None: rho0=np.diag([1,0]).astype(complex)
    def f(t,v):
        Lm=liou(*loop(s*2*np.pi*t/T,**kw),G,gphi); w=Lm@v
        return w-(TR@w)*v
    sol=solve_ivp(f,(0,T),rho0.reshape(-1).astype(complex),method='DOP853',rtol=rtol,atol=1e-13)
    r=sol.y[:,-1].reshape(2,2); return r/np.trace(r)
def chir(T,G,gphi,**kw):
    return tdist(evolve_ode(T,1,G,gphi,**kw),evolve_ode(T,-1,G,gphi,**kw))
