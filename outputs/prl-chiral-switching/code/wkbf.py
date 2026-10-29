import numpy as np
"""Fast Hamilton equations for h(X,rho)=Re Tr[X l0(rho)] + 2V ReA ImA, A=Tr(X n rho)."""
nM=np.array([[0,0],[0,1]],complex); sx=np.array([[0,1],[1,0]],complex); sm=np.array([[0,1],[0,0]],complex); sp_=sm.conj().T
def unpack(y):
    rgg,ree,u,v,pgg,pee,pu,pv=y
    rho=np.array([[rgg,u-1j*v],[u+1j*v,ree]]); X=np.array([[pgg,(pu-1j*pv)/2],[(pu+1j*pv)/2,pee]]); return rho,X
def l0(r,D,Om):
    H=-D*nM+Om/2*sx; return -1j*(H@r-r@H)+sm@r@sp_-0.5*(sp_@sm@r+r@sp_@sm)
def l0d(X,D,Om):
    H=-D*nM+Om/2*sx; return 1j*(H@X-X@H)+sp_@X@sm-0.5*(sp_@sm@X+X@sp_@sm)
def ham(y,D,Om,V):
    r,X=unpack(y); A=np.trace(X@nM@r); return (np.trace(X@l0(r,D,Om))).real+2*V*A.real*A.imag
def flow(t,y,D,Om,V):
    r,X=unpack(y); A=np.trace(X@nM@r)
    G=l0(r,D,Om)+2*V*(A.imag*(nM@r+r@nM)/2+A.real*(nM@r-r@nM)/(2j))
    K=-(l0d(X,D,Om)+2*V*(A.imag*(X@nM+nM@X)/2+A.real*(X@nM-nM@X)/(2j)))
    return np.array([G[0,0].real,G[1,1].real,G[1,0].real,G[1,0].imag,K[0,0].real,K[1,1].real,2*K[1,0].real,2*K[1,0].imag])
def lin(y,D,Om,V,h=1e-7):
    return np.array([(flow(0,y+h*e,D,Om,V)-flow(0,y-h*e,D,Om,V))/(2*h) for e in np.eye(8)]).T
if __name__=="__main__":
    import wkbq, time
    rng=np.random.default_rng(1)
    for _ in range(3):
        y=rng.normal(size=8); a=flow(0,y,3.7,2.0,10.); b=wkbq.flow(0,y,3.7,2.0,10.)
        print('max diff vs sympy:',np.abs(a-b).max(),' h diff:',abs(ham(y,3.7,2.,10.)-wkbq.f_h(y,3.7,2.,10.)))
    t=time.time(); [flow(0,y,3.7,2.,10.) for _ in range(10000)]; print('us/call',(time.time()-t)*100)
