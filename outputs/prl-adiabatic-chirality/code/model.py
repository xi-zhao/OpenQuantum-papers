import numpy as np
from scipy.linalg import expm, eig
I2=np.eye(2)
def liou(d1,O1,G,gphi):
    H=np.array([[d1/2,-O1],[-O1,-d1/2-1j*G]])
    L=np.sqrt(gphi)*np.array([[0,0],[0,1]])
    # row-major vec: vec(A rho B)= kron(A,B.T) vec(rho)
    LL=L.conj().T@L
    return (-1j*np.kron(H,I2)+1j*np.kron(I2,H.conj()) + np.kron(L,L.conj())
            -0.5*np.kron(LL,I2)-0.5*np.kron(I2,LL.T))
def loop(phi,phi0=2*np.pi/3,a=0.5,c=0.5,b=0.5):
    return a*np.sin(phi+phi0), c+b*np.cos(phi+phi0)
def evolve(T,s,G,gphi,rho0=None,nst=4000,**kw):
    if rho0 is None: rho0=np.diag([1,0]).astype(complex)
    v=rho0.reshape(-1); dt=T/nst
    for k in range(nst):
        ph=s*2*np.pi*(k+0.5)/nst
        v=expm(liou(*loop(ph,**kw),G,gphi)*dt)@v
        v/=np.trace(v.reshape(2,2))
    r=v.reshape(2,2); return r/np.trace(r)
def tdist(a,b):
    return 0.5*np.abs(np.linalg.eigvalsh((a-b+(a-b).conj().T)/2)).sum()
def dominant(phi,G,gphi,**kw):
    Lm=liou(*loop(phi,**kw),G,gphi); w,vl,vr=eig(Lm,left=True)
    i=np.argmax(w.real); return w,vl,vr,i
def amp(G,gphi,phi=0.0,h=1e-5,**kw):
    # first-order adiabatic response: delta rho = (dphi/dt) * sum_{k!=0} r_k <l_k|d r0>/(lam0-lam_k)
    def r0(p):
        w,vl,vr,i=dominant(p,G,gphi,**kw); r=vr[:,i].reshape(2,2); return r/np.trace(r)
    w,vl,vr,i=dominant(phi,G,gphi,**kw)
    dr=(r0(phi+h)-r0(phi-h)).reshape(-1)/(2*h)
    # solve (L - lam0) x = dr projected; normalized-state correction
    Lm=liou(*loop(phi,**kw),G,gphi)
    r=vr[:,i]/np.trace(vr[:,i].reshape(2,2)); 
    # x = (L-lam0)^+ dr on complement, then remove component along r via trace
    M=Lm-w[i]*np.eye(4)
    x=np.linalg.lstsq(M,dr,rcond=None)[0]
    x=x-np.trace(x.reshape(2,2))*r
    # chirality at first order: rho_s = r0 + s*(2pi/T)*x  -> C = (2pi/T)*||x||_1
    X=x.reshape(2,2); return 2*np.pi*0.5*np.abs(np.linalg.eigvals(X)).sum(), w
def amp2(G,gphi,phi=0.0,h=1e-5,**kw):
    def r0(p):
        w,vl,vr,i=dominant(p,G,gphi,**kw); r=vr[:,i]; return r/np.trace(r.reshape(2,2))
    Lm=liou(*loop(phi,**kw),G,gphi); w,vr=np.linalg.eig(Lm); vl=np.linalg.inv(vr)
    i=np.argmax(w.real); dr=(r0(phi+h)-r0(phi-h))/(2*h)
    x=sum(vr[:,k]*(vl[k]@dr)/(w[i]-w[k]) for k in range(4) if k!=i)  # (lam0 - L)^-1 on complement... sign irrelevant
    r=r0(phi); x=x-np.trace(x.reshape(2,2))*r
    X=x.reshape(2,2); X=(X+X.conj().T)/2
    return 2*np.pi*np.abs(np.linalg.eigvalsh(X)).sum(), w
