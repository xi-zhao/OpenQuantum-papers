import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from ex import liou, ops, obs_vec
def unvec(N,v):
    Ne,Jx,B,idx,D=ops(N); d=Ne.shape[0]
    full=np.zeros(d*d,complex); full[idx]=v; return full.reshape(d,d,order='F')
def ems(N,D,Om,V):
    """stationary state and two extreme metastable states (Macieszczak et al. PRL 116, 240404)"""
    L=liou(N,D,Om,V); tr,nv=obs_vec(N)
    w,R=sla.eigs(L,k=2,sigma=1e-6); o=np.argsort(-w.real); w=w[o]; R=R[:,o]
    wl,Lv=sla.eigs(L.T.tocsc(),k=2,sigma=1e-6); o=np.argsort(-wl.real); Lv=Lv[:,o]
    rss=R[:,0]/(tr@R[:,0]); r1=R[:,1]; l1=Lv[:,1]; r1=r1/(l1@r1)
    # l1 as operator: Tr(l1^dag rho)= l1(vec)·rho(vec) -> operator X with X_ij such that sum X_ji rho_ij
    X=unvec(N,l1).T
    k=np.argmax(abs(np.diag(X))); ph=np.exp(-1j*np.angle(X[k,k])); X=X*ph; r1=r1/ph
    herm_err=np.abs(X-X.conj().T).max()/np.abs(X).max()
    X=(X+X.conj().T)/2; ev=np.linalg.eigvalsh(X)
    rho_a=rss+ev.max()*r1; rho_b=rss+ev.min()*r1
    na,nb=(nv@rho_a).real/N,(nv@rho_b).real/N
    return w[1].real, rss, (rho_a,rho_b) if na>nb else (rho_b,rho_a)
def evolve(N,rho0,path,T,steps):
    """BDF2 (A-stable, stiffly accurate); first step backward Euler. L evaluated at new time."""
    tr,nv=obs_vec(N); dt=T/steps; v=rho0.copy(); vold=None
    I=sp.identity(len(v),format='csc')
    traj=[]
    for k in range(steps):
        D,Om,V=path((k+1)*dt); L=liou(N,D,Om,V)
        if vold is None:
            vn=sla.spsolve((I-dt*L).tocsc(),v)
        else:
            vn=sla.splu((I-(2/3)*dt*L).tocsc()).solve((4*v-vold)/3)
        vold,v=v,vn
        traj.append((nv@v).real/N)
    return v,np.array(traj)
