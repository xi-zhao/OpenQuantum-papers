import numpy as np
from scipy.optimize import brentq
# Mean-field 3/4-level ladder model of Xie et al. (Sci. Bull. 2025), units 2pi*MHz.
P3=dict(Oc=12.21,g21=27.22,g32=0.29,g31=0.016,g3=20.31,alpha=-273.89)
def ops(d):
    E=np.eye(d); return lambda i,j: np.outer(E[i],E[j])
def liou(H,Ls):
    d=H.shape[0]; I=np.eye(d)
    L=-1j*(np.kron(I,H)-np.kron(H.T,I))
    for A in Ls:
        AdA=A.conj().T@A
        L+=np.kron(A.conj(),A)-0.5*(np.kron(I,AdA)+np.kron(AdA.T,I))
    return L
def ss_lin(H,Ls):
    d=H.shape[0]; L=liou(H,Ls)
    M=np.vstack([L,np.eye(d).flatten(order='F')[None,:]]); b=np.zeros(d*d+1,complex); b[-1]=1
    r=np.linalg.lstsq(M,b,rcond=None)[0]; return r.reshape(d,d,order='F')
def H3(Op,Dt,p): 
    return np.array([[0,Op/2,0],[Op/2,0,p['Oc']/2],[0,p['Oc']/2,-Dt]],complex)
def L3(p):
    s=ops(3); return [np.sqrt(p['g21'])*s(0,1),np.sqrt(p['g32'])*s(1,2),np.sqrt(p['g31'])*s(0,2),np.sqrt(p['g3'])*s(2,2)]
def F3(x,Dc,Op,p=P3):
    r=ss_lin(H3(Op,Dc-p['alpha']*x,p),L3(p)); return r[2,2].real-x
def roots3(Dc,Op,p=P3,xs=np.concatenate([np.linspace(0,0.02,400),np.linspace(0.02,0.5,400)[1:]])):
    f=np.array([F3(x,Dc,Op,p) for x in xs]); out=[]
    for i in np.flatnonzero(np.sign(f[:-1])!=np.sign(f[1:])): out.append(brentq(F3,xs[i],xs[i+1],args=(Dc,Op,p)))
    return np.array(out)
from scipy.integrate import solve_ivp
def mf_rhs3(t,y,path,p=P3):
    r=(y[:9]+1j*y[9:]).reshape(3,3,order='F')
    Dc,Op=path(t); H=H3(Op,Dc-p['alpha']*r[2,2].real,p)
    dr=(liou(H,L3(p))@r.flatten(order='F'))
    return np.concatenate([dr.real,dr.imag])
def ss_branch(Dc,Op,x,p=P3):
    r=ss_lin(H3(Op,Dc-p['alpha']*x,p),L3(p)).flatten(order='F'); return np.concatenate([r.real,r.imag])
def run3(path,T,y0):
    s=solve_ivp(mf_rhs3,[0,T],y0,args=(path,),rtol=1e-8,atol=1e-10,method='LSODA',dense_output=True)
    return s
def cusp3(p=P3):
    # scan Op upward along the wedge until bistable interval closes
    best=None
    for Op in np.linspace(8,40,129):
        Ds=np.linspace(-120,-20,201); b=[D for D in Ds if len(roots3(D,Op,p))>=3]
        if b: best=(Op,min(b),max(b))
    return best
