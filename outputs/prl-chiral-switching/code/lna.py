import numpy as np
from mf import roots
from mfdyn import rhs, steady
V=10.
# single-atom operators: A1=|e><e|, A2=sigma_x/2, A3=-sigma_y/2  (so <A>=(n, Re rho_eg, Im rho_eg))
P=np.array([[0,0],[0,1]],complex)          # basis (g,e): |e><e|
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]],complex)
A=[P,sx/2,sy/2]
sm=np.array([[0,1],[0,0]],complex)         # |g><e| in basis (g,e)
def L1(rho,Dt,Om):
    H=-Dt*P+Om/2*sx
    return -1j*(H@rho-rho@H)+sm@rho@sm.conj().T-0.5*(sm.conj().T@sm@rho+rho@sm.conj().T@sm)
def state(n,c):
    return np.array([[1-n,np.conj(c)],[c,n]],complex)   # rho_eg = <e|rho|g> = rho[1,0]
def D1(D,Om,n,c):
    """single-atom local-noise matrix: dS/dt - M S - S M^T with frozen detuning (see notes)"""
    Dt=D-V*n; rho=state(n,c)
    mean=lambda r: np.array([np.trace(r@a).real for a in A])
    S=lambda r: np.array([[np.trace(r@(a@b+b@a)/2).real for b in A] for a in A])-np.outer(mean(r),mean(r))
    # linear map for means: d<A>/dt = M <A> + b  (frozen Dt): get M by finite differences on the 3 coordinates
    def f(y): r=state(y[0],y[1]+1j*y[2]); return mean(L1(r,Dt,Om))
    y0=mean(rho); h=1e-6; M=np.array([(f(y0+h*e)-f(y0-h*e))/(2*h) for e in np.eye(3)]).T
    dS=(lambda e: (S(rho+e*L1(rho,Dt,Om))-S(rho-e*L1(rho,Dt,Om)))/(2*e))(1e-6)
    return dS-M@S(rho)-S(rho)@M.T, S(rho)
def J(D,Om,y):
    f=lambda y: np.array(rhs(0,y,lambda t:(D,Om,V))); h=1e-7
    return np.array([(f(y+h*e)-f(y-h*e))/(2*h) for e in np.eye(3)]).T
if __name__=="__main__":
    from scipy.linalg import solve_continuous_lyapunov as lyap
    for D,Om,br in [(4.9,2.0,0),(3.0,2.0,0),(2.0,1.0,0)]:
        y=np.array(steady(D,Om,V,br)); Dm,S=D1(D,Om,y[0],y[1]+1j*y[2]); Jm=J(D,Om,y)
        C=lyap(Jm,-Dm); print(D,Om,'n=%.4f'%y[0],'D1 eig',np.round(np.linalg.eigvalsh((Dm+Dm.T)/2),4),' LNA N*Var(n)=%.4f'%C[0,0])
