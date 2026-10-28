"""Large-N WKB (coherent-state / Doi-Peliti) Hamiltonian for the permutation-symmetric Liouvillian
   L = -i[-Delta N_e + Om J_x + (V/2N) N_e(N_e-1), .] + sum_i D[sigma_i^-].
   Vectorised rho in the symmetric sector = N bosons in 4 single-atom Liouville modes; large N gives
     h(X,rho) = Re Tr[X l0(rho)] + 2V Re(A) Im(A),   A = Tr(X n rho),
   rho: single-atom 'position' (Hermitian 2x2), X: conjugate 'momentum' (Hermitian 2x2), pairing Tr(X rho).
   X = I reproduces the mean-field Bloch equations; rare-event exponents are actions  int Tr(X d rho)."""
import numpy as np, sympy as sp
q=sp.symbols('rgg ree u v',real=True); p=sp.symbols('pgg pee pu pv',real=True)
D_,O_,V_=sp.symbols('Delta Omega V',real=True)
rgg,ree,u,v=q; pgg,pee,pu,pv=p
I=sp.I
rho=sp.Matrix([[rgg,u-I*v],[u+I*v,ree]])            # basis (g,e); rho_eg = u+iv
X=sp.Matrix([[pgg,(pu-I*pv)/2],[(pu+I*pv)/2,pee]])   # Tr(X rho) = pgg rgg + pee ree + pu u + pv v
n=sp.Matrix([[0,0],[0,1]]); sx=sp.Matrix([[0,1],[1,0]]); sm=sp.Matrix([[0,1],[0,0]])
H0=-D_*n+O_/2*sx
l0=-I*(H0*rho-rho*H0)+sm*rho*sm.H-(sm.H*sm*rho+rho*sm.H*sm)/2
A=sp.expand((X*n*rho).trace())
h=sp.expand(sp.re(sp.expand((X*l0).trace()))+2*V_*sp.re(A)*sp.im(A))
z=list(q)+list(p)
grad=[sp.diff(h,s) for s in z]
f_h=sp.lambdify([z,D_,O_,V_],h,'numpy')
f_g=sp.lambdify([z,D_,O_,V_],grad,'numpy')
Hess=sp.Matrix(grad).jacobian(z); f_H=sp.lambdify([z,D_,O_,V_],Hess,'numpy')
def flow(t,y,D,Om,V):
    g=np.array(f_g(y,D,Om,V),float)
    return np.concatenate([g[4:],-g[:4]])      # dq/dt = dh/dp, dp/dt = -dh/dq
def lin(y,D,Om,V):
    Hm=np.array(f_H(y,D,Om,V),float); return np.vstack([Hm[4:],-Hm[:4]])
def mf_point(D,Om,V,n0):
    Dt=D-V*n0; c=(-1j*Om/2*(1-2*n0))/(0.5-1j*Dt)
    return np.array([1-n0,n0,c.real,c.imag,1,1,0,0],float)
if __name__=="__main__":
    from mf import roots
    D,Om,V=4.0,2.0,10.
    for n0 in roots(D,Om,V):
        y=mf_point(D,Om,V,n0)
        print('n=%.4f  h=%.2e  |flow|=%.2e  eig='%(n0,f_h(y,D,Om,V),np.abs(flow(0,y,D,Om,V)).max()),np.round(np.sort_complex(np.linalg.eigvals(lin(y,D,Om,V))),4))
