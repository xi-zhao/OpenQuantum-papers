import numpy as np, json
import scipy.sparse as sp, scipy.sparse.linalg as sla
# Gate 2 core test: does the fold normal form, with NO fitted parameter, predict Re(Delta_EP)-Delta_f ~ c N^{-1/2}?
# Normal form (from MF centre manifold, data/normalform.json): dx = (a mu + b x^2) dt + sqrt(Q/N) dW, mu = Delta - Delta_f.
# Fokker-Planck generator on a finite x-interval, with the high branch at x>0 replaced by its finite intrabasin
# relaxation rate kappa: here we only test the *scaling form* of the ghost-escape spectrum, i.e. eigenvalues of
# the FP operator in scaled variables  x = N^{-1/3} X ,  mu = N^{-2/3} M  ->  lambda = N^{-1/3} Lambda(M).
nf=json.load(open('../data/normalform.json')); a,b,Q=nf['a'],nf['b'],nf['Q']
def fp_eigs(M,L=6.0,n=1200,k=4):
    X=np.linspace(-L,L,n); h=X[1]-X[0]; drift=a*M+b*X**2; Dif=Q/2
    # generator for densities: dp/dt = -d(drift p)/dX + Dif d2p/dX2 ; absorbing at -L (escaped), reflecting-ish at +L
    main=-2*Dif/h**2*np.ones(n); up=Dif/h**2-drift[1:]/(2*h); lo=Dif/h**2+drift[:-1]/(2*h)
    G=sp.diags([main,lo,up],[0,-1,1],format='csc')
    w=sla.eigs(G,k=k,sigma=0.0,return_eigenvectors=False); return np.sort(w.real)[::-1]
for M in [-3,-1,0,1,3]:
    print('M=%+d  Lambda_1..3 ='%M,np.round(fp_eigs(M)[:3],4))
