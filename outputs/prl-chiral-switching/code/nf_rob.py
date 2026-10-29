import numpy as np, json, scipy.sparse as sp, scipy.sparse.linalg as sla
nf=json.load(open('../data/normalform.json')); a,b,Q=nf['a'],nf['b'],nf['Q']
def L1(M,L,n,bc):
    X=np.linspace(-L,L,n); h=X[1]-X[0]; dr=a*M+b*X**2; Df=Q/2
    main=-2*Df/h**2*np.ones(n); up=Df/h**2-dr[1:]/(2*h); lo=Df/h**2+dr[:-1]/(2*h)
    G=sp.lil_matrix(sp.diags([main,lo,up],[0,-1,1]))
    if bc=='refl':   # zero-flux at +L
        G[n-1,n-1]=-Df/h**2+dr[n-1]/(2*h)*0-dr[n-1]/h*0
    w=sla.eigs(G.tocsc(),k=3,sigma=0.0,return_eigenvectors=False); return np.sort(w.real)[::-1][0]
for M in [-1,-3,-6]:
    print(M,[round(L1(M,L,n,bc),4) for L,n,bc in [(6,1200,'abs'),(12,2400,'abs'),(6,2400,'abs'),(12,2400,'refl'),(20,4000,'abs')]])
