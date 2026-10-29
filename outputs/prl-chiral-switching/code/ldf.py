import numpy as np, json
# exact quasi-potential along n from P_ss(k): Phi_N(n) = -(ln P(k) - ln P(k_hi))/N ; extrapolate to N->inf with a+b/N
from mf import roots
V=10.;Om=2.0
P={N:json.load(open(f'../data/Pk_N{N}.json')) for N in (30,40,50)}
for D in ['3.7','3.9']:
    r=roots(float(D),Om,V); ng=np.linspace(0.15,0.45,13); rows=[]
    for N in (30,40,50):
        p=np.array(P[N][D]); k=np.arange(N+1)/N; lp=np.log(np.clip(p,1e-300,None))
        ihi=np.argmax(np.where(k>r[1],lp,-np.inf)); phi=-(lp-lp[ihi])/N
        rows.append(np.interp(ng,k,phi))
    rows=np.array(rows); Ns=np.array([30,40,50.])
    ext=[np.polyfit(1/Ns,rows[:,j],1)[1] for j in range(len(ng))]
    print(D,'MF lo,saddle,hi=',np.round(r,3)); print('   n   :',np.round(ng,3)); print('  Phi_inf:',np.round(ext,4))
