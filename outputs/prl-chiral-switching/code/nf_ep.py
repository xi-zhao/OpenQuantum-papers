import numpy as np, json
from scipy.optimize import brentq
from nf_spec import fp_eigs
nf=json.load(open('../data/normalform.json'))
kappa=0.7506   # MF high-branch intrabasin relaxation rate at the fold (|lamJ| computed earlier)
# ghost-escape mode of the normal form: lambda_ghost(N,mu) = N^{-1/3} * Lambda1(M),  M = N^{2/3} mu   (mu<0: outside wedge)
# crossing condition lambda_ghost = -kappa  ->  Lambda1(M*) = -kappa N^{1/3}  ->  mu* = N^{-2/3} M*
Ms=np.linspace(-60,3,253); L1=np.array([fp_eigs(M,L=12.0,n=2400)[0] for M in Ms])
json.dump(dict(Ms=Ms.tolist(),L1=L1.tolist()),open('../data/nf_L1.json','w'))
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
for N in [20,30,40,50,58]:
    tgt=-kappa*N**(1/3)
    i=np.flatnonzero(np.diff(np.sign(L1-tgt)))
    if len(i)==0: print(N,'no crossing in M range'); continue
    Mx=np.interp(tgt,L1[i[0]:i[0]+2][::-1] if L1[i[0]]>L1[i[0]+1] else L1[i[0]:i[0]+2],Ms[i[0]:i[0]+2][::-1] if L1[i[0]]>L1[i[0]+1] else Ms[i[0]:i[0]+2])
    mu=Mx*N**(-2/3); act=EP[N][1]-nf['Df']
    print(N,'normal-form crossing mu*=%.4f   exact ReD_EP-Df=%.4f   ratio=%.2f'%(mu,act,mu/act))
print('large-|M| asymptote check (deterministic ghost: Lambda ~ -c sqrt(-M)?):',np.round(L1[:5],2),np.round(L1[:5]/np.sqrt(-Ms[:5]),3))
