import numpy as np, json
# Test of a mechanism: the avoided-crossing coupling g_N = gap_min/2 should scale like the tunnelling
# amplitude between the ghost of the low branch and the high branch, i.e. g_N ~ sqrt(P_ss(ghost)/P_ss(high)),
# with P_ss the exact steady-state excitation-number distribution at the same real Delta.
av=json.load(open('../data/g2_avoid.json'))
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
N=np.array([r[0] for r in EP.values()],float); Im=np.array([r[2] for r in EP.values()])
for M in [40]:
    P=json.load(open(f'../data/phi_N{M}.json'))
    for D,p in P.items():
        p=np.array(p); k=np.arange(M+1)/M; phi=-np.log(np.clip(p,1e-300,None))/M
        print(M,D,'Phi(n) min at n=%.3f'%k[np.argmin(phi)],' Phi coarse',np.round(phi[::4],3))
