import numpy as np, json
# Gate 2(b): is the exponential in g_N (coupling) an independent, measurable physical rate?
# Candidate: hi->lo switching rate k_hl measured from the real-Delta exact spectrum just inside the wedge.
# Compare d ln g_N / dN with -s_hl(Delta) at Delta = Re Delta_EP(N), and with the EP imaginary-part slope.
g=np.array(json.load(open('../data/gN.json'))); N=g[:,0]; gN=g[:,2]
sl=-np.diff(np.log(gN))/np.diff(N); print('g_N local decay rate per atom:',np.round(sl,4))
EP=dict((int(r[0]),r) for r in json.load(open('../data/cep_all2.json')))
Im=np.array([EP[int(n)][2] for n in N]); print('Im D_EP local slope:',np.round(-np.diff(np.log(Im))/np.diff(N),4))
print('ratio Im/g (should be const = 2/|dslope| if 2x2 picture):',np.round(Im/gN,3))
