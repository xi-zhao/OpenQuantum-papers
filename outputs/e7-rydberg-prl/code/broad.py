import json,numpy as np
E7=np.array([1,1.28558,1.87939,1.96962,2.53209,2.87939,3.70167])
def S(e,w,om,sig): return (w[None,:]*np.exp(-(om[:,None]-e[None,:])**2/(2*sig**2))/np.sqrt(2*np.pi)/sig).sum(1)
for f in ('lscan_th110_t0.2','lscan_th110_t0.3'):
    R=[r for r in json.load(open(f'../data/{f}.json')) if r['L']>=26]
    m1=np.mean([r['ep'][0] for r in R[-3:]])
    om=np.linspace(0.8,4.0,321)*m1
    for sig in (0.03,0.06):
        Sp=np.mean([S(np.array(r['ep']),np.array(r['wp']),om,sig*m1) for r in R],0)
        S0=np.mean([S(np.array(r['e0']),np.array(r['w0']),om,sig*m1) for r in R],0)
        print(f,sig)
        for nm,Sx in (('pi',Sp),('0',S0)):
            # local maxima above 2.2
            x=om/m1; i=[k for k in range(1,len(x)-1) if Sx[k]>Sx[k-1] and Sx[k]>Sx[k+1] and Sx[k]>1e-3*Sx.max()]
            print(' ',nm,' '.join(f"{x[k]:.3f}({Sx[k]/Sx.max():.1e})" for k in i))
