import json,numpy as np,sys
E7={'m2':1.28558,'m3':1.87939,'m4':1.96962}
for L in (28,32):
    R=json.load(open(f'../data/thscan_L{L}.json'))
    print("L",L)
    for r in R:
        ep=np.array(r['ep']);wp=np.array(r['wp']);e0=np.array(r['e0']);w0=np.array(r['w0'])
        m1=ep[0]
        if m1*L<6: continue
        pi=ep[wp>1e-6*wp[0]]/m1; z=e0[w0>1e-6*w0[0]]/m1
        print(f" th={r['th']:4.0f} t={r['t']:.2f} m1={m1:.3f} m1L={m1*L:5.1f}  m2={z[0]:.4f}({z[0]/E7['m2']-1:+.3f}) m3={pi[1]:.4f}({pi[1]/E7['m3']-1:+.3f}) m4={z[1]:.4f}({z[1]/E7['m4']-1:+.3f})  wts pi:{wp[1]/wp[0]:.3f} 0:{w0[1]/w0[0]:.3f}")
