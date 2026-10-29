import numpy as np
from mf import roots
V=10.
C=(4.55,3.0); R=(1.45,1.7)
def pt(th): return C[0]+R[0]*np.cos(th), C[1]+R[1]*np.sin(th)
th=np.linspace(0,2*np.pi,1441)[:-1]
bi=np.array([len(roots(*pt(t),V))==3 for t in th])
if __name__=="__main__":
    i=np.flatnonzero(bi); print(len(i), th[i].min(),th[i].max()); 
    edges=np.flatnonzero(np.diff(bi.astype(int))!=0); print('edges',th[edges])
    tm=th[i].mean(); print('mid',tm,pt(tm),roots(*pt(tm),V))
