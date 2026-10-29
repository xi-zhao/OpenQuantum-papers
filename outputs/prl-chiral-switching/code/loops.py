import numpy as np
from mf import roots
V=10.
TH0=4.089328813438765
P0=(3.7038915919767983,1.6194358071727704)
LOOPS={'A':dict(C=(4.7,3.3),R=(abs(3.7038915919767983-4.7)/np.sqrt(1-((1.6194358071727704-3.3)/2.0)**2),2.0)),      # encloses upper cusp (EP3)
       'B':dict(C=(3.70389,1.2),R=(1.3,0.4194358))} # does not enclose any cusp
def start_angle(L):
    C,R=L['C'],L['R']; return np.arctan2((P0[1]-C[1])/R[1],(P0[0]-C[0])/R[0])
def mkpath(name,T,sgn):
    L=LOOPS[name]; C,R=L['C'],L['R']; t0=start_angle(L)
    def path(t):
        th=t0+sgn*2*np.pi*t/T
        return C[0]+R[0]*np.cos(th), C[1]+R[1]*np.sin(th), V
    return path
if __name__=="__main__":
    for k in LOOPS:
        p=mkpath(k,1.,1); print(k,p(0)[:2],P0)
        s=np.linspace(0,1,2001); nb=np.array([len(roots(*p(x)[:2],V)) for x in s])
        print(''.join('3' if x==3 else '.' for x in nb[::20]))
