"""Exploration (not a test): decompose the avoided-crossing coupling into cross matrix elements.
   L(D) = L0 + D*L1 is linear in Delta.  Adiabatic modes a (hi intrabasin) and b (ghost) at D_off away from the crossing:
   X_ab = <l_a|L1|r_b>, X_ba = <l_b|L1|r_a> (biorthonormal l_i.r_i=1).  2x2 theory: g^2 = (la-lb)^2 X_ab X_ba / (sa-sb)^2.
   Physical normalisation: r scaled so that max_k |r_kk| (excitation-number diagonal) = 1."""
import numpy as np, scipy.sparse.linalg as sla, json, sys
from ex import liou, ops
V=10.;Om=2.0
def diagk(N,vec,kind):
    Ne,Jx,B,idx,D=ops(N); d=Ne.shape[0]; full=np.zeros(d*d,complex); full[idx]=vec
    M=full.reshape(d,d,order='F')
    if kind=='left': M=M.T            # left vector l: Tr(l^T-ish rho) pairing -> operator X with X_ji
    ne=np.round(np.diag(Ne).real).astype(int); P=np.zeros(N+1,complex); A=np.zeros(N+1)
    for i in range(d): P[ne[i]]+=M[i,i]; A[ne[i]]=max(A[ne[i]],abs(M[i,i]))
    return P,A
def modes(N,D,target=-0.75):
    L=liou(N,D,Om,V)
    w,R=sla.eigs(L,k=6,sigma=target,tol=1e-13); wl,Lv=sla.eigs(L.T.tocsc(),k=6,sigma=target,tol=1e-13)
    sel=[i for i in np.argsort(abs(w-target)) if abs(w[i].imag)<1e-6][:2]
    out=[]
    for i in sel:
        j=np.argmin(abs(wl-w[i])); r=R[:,i]; l=Lv[:,j]
        P,_=diagk(N,r,'right'); s=P[np.argmax(abs(P))]; r=r/s; l=l/(l@r)
        out.append((w[i].real,r,l))
    return out
def L1(N):
    a=liou(N,1.0,Om,V); b=liou(N,0.0,Om,V); return (a-b).tocsc()
if __name__=="__main__":
    av=dict((int(r[0]),r) for r in json.load(open('../data/gN.json')))
    res=[]
    for N in [24,32,40,48]:
        Dm=av[N][1]; g=av[N][2]; Lp=L1(N)
        for side in (-1,+1):
            D=Dm+side*max(12*g,0.02)
            m=modes(N,D)
            # identify a (hi intrabasin) vs b (ghost) by where |r| lives: ghost has weight at n<0.25
            def lowfrac(r):
                P,A=diagk(N,r,'right'); k=np.arange(N+1)/N; return np.abs(P[k<0.25]).sum()/np.abs(P).sum()
            m.sort(key=lambda t: lowfrac(t[1]))     # a first (less low weight), b second
            (la,ra,lla),(lb,rb,llb)=m
            Xab=lla@(Lp@rb); Xba=llb@(Lp@ra)
            res.append(dict(N=N,D=D,side=side,la=la,lb=lb,Xab=[Xab.real,Xab.imag],Xba=[Xba.real,Xba.imag],g=g,
                            lowa=lowfrac(ra),lowb=lowfrac(rb)))
            print('N=%d side=%+d D=%.4f la=%.4f lb=%.4f |Xab|=%.3e |Xba|=%.3e |Xab Xba|=%.3e  g^2=%.3e  low(a,b)=(%.2f,%.2f)'%(
                N,side,D,la,lb,abs(Xab),abs(Xba),abs(Xab*Xba),g*g,lowfrac(ra),lowfrac(rb)),flush=True)
    json.dump(res,open('../data/xel.json','w'))
