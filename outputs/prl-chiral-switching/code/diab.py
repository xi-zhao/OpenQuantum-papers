"""Biorthogonal diabatic two-mode reduction at the avoided crossing (exploration, Omega=2).
   L(D)=L0+D*L1 exactly. At D=Dmin: right R (2 cols), left Lt (biorthonormal, Lt^T R = I).
   M0=Lt^T L0 R, M1=Lt^T L1 R ; diabatic basis T diagonalises M1 -> constant coupling V_ab,V_ba.
   Gauge: diabatic right vectors scaled so max |population| = 1 ('density gauge')."""
import numpy as np, scipy.sparse.linalg as sla, json, sys
from ex import liou, ops
V=10.; Om=float(sys.argv[2]) if len(sys.argv)>2 else 2.0
def blocks(N):
    Ne,Jx,B,idx,D=ops(N); d=Ne.shape[0]; ne=np.round(np.diag(Ne).real).astype(int); return d,ne,idx
def profile(N,vec):
    d,ne,idx=blocks(N); full=np.zeros(d*d,complex); full[idx]=vec; M=full.reshape(d,d,order='F')
    P=np.zeros(N+1,complex); A=np.zeros(N+1)
    for i in range(d): P[ne[i]]+=M[i,i]; A[ne[i]]=max(A[ne[i]],abs(M[i,i]))
    # all entries (incl. coherences) with row-excitation k : max abs
    Afull=np.zeros(N+1); absM=np.abs(M)
    for i in range(d): Afull[ne[i]]=max(Afull[ne[i]],absM[i].max(),absM[:,i].max())
    return P,A,Afull
def run(N,Dmin,target=-0.76):
    L0=liou(N,0.0,Om,V); L1=(liou(N,1.0,Om,V)-L0).tocsc(); L=(L0+Dmin*L1).tocsc()
    w,R=sla.eigs(L,k=6,sigma=target,tol=1e-13); wl,Lv=sla.eigs(L.T.tocsc(),k=6,sigma=target,tol=1e-13)
    sel=[i for i in np.argsort(abs(w-target)) if abs(w[i].imag)<1e-6][:2]
    R=R[:,sel]; lam=w[sel]
    jj=[int(np.argmin(abs(wl-l))) for l in lam]; Lv=Lv[:,jj]
    Lt=Lv@np.linalg.inv(R.T@Lv)
    M0=Lt.T@(L0@R); M1=Lt.T@(L1@R)
    s,T=np.linalg.eig(M1)
    Rd=R@T; Ld=Lt@np.linalg.inv(T).T
    # density gauge
    c=np.array([np.max(np.abs(profile(N,Rd[:,i])[0])) for i in range(2)])
    Rd=Rd/c; Ld=Ld*c
    H0=Ld.T@(L0@Rd); H1=Ld.T@(L1@Rd)
    e=np.diag(H0)+Dmin*np.diag(H1)
    # order: a = hi intrabasin (small |slope|), b = ghost (large |slope|)
    o=np.argsort(np.abs(np.diag(H1).real)); P=np.eye(2)[:,o]
    H0=P.T@H0@P; H1=P.T@H1@P; Rd=Rd@P; Ld=Ld@P; e=e[o]
    Vab,Vba=H0[0,1],H0[1,0]
    prof={}
    for nm,vec in [('ra',Rd[:,0]),('rb',Rd[:,1]),('la',Ld[:,0]),('lb',Ld[:,1])]:
        Pp,A,Af=profile(N,vec); prof[nm]=dict(pop=np.abs(Pp).tolist(),maxdiag=A.tolist(),maxall=Af.tolist())
    return dict(N=N,Dmin=Dmin,lam=[x.real for x in lam],slopes=[complex(x).real for x in np.diag(H1)],e=[complex(x).real for x in e],
                Vab=[Vab.real,Vab.imag],Vba=[Vba.real,Vba.imag],offdiagH1=abs(H1[0,1])+abs(H1[1,0]),
                gap2x2=float(2*np.sqrt(abs(Vab*Vba))),prof=prof)
if __name__=="__main__":
    N=int(sys.argv[1])
    gN=dict((int(r[0]),r) for r in json.load(open('../data/gN.json')))
    o=run(N,gN[N][1]); o['g_exact']=gN[N][2]
    json.dump(o,open(f'../data/diab_N{N}_Om{Om}.json','w'))
    print('N=%d Dmin=%.5f lam=%s slopes=%s |Vab|=%.3e |Vba|=%.3e 2x2 half-gap=%.3e exact g=%.3e ratio=%.4f offdiagH1=%.1e'%(
        N,o['Dmin'],np.round(o['lam'],5),np.round(o['slopes'],4),np.hypot(*o['Vab']),np.hypot(*o['Vba']),o['gap2x2']/2,o['g_exact'],o['gap2x2']/2/o['g_exact'],o['offdiagH1']),flush=True)
