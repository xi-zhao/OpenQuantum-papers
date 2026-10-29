import numpy as np, json, os
from figstyle import plt, C, setup, tag
setup()
fig,ax=plt.subplots(1,2,figsize=(7.0,2.4))
# (a) Gate A: switching action hi->saddle, exact-H WKB (Lin) vs gMAM-LNA vs exact s_hl
ex={3.62:0.2500,3.65:0.2400,3.7:0.2222,3.75:0.2040,3.8:0.1860,3.9:0.1519,4.0:0.1209,4.2:0.0669}
lin={3.8:0.18588,3.9:0.15358,4.0:0.12337,4.2:0.06878}
gm={3.7:0.2741,3.8:0.2270,4.0:0.1453,4.2:0.0774}
D=sorted(ex); ax[0].plot(D,[ex[d] for d in D],'-',color=C['ink2'],lw=1.0,label=r'exact $s_{\rm hl}$ ($N=40$–$60$)')
ax[0].plot(sorted(lin),[lin[d] for d in sorted(lin)],'o',color=C['blue'],ms=4.5,label='instanton (this work)')
ax[0].plot(sorted(gm),[gm[d] for d in sorted(gm)],'^',color=C['orange'],ms=4,mfc='none',label='Gaussian-noise gMAM')
ax[0].plot([3.62],[0.0517],'v',color=C['aqua'],ms=4,mfc='none',label='incoherent rate eq.')
ax[0].set_xlabel(r'$\Delta/\gamma$'); ax[0].set_ylabel('switching exponent'); ax[0].legend(loc='upper right',fontsize=6.2); tag(ax[0],'(a)',x=0.05,y=0.08)
# (b) Gate A': quasipotential Phi(n) at Delta=3.7: WKB vs exact finite-N
P={N:json.load(open(f'../data/Pk_N{N}.json'))['3.7'] for N in (30,40,50)}
cl={30:C['aqua'],40:C['blue'],50:C['orange']}
for N in (30,40,50):
    p=np.array(P[N]); lp=np.log(np.clip(p,1e-300,None)); k=np.arange(N+1)/N
    m=k>0.17; ax[1].plot(k[m],(-(lp-lp.max())/N)[m],'-',color=cl[N],lw=1.0,label=r'exact, $N=%d$'%N)
w=[(0.25,0.06224),(0.30,0.02918),(0.35,0.00882)]; ax[1].plot([a for a,b in w],[b for a,b in w],'o',color=C['ink'],ms=4.5,label=r'instanton $\Phi(n)$, a priori')
ax[1].set_xlim(0.17,0.47); ax[1].set_ylim(-0.002,0.13); ax[1].set_xlabel(r'$n$'); ax[1].set_ylabel(r'$\Phi(n)$  at $\Delta=3.7$'); ax[1].legend(loc='upper right',fontsize=6.2); tag(ax[1],'(b)',x=0.05,y=0.08)
plt.tight_layout(); plt.savefig('../figs/rev_fig3.pdf'); plt.savefig('../figs/rev_fig3.png',dpi=200)
