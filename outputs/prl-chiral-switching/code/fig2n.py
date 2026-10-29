import numpy as np, json
from figstyle import plt, C, setup, tag
setup()
fig,ax=plt.subplots(1,2,figsize=(7.0,2.4))
# (a) mode profiles in excitation number at N=48 near the crossing (log scale)
pr=json.load(open('../data/fig_profile.json'))['48']; N=48; k=np.arange(N+1)/N
ax[0].semilogy(k,np.clip(pr['Pss'],1e-14,None),color=C['ink2'],lw=1.0,label=r'$P_{\rm ss}(n)$')
labs=['high-branch relaxation','ghost escape']; cl=[C['blue'],C['orange']]
m=sorted(pr['modes'],key=lambda x: -np.array(x['A'])[:int(0.25*N)].sum())[::-1]
for mm,lb,c in zip(m,labs,cl): ax[0].semilogy(k,np.clip(mm['A'],1e-14,None),color=c,lw=1.4,label=lb+r' ($\lambda=%.2f$)'%mm['lam'])
nf=json.load(open('../data/normalform.json'))['nf']; ax[0].axvline(nf,color=C['ink2'],lw=0.6,ls='--'); ax[0].text(nf+0.01,2e-12,r'$n_f$',fontsize=7,color=C['ink2'])
ax[0].set_xlim(0,0.7); ax[0].set_ylim(1e-12,3); ax[0].set_xlabel(r'$n=k/N$'); ax[0].set_ylabel('weight (normalised)'); ax[0].legend(loc='lower right',fontsize=6.2); tag(ax[0],'(a)')
# (b) EP location: Re shift vs fold (power) and Im (exponential)
ep=np.array(json.load(open('../data/cep_all2.json'))); Nn=ep[:,0]; Df=3.59963689
ax2=ax[1]
ax2.semilogy(Nn,ep[:,2],'o',color=C['orange'],ms=3.2,label=r'Im$\,\Delta_{\rm EP}$')
ax2.semilogy(Nn,Df-ep[:,1],'s',color=C['blue'],ms=3.0,label=r'$\Delta_f-$Re$\,\Delta_{\rm EP}$')
x=np.linspace(8,60,50); s=Nn>=30
k1=np.polyfit(Nn[s],np.log(ep[s,2]),1); ax2.semilogy(x,np.exp(np.polyval(k1,x)),'-',color=C['orange'],lw=0.8)
q=np.polyfit(np.log(Nn[Nn>=24]),np.log(Df-ep[Nn>=24,1]),1); ax2.semilogy(x,np.exp(np.polyval(q,np.log(x))),'-',color=C['blue'],lw=0.8)
ax2.text(36,5e-3,r'$\propto e^{-%.3fN}$'%(-k1[0]),fontsize=7,color=C['orange']); ax2.text(36,0.42,r'$\propto N^{%.2f}$'%q[0],fontsize=7,color=C['blue'])
ax2.set_xlabel('$N$'); ax2.set_ylabel(r'distance$/\gamma$'); ax2.legend(loc='lower left'); tag(ax2,'(b)',x=0.9)
plt.tight_layout(); plt.savefig('../figs/rev_fig2.pdf'); plt.savefig('../figs/rev_fig2.png',dpi=200)
