import numpy as np, json
from figstyle import plt, C, setup, tag
from folds import fold_curve
from mf import roots, cusp
setup(); V=10.
fig,ax=plt.subplots(1,2,figsize=(7.0,2.45),gridspec_kw=dict(width_ratios=[1,1.15]))
# (a) MF wedge, folds, cusp, EP locus (real part) at Omega=2
Dg,Og=np.meshgrid(np.linspace(1.5,7,281),np.linspace(0.05,5,250))
nb=np.vectorize(lambda d,o: len(roots(d,o,V)))(Dg,Og)
ax[0].contourf(Dg,Og,nb==3,levels=[0.5,1.5],colors=[C['fill']])
fc=fold_curve(); s=fc[:,1]>1.5; ax[0].plot(fc[s,1],fc[s,2],'.',ms=0.5,color=C['blue'])
for n0,Dc,Oc in cusp(V):
    if Dc>1.5: ax[0].plot(Dc,Oc,'*',color=C['yellow'],ms=8,mec=C['ink'],mew=0.4,zorder=5)
ep=np.array(json.load(open('../data/cep_all2.json')))
ax[0].plot(ep[ep[:,0]>=20,1],2.0+0*ep[ep[:,0]>=20,1],'|',color=C['orange'],ms=5,mew=0.8)
ax[0].annotate(r'Re$\,\Delta_{\rm EP}(N)$',xy=(3.08,2.0),xytext=(1.65,2.75),fontsize=7,color=C['ink'],arrowprops=dict(arrowstyle='->',lw=0.6,color=C['ink2']))
ax[0].text(4.15,2.45,'bistable',fontsize=7,color=C['ink2']); ax[0].text(5.6,4.55,'cusp',fontsize=7,color=C['ink2'])
ax[0].text(3.15,1.25,'fold',fontsize=7,color=C['blue'],rotation=35)
ax[0].set_xlim(1.5,7); ax[0].set_ylim(0,5); ax[0].set_xlabel(r'$\Delta/\gamma$'); ax[0].set_ylabel(r'$\Omega/\gamma$'); tag(ax[0],'(a)')
# (b) exact spectrum near the left fold at Omega=2 (dense scans N=30,40): the two real modes near -0.75 and the oscillatory pair
cols={30:C['aqua'],40:C['blue']}
for N in (30,40):
    d=json.load(open(f'../data/dense_N{N}.json')); D=np.array([r[0] for r in d])
    W=[[complex(a,b) for a,b in r[1]] for r in d]
    re=[sorted([-x.real for x in w if abs(x.imag)<1e-6])[:2] for w in W]
    osc=np.array([min([-x.real for x in w if abs(x.imag)>1e-3],default=np.nan) for w in W])
    l1=np.array([x[0] for x in re]); l2=np.array([x[1] if len(x)>1 else np.nan for x in re])
    ax[1].plot(D,l1,'-',color=cols[N],lw=1.4,label=r'$N=%d$'%N); ax[1].plot(D,l2,'-',color=cols[N],lw=1.4)
    ax[1].plot(D,osc,':',color=cols[N],lw=1.0)
ep=np.array(json.load(open('../data/cep_all2.json'))); e=dict((int(r[0]),r) for r in ep)
for N in (30,40): ax[1].plot(e[N][1],-e[N][3],'x',color=cols[N],ms=5,mew=1.2)
ax[1].axhline(0.7506,color=C['ink2'],lw=0.6,ls='--'); ax[1].text(3.50,0.79,r'MF $\kappa_{\rm hi}$',fontsize=6.5,color=C['ink2'])
ax[1].axvline(3.5996,color=C['blue'],lw=0.6); ax[1].text(3.607,1.15,'fold',fontsize=6.5,color=C['blue'])
ax[1].text(3.43,0.33,'ghost\nescape',fontsize=6.5,color=C['ink2'])
ax[1].text(3.13,1.06,r'oscillatory pair ($\cdots$)',fontsize=6.3,color=C['ink2'])
ax[1].set_xlabel(r'$\Delta/\gamma$'); ax[1].set_ylabel(r'$-{\rm Re}\,\lambda/\gamma$'); ax[1].set_xlim(2.9,3.7); ax[1].set_ylim(0,1.3)
ax[1].legend(loc='lower left'); tag(ax[1],'(b)',x=0.9)
plt.tight_layout(); plt.savefig('../figs/rev_fig1.pdf'); plt.savefig('../figs/rev_fig1.png',dpi=200)
