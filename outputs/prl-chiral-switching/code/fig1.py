import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from folds import fold_curve
from mf import roots, cusp
from mfdyn import rhs, steady
from loops import LOOPS, mkpath, P0, V, start_angle
plt.rcParams.update({'font.size':8,'font.family':'serif','mathtext.fontset':'cm'})
fig,ax=plt.subplots(1,2,figsize=(7.0,2.6),gridspec_kw=dict(width_ratios=[1,1.25]))
# (a) phase diagram
Dg,Og=np.meshgrid(np.linspace(0,7.5,301),np.linspace(0.05,5,250))
nb=np.vectorize(lambda d,o: len(roots(d,o,V)))(Dg,Og)
ax[0].contourf(Dg,Og,nb==3,levels=[0.5,1.5],colors=['#c9d9ef'])
fc=fold_curve(); ax[0].plot(fc[:,1],fc[:,2],'.',ms=0.4,color='#2b5d9b')
for n0,Dc,Oc in cusp(V): ax[0].plot(Dc,Oc,'*',color='#d4a017',ms=8,mec='k',mew=0.4,zorder=5)
cols={'A':'#b2182b','B':'#1b7837'}
for k in 'AB':
    p=mkpath(k,1,1); s=np.linspace(0,1,400); xy=np.array([p(t)[:2] for t in s])
    ax[0].plot(xy[:,0],xy[:,1],color=cols[k],lw=1.2)
    i=100 if k=='A' else 140; ax[0].annotate('',xy=xy[i+3],xytext=xy[i],arrowprops=dict(arrowstyle='->',color=cols[k],lw=1.2))
    ax[0].text(*(xy[200]+np.array([0.12,0.12 if k=='A' else -0.25])),k,color=cols[k],fontsize=9)
ax[0].plot(*P0,'ko',ms=3.5)
ax[0].set_xlim(0,7.5); ax[0].set_ylim(0,5); ax[0].set_xlabel(r'$\Delta/\gamma$'); ax[0].set_ylabel(r'$\Omega/\gamma$')
ax[0].text(0.03,0.93,'(a)',transform=ax[0].transAxes)
# (b) MF dynamics along loops: n vs phase
T=1000.
for k,ls in zip('AB',['-','--']):
    for sgn in (1,-1):
        s=solve_ivp(rhs,[0,T],steady(*P0,V,2),args=(mkpath(k,T,sgn),),rtol=1e-8,atol=1e-10,dense_output=True,method='LSODA')
        t=np.linspace(0,T,2000); n=s.sol(t)[0]
        ax[1].plot(t/T,n,ls,color=cols[k] if sgn>0 else ('#ef8a62' if k=='A' else '#7fbf7b'),lw=1.1,
                   label=f'{k}, {"CCW" if sgn>0 else "CW"}')
ax[1].set_xlabel(r'$t/T$'); ax[1].set_ylabel(r'$n=\langle \hat N_e\rangle/N$'); ax[1].legend(frameon=False,fontsize=7,ncol=2,loc='center right')
r=roots(*P0,V); [ax[1].axhline(x,color='gray',lw=0.4,ls=':') for x in (r[0],r[2])]
ax[1].text(0.03,0.93,'(b)',transform=ax[1].transAxes); ax[1].set_xlim(0,1)
plt.tight_layout(); plt.savefig('../figs/fig1.pdf'); plt.savefig('../figs/fig1.png',dpi=200)
