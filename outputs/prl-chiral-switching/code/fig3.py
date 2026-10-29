import numpy as np, json, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from predict import predict, tables, segs
plt.rcParams.update({'font.size':8,'font.family':'serif','mathtext.fontset':'cm'})
R=[json.loads(l) for l in open('../data/batch.jsonl')]
fig,ax=plt.subplots(1,2,figsize=(7.0,2.5))
cols={10:'#fdae61',20:'#d7191c',30:'#7b3294',40:'#2c7bb6'}
Tg=np.logspace(0.5,6,30)
for loop,a,mk in [('A',ax[0],'o'),('B',ax[0],'s')]:
    for N in [10,20,30,40]:
        th,lam,pss,nh,nl=tables(N,loop); t0=segs[loop][2]; dn=float(nh(t0)-nl(t0))
        Cp=[(predict(N,loop,T,1)[0]-predict(N,loop,T,-1)[0])/dn for T in Tg]
        a.semilogx(Tg,Cp,'-' if loop=='A' else '--',color=cols[N],lw=1,label=f'$N={N}$' if loop=='A' else None)
        rs={(round(r['T'],3),r['sgn']):r['nf'] for r in R if r['N']==N and r['loop']==loop}
        Ts=sorted(set(k[0] for k in rs if (k[0],1) in rs and (k[0],-1) in rs))
        a.semilogx(Ts,[(rs[(T,1)]-rs[(T,-1)])/dn for T in Ts],mk,ms=3,mfc='none' if loop=='B' else cols[N],color=cols[N],lw=0)
ax[0].axhline(1,color='gray',lw=0.5,ls=':'); ax[0].axhline(-1,color='gray',lw=0.5,ls=':'); ax[0].axhline(0,color='k',lw=0.4)
ax[0].set_xlabel(r'$\gamma T$'); ax[0].set_ylabel(r'chirality $\mathcal{C}$'); ax[0].legend(frameon=False,fontsize=6.5,loc='center right')
ax[0].text(0.02,0.80,'(a)  loop A',transform=ax[0].transAxes,fontsize=7); ax[0].text(0.02,0.14,'loop B',transform=ax[0].transAxes,fontsize=7)
# (b) collapse with lambda_min along the switching arc
for loop,mk in [('A','o'),('B','s')]:
    sw=1 if loop=='A' else -1
    for N in [10,20,30,40]:
        th,lam,pss,nh,nl=tables(N,loop); t0=segs[loop][2]; dn=float(nh(t0)-nl(t0))
        x=np.linspace(th[0],t0,400) if sw==1 else np.linspace(t0,th[-1],400); lm=lam(x).min()
        Cp=np.array([(predict(N,loop,T,1)[0]-predict(N,loop,T,-1)[0])/dn for T in Tg])
        ax[1].semilogx(Tg*lm,np.abs(Cp),'-' if loop=='A' else '--',color=cols[N],lw=1)
ax[1].set_xlabel(r'$T\,\lambda_{\min}(N)$'); ax[1].set_ylabel(r'$|\mathcal{C}|$')
ax[1].set_xlim(1e-1,1e3); ax[1].text(0.85,0.9,'(b)',transform=ax[1].transAxes)
plt.tight_layout(); plt.savefig('../figs/fig3.pdf'); plt.savefig('../figs/fig3.png',dpi=200)
