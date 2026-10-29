import numpy as np, json, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8,'font.family':'serif','mathtext.fontset':'cm'})
fig,ax=plt.subplots(1,3,figsize=(7.0,2.3))
# (a) slow spectrum vs Delta at Om=2, N=40 (scanN)
d=json.load(open('../data/scanN_Om2.json'))
from mf import roots
for N,c in [('20','#92c5de'),('30','#4393c3'),('50','#053061')]:
    rows=d[N]; D=np.array([r[0] for r in rows]); l=np.array([[-x[0] for x in r[1]] for r in rows])
    ax[0].semilogy(D,l[:,0],'-o',ms=2,color=c,lw=1,label=f'$N={N}$')
    ax[0].semilogy(D,l[:,1],':',color=c,lw=0.8)
ax[0].axvspan(3.5996,4.5750,color='#c9d9ef',zorder=0)
ax[0].set_xlabel(r'$\Delta/\gamma$'); ax[0].set_ylabel(r'$-\mathrm{Re}\,\lambda_{1,2}/\gamma$'); ax[0].legend(frameon=False,fontsize=6.5,loc='lower left')
ax[0].set_ylim(3e-4,2); ax[0].text(0.90,0.9,'(a)',transform=ax[0].transAxes)
# (b) exponential scaling of lambda1 inside wedge
r=json.load(open('../data/sfit.json')); Ns=np.array([20,30,40,50,60])
for D,Om,l,nb in r:
    if abs(Om-2.0)<1e-6 and abs(D-4.0)<1e-3: ax[1].semilogy(Ns,l,'o-',ms=3,color='#b2182b',label=r'$\Delta=4.0$')
    if abs(Om-2.0)<1e-6 and abs(D-4.3)<1e-3: ax[1].semilogy(Ns,l,'s-',ms=3,color='#ef8a62',label=r'$\Delta=4.3$')
    if abs(Om-2.0)<1e-6 and abs(D-3.8)<1e-3: ax[1].semilogy(Ns,l,'^-',ms=3,color='#67a9cf',label=r'$\Delta=3.8$')
ax[1].set_xlabel('$N$'); ax[1].set_ylabel(r'$\lambda_{\rm sw}=-\lambda_1$'); ax[1].legend(frameon=False,fontsize=6.5); ax[1].text(0.85,0.9,'(b)',transform=ax[1].transAxes)
# (c) EP of exact L at complex Delta: Im Delta_EP vs N
c=json.load(open('../data/cep_left.json')) if __import__('os').path.exists('../data/cep_left.json') else None
c=np.array(json.load(open('../data/cep_all.json')))
ax[2].semilogy(c[:,0],c[:,2],'o',ms=2.5,color='#5e3c99',label=r'$\mathrm{Im}\,\Delta_{\rm EP}$')
lam=[]
for D,Om,l,nb in r:
    if abs(Om-2.0)<1e-6 and abs(D-3.65)<1e-3: pass
ax[2].set_xlabel('$N$'); ax[2].set_ylabel(r'$\mathrm{Im}\,\Delta_{\rm EP}/\gamma$')
k=np.polyfit(c[c[:,0]>=30,0],np.log(c[c[:,0]>=30,2]),1); x=np.linspace(8,52,50)
ax[2].semilogy(x,np.exp(np.polyval(k,x)),'-',lw=0.7,color='gray')
ax[2].text(0.35,0.75,r'$\propto e^{-%.3f N}$'%(-k[0]),transform=ax[2].transAxes,fontsize=7)
ax[2].text(0.85,0.9,'(c)',transform=ax[2].transAxes)
plt.tight_layout(); plt.savefig('../figs/fig2.pdf'); plt.savefig('../figs/fig2.png',dpi=200)
