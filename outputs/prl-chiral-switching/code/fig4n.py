import numpy as np, json
from figstyle import plt, C, setup, tag
setup()
fig,ax=plt.subplots(1,2,figsize=(7.0,2.4))
# (a) v4: ln g_N vs N*Phi_WKB at unseen Omega=2.25, 2.75 (+ v3 points with Phi_exact at 2.5, 3.0)
p4=json.load(open('../data/v4_prediction_frozen.json')); S4={(o['Om'],o['N']):o['g'] for o in map(json.loads,open('../data/v4_SEALED.jsonl'))}
p3=json.load(open('../data/v3_prediction_frozen.json')); S3={(o['Om'],o['N']):o['g'] for o in map(json.loads,open('../data/v3_g_SEALED.jsonl'))}
sets=[('2.25',p4,S4,C['blue'],'o',r'$\Omega=2.25$ (instanton)'),('2.75',p4,S4,C['orange'],'s',r'$\Omega=2.75$ (instanton)'),
      ('2.5',p3,S3,C['aqua'],'^',r'$\Omega=2.5$ (exact $\Phi$)'),('3.0',p3,S3,C['yellow'],'v',r'$\Omega=3.0$ (exact $\Phi$)')]
for Om,p,S,c,mk,lb in sets:
    X=np.array(p[Om]['NPhi']); Y=np.log([S[(float(Om),N)] for N in p[Om]['N']])
    ax[0].plot(X,Y-Y[0]+X[0]*0,mk,color=c,ms=4.5,label=lb,mfc=c if 'instanton' in lb else 'none')
xx=np.linspace(2.8,6.8,10); ax[0].plot(xx,-(xx-3.25),'-',color=C['ink2'],lw=0.8)
ax[0].text(5.0,-1.0,'slope $-1$',fontsize=7,color=C['ink2'])
ax[0].set_xlabel(r'$N\,\Phi(n_f;\Delta_{\min})$'); ax[0].set_ylabel(r'$\ln g_N-\ln g_{24}$'); ax[0].legend(loc='lower left',fontsize=6.2); tag(ax[0],'(a)',x=0.9)
# (b) summary of pre-registered tests: observed/predicted local exponent and slope
labs=['v1*\nfold state','v2\nbottleneck','v3\n$\\Omega$=2.5','v3\n$\\Omega$=3.0','v4\n$\\Omega$=2.25','v4\n$\\Omega$=2.75']
slope=[None,-0.05,0.957,0.878,0.958,0.929]; loc=[0.14/0.52,0.14/0.33,0.944,0.878,0.963,0.934]
x=np.arange(6); ax[1].axhspan(0.85,1.15,color=C['fill'],lw=0)
ax[1].plot(x[2:],np.array(slope[2:]),'o',color=C['blue'],ms=5,label='fit slope')
ax[1].plot(x,loc,'s',color=C['orange'],ms=4.5,mfc='none',mew=1.2,label='local exponent obs/pred')
ax[1].plot([1],[slope[1]],'o',color=C['blue'],ms=5)
ax[1].axhline(1,color=C['ink2'],lw=0.6); ax[1].set_xticks(x); ax[1].set_xticklabels(labs,fontsize=6.2); ax[1].set_ylim(-0.15,1.25)
ax[1].text(0.15,1.08,'fail',fontsize=7,color=C['ink2']); ax[1].text(3.55,1.18,'pass (pre-registered band)',fontsize=6.5,color=C['ink2'])
ax[1].axhline(1.016,xmin=0.62,xmax=0.66,color=C['ink2'],lw=0)
null=[None,None,0.999,0.842,1.016,0.941]
ax[1].plot(x[2:],np.array(null[2:]),'x',color=C['ink2'],ms=5,mew=1.0,label='constant-$\\sigma$ null (local)')
ax[1].set_ylabel('ratio'); ax[1].legend(loc='lower right',fontsize=6.2); tag(ax[1],'(b)',x=0.03,y=0.08)
plt.tight_layout(); plt.savefig('../figs/rev_fig4.pdf'); plt.savefig('../figs/rev_fig4.png',dpi=200)
