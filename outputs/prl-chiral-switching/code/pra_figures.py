"""Publication figures for the finite-size spectral-structure manuscript.

All plotted values are drawn from stored audited numerical outputs. Tables
and source records accompany the manuscript; line style/marker complements
color. The two categorical colors pass the dataviz light-palette validator.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mf import roots

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/pra_spectral'
FIG = ROOT/'figs'
blue, orange = '#2a78d6', '#eb6834'
ink, muted, surface = '#111111', '#595959', '#ffffff'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,
  'text.color':ink,'axes.labelcolor':ink,'xtick.color':ink,'ytick.color':ink,
  'axes.edgecolor':'#9a9a9a','axes.spines.top':False,'axes.spines.right':False,
  'axes.facecolor':surface,'figure.facecolor':surface,'savefig.facecolor':surface,
  'grid.color':'#e3e3e3','grid.linewidth':.5,'lines.linewidth':1.25,
  'pdf.fonttype':42,'ps.fonttype':42})


def load(path):return json.loads((ROOT/path).read_text())


def finish(fig, name):
  fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight')
  fig.savefig(FIG/f'{name}.png',dpi=180,bbox_inches='tight')
  plt.close(fig)


def line(ax, x, y, color, label, marker='o', style='-'):
  ax.plot(x,y,color=color,label=label,marker=marker,ls=style,ms=4.5,
          markeredgecolor=surface,markeredgewidth=.65)


def spectrum():
  modes=load('data/pra_spectral/modes_N32.json')['records']
  original=load('data/biorthogonal_audit/analysis.json')['records']
  fig,axs=plt.subplots(1,3,figsize=(7.05,2.45),layout='constrained')
  ax=axs[0]
  for d in np.linspace(2.95,4.85,500):
    r=roots(d,2.,10.)
    for i,n in enumerate(r):
      ax.plot(d,n,'.',color=muted if len(r)==3 and i==1 else blue,ms=1.2)
  dstar=np.array([r['detuning'] for r in original])
  ax.plot(dstar,np.full(len(dstar),.025),'|',color=orange,ms=6,label='Crossings')
  ax.axvline(3.599636890457168,color=muted,lw=.7,ls=':')
  ax.text(3.67,.13,'Left fold',fontsize=7)
  ax.text(4.12,.465,'Stable',fontsize=7)
  ax.text(4.07,.25,'Unstable',fontsize=7)
  ax.set(xlabel=r'$\Delta/\gamma$',ylabel=r'MF population $n$',ylim=(0,.5))
  ax.set_title('(a) Mean-field stationary branches',loc='left',fontsize=8)
  ax.legend(frameon=False,fontsize=7,loc='upper left')
  ax=axs[1];x=np.array([r['offset'] for r in modes]);lam=np.array([r['eigenvalues'] for r in modes])
  line(ax,x,lam[:,0],blue,r'$\lambda_-$','o')
  line(ax,x,lam[:,1],orange,r'$\lambda_+$','s')
  ax.plot(x,[r['high_meanfield_rate'] for r in modes],':',color=muted,label='High MF')
  ax.text(.025,-.72,r'$\lambda_-$',fontsize=8)
  ax.text(.025,-.64,r'$\lambda_+$',fontsize=8)
  ax.set(xlabel=r'$(\Delta-\Delta_*)/\gamma$',ylabel=r'$\lambda/\gamma$',ylim=(-.99,-.56))
  ax.set_title('(b) Real decay modes, $N=32$',loc='left',fontsize=8)
  ax.legend(frameon=False,fontsize=6.4,loc='lower left')
  ax=axs[2]
  N=np.array([r['N'] for r in original]);g=np.array([r['g'] for r in original])
  line(ax,N,g,blue,'Half-gap','o')
  ax.set_yscale('log')
  ax.text(.07,.13,r'$\sigma_{24:56}=0.14050$'+'\nFinite-window fit',transform=ax.transAxes,fontsize=7)
  ax.set(xlabel='Atom number $N$',ylabel=r'$g_N/\gamma$',xticks=[24,32,40,48,56])
  ax.set_title('(c) Small finite-size half-gap',loc='left',fontsize=8)
  finish(fig,'pra_spectrum')


def residues():
  source=load('data/pra_spectral/modes_N32.json');rows=source['records']
  fig,axs=plt.subplots(1,3,figsize=(7.05,2.45),layout='constrained')
  x=np.array([r['offset'] for r in rows]);a=np.array([r['residues'] for r in rows])
  for column,ax in enumerate(axs[:2]):
    scale=1e3 if column==0 else 1.
    line(ax,x,scale*a[:,0,column],blue,r'$\lambda_-$','o')
    line(ax,x,scale*a[:,1,column],orange,r'$\lambda_+$','s')
    ax.axhline(0,color=muted,lw=.5)
    ax.set_xlabel(r'$(\Delta-\Delta_*)/\gamma$')
    ax.legend(frameon=False,fontsize=7,loc='center left' if column==0 else 'lower right')
    if column==0:
      ax.text(.69,.83,r'$\lambda_-$',transform=ax.transAxes,fontsize=7)
      ax.text(.69,.10,r'$\lambda_+$',transform=ax.transAxes,fontsize=7)
    else:
      ax.text(.65,.80,r'$\lambda_+$',transform=ax.transAxes,fontsize=7)
      ax.text(.20,.18,r'$\lambda_-$',transform=ax.transAxes,fontsize=7)
  axs[0].set_title('(a) High-branch perturbation',loc='left',fontsize=8)
  axs[0].set_ylabel(r'Signed density residue ($10^{-3}$)')
  axs[1].set_title('(b) Fold-state preparation',loc='left',fontsize=8)
  axs[1].set_ylabel('Signed density residue')
  ax=axs[2]
  for offset,color,marker in [(-.06,blue,'o'),(.06,orange,'s')]:
    r=next(r for r in rows if abs(r['offset']-offset)<1e-9)
    values=[p['low_region_absolute_fraction'] for p in r['eigenmode_population_location']]
    line(ax,[0,1],values,color,rf'$\delta={offset:+.2f}$',marker)
  ax.set_xticks([0,1],[r'$\lambda_-$',r'$\lambda_+$'])
  ax.set(ylabel=r'Low-density fraction of $|P_r(k)|$',ylim=(0,.28))
  ax.set_title('(c) Profile character exchanges',loc='left',fontsize=8)
  ax.legend(frameon=False,fontsize=7,loc='upper center')
  ax.text(.05,.68,r'$\delta=-0.06$',transform=ax.transAxes,fontsize=7)
  ax.text(.53,.36,r'$\delta=+0.06$',transform=ax.transAxes,fontsize=7)
  finish(fig,'pra_residues')


def couplings():
  a=load('data/biorthogonal_audit/analysis.json');r=a['records'];N=np.array([x['N'] for x in r])
  fig,axs=plt.subplots(1,2,figsize=(7.05,2.45),layout='constrained')
  ax=axs[0]
  line(ax,N,[x['Vab'] for x in r],blue,r'$|V_{ab}|$','o')
  line(ax,N,[x['Vba'] for x in r],orange,r'$|V_{ba}|$','s')
  line(ax,N,[x['g'] for x in r],muted,r'$g_N$','^','--')
  ax.set_yscale('log');ax.set(xlabel='Atom number $N$',ylabel=r'Rate in units of $\gamma$',xlim=(22,64))
  for key,label in [('Vab',r'$V_{ab}$'),('Vba',r'$V_{ba}$'),('g',r'$g_N$')]:
    ax.text(57,r[-1][key],label,va='center',fontsize=8)
  ax.legend(frameon=False,fontsize=7,loc='lower left')
  ax.set_title('(a) Diabatic couplings in the population gauge',loc='left',fontsize=8)
  ax=axs[1]
  for key,color,marker,label in [('Vab',blue,'o',r'$\alpha_{ab}$'),('Vba',orange,'s',r'$\alpha_{ba}$')]:
    fit=a['fits']['24-56']; vals=[fit[key]['alpha'],fit['trace_'+key]['alpha'],fit['balanced_'+key]['alpha']]
    line(ax,[0,1,2],vals,color,label,marker)
  ax.axhline(a['fits']['24-56']['g']['alpha'],color=muted,ls='--',lw=1,label=r'$\sigma$')
  ax.set_xticks([0,1,2],['Population peak','Trace norm','Balanced'])
  ax.set(ylabel='Finite-window exponent',ylim=(-.02,.32))
  ax.set_title('(b) Individual exponents depend on normalization',loc='left',fontsize=8)
  ax.legend(frameon=False,fontsize=7,loc='center left')
  ax.text(1.15,.31,r'$\alpha_{ba}$',fontsize=8)
  ax.text(1.15,.03,r'$\alpha_{ab}$',fontsize=8)
  finish(fig,'pra_couplings')


def coherence():
  source=load('data/coherence_mechanism/fraction_N48_q48.json')
  r=[]
  for N in (24,32,40,48):
    d=load(f'data/coherence_mechanism/fraction_derivatives_N{N}.json')
    r.append(dict(N=N,zc=d['z'],quadratic_g=d['quadratic_g'],epsilon=-d['f'],
                  quadratic_relative_error=d['relative_prediction_error']))
  last=r[-1]
  fig,axs=plt.subplots(1,3,figsize=(7.05,2.45),layout='constrained')
  ax=axs[0];x=np.linspace(-1.16,1.16,250)
  ax.plot(x,x*x-1,'--',color=orange,label='Quadratic')
  points={}
  for z,f in source['samples']:
    xx=(z-last['zc'])/last['quadratic_g']
    if abs(xx)<=1.16:points.setdefault(round(xx,4),(xx,f/last['epsilon']))
  p=np.array(sorted(points.values()));ax.scatter(p[:,0],p[:,1],s=15,color=blue,edgecolors=surface,linewidths=.5,label='Exact')
  ax.axhline(0,color=muted,lw=.6)
  ax.set(xlabel=r'$(z-z_c)/g_{\rm quad}$',ylabel=r'$f_N(z)/\epsilon_N$',ylim=(-1.1,.5))
  ax.set_title('(a) Two self-consistent roots, $N=48$',loc='left',fontsize=7.6)
  ax.legend(frameon=False,fontsize=7,loc='upper center')
  ax.text(.30,.58,'Two roots',transform=ax.transAxes,fontsize=7)
  N=np.array([v['N'] for v in r]);ax=axs[1]
  line(ax,N,[v['epsilon'] for v in r],blue,'Depth','o');ax.set_yscale('log')
  ax.set(xlabel='Atom number $N$',ylabel=r'$\epsilon_N/\gamma$',xticks=N)
  ax.set_title('(b) Minimum depth',loc='left',fontsize=8)
  ax.text(.08,.12,r'$\alpha_\epsilon=0.26764$',transform=ax.transAxes,fontsize=7)
  ax=axs[2];error=100*abs(np.array([v['quadratic_relative_error'] for v in r]))
  line(ax,N,error,orange,'Error','s');ax.set_yscale('log')
  ax.set(xlabel='Atom number $N$',ylabel='Relative half-gap error (%)',xticks=N)
  ax.set_title('(c) Quadratic approximation',loc='left',fontsize=8)
  ax.text(.06,.14,r'$g_{\rm quad}=\sqrt{2\epsilon_N/c_N}$',transform=ax.transAxes,fontsize=7)
  finish(fig,'pra_coherence')


def dynamics():
  a=load('data/pra_spectral/modes_N32.json');r=next(r for r in a['records'] if abs(r['offset'])<1e-10)
  d=r['dynamics'];t=np.array(d['times']);start=np.array(r['initial_signal']);factor=np.exp(-np.mean(r['eigenvalues'])*t)
  fig,axs=plt.subplots(1,2,figsize=(7.05,2.45),layout='constrained')
  for i,ax in enumerate(axs):
    for key,label,color,ls in [('exact','Full evolution',ink,'-'),('two_modes','Two real modes',orange,'--'),('four_modes','Four modes',blue,':')]:
      y=np.array(d[key])[:,i]/start[i]*factor
      ax.plot(t,y,color=color,ls=ls,label=label)
    ax.set(xlabel=r'$t\gamma$',ylabel=r'$e^{-\bar\lambda t}\delta n(t)/\delta n(0)$',xlim=(0,12))
    ax.set_title('(a) High-branch perturbation' if i==0 else '(b) Fold-state preparation',loc='left',fontsize=8)
    ax.legend(frameon=False,fontsize=6.5,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=3)
    ax.text(.04,.93,'Envelope removed',transform=ax.transAxes,fontsize=6.5)
    low,high=ax.get_ylim();ax.set_ylim(low,high+.16*(high-low))
    ax.set_title('(a) High-branch perturbation' if i==0 else '(b) Fold-state preparation',loc='left',fontsize=8,pad=27)
  finish(fig,'pra_dynamics')


def ep_loop():
  source=load('data/pra_spectral/ep_N32_loop_s17_steps24.json')
  theta=np.array([p['theta'] for p in source['points']])/(2*np.pi)
  values=np.array([p['full_operator_branches'] for p in source['points']])
  fig,axs=plt.subplots(1,2,figsize=(7.05,2.45),layout='constrained')
  for i,ax in enumerate(axs):
    for j,color,style in [(0,blue,'-'),(1,orange,'--')]:
      ax.plot(theta,values[:,j,i],color=color,ls=style,label=f'Branch {j+1}')
    ax.axvline(1,color=muted,lw=.6,ls=':')
    ax.set_xlabel(r'Loop turns $\theta/(2\pi)$')
    ax.set_ylabel(r'$\mathrm{Re}\,\lambda/\gamma$' if i==0 else r'$\mathrm{Im}\,\lambda/\gamma$')
    ax.set_title('(a) Exchange after one circuit' if i==0 else '(b) Return after two circuits',loc='left',fontsize=8)
    ax.legend(frameon=True,fontsize=7,loc='upper right',facecolor=surface,edgecolor='none',framealpha=.95,fancybox=False)
  finish(fig,'pra_ep_loop')


if __name__=='__main__':
  FIG.mkdir(exist_ok=True)
  spectrum();residues();couplings();coherence();dynamics();ep_loop()
