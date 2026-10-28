"""Static research figure; source numbers are in the accompanying report/JSON."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
x=json.loads((ROOT/'data/biorthogonal_audit/analysis.json').read_text())
r=x['records']; N=np.array([d['N'] for d in r])
blue,orange='#2a78d6','#eb6834'
ink,muted,surface='#0b0b0b','#52514e','#fcfcfb'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':ink,
  'axes.labelcolor':ink,'xtick.color':muted,'ytick.color':muted,
  'axes.edgecolor':'#c3c2b7','axes.spines.top':False,'axes.spines.right':False,
  'axes.facecolor':surface,'figure.facecolor':surface,'savefig.facecolor':surface,
  'grid.color':'#e1e0d9','grid.linewidth':.6,'lines.linewidth':1.5})
fig,axs=plt.subplots(1,3,figsize=(13.4,4.3),layout='constrained')

def line(ax,y,label,color,marker,style='-'):
  ax.semilogy(N,y,label=label,color=color,marker=marker,linestyle=style,
              markersize=6.5,markeredgecolor=surface,markeredgewidth=1.0)

ax=axs[0]
line(ax,[d['Vab'] for d in r],r'$|V_{ab}|$',blue,'o')
line(ax,[d['Vba'] for d in r],r'$|V_{ba}|$',orange,'s')
line(ax,[d['g'] for d in r],r'$g$',muted,'^','--')
ax.set_title('(a) Couplings in the population gauge',loc='left',fontsize=10.5)
ax.set_ylabel('Magnitude (decay rate = 1)')
ax.text(57,r[-1]['Vab'],r'$V_{ab}$',va='center')
ax.text(57,r[-1]['Vba'],r'$V_{ba}$',va='center')
ax.text(57,r[-1]['g'],r'$g$',va='center')
ax.set_xlim(22,63);ax.legend(loc='lower left',frameon=False,fontsize=9)
ax.text(.43,.71,r'$\alpha_{ab}=0.00964$'+'\n'+r'$\alpha_{ba}=0.27137$',transform=ax.transAxes,fontsize=9)

ax=axs[1]
line(ax,[d['contractions']['ba']['raw_abs_sum'] for d in r],r'$\sum_i |w_{b,i}(Lr_a)_i|$',blue,'o')
line(ax,[d['Vba'] for d in r],r'$|w_b^T Lr_a|$',orange,'s')
ax.set_title('(b) Cancellation in the direct contraction',loc='left',fontsize=10.5)
ax.set_ylabel('Magnitude')
ax.text(57,r[-1]['contractions']['ba']['raw_abs_sum'],'Envelope',va='center',fontsize=9)
ax.text(57,r[-1]['Vba'],'Net',va='center',fontsize=9)
ax.set_xlim(22,66);ax.legend(loc='lower left',frameon=False,fontsize=9)
ax.text(.04,.70,'Envelope exponent: 0.00468\nCancellation exponent: 0.26669',transform=ax.transAxes,fontsize=9)

ax=axs[2]
for key,label,color,marker in [('Vab',r'$V_{ab}$',blue,'o'),('Vba',r'$V_{ba}$',orange,'s')]:
  axis=0 if key=='Vab' else 1
  alpha=[x['fits']['24-56'][key]['alpha'],x['fits']['24-56']['trace_'+key]['alpha'],
         x['fits']['24-56']['balanced_'+key]['alpha']]
  ax.plot([0,1,2],alpha,color=color,marker=marker,label=label,markersize=6.5,
          markeredgecolor=surface,markeredgewidth=1.)
ax.axhline(x['fits']['24-56']['g']['alpha'],color=muted,ls='--',lw=1,label=r'$\sigma=0.14050$')
ax.set_xticks([0,1,2],['Population\npeak','Physical\ntrace norm','Balanced'])
ax.set_title('(c) Individual exponents depend on gauge',loc='left',fontsize=10.5)
ax.set_ylabel('Fitted exponent')
ax.set_ylim(-.02,.32);ax.legend(loc='lower right',frameon=False,fontsize=9)
ax.annotate('Both equal $\\sigma$ here',xy=(2,.1405),xytext=(.62,.20),
            arrowprops={'arrowstyle':'-','color':muted,'lw':.7},fontsize=9)
for a in axs:
  a.grid(axis='y',which='major')
  a.set_axisbelow(True)
for a in axs[:2]: a.set_xlabel('Number of atoms, N')
fig.suptitle('Biorthogonal reduction: a reliable gap decomposition, not yet an action prediction',fontsize=13)
(ROOT/'figs').mkdir(exist_ok=True)
fig.savefig(ROOT/'figs/biorthogonal_audit.png',dpi=170)
svg_path=ROOT/'figs/biorthogonal_audit.svg'
fig.savefig(svg_path)
svg_path.write_text('\n'.join(line.rstrip() for line in svg_path.read_text().splitlines())+'\n')
