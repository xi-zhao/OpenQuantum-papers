"""Observable-focused figures; use only saved validated dynamics and spectra."""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from make_figures import COLORS, MARKERS, save
from spectral_analysis import constants

ROOT = Path(__file__).resolve().parents[1]


def observable_dynamics():
  data = json.loads((ROOT/'data/verified_observables.json').read_text())
  row = next(r for r in data['curves'] if r['n']==12)
  t = np.array(row['time'])
  theory = row['theory']
  amplitudes = row['modal_amplitudes']
  fig, axes = plt.subplots(1,3,figsize=(7.05,2.55),layout='constrained')
  for name,key,amplitude,color,marker,label in [
    ('density','density_deviation','density_asymptotic_amplitude',COLORS[0],'o',r'$N/n-\nu_*$'),
    ('imbalance','imbalance','imbalance_asymptotic_amplitude',COLORS[1],'s',r'$D=n_1-n_2$')]:
    y=np.array(row[key]); norm=y[0]
    axes[0].semilogy(t,y/norm,color=color,lw=1.4,label=label)
    points=np.arange(0,len(t),20)
    axes[0].scatter(t[points],y[points]/norm,s=15,marker=marker,color=color,edgecolors='white',linewidths=.4,zorder=3)
    asym=amplitudes[amplitude]*np.exp(-theory[name]*t)/norm
    mask=t>=3
    axes[0].semilogy(t[mask],asym[mask],color='#333333',ls='--',lw=.65)
  axes[0].semilogy(t,np.exp(-theory['gap_all']*t),color='#888888',lw=.9,ls=':',label=r'$e^{-\Delta_{\rm all}t}$ (reference)')
  axes[0].set(title='(a) Two occupation observables',xlabel=r'$\gamma t$',ylabel='Deviation / initial deviation',xlim=(0,8),ylim=(1e-6,1.5))
  axes[0].legend(frameon=False,loc='lower left',fontsize=6.5,handlelength=1.5)
  axes[0].text(5.2,.12,r'$\Delta_{\rm all}=0.572\gamma$',fontsize=6.5)
  axes[0].text(5.25,.003,r'$\Gamma_D=1.118\gamma$',fontsize=6.5)
  axes[0].text(5.15,1.4e-4,r'$\Gamma_N=1.679\gamma$',fontsize=6.5)

  for idx,eta in enumerate([0.,.1,1.]):
    effective=np.array(row['readouts'][str(eta)]['effective_rate'])
    axes[1].plot(t,effective,color=COLORS[idx],marker=MARKERS[idx],markevery=35,
                 markersize=4,markeredgecolor='white',markeredgewidth=.4,label=rf'$\eta={eta:g}$')
  axes[1].axhline(theory['density'],lw=.7,ls='--',color='#555555')
  axes[1].axhline(theory['imbalance'],lw=.7,ls='--',color='#555555')
  axes[1].text(10,theory['density']+.018,r'$\Gamma_N$',fontsize=7)
  axes[1].text(10,theory['imbalance']+.018,r'$\Gamma_D$',fontsize=7)
  axes[1].set(title='(b) Unequal readout weights',xlabel=r'$\gamma t$',ylabel=r'$-\partial_t\ln|\langle O_\eta\rangle-\nu_*|/\gamma$',xlim=(0,12),ylim=(.75,1.79))
  axes[1].legend(frameon=False,loc='lower left',fontsize=6.5,ncol=1)

  axes[2].semilogy(t,row['survival'],color=COLORS[0],lw=1.4)
  for time in [2,5,8]:
    j=np.argmin(abs(t-time));value=row['survival'][j]
    axes[2].scatter([time],[value],s=20,color=COLORS[0],edgecolors='white',linewidths=.4,zorder=3)
    axes[2].annotate(f'{value:.1e}',(time,value),xytext=(-35,7),textcoords='offset points',fontsize=6.5)
  axes[2].set(title='(c) Conditional survival cost',xlabel=r'$\gamma t$',ylabel=r'$P_{\rm surv}(t)$',xlim=(0,8.4),ylim=(1e-27,1.1))
  axes[2].text(.3,1e-21,r'$n=12,\ d=3,\ z=1/144$',fontsize=7)
  save(fig,'observable_dynamics')
  with (ROOT/'data/observable_dynamics_table.csv').open('w') as f:
    writer=csv.writer(f,lineterminator='\n')
    writer.writerow(['time','density_deviation','imbalance','survival','rate_eta0','rate_eta0.1','rate_eta1'])
    for i,time in enumerate(t):
      writer.writerow([time,row['density_deviation'][i],row['imbalance'][i],row['survival'][i],
                       *[row['readouts'][str(eta)]['effective_rate'][i] for eta in [0.,.1,1.]]])


def rate_hierarchy():
  perron=json.loads((ROOT/'data/verified_perron.json').read_text())
  precise=json.loads((ROOT/'data/verified_precision.json').read_text())
  lookup={(r['n'],r['d']):r for r in perron['points']}
  fig,axes=plt.subplots(1,2,figsize=(3.35,2.1),layout='constrained')
  rows=[]
  for idx,d in enumerate([2,3,4]):
    ps=[r for r in precise['points'] if r['d']==d]
    ns=np.array([r['n'] for r in ps]);ratios_all=[];ratios_obs=[]
    for pr in ps:
      r=lookup[pr['n'],d]
      all_gap=min(pr['gap'],r['gap_one'])
      ratios_all.append(r['gap_two']/all_gap)
      ratios_obs.append(pr['gap']/r['gap_two'])
      rows.append([d,pr['n'],all_gap,r['gap_two'],pr['gap'],ratios_all[-1],ratios_obs[-1]])
    style=dict(color=COLORS[idx],marker=MARKERS[idx],markersize=4,markeredgecolor='white',markeredgewidth=.4,label=rf'$d={d}$')
    axes[0].plot(1/ns,ratios_all,**style)
    axes[1].plot(1/ns,ratios_obs,**style)
    limit=constants(d)['omega']
    axes[1].plot([0,.064],[limit,limit],color=COLORS[idx],ls='--',lw=.7)
    axes[1].text(.045,ratios_obs[0]+.012,rf'$d={d}$',fontsize=6)
  axes[0].axhline(2,color='#777777',ls='--',lw=.7)
  axes[0].set(title='(a) Parity restriction',ylabel=r'$\Gamma_D/\Delta_{\rm all}$',ylim=(1.91,2.015))
  axes[1].set(title='(b) Observable selection',ylabel=r'$\Gamma_N/\Gamma_D$',ylim=(1.20,1.76))
  for ax in axes:
    ax.set(xlabel=r'$1/n$',xlim=(-.001,.068));ax.set_xticks([0,.03,.06])
  axes[0].legend(frameon=False,fontsize=6,loc='lower left',handlelength=1)
  save(fig,'rate_hierarchy')
  with (ROOT/'data/rate_hierarchy_table.csv').open('w') as f:
    writer=csv.writer(f,lineterminator='\n')
    writer.writerow(['d','n','gap_all_z1','Gamma_D_z1','Gamma_N_z1','D_over_all','N_over_D'])
    writer.writerows(rows)


if __name__=='__main__':
  observable_dynamics()
  rate_hierarchy()
