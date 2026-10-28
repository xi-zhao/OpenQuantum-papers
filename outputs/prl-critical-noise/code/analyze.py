"""Build figures and verification records from deterministic no-click data.

Run from any directory. No downloads, random trajectories, or external services.
All plotted values originate in the saved numerical runs. The old-paper formula
is transcribed independently for an explicit normalization check.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('first_click', ROOT/'code/first_click.py')
soft = importlib.util.module_from_spec(spec)
spec.loader.exec_module(soft)
DATA = ROOT/'data'
FIG = ROOT/'figures'
BLUE = '#2a78d6'
ORANGE = '#eb6834'
INK = '#242424'


def load(name):
  return json.loads((DATA/name).read_text())


def old_paper_survival(omega, gamma, recoil, t):
  """Fedorov et al. Eq.84, resonance alpha=pi/2, density 1/pi."""
  alpha = np.pi/2
  tau = np.sqrt(2)*abs(omega)*t*np.sqrt(recoil*np.sin(alpha)/gamma)
  numerator = np.sqrt(2/np.pi)*(2*recoil*gamma/(abs(omega)**2*np.sin(alpha)))**.25
  denominator = np.sqrt(
    np.sin(alpha/2)*np.sin(2*tau*np.cos(alpha/2))
    + np.cos(alpha/2)*np.sinh(2*tau*np.sin(alpha/2)))
  return numerator/denominator


def exact_gaussian_survival(d, t, omega=1.):
  """H=(p^2+omega^2*q^2)/2-id*q^2/2, initial oscillator ground state."""
  z = np.sqrt(omega**2-1j*d)
  denominator = np.cos(z*t)+1j*omega/z*np.sin(z*t)
  width = (omega*np.cos(z*t)+1j*z*np.sin(z*t))/denominator
  return np.sqrt(omega/width.real)/abs(denominator)


def draw():
  FIG.mkdir(exist_ok=True)
  plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 9,
    'axes.labelsize': 9, 'axes.titlesize': 10, 'legend.fontsize': 8,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.edgecolor': '#898781', 'xtick.color': INK, 'ytick.color': INK,
    'axes.labelcolor': INK, 'text.color': INK,
    'lines.linewidth': 1.5, 'savefig.dpi': 180,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
  })
  base = load('soft_mode.json')
  moments = load('soft_mode_moments.json')
  newest = {r['d']:r for r in base['results']}
  newest.update({r['d']:r for r in moments['results']})
  rows = [newest[k] for k in sorted(newest)]
  constants = moments['constants']
  fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8), layout='constrained')
  d = np.array([r['d'] for r in rows])
  ax[0].loglog(d, [r['median'] for r in rows], 'o-', color=BLUE,
              ms=5, label='Median')
  ax[0].loglog(d, [r['mean'] for r in rows], 's-', color=ORANGE,
              ms=5, label='Mean')
  strong = np.geomspace(30, 1e5, 100)
  ax[0].loglog(strong, constants['median_coefficient']/strong,
              '--', color=BLUE, lw=1)
  ax[0].loglog(strong, constants['mean_coefficient']*strong**(-.75),
              '--', color=ORANGE, lw=1)
  ax[0].text(500, 7e-5, r'Median $\propto d^{-1}$', fontsize=8)
  ax[0].text(300, .12, r'Mean $\propto d^{-3/4}$', fontsize=8)
  ax[0].set(xlabel=r'Monitoring strength $d$', ylabel=r'First-click time $\tau$',
            title='(a) Distinct quantile and moment scales')
  ax[0].legend(frameon=False, loc='lower left')
  r = newest[100000.]
  u = np.sqrt(r['d'])*np.array(r['tau'])
  mask = (u>.03)&(u<7)
  ax[1].semilogy(u[mask],np.array(r['survival'])[mask]*r['d']**.25/
                constants['density_at_zero'], color=BLUE,
                label=r'Quartic, $d=10^5$')
  uu=np.geomspace(.03,7,500)
  ax[1].semilogy(uu,soft.boundary_shape(uu),'--',color=INK,
                label='Known boundary kernel')
  gauss=exact_gaussian_survival(1e5,uu/np.sqrt(1e5))*1e5**.25*np.sqrt(np.pi)
  ax[1].semilogy(uu,gauss,':',color=ORANGE,label='Gaussian control')
  ax[1].set(xlabel=r'$u=\sqrt{d}\,\tau$',ylabel=r'$d^{1/4}S_d/|\psi_0(0)|^2$',
            title='(b) Kernel already contained in Ref. [1]')
  ax[1].legend(frameon=False, loc='upper right')
  fig.savefig(FIG/'first_click.pdf')
  fig.savefig(FIG/'first_click.png')
  plt.close(fig)

  micro=load('dicke.json')+load('dicke_large.json')
  fig, ax=plt.subplots(1,2,figsize=(7.0,2.6),layout='constrained')
  for n, style, color in [(16,':','#898781'),(128,'--',ORANGE),(512,'-',BLUE)]:
    r=next(r for r in micro if r['N']==n and r['d']==1.)
    ax[0].plot(r['tau'],r['survival'],style,color=color,label=fr'$N={n}$')
  r=next(r for r in base['results'] if r['d']==1.)
  ax[0].plot(r['tau'],r['survival'],'-.',color=INK,label='Soft mode')
  ax[0].set(xlim=(0,7),ylim=(0,1),xlabel=r'$\tau=t/(4N)^{1/3}$',ylabel='Survival probability',
            title='(a) Microscopic approach at fixed $d=1$')
  ax[0].legend(frameon=False)
  differences=[]
  for r in micro:
    if r['d'] != 1.: continue
    s=next(s for s in base['results'] if s['d']==1.)
    error=float(np.max(np.abs(np.array(r['survival'])-np.interp(r['tau'],s['tau'],s['survival']))))
    differences.append((r['N'],error))
  differences.sort()
  sizes,errors=np.array(differences).T
  ax[1].loglog(sizes,errors,'o-',ms=5,color=BLUE,label='Measured discrepancy')
  ax[1].loglog(sizes,errors[-1]*(sizes/sizes[-1])**(-1/3),'--',color=INK,label=r'$N^{-1/3}$ guide')
  ax[1].set(xlabel=r'Atom number $N$',ylabel=r'$\max_\tau|S_N-S_{\rm soft}|$',
            title='(b) Finite-size error remains visible')
  ax[1].legend(frameon=False)
  fig.savefig(FIG/'microscopic.pdf')
  fig.savefig(FIG/'microscopic.png')
  plt.close(fig)
  return rows


def verify(rows):
  fine=load('soft_mode_moments.json')['results']
  refined=load('soft_mode_refined.json')['results']
  micro=load('dicke_large.json')
  cutoff=load('dicke_cutoff_check.json')
  space=[]
  for r in fine:
    same=next((s for s in refined if s['d']==r['d']),None)
    if same:
      space.append(dict(d=r['d'], **{key:r[key]/same[key]-1 for key in ['mean','median','second_moment']}))
  cut=[]
  for r in cutoff:
    s=next(s for s in micro if s['N']==r['N'] and s['d']==r['d'])
    cut.append(dict(N=r['N'],d=r['d'],max_survival_difference=float(
      np.max(np.abs(np.array(r['survival'])-s['survival']))),
      relative_gap_difference=r['gap_scaled']/s['gap_scaled']-1))
  u=np.geomspace(.001,10,400)
  omega,gamma,recoil=1.,10.,1e-4
  d=omega**2/(2*gamma*recoil)
  t=u/(2*recoil*np.sqrt(d))
  old=old_paper_survival(omega,gamma,recoil,t)
  ours=d**(-.25)*soft.boundary_shape(u)/np.pi
  equality=float(np.max(np.abs(old/ours-1)))
  independent=load('soft_mode_converged.json')['results']
  rerun=max(abs(r['mean']/next(s for s in independent if s['d']==r['d'])['mean']-1)
            for r in fine if any(s['d']==r['d'] for s in independent))
  report=dict(
    interpretation='Candidate rejected for PRL novelty; numerical success is separate.',
    software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__),
    old_paper_eq84_max_relative_difference=equality,
    repeated_ground_run_max_mean_relative_difference=rerun,
    soft_mode_grid_step_comparison=space,
    microscopic_cutoff_48_to_64=cut,
    maximum_emission_normalization_error=max(abs(r['emission_probability']-1) for r in fine),
    maximum_mean_flux_relative_difference=max(abs(r['mean_flux_relative_difference']) for r in fine),
    maximum_survival_end=max(r['survival_end'] for r in fine),
    maximum_survival_increase=max(r['maximum_survival_increase'] for r in fine),
    limitations=[
      'No uniform microscopic long-time moment convergence theorem.',
      'Microscopic means are truncated; not used as full means.',
      'Strong-d finite-size behavior is not extrapolated from fixed N.',
      'No full Lindblad purity evolution was performed.',
      'The highest-d data retain small discretization errors.',
      'Numerical checks do not establish novelty or PRL significance.'
    ])
  assert equality<1e-12
  assert max(r['max_survival_difference'] for r in cut)<2e-9
  assert report['maximum_emission_normalization_error']<1e-5
  assert report['maximum_mean_flux_relative_difference']<3e-6
  assert report['maximum_survival_increase']<1e-10
  assert rerun<1e-8
  (DATA/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
  # A plain-text table accompanies the static figures.
  table=['d\tmedian\tmean\tsecond_moment']
  table += ['\t'.join(f'{r[k]:.12g}' for k in ['d','median','mean','second_moment']) for r in rows]
  (DATA/'figure_table.tsv').write_text('\n'.join(table)+'\n')
  print(json.dumps(report,indent=2))


if __name__=='__main__':
  verify(draw())
