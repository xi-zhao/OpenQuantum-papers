"""Static scientific figure; the accompanying report provides the data table."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/coherence_mechanism'
summary = json.loads((DATA/'analysis.json').read_text())
rows = summary['rows']
source = json.loads((DATA/'fraction_N48_q48.json').read_text())
r48 = rows[-1]
blue, orange = '#2a78d6', '#eb6834'
ink, muted, surface = '#0b0b0b', '#52514e', '#fcfcfb'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
  'text.color': ink, 'axes.labelcolor': ink, 'xtick.color': muted,
  'ytick.color': muted, 'axes.edgecolor': '#c3c2b7',
  'axes.spines.top': False, 'axes.spines.right': False,
  'axes.facecolor': surface, 'figure.facecolor': surface,
  'savefig.facecolor': surface, 'grid.color': '#e1e0d9',
  'grid.linewidth': .6, 'lines.linewidth': 1.5})
fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.3), layout='constrained')

ax = axes[0]
x = np.linspace(-1.18, 1.18, 250)
ax.plot(x, x*x-1, '--', color=orange, label='Quadratic expansion')
points = {}
for z, f in source['samples']:
  sx = (z-r48['zc'])/r48['quadratic_g']
  if abs(sx) <= 1.18:
    points.setdefault(round(sx, 4), (sx, f/r48['epsilon']))
samples = np.array(sorted(points.values()))
ax.scatter(samples[:, 0], samples[:, 1], s=36, color=blue,
  edgecolors=surface, linewidths=.8, label='Exact elimination', zorder=3)
ax.axhline(0, color=muted, lw=.8)
for z in source['roots']:
  sx = (z-r48['zc'])/r48['quadratic_g']
  ax.plot(sx, 0, 'o', color=blue, markersize=7, markeredgecolor=surface)
ax.text(.04, .90, '$N=48$', transform=ax.transAxes)
ax.text(.30, .56, 'Two self-consistent roots', transform=ax.transAxes, fontsize=9)
ax.text(.52, .31, 'Quadratic', transform=ax.transAxes, fontsize=9)
ax.set(xlabel=r'$(z-z_c)/g_{\rm quad}$', ylabel=r'$f_N(z)/\epsilon_N$',
  xlim=(-1.22, 1.22), ylim=(-1.12, .60))
ax.set_title('(a) Coherence feedback retains spectral memory', loc='left', fontsize=10)
ax.legend(loc='upper center', bbox_to_anchor=(.5, 1.02), frameon=False, fontsize=8)

N = np.array([r['N'] for r in rows])
ax = axes[1]
epsilon = np.array([r['epsilon'] for r in rows])
ax.semilogy(N, epsilon, '-o', color=blue, markersize=6.5,
  markeredgecolor=surface, markeredgewidth=1.)
alpha = summary['fits_24_48']['epsilon']['alpha']
ax.text(.06, .16, rf'Finite-window exponent: {alpha:.5f}'+'\n'+
  'Action origin remains unproved', transform=ax.transAxes, fontsize=9)
ax.text(49, epsilon[-1], r'$\epsilon_N$', va='center', fontsize=10)
ax.set(xlabel='Number of atoms, N', ylabel=r'$\epsilon_N=-f_N(z_c)$',
  xlim=(22, 54), ylim=(7e-6, 2e-2), xticks=N)
ax.set_title('(b) The unresolved exponentially small depth', loc='left', fontsize=10)

ax = axes[2]
error = 100*np.abs([r['quadratic_relative_error'] for r in rows])
ax.semilogy(N, error, '-s', color=orange, markersize=6.5,
  markeredgecolor=surface, markeredgewidth=1.)
ax.text(.42, .86, r'$g_{\rm quad}=\sqrt{2\epsilon_N/c_N}$', transform=ax.transAxes)
ax.annotate(f'{error[-1]:.6f}%', xy=(48, error[-1]), xytext=(32, .002),
  arrowprops={'arrowstyle': '-', 'color': muted, 'lw': .7}, fontsize=9)
ax.set(xlabel='Number of atoms, N', ylabel='Absolute relative gap error (%)',
  xlim=(22, 51), ylim=(3e-4, 2.5), xticks=N)
ax.set_title('(c) Local quadratic description improves', loc='left', fontsize=10)
for ax in axes:
  ax.grid(axis='y', which='major')
  ax.set_axisbelow(True)
fig.suptitle('An exact coherence reduction narrows the exponent problem; it does not yet solve it', fontsize=12)
fig.savefig(ROOT/'figs/coherence_mechanism.png', dpi=170)
svg = ROOT/'figs/coherence_mechanism.svg'
fig.savefig(svg)
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
