"""Plot stored spectra with one common periodic-gap reference; no new ED."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'axes.linewidth': 0.6,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6})
BL, OR, INK, MU = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e'
E7 = np.array([1, 1.28558, 1.87939, 1.96962, 2.53209, 2.87939, 3.70167])
ODD = {0, 2, 5}


def spectrum(e, w, omega, sigma):
  return (w[None, :] * np.exp(-0.5 * ((omega[:, None] - e[None, :]) / sigma)**2)
          / (np.sqrt(2 * np.pi) * sigma)).sum(1)


periodic = next(r for r in json.loads((ROOT / 'data/lscan_th110_t0.2.json').read_text())
                if r['L'] == 36)
levels = next(r for r in json.loads((ROOT / 'data/line_L36.json').read_text())
              if r['th'] == 110 and r['t'] == 0.2)
opened = json.loads((ROOT / 'data/open2_L27_th110_t0.2.json').read_text())
m_ref = levels['ep'][0]
odd_edge = 1 + levels['e0'][0] / m_ref
omega = np.linspace(0.7, 4.0, 1200)
curves = []
for data, ek, wk, sigma in [(periodic, 'ep', 'wp', 0.03),
                            (periodic, 'e0', 'w0', 0.03),
                            (opened['s_taper'], 'e', 'w', 0.05),
                            (opened['u_taper'], 'e', 'w', 0.05)]:
  curves.append(spectrum(np.array(data[ek]) / m_ref, np.array(data[wk]), omega, sigma))

fig, axes = plt.subplots(2, 1, figsize=(3.4, 3.7), sharex=True,
                         gridspec_kw={'hspace': 0.10})
for ax, (staggered, uniform, title) in zip(axes, [
    (*curves[:2], r'periodic, $L=36$'),
    (*curves[2:], r'open, $L=27$, tapered drive')]):
  scale = staggered.max()
  ax.semilogy(omega, staggered / scale, color=OR, lw=1.2, label='staggered')
  ax.semilogy(omega, uniform / scale, color=BL, lw=1.2, ls='--', label='uniform')
  for i, mass in enumerate(E7):
    ax.axvline(mass, color=OR if i in ODD else BL, lw=0.5, ls=':', alpha=0.8)
  ax.axvline(2, color=MU, lw=0.65, ls='--')
  ax.axvline(odd_edge, color=MU, lw=0.65, ls='-.')
  ax.set_ylim(3e-4, 2)
  ax.set_xlim(0.7, 4.0)
  ax.text(0.98, 0.94, title, transform=ax.transAxes, ha='right', va='top', color=INK)
  ax.set_ylabel('response (norm.)')
  for spine in ('top', 'right'):
    ax.spines[spine].set_visible(False)

for i, mass in enumerate(E7):
  axes[0].text(mass + {2: -0.05, 3: 0.05}.get(i, 0), 2.3, f'$A_{i+1}$',
               ha='center', va='bottom', fontsize=7, color=INK)
axes[0].legend(loc='upper right', frameon=False, fontsize=7, bbox_to_anchor=(1, 0.83))
axes[1].set_xlabel(r'$\omega/m_{\rm ref}$')
axes[0].text(-0.17, 1, '(a)', transform=axes[0].transAxes, fontsize=9)
axes[1].text(-0.17, 1, '(b)', transform=axes[1].transAxes, fontsize=9)
fig.savefig(ROOT / 'figs/fig2_spectra.pdf', bbox_inches='tight')
fig.savefig(ROOT / 'figs/fig2_spectra.png', dpi=200, bbox_inches='tight')
with (ROOT / 'figs/fig2_spectra.csv').open('w', newline='') as stream:
  writer = csv.writer(stream)
  writer.writerow(['omega_over_m_ref', 'periodic_staggered', 'periodic_uniform',
                   'open_staggered', 'open_uniform'])
  writer.writerows(zip(omega, curves[0] / curves[0].max(), curves[1] / curves[0].max(),
                      curves[2] / curves[2].max(), curves[3] / curves[2].max()))
print(f'm_ref={m_ref:.10f}, periodic gap-sum references: 2, {odd_edge:.8f}')
