"""Size dependence of the A4 binding and of spectral weights from stored periodic spectra; no new ED."""
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
INK, MU = '#0b0b0b', '#52514e'
BL, OR, GR = '#2a78d6', '#eb6834', '#1baf7a'
B_E7 = 2 - 2 * np.cos(np.pi / 18)


def lines(row):
  m1 = row['ep'][0]
  even = [(e / m1, w) for e, w in zip(row['e0'], row['w0']) if e / m1 > 0.3 and w > 1e-8]
  odd = [(e / m1, w) for e, w in zip(row['ep'], row['wp']) if w > 1e-8]
  return m1, even, odd


fig, (a, b) = plt.subplots(2, 1, figsize=(3.4, 2.9), sharex=True, gridspec_kw={'hspace': 0.12})
from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
table = []
for t, color, marker in ((0.2, BL, 'o'), (0.3, OR, 's')):
  rows = json.loads((ROOT / f'data/lscan_th110_t{t}.json').read_text())
  sizes, binding = [], []
  for row in rows:
    m1, even, odd = lines(row)
    sizes.append(row['L'])
    binding.append(2 - even[1][0])
    if t == 0.2:
      table.append([row['L'], m1, 2 - even[1][0], odd[0][1], even[0][1], odd[1][1], even[1][1],
                    even[2][0], even[2][1], odd[2][0], odd[2][1]])
  a.plot(sizes, binding, marker=marker, ms=3.5, lw=1.0, color=color, mec='white', mew=0.4,
         label=f'$t={t}$')
a.axhline(B_E7, color=MU, lw=0.7, ls='--')
a.text(18.3, B_E7 * 1.06, r'$E_7$: $2-2\cos(\pi/18)$', fontsize=7, color=MU, ha='left', va='bottom')
a.set_yscale('log')
a.set_ylim(0.025, 0.2)
a.yaxis.set_major_locator(FixedLocator([0.03, 0.05, 0.1, 0.2]))
a.yaxis.set_minor_locator(NullLocator())
a.yaxis.set_major_formatter(ScalarFormatter())
a.set_ylabel(r'$(2\widehat m_1-\widehat m_4)/\widehat m_1$')
a.legend(frameon=False, fontsize=7, loc='upper right')
a.text(-0.2, 1.0, '(a)', transform=a.transAxes, fontsize=9)

data = np.array(table)
ref = data[-1]
series = [(3, '$A_1$', OR, '-', 'o'), (5, '$A_3$', OR, '--', '^'),
          (4, '$A_2$', BL, '-', 's'), (6, '$A_4$', BL, '--', 'v'),
          (8, 'even scatt.', GR, ':', 'D'), (10, 'odd scatt.', MU, ':', 'P')]
for col, label, color, ls, marker in series:
  keep = data[:, 0] >= (24 if col >= 8 else 0)
  b.plot(data[keep, 0], data[keep, col] / ref[col], marker=marker, ms=3, lw=1.0, ls=ls,
         color=color, mec='white', mew=0.3, label=label)
lref = np.arange(24, 37)
b.plot(lref, 36 / lref, color=MU, lw=0.6, ls='-', alpha=0.6)
b.text(24.0, 36 / 24 * 0.86, r'$\propto 1/L$', fontsize=6.5, color=MU)
b.set_yscale('log')
b.set_ylim(0.45, 3.0)
b.yaxis.set_major_locator(FixedLocator([0.5, 0.8, 1, 1.5, 2, 3]))
b.yaxis.set_minor_locator(NullLocator())
b.yaxis.set_major_formatter(ScalarFormatter())
b.set_xlabel('$L$')
a.set_xlim(17, 37)
b.set_ylabel(r'$w(L)/w(36)$')
b.legend(frameon=False, fontsize=6.5, loc='lower left', ncol=3, handlelength=1.8, columnspacing=0.8, borderaxespad=0.2)
b.text(-0.2, 1.0, '(b)', transform=b.transAxes, fontsize=9)
for ax in (a, b):
  for spine in ('top', 'right'):
    ax.spines[spine].set_visible(False)
fig.savefig(ROOT / 'figs/fig3_size.pdf', bbox_inches='tight')
fig.savefig(ROOT / 'figs/fig3_size.png', dpi=200, bbox_inches='tight')
with (ROOT / 'figs/fig3_size.csv').open('w', newline='') as stream:
  writer = csv.writer(stream, lineterminator='\n')
  writer.writerow(['L', 'gap1', 'B4_over_m1', 'w_A1', 'w_A2', 'w_A3', 'w_A4',
                   'even_scatt_E_over_m1', 'even_scatt_w', 'odd_scatt_E_over_m1', 'odd_scatt_w'])
  writer.writerows(table)
