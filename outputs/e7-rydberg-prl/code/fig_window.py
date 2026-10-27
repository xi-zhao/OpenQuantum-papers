"""Display the stored finite-size gap-ratio scan; no phase or error-budget inference."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size': 8, 'font.family': 'serif', 'axes.linewidth': 0.6})
INK, MU = '#0b0b0b', '#52514e'
COLORS = {28: '#2a78d6', 32: '#eb6834', 36: '#1baf7a'}
STYLES = {28: ':', 32: '--', 36: '-'}
E7 = np.array([2 * np.cos(5 * np.pi / 18), 2 * np.cos(np.pi / 9), 2 * np.cos(np.pi / 18)])
DC, VC = -1.5011, -1.6260
rows = [r for r in json.loads((ROOT / 'data/grid_L32.json').read_text()) if r]
deltas = np.array(sorted({r['D'] for r in rows}))
interactions = np.array(sorted({r['V'] for r in rows}))
errors = np.full((len(interactions), len(deltas)), np.nan)
grid_table = []
for row in rows:
  m1 = row['ep'][0]
  ratios = np.array([row['e0'][0], row['ep'][1], row['e0'][1]]) / m1
  error = 100 * np.max(np.abs(ratios / E7 - 1))
  shown = m1 * 32 >= 18 and row['e0'][0] >= 0.5 * m1
  grid_table.append([row['D'], row['V'], m1, *ratios, error, int(shown)])
  if shown:
    errors[np.searchsorted(interactions, row['V']), np.searchsorted(deltas, row['D'])] = error

fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.5),
                         gridspec_kw={'wspace': 0.55, 'width_ratios': [1, 1.1]})
a, b = axes
cmap = LinearSegmentedColormap.from_list('deviation', ['#cde2fb', '#6da7ec', '#256abf', '#0d366b'])
mesh = a.pcolormesh(deltas - DC, interactions - VC, errors, cmap=cmap,
                    vmin=0, vmax=12, shading='nearest')
a.contour(deltas - DC, interactions - VC, errors, levels=[3], colors=INK, linewidths=0.8)
a.plot(0, 0, marker='*', color='#eb6834', ms=8, mec=INK, mew=0.5, ls='none')
a.text(0.03, 0.03, 'TCI', color=INK, fontsize=7)
angle = np.deg2rad(110)
a.annotate('', xy=(0.4 * np.cos(angle), 0.4 * np.sin(angle)), xytext=(0, 0),
            arrowprops={'arrowstyle': '->', 'color': INK, 'lw': 0.8})
a.text(0.4 * np.cos(angle) - 0.02, 0.4 * np.sin(angle) + 0.02,
        r'$\theta=110^\circ$', fontsize=7, ha='right')
a.set_xlabel(r'$(\Delta-\Delta_c)/\Omega$')
a.set_ylabel(r'$(V_2-V_{2c})/\Omega$')
cb = fig.colorbar(mesh, ax=a, pad=0.02, extend='max')
cb.set_label(r'max gap-ratio deviation $\delta_{32}$ (%)', fontsize=7)
cb.outline.set_linewidth(0.4)
a.text(-0.25, 1, '(a)', transform=a.transAxes, fontsize=9)

markers = ('o', 's', '^')
names = (r'$\widehat m_2$', r'$\widehat m_3$', r'$\widehat m_4$')
cut_table = []
for length, filename in [(28, 'thscan_L28'), (32, 'thscan_L32'), (36, 'line_L36')]:
  points = sorted([r for r in json.loads((ROOT / f'data/{filename}.json').read_text())
                   if r['th'] == 110], key=lambda r: r['t'])
  m1 = np.array([r['ep'][0] for r in points])
  ratios = np.array([[r['e0'][0], r['ep'][1], r['e0'][1]] for r in points]) / m1[:, None]
  deviations = 100 * (ratios / E7 - 1)
  for row, gap, ratios_row, deviations_row in zip(points, m1, ratios, deviations):
    cut_table.append([length, row['t'], gap, *ratios_row, *deviations_row])
  for i in range(3):
    b.plot(m1, deviations[:, i], marker=markers[i], ms=3.5, lw=1.0,
            ls=STYLES[length], color=COLORS[length], mec='white', mew=0.4,
            label=f'$L={length}$' if i == 0 else None)
b.axhspan(-3, 3, color='#f0efec', zorder=0)
b.axhline(0, color=MU, lw=0.6, ls='--')
b.set_xlabel(r'$\widehat m_1/\Omega$ (along $\theta=110^\circ$)')
b.set_ylabel('gap-ratio deviation from $E_7$ (%)')
b.set_ylim(-11, 6)
b.set_xlim(0.4, 1.35)
b.text(0.42, 3.3, r'$\pm3\%$', fontsize=7, color=MU, ha='left')
h1, labels = b.get_legend_handles_labels()
h2 = [Line2D([], [], marker=marker, color=MU, ls='none', ms=3.5) for marker in markers]
b.legend(h1 + h2, labels + list(names), frameon=False, fontsize=7, loc='lower right', ncol=2)
for spine in ('top', 'right'):
  b.spines[spine].set_visible(False)
b.text(-0.2, 1, '(b)', transform=b.transAxes, fontsize=9)
fig.savefig(ROOT / 'figs/fig1_window.pdf', bbox_inches='tight')
fig.savefig(ROOT / 'figs/fig1_window.png', dpi=200, bbox_inches='tight')
for filename, header, table in [
    ('fig1_grid.csv', ['Delta', 'V2', 'gap1', 'ratio2', 'ratio3', 'ratio4', 'max_error_percent', 'shown'], grid_table),
    ('fig1_cut.csv', ['L', 't', 'gap1', 'ratio2', 'ratio3', 'ratio4', 'error2_percent', 'error3_percent', 'error4_percent'], cut_table)]:
  with (ROOT / 'figs' / filename).open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(header)
    writer.writerows(table)
