"""Small static research figure from the frozen physical-effect audit.

Raw density units throughout. Does not regenerate the manuscript's six figures.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/physical_effects'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
INK, GRID, SURFACE = '#0b0b0b', '#e1e0d9', '#fcfcfb'


def main(args):
  data = json.loads((DATA/args.input).read_text())
  output = ROOT/'figs'/args.output_stem
  targets = [output.with_suffix(s) for s in ('.pdf', '.png', '.csv')]
  if any(p.exists() for p in targets):
    raise FileExistsError('Preserve previous figure; use a new output name in the script')
  times = np.asarray(data['times'])
  scan32 = next(s for s in data['scans'] if s['N'] == 32)
  scan48 = next(s for s in data['scans'] if s['N'] == 48)
  center32 = next(r for r in scan32['records'] if r['offset'] == 0.)
  plt.rcParams.update({'font.size': 10, 'font.family': 'DejaVu Sans', 'text.color': INK,
    'axes.labelcolor': INK, 'xtick.color': '#52514e', 'ytick.color': '#52514e',
    'axes.facecolor': SURFACE, 'figure.facecolor': SURFACE, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.edgecolor': '#c3c2b7', 'pdf.fonttype': 42,
    'lines.linewidth': 1.8, 'legend.frameon': False})
  fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4), layout='constrained')
  ax = axes[0, 0]
  series = [('Full PI, N=32', 'full', BLUE, '-'),
    ('Independent atoms', 'noninteracting', ORANGE, '--'),
    ('Phase-ensemble MF', 'meanfield_phase_ensemble', AQUA, ':')]
  for label, key, color, style in series:
    ax.plot(times, np.asarray(center32[key])[:, 0], color=color, linestyle=style, label=label)
  ax.set(xlim=(0, 6), ylim=(-.18, .20), xlabel=r'$\gamma t$', ylabel=r'$\delta n(t)$',
    title='(a) Coherence-to-population response')
  ax.legend(loc='lower right', fontsize=9)
  ax.text(.98, .95, 'Same initial populations; different coherences', ha='right', va='top',
    transform=ax.transAxes, fontsize=9)
  ax = axes[0, 1]
  control = data['fixed_detuning_control']
  sizes = [r['N'] for r in control['records']]
  peaks = [r['peak'][0] for r in control['records']]
  ax.plot(sizes, peaks, '-o', color=BLUE, markersize=7, markeredgecolor=SURFACE,
    label='Interacting, fixed detuning')
  free_peak = float(np.max(abs(np.asarray(control['noninteracting'])[:, 0])))
  ax.axhline(free_peak, color=ORANGE, linestyle='--', label='Independent atoms')
  ax.set(xlim=(0, 51), ylim=(0, .19), xlabel='$N$', ylabel=r'$\max_t |\delta n(t)|$',
    title='(b) Fixed preparation and detuning')
  ax.legend(loc='lower right', fontsize=9)
  ax.annotate(f'{peaks[-1]:.4f}', (sizes[-1], peaks[-1]), xytext=(-5, -17),
    textcoords='offset points', ha='right', fontsize=9)
  ax.text(48, free_peak+.006, f'{free_peak:.4f}', ha='right', fontsize=9)
  ax = axes[1, 0]
  width = scan48['mixing_width']
  near = [r for r in scan48['records'] if abs(r['offset']) <= 4.01*width]
  x = np.array([r['offset']/width for r in near])
  residues = np.array([r['modes']['residues'] for r in near])[:, :, 1]
  ax.plot(x, residues[:, 0], '-o', color=BLUE, markersize=7, label=r'$a_-$')
  ax.plot(x, residues[:, 1], '--s', color=ORANGE, markersize=7, label=r'$a_+$')
  ax.set(xlabel=r'$(\Delta-\Delta_*)/w_{48}$', ylabel='Signed density residue',
    title='(c) Strong residue exchange, N=48')
  ax.legend(loc='center right', fontsize=10)
  ax.text(x[-1], residues[-1, 0]+.000025, r'$a_-$', ha='right')
  ax.text(x[-1], residues[-1, 1]+.000025, r'$a_+$', ha='right')
  ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
  ax = axes[1, 1]
  hi = min(scan48['records'], key=lambda r: abs(r['offset']-width))
  lo = min(scan48['records'], key=lambda r: abs(r['offset']+width))
  for label, key, color, style in [('Full PI', 'full', BLUE, '-'),
      ('Mean field', 'meanfield_phase_ensemble', ORANGE, '--')]:
    difference = np.asarray(hi[key])[:, 1]-np.asarray(lo[key])[:, 1]
    ax.plot(times, difference, color=color, linestyle=style, label=label)
  ax.set(xlim=(0, 12), xlabel=r'$\gamma t$', ylabel=r'$s(\Delta_*+w)-s(\Delta_*-w)$',
    title='(d) Small full-signal change, same preparation')
  ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
  ax.legend(loc='upper right', fontsize=9)
  ax.text(.98, .82, r'Peak $1.89\times10^{-6}$', transform=ax.transAxes, ha='right', fontsize=10)
  for ax in axes.flat:
    ax.grid(axis='y', color=GRID, linewidth=.6)
    ax.set_axisbelow(True)
  fig.suptitle('Physical-effect test: mode exchange versus measured density', fontsize=14)
  fig.savefig(targets[0], facecolor=SURFACE)
  fig.savefig(targets[1], dpi=170, facecolor=SURFACE)
  rows = []
  for label, key, _, _ in series:
    rows += [('a', label, float(t), float(y)) for t, y in zip(times, np.asarray(center32[key])[:, 0])]
  rows += [('b', 'Interacting, fixed detuning', n, y) for n, y in zip(sizes, peaks)]
  rows += [('b', 'Independent atoms', n, free_peak) for n in sizes]
  for i in range(2):
    rows += [('c', ('a_minus', 'a_plus')[i], float(xx), float(y)) for xx, y in zip(x, residues[:, i])]
  for label, key in [('Full PI', 'full'), ('Mean field', 'meanfield_phase_ensemble')]:
    diff = np.asarray(hi[key])[:, 1]-np.asarray(lo[key])[:, 1]
    rows += [('d', label, float(t), float(y)) for t, y in zip(times, diff)]
  with targets[2].open('w', newline='') as stream:
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(['panel', 'series', 'x', 'y'])
    writer.writerows(rows)
  print(json.dumps(dict(outputs=[str(p) for p in targets], table_rows=len(rows))))


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--input', default='effects_2026-11-01.json')
  parser.add_argument('--output-stem', default='physical_effects_2026-11-01_v2')
  main(parser.parse_args())
