"""Summarize finite-N coherence results without fitting an instanton action."""
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import scipy
import qutip

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/coherence_mechanism'


def fitted_alpha(N, y):
  slope, intercept = np.polyfit(N, np.log(y), 1)
  return dict(alpha=float(-slope), intercept=float(intercept),
    maximum_log_residual=float(np.max(abs(np.log(y)-(slope*N+intercept)))))


def main():
  rows = []
  for N in (24, 32, 40, 48):
    fraction = json.loads((DATA/f'fraction_N{N}_q{N}.json').read_text())
    deriv = json.loads((DATA/f'fraction_derivatives_N{N}.json').read_text())
    rows.append(dict(N=N, detuning=fraction['detuning'],
      full_g=deriv['reference_g'], reconstructed_g=fraction['g'],
      root_reconstruction_relative_error=fraction['reference_relative_error'],
      epsilon=-deriv['f'], curvature=deriv['f_second'], zc=deriv['z'],
      quadratic_g=deriv['quadratic_g'],
      quadratic_relative_error=deriv['relative_prediction_error'],
      stationarity_residual=abs(deriv['f_prime'])))
  N = np.array([r['N'] for r in rows])
  fits = {key: fitted_alpha(N, np.array([r[key] for r in rows]))
          for key in ('full_g', 'epsilon', 'curvature', 'quadratic_g')}
  fits['quadratic_alpha_identity'] = (fits['epsilon']['alpha']-fits['curvature']['alpha'])/2
  local = []
  for a, b in zip(rows[:-1], rows[1:]):
    row = dict(Nmin=a['N'], Nmax=b['N'])
    for key in ('full_g', 'epsilon', 'curvature', 'quadratic_g'):
      row[key] = float(-np.log(b[key]/a[key])/(b['N']-a['N']))
    local.append(row)
  truncations = []
  for size in (40, 56):
    source = json.loads((DATA/f'coherence_N{size}.json').read_text())
    for item in source['records']:
      if item['qmax'] not in (6, 8, 10, 12, 14, 16, 20, 24, size):
        continue
      truncations.append(dict(N=size, qmax=item['qmax'],
        dimension=item['dimension'], pair=item.get('pair'),
        halfgap=item.get('halfgap'), relative_gap_error=item.get('relative_gap_error'),
        mean_eigenvalue_error=None if not item.get('pair') else float(
          np.mean(item['pair'])-np.mean(source['reference_pair']))))
  attempts = []
  for path in sorted(DATA.glob('fold_*.json')):
    for row in json.loads(path.read_text()):
      attempts.append(dict(source=path.name, **row))
  inputs = sorted(DATA.glob('*.json'))
  support = [ROOT/'code'/name for name in ('mf.py', 'wkbr.py', 'wkbq.py')]
  hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in inputs+support if p.name != 'analysis.json'}
  result = dict(scope='Finite-N spectral reduction at supplied exact crossing detunings; not a large-N action prediction.',
    provenance_notes=[
      'fraction_N24_q24.json predates local geometric root bracketing and used fixed endpoints -0.81 and -0.71; other fraction files use local bracketing.',
      'Original fraction and derivative outputs predate final root-residual and strict-real-branch guards; their eight roots were back-substituted and revalidated by the final validation script.',
      'Instanton attempts are historical exploratory runs, not successful fold-action estimates; older BVP filenames predate explicit tolerance naming.'
    ],
    rows=rows, fits_24_48=fits, local_exponents=local, truncations=truncations,
    instanton_attempts=attempts,
    instanton_conclusion='No validated path to the fold; no S_fold estimate.',
    versions=dict(python=platform.python_version(), numpy=np.__version__,
      scipy=scipy.__version__, qutip=qutip.__version__), input_sha256=hashes)
  (DATA/'analysis.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
  print(json.dumps(dict(rows=rows, fits=fits), indent=2), flush=True)


if __name__ == '__main__':
  main()
