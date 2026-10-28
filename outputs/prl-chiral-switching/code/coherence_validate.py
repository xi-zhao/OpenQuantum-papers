"""Validate coherence elimination by reconstructing full right eigenmodes.

No full-L eigenvectors are inputs. Root locations come from the continued
fraction; the full generator is used only for a posteriori residuals. The
reciprocal condition estimates diagnose elimination denominators locally,
not poles elsewhere in the spectral plane.
"""
from pathlib import Path
import argparse
import hashlib
import json
import warnings
import numpy as np
from scipy.linalg import eig, eigvals, lu_factor, lu_solve
from scipy.linalg.lapack import get_lapack_funcs
from coherence_fraction import Fraction, DATA, checked_root, nearest_real_index
from coherence_fraction_derivatives import derivatives
from biorthogonal_generator import generators
from biorthogonal_audit import mm


def reconstruct(N, detuning, z):
  model = Fraction(N, detuning)
  sigma = np.zeros_like(model.A[-1])
  transfers = [None] * N
  conditioning = []
  for k in range(N, 0, -1):
    X = z * np.eye(len(model.A[k])) - model.A[k] - sigma
    lu, piv = lu_factor(X)
    gecon = get_lapack_funcs('gecon', (X,))
    rcond, info = gecon(lu, np.linalg.norm(X, 1))
    if info != 0:
      raise RuntimeError(f'Condition estimate failed at q={k}: {info}')
    conditioning.append(dict(q=k, reciprocal_condition_1=float(rcond)))
    transfers[k-1] = lu_solve((lu, piv), model.down[k-1])
    sigma = mm(model.up[k-1], transfers[k-1])
  K = model.A[0].real + 2 * sigma.real
  values, vectors = eig(K)
  selected = nearest_real_index(values, z)
  r0 = vectors[:, selected].real
  positive = [r0.astype(complex)]
  for transfer in transfers:
    positive.append(mm(transfer, positive[-1][:, None])[:, 0])

  L0, L1, space = generators(N)
  L = L0 + detuning * L1
  q = space['ne'][space['row']] - space['ne'][space['col']]
  coo = L.tocoo()
  max_step = int(np.max(abs(q[coo.row]-q[coo.col])))
  if max_step != 1:
    raise RuntimeError('Generator is not tridiagonal in coherence order')
  full = np.zeros(L.shape[0], complex)
  for k, block in enumerate(positive):
    ids = np.flatnonzero(q == k)
    full[ids] = block
    if k:
      full[space['swap'][ids]] = block.conj()
  norm = np.linalg.norm(full)
  residual = float(np.linalg.norm(L @ full-z*full)/norm)
  trace = float(abs(np.dot(space['trace'], full))/norm)
  hermitian_error = float(np.linalg.norm(full-full[space['swap']].conj())/norm)
  result = dict(N=N, detuning=detuning, root=z,
    scalar_residual=float(values[selected].real-z),
    full_mode_residual=residual, normalized_trace=trace,
    hermiticity_error=hermitian_error,
    q0_norm_fraction=float(np.linalg.norm(r0)/norm),
    nearest_other_effective_eigenvalue_distance=float(np.min(
      abs(np.delete(values, selected)-values[selected]))),
    max_coherence_step=max_step,
    minimum_reciprocal_condition_1=min(r['reciprocal_condition_1'] for r in conditioning),
    denominator_conditioning=conditioning)
  if residual > 1e-9 or trace > 1e-8 or hermitian_error > 1e-12:
    raise RuntimeError(f'Back-substitution failed: {result}')
  return result


def finite_difference():
  source = json.loads((DATA/'fraction_N24_q24.json').read_text())
  model = Fraction(24, source['detuning'])
  z = source['center']
  K, K1, K2 = derivatives(model, z)
  rows = []
  for dz in (1e-3, 3e-4, 1e-4):
    plus, minus = model.operator(z+dz), model.operator(z-dz)
    rows.append(dict(step=dz,
      first_relative_error=float(np.linalg.norm((plus-minus)/(2*dz)-K1)/np.linalg.norm(K1)),
      second_relative_error=float(np.linalg.norm((plus+minus-2*K)/dz**2-K2)/np.linalg.norm(K2))))
  if rows[-1]['first_relative_error'] > 1e-4 or rows[-1]['second_relative_error'] > 1e-3:
    raise RuntimeError(f'Analytic derivative check failed: {rows}')
  return rows


def failure_checks():
  outcomes = []

  def expect_error(label, action, kind):
    try:
      action()
    except kind as error:
      outcomes.append(dict(case=label, rejected=True, reason=str(error)))
    else:
      raise AssertionError(f'Invalid input accepted: {label}')

  spectrum = np.array([-.80, -.72])
  def discontinuous(z):
    return float(spectrum[np.argmin(abs(spectrum-z))]-z)
  expect_error('nearest-branch jump with no root',
    lambda: checked_root(discontinuous, -.78, -.74), RuntimeError)
  tiny_pair = eigvals(np.array([[-.75, -1e-9], [1e-9, -.75]]))
  expect_error('complex pair with small imaginary parts',
    lambda: nearest_real_index(tiny_pair, -.75), RuntimeError)
  expect_error('complex detuning', lambda: Fraction(2, 3.2+1e-5j), ValueError)
  model = Fraction(2, 3.2)
  expect_error('complex spectral parameter', lambda: model.operator(-.75+1e-5j), ValueError)
  expect_error('negative coherence cutoff', lambda: model.operator(-.75, -1), ValueError)
  root = checked_root(lambda z: (z+.8)*(z+.72), -.82, -.76)
  if abs(root+.8) > 1e-11:
    raise AssertionError('Valid smooth root rejected or displaced')
  outcomes.append(dict(case='valid smooth root', accepted=True, root=root))
  (DATA/'failure_path_validation.json').write_text(json.dumps(outcomes, indent=2)+'\n')
  return outcomes


def main(failure_only=False):
  failures = failure_checks()
  if failure_only:
    print(json.dumps(failures, indent=2), flush=True)
    return
  with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    validation = dict(derivatives=finite_difference(), failure_paths=failures, reconstructed_modes=[])
    for N in (24, 32, 40, 48):
      source = json.loads((DATA/f'fraction_N{N}_q{N}.json').read_text())
      for z in source['roots']:
        row = reconstruct(N, source['detuning'], z)
        validation['reconstructed_modes'].append(row)
        print(json.dumps({k: v for k, v in row.items() if k != 'denominator_conditioning'}), flush=True)
  validation['warnings'] = sorted(set(str(w.message) for w in caught))
  code = Path(__file__).resolve().parent
  validation['source_sha256'] = {name: hashlib.sha256((code/name).read_bytes()).hexdigest()
    for name in ('coherence_fraction.py', 'coherence_fraction_derivatives.py', 'coherence_validate.py')}
  (DATA/'validation.json').write_text(json.dumps(validation, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--failure-only', action='store_true')
  main(parser.parse_args().failure_only)
