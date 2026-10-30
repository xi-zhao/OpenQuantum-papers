"""Observable residues including the source carried by initial coherences.

Real physical parameters and real spectral roots only. The source and K'(z)
are eliminated by the same coherence-order recurrence as K(z). Frozen roots
are inputs; archived full-generator modes are used only for validation.
Outputs are new JSON records, never replacements for the original samples.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from scipy.linalg import eig, lu_factor, lu_solve, solve
from scipy.linalg.lapack import get_lapack_funcs
from scipy.sparse.linalg import expm_multiply
from biorthogonal_audit import mm
from biorthogonal_generator import generators
from coherence_fraction import Fraction, nearest_real_index
from pra_modes import product_state, fixed, jacobian, DF, NF, OMEGA, V
from mf import roots

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/coherence_response'
PREPARATIONS = ('high_difference', 'fold_state', 'fold_coherence_difference')


def real_checked(value, label, tolerance=1e-8):
  value = np.asarray(value)
  if np.max(abs(value.imag)) > tolerance:
    raise RuntimeError(f'{label} is not real: {np.max(abs(value.imag))}')
  return value.real


def inputs(N, detuning, space):
  qhi = fixed(detuning, roots(detuning, OMEGA, V)[-1])
  values, vectors = eig(jacobian(detuning, qhi))
  realids = np.flatnonzero(abs(values.imag) < 1e-10)
  selected = realids[np.argmax(values[realids].real)]
  direction = vectors[:, selected].real
  direction *= np.sign(direction[0])/np.linalg.norm(direction)
  high = product_state(qhi+.005*direction, space)-product_state(qhi-.005*direction, space)
  folded = product_state(fixed(DF, NF), space)
  diagonal = space['trace'].astype(bool)
  dephased = np.where(diagonal, folded, 0.)
  # This is global excitation-number dephasing, not local atomwise dephasing.
  coherence = folded-dephased
  states = np.column_stack([high, folded, coherence])
  return states, dict(high_state=qhi.tolist(), high_direction=direction.tolist(),
    high_step=.005, fold_state=fixed(DF, NF).tolist(),
    trace=real_checked(space['trace']@states, 'input traces').tolist(),
    identical_q0_error=float(np.max(abs(coherence[diagonal]))),
    dephasing='sum_k Pi_k rho Pi_k; preserves all q=0 Dicke populations')


def eliminate(model, z, sources, ids):
  if not np.isreal(z):
    raise ValueError('Conjugate-sector source reduction requires real z')
  sigma = np.zeros_like(model.A[-1])
  derivative = sigma.copy()
  source = sources[ids[-1]].copy()
  transfers, factors, conditioning = {}, {}, []
  for k in range(model.qmax, 0, -1):
    X = z*np.eye(len(model.A[k]))-model.A[k]-sigma
    lu, piv = lu_factor(X)
    rcond, info = get_lapack_funcs('gecon', (X,))(lu, np.linalg.norm(X, 1))
    if info:
      raise RuntimeError(f'Condition estimate failed at q={k}: {info}')
    conditioning.append(float(rcond))
    transfer = lu_solve((lu, piv), model.down[k-1])
    propagated = lu_solve((lu, piv), source)
    derivative = -mm(model.up[k-1], lu_solve((lu, piv),
      mm(np.eye(len(X))-derivative, transfer)))
    sigma = mm(model.up[k-1], transfer)
    if k == 1:
      effective_source = sources[ids[0]].real+2*mm(model.up[0], propagated).real
    else:
      source = sources[ids[k-1]]+mm(model.up[k-1], propagated)
    transfers[k] = transfer
    factors[k] = (lu, piv)
  K = model.A[0].real+2*sigma.real
  K1 = 2*derivative.real
  return K, K1, effective_source, transfers, factors, min(conditioning)


def residue(model, z, sources, L, space, ids, observable):
  hermitian_error = np.max(abs(sources-sources[space['swap']].conj()))
  if hermitian_error > 1e-12:
    raise ValueError('Real-axis conjugate-source recurrence needs Hermitian inputs')
  K, K1, source, transfers, factors, rcond = eliminate(model, z, sources, ids)
  values, left, right = eig(K, left=True, right=True)
  index = nearest_real_index(values, z)
  mu = float(values[index].real)
  if abs(mu-z) > 1e-9:
    raise RuntimeError(f'Input is not a self-consistent root: {mu-z}')
  r, l = right[:, index].real, left[:, index].real
  l /= np.dot(l, r)
  denominator = float(1-np.dot(l, mm(K1, r[:, None])[:, 0]))
  if abs(denominator) < 1e-10:
    raise RuntimeError('Unresolved multiple root: derivative pairing is too small')
  visibility = float(np.dot(observable[ids[0]], r))
  numerator = visibility*np.dot(l, source)
  population_numerator = visibility*np.dot(l, sources[ids[0]].real)
  a = numerator/denominator

  R, W = np.zeros(L.shape[0], complex), np.zeros(L.shape[0], complex)
  R[ids[0]], W[ids[0]] = r, l
  rp, lp = r.astype(complex), l.astype(complex)
  for k in range(1, model.qmax+1):
    rp = mm(transfers[k], rp[:, None])[:, 0]
    lp = lu_solve(factors[k], np.dot(lp, model.up[k-1]), trans=1)
    R[ids[k]], W[ids[k]] = rp, lp
    R[space['swap'][ids[k]]] = rp.conj()
    W[space['swap'][ids[k]]] = lp.conj()
  pairing = np.dot(W, R)
  direct = np.dot(observable, R)*np.dot(W, sources)/pairing
  diagnostics = dict(
    scalar_root_residual=mu-z,
    right_residual=float(np.linalg.norm(L@R-z*R)/np.linalg.norm(R)),
    left_residual=float(np.linalg.norm(L.T@W-z*W)/np.linalg.norm(W)),
    pairing_relative_error=float(abs(pairing-denominator)/abs(denominator)),
    source_contraction_relative_error=float(np.linalg.norm(np.dot(W, sources)-np.dot(l, source))/
      max(np.linalg.norm(np.dot(W, sources)), 1e-30)),
    back_substituted_residue_relative_error=float(np.linalg.norm(direct-a)/max(np.linalg.norm(a), 1e-30)),
    minimum_denominator_reciprocal_condition_1=rcond,
    other_effective_eigenvalue_distance=float(np.min(abs(np.delete(values, index)-mu))))
  if max(diagnostics[k] for k in ('right_residual', 'left_residual')) > 1e-8:
    raise RuntimeError(f'Full-mode reconstruction failed: {diagnostics}')
  if diagnostics['pairing_relative_error'] > 1e-8:
    raise RuntimeError(f'Derivative pairing failed: {diagnostics}')
  return dict(root=z, derivative_pairing=denominator, numerator=numerator.tolist(),
    residues=a.tolist(), population_source_only=(population_numerator/denominator).tolist(),
    without_derivative_normalization=numerator.tolist(), diagnostics=diagnostics)


def archived_reference(N, detuning, sources, observable):
  path = ROOT/f'data/biorthogonal_audit/N{N}_float64_s17_t-0.76.npz'
  record = json.loads(path.with_suffix('.json').read_text())
  with np.load(path, allow_pickle=False) as data:
    if abs(float(data['detuning'])-detuning) > 1e-12:
      raise RuntimeError('Archived projection is at a different detuning')
    R, W = data['R'], data['W']
  values, transform = eig(np.asarray(record['H']))
  order = np.argsort(values.real)
  transform = transform[:, order]
  eigenR = mm(R, transform)
  eigenW = mm(W, np.linalg.inv(transform).T)
  a = np.dot(observable, eigenR)[:, None]*mm(eigenW.T, sources)
  return real_checked(values[order], 'reference eigenvalues'), real_checked(a, 'reference residues')


def resolvent_control():
  # The dense solve does not use a coherence recurrence or archived eigenmodes.
  N, detuning, z = 4, 3.2, .37
  L0, L1, space = generators(N)
  L = L0+detuning*L1
  q = space['ne'][space['row']]-space['ne'][space['col']]
  ids = [np.flatnonzero(q == k) for k in range(N+1)]
  source, _ = inputs(N, detuning, space)
  model = Fraction(N, detuning)
  K, _, effective, _, _, _ = eliminate(model, z, source, ids)
  reduced = solve(z*np.eye(len(K))-K, effective)
  full = solve(z*np.eye(L.shape[0])-L.toarray(), source)
  error = float(np.linalg.norm(reduced-full[ids[0]])/np.linalg.norm(full[ids[0]]))
  if error > 1e-11:
    raise RuntimeError(f'Small-N resolvent control failed: {error}')
  return dict(N=N, detuning=detuning, z=z, relative_error=error)


def run(N, dynamics=False):
  start = time.monotonic()
  frozen_path = ROOT/f'data/coherence_mechanism/fraction_N{N}_q{N}.json'
  frozen = json.loads(frozen_path.read_text())
  detuning = frozen['detuning']
  L0, L1, space = generators(N)
  L = L0+detuning*L1
  observable = space['trace']*space['ne'][space['row']]/N
  q = space['ne'][space['row']]-space['ne'][space['col']]
  ids = [np.flatnonzero(q == k) for k in range(N+1)]
  source, preparation = inputs(N, detuning, space)
  model = Fraction(N, detuning)
  modes = [residue(model, z, source, L, space, ids, observable) for z in frozen['roots']]
  computed = np.array([r['residues'] for r in modes])
  eigenvalues, reference = archived_reference(N, detuning, source, observable)
  relative_error = np.linalg.norm(computed-reference, axis=0)/np.linalg.norm(reference, axis=0)
  if max(relative_error) > 1e-6:
    raise RuntimeError(f'Independent full-mode residue comparison failed: {relative_error}')
  population_only = np.array([r['population_source_only'] for r in modes])
  without_derivative = np.array([r['without_derivative_normalization'] for r in modes])
  g = np.diff(frozen['roots'])[0]/2
  result = dict(N=N, detuning=detuning, preparations=list(PREPARATIONS),
    preparation=preparation, roots=frozen['roots'], modes=modes,
    reference_eigenvalues=eigenvalues.tolist(), reference_residues=reference.tolist(),
    residue_relative_error=relative_error.tolist(),
    omitted_source_relative_error=(np.linalg.norm(population_only-reference, axis=0)/
      np.linalg.norm(reference, axis=0)).tolist(),
    omitted_derivative_relative_error=(np.linalg.norm(without_derivative-reference, axis=0)/
      np.linalg.norm(reference, axis=0)).tolist(),
    pair_combinations=dict(A=np.sum(computed, axis=0).tolist(),
      B=(g*(computed[1]-computed[0])).tolist(), mean_rate=float(np.mean(frozen['roots'])), g=float(g)))
  if N in (32, 40):
    name = f'modes_N{N}{"_quick" if N == 40 else ""}.json'
    old = json.loads((ROOT/'data/pra_spectral'/name).read_text())
    central = min(old['records'], key=lambda r: abs(r['offset']))
    old_residues = np.array(central['residues'])
    error = np.linalg.norm(computed[:, :2]-old_residues, axis=0)/np.linalg.norm(old_residues, axis=0)
    if max(error) > 1e-6:
      raise RuntimeError(f'Archived physical residue comparison failed: {error}')
    result['archived_physical_relative_error'] = error.tolist()
  if dynamics:
    times = np.linspace(0., 12., 121)
    evolved = expm_multiply(L, source[:, 2], start=0., stop=12., num=len(times),
      traceA=L.diagonal().sum())
    signal = real_checked(np.dot(evolved, observable), 'coherence difference signal')
    pair = np.dot(np.exp(times[:, None]*np.array(frozen['roots'])[None, :]), computed[:, 2])
    slope = real_checked(np.dot(observable, L@source[:, 2]), 'initial slope').item()
    if abs(slope-NF) > 1e-10 or abs(signal[0]) > 1e-12:
      raise RuntimeError('Identical-population preparation check failed')
    result['dynamics'] = dict(times=times.tolist(), full=signal.tolist(), two_real_modes=pair.tolist(),
      initial_slope=slope, exact_initial_slope=NF,
      peak_absolute_signal=float(np.max(abs(signal))),
      peak_time=float(times[np.argmax(abs(signal))]),
      maximum_trace_error=float(np.max(abs(np.dot(evolved, space['trace'])))),
      pair_error_6_to_12=float(np.linalg.norm(pair[times >= 6]-signal[times >= 6])/
        np.linalg.norm(signal[times >= 6])))
  result['seconds'] = time.monotonic()-start
  return result


def main(args):
  DATA.mkdir(exist_ok=True)
  output = DATA/args.output
  if output.exists():
    raise FileExistsError(f'Preserve previous evidence; choose a new --output: {output}')
  control = resolvent_control()
  records = []
  for N in args.sizes:
    record = run(N, dynamics=N == args.dynamics_size)
    records.append(record)
    print(json.dumps({k: record[k] for k in ('N', 'residue_relative_error',
      'omitted_source_relative_error', 'omitted_derivative_relative_error', 'seconds')}), flush=True)
  source_paths = [Path(__file__), *(ROOT/'code'/f'{name}.py' for name in
    ('coherence_fraction', 'biorthogonal_generator', 'biorthogonal_audit', 'pra_modes', 'mf'))]
  for N in args.sizes:
    source_paths.extend([ROOT/f'data/coherence_mechanism/fraction_N{N}_q{N}.json',
      ROOT/f'data/biorthogonal_audit/N{N}_float64_s17_t-0.76.json',
      ROOT/f'data/biorthogonal_audit/N{N}_float64_s17_t-0.76.npz'])
    if N in (32, 40):
      source_paths.append(ROOT/f'data/pra_spectral/modes_N{N}{"_quick" if N == 40 else ""}.json')
  result = dict(created='2026-10-30', small_N_resolvent_control=control, records=records,
    source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
    environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
    interpretation='Finite-size Schur response with dressed source and derivative normalization; '
      'no new general residue theorem, asymptotic 1/g law, or action prediction.')
  output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
  print(json.dumps(dict(output=str(output), small_N_resolvent_control=control)), flush=True)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--sizes', type=int, nargs='+', default=[24, 32, 40, 48])
  parser.add_argument('--dynamics-size', type=int, default=32)
  parser.add_argument('--output', default='response_2026-10-30.json')
  main(parser.parse_args())
