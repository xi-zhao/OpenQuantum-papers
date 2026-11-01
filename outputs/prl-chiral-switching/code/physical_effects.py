"""Fixed-preparation tests of finite-rate mode mixing in density signals.

Archived crossings are inputs. Full transients use the original PI generator;
no rates or amplitudes are fitted. Mean-field global-dephasing controls evolve
an ensemble of phases, not its initially averaged Bloch vector. New output
files only; previous scientific evidence is never overwritten.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import eig, eigvalsh
from scipy.sparse.linalg import expm_multiply
from biorthogonal_audit import mm, project
from biorthogonal_generator import generators
from coherence_response import real_checked
from pra_modes import product_state, fixed, jacobian, DF, NF, OMEGA, V
from mf import roots

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/physical_effects'
PREPARATIONS = ('fold_coherence_difference', 'fixed_high_difference')
TIMES = np.linspace(0., 12., 121)


def crossing(N):
  return json.loads((ROOT/f'data/biorthogonal_audit/N{N}_float64_s17_t-0.76.json').read_text())


def preparation_parameters(detuning):
  high = fixed(detuning, roots(detuning, OMEGA, V)[-1])
  values, vectors = eig(jacobian(detuning, high))
  ids = np.flatnonzero(abs(values.imag) < 1e-10)
  direction = vectors[:, ids[np.argmax(values[ids].real)]].real
  direction *= np.sign(direction[0])/np.linalg.norm(direction)
  return dict(fold=fixed(DF, NF).tolist(), high=high.tolist(),
    high_plus=(high+.005*direction).tolist(), high_minus=(high-.005*direction).tolist(),
    high_direction=direction.tolist(), high_step=.005,
    seed_detuning=float(detuning))


def trace_distance(vector, space):
  full = np.zeros((space['d'], space['d']), complex)
  full[space['row'], space['col']] = vector
  norm = 0.
  for tj in range(space['N'], -1, -2):
    ids = np.flatnonzero(space['twoj'] == tj)
    norm += np.sum(abs(eigvalsh(full[np.ix_(ids, ids)])))
  # Right blocks already sum multiplicities: no extra representation weight.
  return float(norm/2)


def preparations(parameters, space):
  folded = product_state(parameters['fold'], space)
  dephased = np.where(space['trace'].astype(bool), folded, 0.)
  high = product_state(parameters['high_plus'], space)-product_state(parameters['high_minus'], space)
  source = np.column_stack([folded-dephased, high])
  if np.max(abs(np.dot(space['trace'], source))) > 1e-11:
    raise RuntimeError('Initial difference is not traceless')
  if np.max(abs(source-source[space['swap']].conj())) > 1e-12:
    raise RuntimeError('Initial difference is not Hermitian')
  return source, [trace_distance(source[:, i], space) for i in range(2)]


def classical_signals(parameters, detuning, coupling, phases=64):
  """Independent Bloch trajectories; coupling=V*(N-1)/N for factorization.

  Uniform phase mixing represents the global-dephased product. Nonlinear
  evolution is applied before averaging. coupling=0 gives an exact one-body
  reference for the measured density, even when the input has correlations.
  """
  fold = np.asarray(parameters['fold'])
  angles = 2*np.pi*np.arange(phases)/phases
  coherence = (fold[1]+1j*fold[2])*np.exp(-1j*angles)
  ensemble = np.column_stack([np.full(phases, fold[0]), coherence.real, coherence.imag])
  initial = np.vstack([fold, ensemble, parameters['high_plus'], parameters['high_minus']])
  def rhs(t, vector):
    n, u, v = vector.reshape(-1, 3).T
    d = detuning-coupling*n
    return np.column_stack([-n-OMEGA*v, -d*v-u/2,
      d*u-OMEGA*(1-2*n)/2-v/2]).ravel()
  ode = solve_ivp(rhs, (0., TIMES[-1]), initial.ravel(), t_eval=TIMES,
    method='DOP853', rtol=2e-11, atol=2e-13)
  if not ode.success:
    raise RuntimeError(ode.message)
  states = ode.y.reshape(len(initial), 3, len(TIMES))
  signal = np.column_stack([states[0, 0]-np.mean(states[1:phases+1, 0], axis=0),
    states[-2, 0]-states[-1, 0]])
  return signal, dict(evaluations=ode.nfev, phases=phases, coupling=coupling,
    maximum_bloch_radius=float(np.max(np.sqrt((2*states[:, 0]-1)**2+
      4*states[:, 1]**2+4*states[:, 2]**2))))


def raw_dynamics(L, source, observable, space):
  evolved = expm_multiply(L, source, start=0., stop=TIMES[-1], num=len(TIMES),
    traceA=L.diagonal().sum())
  signal = real_checked(np.einsum('i,tij->tj', observable, evolved, optimize=False), 'density signal')
  trace_error = float(np.max(abs(np.einsum('i,tij->tj', space['trace'], evolved, optimize=False))))
  if trace_error > 1e-9:
    raise RuntimeError(f'Trace conservation failed: {trace_error}')
  return signal, trace_error


def modal_signal(audit, R, W, source, observable):
  values, T = eig(np.asarray(audit['H']))
  order = np.argsort(values.real)
  values = real_checked(values[order], 'real pair roots')
  T = T[:, order]
  right = mm(R, T)
  left = mm(W, np.linalg.inv(T).T)
  residues = real_checked(np.dot(observable, right)[:, None]*mm(left.T, source), 'real pair residues')
  mean_rate, g = np.mean(values), np.diff(values)[0]/2
  if g <= 0 or values[-1] >= 0:
    raise RuntimeError('Pair is not distinct and decaying')
  A, B = residues.sum(axis=0), g*(residues[1]-residues[0])
  signal = np.exp(mean_rate*TIMES[:, None])*(
    np.cosh(g*TIMES[:, None])*A+np.sinh(g*TIMES[:, None])/g*B)
  collapsed = np.exp(mean_rate*TIMES[:, None])*(A+TIMES[:, None]*B)
  direct = mm(np.exp(TIMES[:, None]*values[None, :]), residues)
  if np.max(abs(signal-direct)) > 1e-9:
    raise RuntimeError('Hyperbolic pair identity failed')
  return signal, collapsed, dict(eigenvalues=values.tolist(), residues=residues.tolist(),
    mean_rate=float(mean_rate), g=float(g), A=A.tolist(), B=B.tolist(),
    g_over_decay=float(g/abs(mean_rate)),
    log10_common_envelope_at_inverse_g=float(mean_rate/g/np.log(10)),
    pair_identity_error=float(np.max(abs(signal-direct))))


def metrics(signal, reference):
  error = signal-reference
  result = dict(maximum_absolute_error=np.max(abs(error), axis=0).tolist())
  for lower in (0., 2., 6.):
    mask = TIMES >= lower
    denominator = np.linalg.norm(reference[mask], axis=0)
    if np.any(denominator == 0):
      raise RuntimeError('Zero reference norm')
    result[f'relative_l2_{lower:g}_to_12'] = (np.linalg.norm(error[mask], axis=0)/denominator).tolist()
  return result


def row(N, offset, frozen, parameters, source, observable, L0, L1, space, phases):
  detuning = frozen['detuning']+offset
  L = (L0+detuning*L1).tocsc()
  signal, trace_error = raw_dynamics(L, source, observable, space)
  if offset == 0.:
    audit = frozen
    with np.load(ROOT/f'data/biorthogonal_audit/N{N}_float64_s17_t-0.76.npz', allow_pickle=False) as data:
      R, W = data['R'], data['W']
    modal_origin = 'archived exact two-dimensional subspace'
  else:
    audit, R, W, _ = project(L0, L1, space, detuning, 4100+N, -.76)
    modal_origin = 'fresh full-generator two-dimensional subspace'
  residuals = audit['diagnostics']['residual_right']+audit['diagnostics']['residual_left']
  if max(residuals) > 1e-8 or audit['diagnostics']['biorthogonality_error'] > 1e-8:
    raise RuntimeError('Two-mode projection validation failed')
  pair, collapsed, modes = modal_signal(audit, R, W, source, observable)
  factorized, mfcheck = classical_signals(parameters, detuning, V*(N-1)/N, phases)
  free, freecheck = classical_signals(parameters, detuning, 0., 4)
  matched = np.zeros_like(free)
  for i, state in enumerate(('fold', 'high')):
    shifted = detuning-V*(N-1)/N*parameters[state][0]
    matched[:, i] = classical_signals(parameters, shifted, 0., 4)[0][:, i]
  return dict(N=N, detuning=float(detuning), offset=float(offset),
    full=signal.tolist(), pair=pair.tolist(), coalesced_shape=collapsed.tolist(),
    meanfield_phase_ensemble=factorized.tolist(), noninteracting=free.tolist(),
    noninteracting_frozen_hartree=matched.tolist(), modes=modes,
    full_peak=np.max(abs(signal), axis=0).tolist(),
    full_peak_time=TIMES[np.argmax(abs(signal), axis=0)].tolist(),
    pair_error=metrics(pair, signal), meanfield_error=metrics(factorized, signal),
    noninteracting_error=metrics(free, signal),
    finite_splitting_effect=dict(
      maximum_absolute_change=np.max(abs(pair-collapsed), axis=0).tolist(),
      change_over_full_l2=(np.linalg.norm(pair-collapsed, axis=0)/np.linalg.norm(signal, axis=0)).tolist(),
      change_over_pair_l2=(np.linalg.norm(pair-collapsed, axis=0)/np.linalg.norm(pair, axis=0)).tolist()),
    trace_error=trace_error, modal_origin=modal_origin,
    spectral_maximum_residual=float(max(residuals)),
    meanfield_diagnostics=mfcheck, noninteracting_diagnostics=freecheck)


def scan(N, phases):
  start = time.monotonic()
  frozen = crossing(N)
  width = 2*frozen['g']/abs(frozen['H1'][0][0]-frozen['H1'][1][1])
  parameters = preparation_parameters(frozen['detuning'])
  L0, L1, space = generators(N)
  observable = space['trace']*space['ne'][space['row']]/N
  source, distances = preparations(parameters, space)
  offsets = sorted(set([-0.06, -4*width, -width, 0., width, 4*width, .06]))
  records = []
  for offset in offsets:
    record = row(N, offset, frozen, parameters, source, observable, L0, L1, space, phases)
    records.append(record)
    print(json.dumps(dict(N=N, offset=offset, full_peak=record['full_peak'],
      g=record['modes']['g'], seconds=time.monotonic()-start)), flush=True)
  def at(offset):
    return min(records, key=lambda r: abs(r['offset']-offset))
  center, lo, hi = at(0.), at(-width), at(width)
  narrow = {}
  for key in ('full', 'pair', 'meanfield_phase_ensemble', 'noninteracting', 'noninteracting_frozen_hartree'):
    c, a, b = (np.asarray(r[key]) for r in (center, lo, hi))
    denominator = np.linalg.norm(c, axis=0)
    narrow[key] = dict(endpoint_difference_max=np.max(abs(b-a), axis=0).tolist(),
      endpoint_difference_l2_over_center=(np.linalg.norm(b-a, axis=0)/denominator).tolist(),
      central_curvature_max=np.max(abs(b+a-2*c), axis=0).tolist(),
      central_curvature_l2_over_center=(np.linalg.norm(b+a-2*c, axis=0)/denominator).tolist())
  return dict(N=N, crossing=frozen['detuning'], mixing_width=width,
    preparations=parameters, initial_trace_distance=distances,
    initial_slope=real_checked(np.dot(observable, (L0+frozen['detuning']*L1)@source), 'initial slope').tolist(),
    narrow_comparison=narrow, records=records, seconds=time.monotonic()-start)


def size_control(sizes, detuning):
  parameters = preparation_parameters(detuning)
  records = []
  for N in sizes:
    L0, L1, space = generators(N)
    source, distance = preparations(parameters, space)
    observable = space['trace']*space['ne'][space['row']]/N
    signal, trace_error = raw_dynamics(L0+detuning*L1, source, observable, space)
    records.append(dict(N=N, full=signal.tolist(), peak=np.max(abs(signal), axis=0).tolist(),
      peak_time=TIMES[np.argmax(abs(signal), axis=0)].tolist(),
      initial_trace_distance=distance, trace_error=trace_error))
    print(json.dumps(dict(control='fixed_detuning', N=N, peak=records[-1]['peak'])), flush=True)
  return dict(detuning=detuning, preparations=parameters, records=records,
    noninteracting=classical_signals(parameters, detuning, 0., 4)[0].tolist())


def main(args):
  DATA.mkdir(exist_ok=True)
  output = DATA/args.output
  if output.exists():
    raise FileExistsError(f'Choose a new output; preserve evidence: {output}')
  records = [scan(N, args.phases) for N in args.sizes]
  fixed_control = size_control(args.control_sizes, crossing(32)['detuning']) if args.control_sizes else None
  paths = [Path(__file__), ROOT/'PHYSICAL_EFFECT_TEST_2026-11-01.md']
  paths += [ROOT/'code'/f'{name}.py' for name in ('pra_modes', 'mf',
    'biorthogonal_audit', 'biorthogonal_generator', 'coherence_response')]
  for N in set(args.sizes)|{32}:
    paths.extend(ROOT/f'data/biorthogonal_audit/N{N}_float64_s17_t-0.76.{suffix}' for suffix in ('json','npz'))
  result = dict(created='2026-11-01', preparations=list(PREPARATIONS), times=TIMES.tolist(),
    scans=records, fixed_detuning_control=fixed_control,
    environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
    source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    interpretation='Exploratory finite-size fixed-preparation comparison. The coalesced shape is '
      'a two-pole algebraic diagnostic, not a CPTP alternative or a universal indistinguishability bound.')
  output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
  print(json.dumps(dict(output=str(output), bytes=output.stat().st_size)), flush=True)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--sizes', nargs='+', type=int, default=[24, 32, 40, 48])
  parser.add_argument('--control-sizes', nargs='*', type=int, default=[1, 2, 4, 8, 16, 24, 32, 40, 48])
  parser.add_argument('--phases', type=int, default=64)
  parser.add_argument('--output', default='effects_2026-11-01.json')
  main(parser.parse_args())
