"""Reproduce the three small independent physical-effects review checks.

Run with PYTHONDONTWRITEBYTECODE=1, VECLIB_MAXIMUM_THREADS=1,
OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and python3 -W error.
Only the selected output JSON is written; existing outputs are refused.
No full transient scan or large-N sensitivity calculation is repeated.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import qr, svdvals
from scipy.sparse import bmat
from scipy.sparse.linalg import eigs, expm_multiply

from biorthogonal_audit import project
from biorthogonal_generator import generators
from physical_effects import preparations

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/physical_effects'
THREAD_VARIABLES = ('VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS')


def sha256(path):
  return hashlib.sha256(path.read_bytes()).hexdigest()


def complex_pairs(values):
  return [[float(z.real), float(z.imag)] for z in values]


def subspace_overlap(X, Y):
  # Principal-angle cosines in the common stored coordinate convention.
  # This uses a conjugate Euclidean inner product for subspace geometry;
  # the spectral left vectors themselves remain the algebraic transpose dual.
  qx = qr(X, mode='economic')[0]
  qy = qr(Y, mode='economic')[0]
  gram = np.einsum('ji,jk->ik', qx.conj(), qy, optimize=False)
  return float(min(svdvals(gram)))


def subspace_tracking(scan, archive):
  N, star, width = scan['N'], scan['crossing'], scan['mixing_width']
  if N not in (24, 32, 40, 48):
    raise ValueError('This review is limited to N=24,32,40,48 subspace checks')
  L0, L1, space = generators(N)
  with np.load(archive, allow_pickle=False) as data:
    R0, W0 = data['R'], data['W']
    if abs(float(data['detuning'])-star) > 1e-12:
      raise RuntimeError('Archived center does not match the scan')
  rows = []
  for sign in (-1, 1):
    previous_right, previous_left = R0, W0
    for step in (1, 2, 3, 4):
      offset = sign*step*width
      projection_seed = 72700+N+10*sign+step
      candidate_seed = 899+N+step
      audit, R, W, L = project(L0, L1, space, star+offset, projection_seed, -.76)
      attempts = []
      # A shift nearly on the target doublet can leave distant Ritz vectors
      # inaccurate even after ARPACK reports convergence. Preserve rejected
      # residuals and retry away from that doublet, without relaxing the bound.
      for candidate_shift, tolerance, ncv in ((-.76, 2e-13, None),
                                             (-.76, 2e-14, 40),
                                             (-.9, 2e-13, None)):
        values, vectors = eigs(L, k=8, sigma=candidate_shift, tol=tolerance, ncv=ncv,
          v0=np.random.default_rng(candidate_seed).normal(size=L.shape[0]))
        residuals = np.linalg.norm(L@vectors-vectors*values, axis=0)/np.linalg.norm(vectors, axis=0)
        accepted = bool(max(residuals) < 1e-9)
        attempts.append(dict(shift=candidate_shift, tolerance=tolerance, ncv=ncv,
          eigenvalues=complex_pairs(values), residuals=residuals.tolist(), accepted=accepted))
        if accepted:
          break
        print(json.dumps(dict(check='candidate_residual_rejected', N=N, step=sign*step,
          shift=candidate_shift, tolerance=tolerance, ncv=ncv,
          maximum_residual=float(max(residuals)))), flush=True)
      else:
        raise RuntimeError('Returned candidate eigenmode failed all residual checks')
      target = np.asarray(audit['lam'])
      matched = np.array([int(np.argmin(abs(values-lam))) for lam in target])
      if len(set(matched.tolist())) != 2:
        raise RuntimeError('Candidate matching is not one-to-one')
      others = np.delete(values, matched)
      outside_distance = float(min(abs(lam-other) for lam in target for other in others))
      real_ids = [i for i in range(len(values)) if abs(values[i].imag) < 1e-7]
      candidates = [(subspace_overlap(previous_right, vectors[:, ids]), list(ids))
        for ids in itertools.combinations(real_ids, 2)]
      candidates.sort(reverse=True)
      best, runner_up = candidates[0], candidates[1]
      if set(best[1]) != set(matched.tolist()) or best[0] <= runner_up[0]:
        raise RuntimeError('The target pair is not the unique best continuation')
      saved = min(scan['records'], key=lambda r: abs(r['offset']-offset))
      saved_error = None
      if abs(saved['offset']-offset) < 1e-12:
        saved_error = float(max(abs(target-np.asarray(saved['modes']['eigenvalues']))))
        if saved_error >= 1e-9:
          raise RuntimeError('Repeated eigenvalues disagree with the saved scan')
      rows.append(dict(step=sign*step, offset=float(offset), detuning=float(star+offset),
        eigenvalues=target.tolist(), projection_seed=projection_seed, candidate_seed=candidate_seed,
        neighbor_right_min_cosine=subspace_overlap(previous_right, R),
        neighbor_left_min_cosine=subspace_overlap(previous_left, W),
        center_right_min_cosine=subspace_overlap(R0, R),
        center_left_min_cosine=subspace_overlap(W0, W),
        candidate_eigenvalues=complex_pairs(values), candidate_residuals=residuals.tolist(),
        candidate_attempts=attempts,
        target_candidate_indices=matched.tolist(), real_candidate_indices=real_ids,
        candidate_pair_count=len(candidates), best_pair_indices=best[1], best_pair_overlap=best[0],
        runner_up_pair_indices=runner_up[1], runner_up_overlap=runner_up[0],
        best_pair_is_target=True, unique_best=True,
        target_to_other_returned_candidate_min_distance=outside_distance,
        invariant_right_error=audit['diagnostics']['invariant_right_error'],
        invariant_left_error=audit['diagnostics']['invariant_left_error'],
        saved_eigenvalue_max_error=saved_error))
      previous_right, previous_left = R, W
      print(json.dumps(dict(check='subspace_tracking', N=N, step=sign*step,
        neighbor_right_min_cosine=rows[-1]['neighbor_right_min_cosine'],
        neighbor_left_min_cosine=rows[-1]['neighbor_left_min_cosine'],
        best_pair_is_target=True, saved_eigenvalue_max_error=saved_error)), flush=True)
  return dict(N=N, dimension=L0.shape[0], crossing=star, mixing_width=width,
    projection_shift=-.76, candidates_per_point=8, candidate_residual_threshold=1e-9,
    candidate_retry_policy='Try shift -0.76 with tolerance 2e-13; on residual failure '
      'try tolerance 2e-14 with ncv=40; if still failing, try shift -0.9 with tolerance '
      '2e-13. Record every attempt; only residual-accepted candidates determine identity.',
    steps=rows, summary=dict(
      rejected_candidate_attempts=sum(not a['accepted'] for r in rows for a in r['candidate_attempts']),
      minimum_neighbor_right_cosine=min(r['neighbor_right_min_cosine'] for r in rows),
      minimum_neighbor_left_cosine=min(r['neighbor_left_min_cosine'] for r in rows),
      maximum_runner_up_overlap=max(r['runner_up_overlap'] for r in rows),
      minimum_target_to_returned_others=min(r['target_to_other_returned_candidate_min_distance'] for r in rows),
      maximum_candidate_residual=max(max(r['candidate_residuals']) for r in rows),
      maximum_saved_eigenvalue_error=max(r['saved_eigenvalue_max_error'] for r in rows
        if r['saved_eigenvalue_max_error'] is not None)),
    scope='Numerical continuation among eight returned candidates at each sampled point; '
      'not a proof of isolation from the entire spectrum. Cosines use stored-coordinate '
      'Euclidean subspace geometry, not physical-state fidelity.')


def parameter_differences(scan):
  N, detuning, parameters = 2, scan['crossing'], scan['preparations']
  L0, derivative, space = generators(N)
  L = L0+detuning*derivative
  source, _ = preparations(parameters, space)
  dim = len(source)
  observable = space['trace']*space['ne'][space['row']]/N
  augmented = bmat([[L, None, None], [derivative, L, None],
    [None, 2*derivative, L]], format='csc')
  initial = np.vstack([source, np.zeros_like(source), np.zeros_like(source)])

  def density(delta, t):
    matrix = (L0+delta*derivative).toarray()
    def rhs(time_value, vector):
      return np.einsum('ij,jk->ik', matrix, vector.reshape(dim, 2), optimize=False).ravel()
    ode = solve_ivp(rhs, (0., t), source.ravel(), method='DOP853',
      rtol=5e-14, atol=1e-15, t_eval=[t])
    if not ode.success:
      raise RuntimeError(ode.message)
    signal = np.einsum('i,ij->j', observable, ode.y[:, -1].reshape(dim, 2), optimize=False)
    if max(abs(signal.imag)) > 1e-11:
      raise RuntimeError('Independent ODE density is not real')
    return signal.real, ode.nfev

  rows = []
  for t in (.1, 1., 3.):
    tangent = expm_multiply(t*augmented, initial, traceA=t*3*L.diagonal().sum())
    signals = [np.einsum('i,ij->j', observable, tangent[i*dim:(i+1)*dim],
      optimize=False) for i in range(3)]
    if max(float(max(abs(s.imag))) for s in signals) > 1e-11:
      raise RuntimeError('Augmented density derivatives are not real')
    signals = [s.real for s in signals]
    for h in (.004, .002):
      samples = [density(detuning+k*h, t) for k in (-2, -1, 0, 1, 2)]
      fm2, fm, f0, fp, fp2 = [s[0] for s in samples]
      first = (fm2-8*fm+8*fp-fp2)/(12*h)
      second = (-fm2+16*fm-30*f0+16*fp-fp2)/(12*h*h)
      first_error = float(max(abs(first-signals[1])))
      second_error = float(max(abs(second-signals[2])))
      if max(first_error, second_error) >= 1e-8:
        raise RuntimeError('Independent five-point difference disagrees with tangents')
      rows.append(dict(time=t, step=h, augmented_density=signals[0].tolist(),
        augmented_first_derivative=signals[1].tolist(), augmented_second_derivative=signals[2].tolist(),
        finite_difference_first=first.tolist(), finite_difference_second=second.tolist(),
        shifted_densities=[s[0].tolist() for s in samples], ode_evaluations=[s[1] for s in samples],
        first_max_error=first_error, second_max_error=second_error))
  fine = [r for r in rows if r['step'] == .002]
  return dict(N=N, detuning=detuning, preparation_reference_N=scan['N'],
    preparations=parameters, dimension=dim, augmented_dimension=3*dim,
    derivative_convention='y_j = partial_Delta^j rho; y_1(0)=y_2(0)=0; '
      'y_2_dot=L y_2+2 L1 y_1; Taylor coefficient is h^2/2.',
    independent_solver=dict(method='DOP853', rtol=5e-14, atol=1e-15,
      finite_difference_order=4, shift_multipliers=[-2, -1, 0, 1, 2]),
    checks=rows, summary=dict(
      maximum_first_error=max(r['first_max_error'] for r in rows),
      maximum_second_error=max(r['second_max_error'] for r in rows),
      fine_step_maximum_first_error=max(r['first_max_error'] for r in fine),
      fine_step_maximum_second_error=max(r['second_max_error'] for r in fine)))


def initial_second_derivatives(scan):
  detuning, parameters = scan['crossing'], scan['preparations']
  rows = []
  for N in (1, 2, 4):
    L0, derivative, space = generators(N)
    L = L0+detuning*derivative
    source, _ = preparations(parameters, space)
    observable = space['trace']*space['ne'][space['row']]/N
    n, u, v = parameters['fold']
    analytic = 3*v-2*(detuning-10*(N-1)/N*n)*u
    value = np.dot(observable, L@(L@source[:, 0]))
    if abs(value.imag) > 1e-11 or abs(value.real-analytic) >= 1e-10:
      raise RuntimeError('Initial second time derivative disagrees with its identity')
    rows.append(dict(N=N, analytic=float(analytic), full_generator=float(value.real),
      absolute_difference=float(abs(value.real-analytic))))
  return dict(detuning=detuning, fold_state=parameters['fold'], gamma=1., omega=2., interaction=10.,
    formula='s_second(0) = 3*gamma*Omega*v_f/2 - '
      'Omega*(Delta - V*(N-1)*n_f/N)*u_f',
    preparation='rho_fold minus its global excitation-number dephasing', checks=rows,
    maximum_absolute_difference=max(r['absolute_difference'] for r in rows))


def main(args):
  for name in THREAD_VARIABLES:
    if os.environ.get(name) != '1':
      raise RuntimeError(f'Set {name}=1 before starting Python')
  if not sys.dont_write_bytecode:
    raise RuntimeError('Set PYTHONDONTWRITEBYTECODE=1 before starting Python')
  if 'error' not in sys.warnoptions:
    raise RuntimeError('Run Python with -W error; warnings must not be suppressed')
  output = DATA/args.output
  if output.exists():
    raise FileExistsError(f'Preserve evidence; choose a new output: {output}')
  effects_path = DATA/'effects_2026-11-01.json'
  sensitivity_path = DATA/'sensitivity_2026-11-01.json'
  archives = {N: ROOT/f'data/biorthogonal_audit/N{N}_float64_s17_t-0.76.npz'
    for N in (24, 32, 40, 48)}
  dependencies = [ROOT/'code'/f'{name}.py' for name in (
    'physical_effects', 'physical_effects_sensitivity', 'biorthogonal_audit',
    'biorthogonal_generator', 'pra_modes', 'mf', 'coherence_response', 'coherence_fraction')]
  paths = [Path(__file__).resolve(), effects_path, sensitivity_path, *archives.values(), *dependencies]
  hashes = {str(p.relative_to(ROOT)): sha256(p) for p in paths}
  effects = json.loads(effects_path.read_text())
  sensitivity = json.loads(sensitivity_path.read_text())
  if sensitivity['input_sha256'] != hashes[str(effects_path.relative_to(ROOT))]:
    raise RuntimeError('Sensitivity input hash does not identify the supplied effects record')
  scans = {N: next(s for s in effects['scans'] if s['N'] == N) for N in archives}
  scan32 = scans[32]
  start = time.monotonic()
  tracking = [subspace_tracking(scans[N], archives[N]) for N in archives]
  differences = parameter_differences(scan32)
  initial = initial_second_derivatives(scan32)
  for p in paths:
    if sha256(p) != hashes[str(p.relative_to(ROOT))]:
      raise RuntimeError(f'Input or source changed during the review: {p}')
  result = dict(created='2026-11-01', status='passed', seconds=time.monotonic()-start,
    script=str(Path(__file__).resolve().relative_to(ROOT)),
    script_sha256=hashes[str(Path(__file__).resolve().relative_to(ROOT))],
    inputs=dict(effects=dict(path=str(effects_path.relative_to(ROOT)), sha256=sha256(effects_path)),
      sensitivity=dict(path=str(sensitivity_path.relative_to(ROOT)), sha256=sha256(sensitivity_path)),
      sensitivity_role='Frozen reviewed input and provenance; no large-N sensitivity evolution rerun.'),
    source_sha256=hashes,
    environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
      platform=platform.platform(), threads={name: os.environ[name] for name in THREAD_VARIABLES},
      warning_options=sys.warnoptions, bytecode_disabled=sys.dont_write_bytecode),
    earlier_attempt_note=dict(
      scope='Historical execution note from the preceding interactive review, not new numerical data.',
      observed='Direct NumPy dense matmul in subspace-overlap evaluation emitted divide-by-zero, '
        'overflow and invalid-value RuntimeWarnings. A subsequent dense scipy.linalg.expm '
        'attempt stopped at its internal matrix-squaring matmul under -W error on this '
        'Accelerate-backed environment.',
      resolution='Warnings were not suppressed. Dense overlap contractions were replaced by '
        'einsum(optimize=False); the independent finite-difference evolution uses DOP853 '
        'with einsum, and augmented tangents use sparse expm_multiply. The successful '
        'versions are rerun here under -W error. Failed-attempt values are not used.'),
    subspace_tracking=tracking, parameter_finite_difference=differences,
    initial_second_time_derivative=initial,
    interpretation='The N24 local subspace continuation is repeated and extended by the same '
      'procedure to N32/40/48. The previously executed N2 parameter finite differences and '
      'N1/2/4 initial second-time-derivative checks are repeated. No main transient scan is '
      'repeated. This is numerical evidence, not full-spectrum isolation or an asymptotic proof.')
  with output.open('x') as stream:
    stream.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
  print(json.dumps(dict(output=str(output), bytes=output.stat().st_size, status=result['status'],
    subspace_summary=[dict(N=r['N'], **r['summary']) for r in tracking],
    finite_difference_summary=differences['summary'],
    initial_second_derivative_max_error=initial['maximum_absolute_difference'],
    seconds=result['seconds']), indent=2), flush=True)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--output', default='review_2026-11-01.json')
  main(parser.parse_args())
