"""Local EP2 verification for the analytic float64 family L(z) = L0 + z L1.

This is a complex commutator continuation, not H rho - rho H.dag(), and is
not a physical CPTP generator off the real axis. No coherence elimination or
2 Re expression is used. The full permutation-invariant Dicke-block operator
is used at every point; only eight nearby eigenpairs are sampled, not its
entire spectrum.

A right and a transpose-left two-dimensional invariant subspace are obtained
by block inverse iteration with QR. With W.T R = I, H = W.T L R and
B = W.T L1 R, D = ((H00-H11)/2)^2 + H01 H10. The analytic derivative D' is
computed with B; the basis-derivative commutator drops out of this similarity
invariant. Newton iterations therefore update the subspaces, rather than
finding a zero of one frozen real-axis 2x2 model. Orthonormal right bases stay
well conditioned at EP2 even though individual eigenvectors coalesce.

All complex numbers in JSON are [real, imaginary]. The 'checks_passed' fields
are explicit numerical diagnostics, not central scientific Acceptance.

Reproduce (existing Python/SciPy/QuTiP only):
  VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  PYTHONDONTWRITEBYTECODE=1 python3 code/pra_ep_verify.py \
  --sizes 24 32 40 --loop-size 32 --loop-steps 24 --stability
"""
import os
for _name in ('VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
  os.environ[_name] = '1'

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time
import traceback
import warnings

import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse as sp
from scipy.optimize import linear_sum_assignment
from scipy.sparse.linalg import LinearOperator, eigs, splu

from biorthogonal_audit import mm, project
from biorthogonal_generator import generators

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/pra_spectral'
LIMITS = dict(subspace_residual=2e-10, mode_residual=2e-9,
              discriminant_over_g2=1e-9, newton_step=2e-11,
              nilpotent_squared_ratio=1e-8, nonzero_nilpotent=1e-4,
              external_separation=0.03, circle_step_over_gap=0.4,
              closure_error=2e-8, root_repeat_difference=2e-9)


def serial(x):
  if isinstance(x, dict):
    return {str(k): serial(v) for k, v in x.items()}
  if isinstance(x, (list, tuple)):
    return [serial(v) for v in x]
  if isinstance(x, np.ndarray):
    return serial(x.tolist())
  if isinstance(x, (complex, np.complexfloating)):
    return [float(x.real), float(x.imag)]
  if isinstance(x, np.generic):
    return x.item()
  return x


def save(name, data):
  if not name.startswith('ep_'):
    raise ValueError('Only ep_* output files are allowed')
  OUT.mkdir(parents=True, exist_ok=True)
  (OUT / name).write_text(json.dumps(serial(data), indent=2, allow_nan=False) + '\n')


def sha(path):
  return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def metadata():
  return dict(python=sys.version, executable=sys.executable,
              numpy=np.__version__, scipy=scipy.__version__,
              threads={n: os.environ[n] for n in
                       ('VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS')},
              sources={str(p.relative_to(ROOT)): sha(p) for p in
                       (Path(__file__), ROOT/'code/biorthogonal_generator.py',
                        ROOT/'code/biorthogonal_audit.py')},
              limits=LIMITS,
              interpretation='Analytic commutator continuation; non-CPTP for complex z. '
                'Numerical local EP evidence, not a global enumeration of EPs or a '
                'central scientific Acceptance result.')


def qr(x):
  return la.qr(x, mode='economic', check_finite=True)[0]


def norm(x):
  return float(np.linalg.norm(x))


def disc(h):
  return ((h[0, 0] - h[1, 1]) / 2)**2 + h[0, 1] * h[1, 0]


def disc_derivative(h, b):
  return ((h[0, 0] - h[1, 1]) * (b[0, 0] - b[1, 1]) / 2
          + b[0, 1] * h[1, 0] + h[0, 1] * b[1, 0])


def residuals(a, values, vectors):
  return [norm(a @ vectors[:, i] - value * vectors[:, i]) / norm(vectors[:, i])
          for i, value in enumerate(values)]


class LocalPair:
  def __init__(self, n, seed, shift_offset):
    self.n, self.seed, self.shift_offset = n, seed, shift_offset
    self.rng = np.random.default_rng(seed)
    self.l0, self.l1, self.space = generators(n, backend='float64')
    self.eye = sp.eye(self.l0.shape[0], format='csc')
    self.input_path = ROOT / f'data/biorthogonal_audit/N{n}_float64_s17_t-0.76.json'
    self.input = json.loads(self.input_path.read_text())
    self.g = self.input['g']
    slopes = np.diag(self.input['H1'])
    self.y0 = 2 * self.g / abs(slopes[0] - slopes[1])
    self.z0 = self.input['detuning'] + 1j * self.y0
    audit, r, v, _ = project(self.l0, self.l1, self.space,
                            self.input['detuning'], seed, -0.76)
    self.initial = dict(saved_input_sha256=sha(self.input_path),
                        real_axis_recalculation=audit, linear_prediction=self.z0,
                        full_dimension=self.l0.shape[0], nnz=self.l0.nnz)
    self.r, self.v = qr(r), qr(v)
    self.centre = complex(np.trace(self.input['H']) / 2
                          + 1j * self.y0 * np.trace(self.input['H1']) / 2)

  def evaluate(self, z, r, v, centre, candidates=False):
    started = time.monotonic()
    l = (self.l0 + z * self.l1).tocsc()
    shift = centre + self.shift_offset
    lu = splu(l - shift * self.eye)
    hist = []
    for iteration in range(80):
      r = qr(lu.solve(r))
      v = qr(lu.solve(v, trans='T'))
      gram = mm(v.T, r)
      w = mm(v, np.linalg.inv(gram).T)
      lr, ltw = l @ r, l.T @ w
      h = mm(w.T, lr)
      rr = norm(lr - mm(r, h)) / norm(r)
      rl = norm(ltw - mm(w, h.T)) / norm(w)
      hist.append([rr, rl])
      if max(rr, rl) < LIMITS['subspace_residual'] / 20:
        break
    converged = max(rr, rl) < LIMITS['subspace_residual'] / 20
    b = mm(w.T, self.l1 @ r)
    d, dp = disc(h), disc_derivative(h, b)
    k = h - np.trace(h) / 2 * np.eye(2)
    singular = la.svdvals(k)
    obs = dict(z=z, shift=shift, H=h, B=b, discriminant=d,
               discriminant_derivative=dp, discriminant_over_g2=abs(d)/self.g**2,
               proposed_newton_step=-d/dp,
               mean_eigenvalue=np.trace(h)/2,
               projected_eigenvalues=la.eigvals(h),
               invariant_right_residual=rr, invariant_transpose_left_residual=rl,
               inverse_iteration_converged=converged,
               inverse_iteration_residual_history=hist,
               biorthogonality_error=norm(mm(w.T, r)-np.eye(2)),
               right_basis_orthogonality_error=norm(mm(r.conj().T, r)-np.eye(2)),
               dual_gram_singular_values=la.svdvals(gram),
               dual_gram_condition=float(np.linalg.cond(gram)),
               spectral_projector_norm=float(la.svdvals(w)[0]),
               nilpotent_part=k, nilpotent_norm=norm(k),
               nilpotent_singular_values=singular,
               nilpotent_squared_ratio=norm(mm(k,k))/norm(k)**2,
               trace_preservation_defect=float(np.max(abs(self.space['trace'] @ l))))
    state = dict(r=r, v=v, w=w, h=h, l=l, lu=lu, shift=shift)
    if candidates:
      obs['spectrum'] = self.spectrum(state)
    obs['seconds'] = time.monotonic() - started
    return obs, state

  def spectrum(self, state):
    l, lu, shift, r, h = (state[k] for k in ('l', 'lu', 'shift', 'r', 'h'))
    opts = dict(k=8, sigma=shift, tol=2e-13, ncv=28, maxiter=1200)
    op = LinearOperator(l.shape, matvec=lu.solve, dtype=complex)
    opt = LinearOperator(l.shape, matvec=lambda x: lu.solve(x, trans='T'), dtype=complex)
    def random_start():
      return self.rng.normal(size=l.shape[0]) + 1j*self.rng.normal(size=l.shape[0])
    vals, rv = eigs(l, OPinv=op, v0=random_start(), **opts)
    leftvals, lv = eigs(l.T, OPinv=opt, v0=random_start(), **opts)
    order = np.argsort(abs(vals - np.trace(h)/2))
    vals, rv = vals[order], rv[:, order]
    il, jl = linear_sum_assignment(abs(vals[:, None] - leftvals[None, :]))
    left_order = jl[np.argsort(il)]
    leftvals, lv = leftvals[left_order], lv[:, left_order]
    projected = la.eigvals(h)
    ip, jp = linear_sum_assignment(abs(projected[:, None] - vals[None, :]))
    selected = jp[np.argsort(ip)]
    outside = [j for j in range(len(vals)) if j not in selected]
    separation = float(np.min(abs(vals[selected, None] - vals[None, outside])))
    # This Sylvester separation concerns only the sampled eight-mode subspace.
    # It is not a certified lower bound against every uncomputed eigenmode.
    q = qr(np.column_stack((r, rv[:, outside])))
    compressed = mm(q.conj().T, l @ q)
    c = compressed[2:, 2:]
    sylvester = np.kron(np.eye(c.shape[0]), h) - np.kron(c.T, np.eye(2))
    dr, dl = residuals(l, vals, rv), residuals(l.T, leftvals, lv)
    pairvals = vals[selected]
    state.update(candidate_values=vals, candidate_right=rv,
                 candidate_left_values=leftvals, candidate_left=lv,
                 selected=selected)
    return dict(candidate_values=vals, candidate_transpose_left_values=leftvals,
                right_residuals=dr, transpose_left_residuals=dl,
                transpose_left_matching=left_order,
                left_right_eigenvalue_mismatches=abs(vals-leftvals),
                selected_by_projected_pair=selected,
                selected_values=pairvals,
                selected_gap=abs(pairvals[1]-pairvals[0]),
                selected_discriminant=(pairvals[1]-pairvals[0])**2/4,
                projected_vs_candidate_discriminant_error=abs(
                  disc(h)-(pairvals[1]-pairvals[0])**2/4),
                nearest_sampled_external_mode_distance=separation,
                external_to_pair_gap_ratio=separation/max(abs(pairvals[1]-pairvals[0]),1e-300),
                sampled_sylvester_separation=float(la.svdvals(sylvester)[-1]),
                sampled_complement_invariance_residual=norm(
                  l @ q - mm(q, compressed))/norm(q),
                scope='Eight eigenpairs of the full operator near the specified shift; '
                  'spectral and Sylvester separations refer to sampled neighbors only.')


def root_checks(obs, g):
  spectrum = obs['spectrum']
  selected = spectrum['selected_by_projected_pair']
  selected_residuals = [spectrum[key][j] for key in
                        ('right_residuals', 'transpose_left_residuals') for j in selected]
  return dict(subspace_converged=obs['inverse_iteration_converged'],
              small_discriminant=obs['discriminant_over_g2'] < LIMITS['discriminant_over_g2'],
              small_newton_step=abs(obs['proposed_newton_step']) < LIMITS['newton_step'],
              right_left_subspace_residual=max(obs['invariant_right_residual'],
                obs['invariant_transpose_left_residual']) < LIMITS['subspace_residual'],
              full_operator_modes=max(selected_residuals) < LIMITS['mode_residual'],
              biorthogonal=obs['biorthogonality_error'] < 1e-9,
              nonzero_nilpotent=obs['nilpotent_norm'] > LIMITS['nonzero_nilpotent'],
              nilpotent_square=obs['nilpotent_squared_ratio'] < LIMITS['nilpotent_squared_ratio'],
              nilpotent_rank_one=float(obs['nilpotent_singular_values'][1] /
                obs['nilpotent_singular_values'][0]) < 1e-8,
              simple_discriminant_zero=abs(obs['discriminant_derivative']) > 0.1*g,
              isolated_sampled_pair=spectrum['nearest_sampled_external_mode_distance'] >
                LIMITS['external_separation'],
              sampled_sylvester_separated=spectrum['sampled_sylvester_separation'] > 1e-4,
              independent_discriminant_agreement=
                spectrum['projected_vs_candidate_discriminant_error']/g**2 < 1e-7)


def solve_root(n, seed=17, shift_offset=0.04+0.013j, perturb=0j, tag=''):
  start = time.monotonic()
  name = f'ep_N{n}{tag}'
  out = dict(N=n, seed=seed, shift_offset=shift_offset, initial_perturbation=perturb,
             provenance=metadata(), iterations=[], checks_passed=False)
  try:
    engine = LocalPair(n, seed, shift_offset)
    out['initial'] = engine.initial
    z = engine.z0 + perturb*engine.y0
    r, v, centre = engine.r, engine.v, engine.centre
    for it in range(12):
      obs, state = engine.evaluate(z, r, v, centre)
      out['iterations'].append(obs)
      print(f'N={n}{tag} Newton {it}: z={z:.14g}, |D|/g^2='
            f'{obs["discriminant_over_g2"]:.3g}, |dz|='
            f'{abs(obs["proposed_newton_step"]):.3g}', flush=True)
      if not obs['inverse_iteration_converged']:
        raise RuntimeError('Block inverse iteration failed to meet its residual threshold')
      if (obs['discriminant_over_g2'] < LIMITS['discriminant_over_g2'] and
          abs(obs['proposed_newton_step']) < LIMITS['newton_step']):
        break
      step = obs['proposed_newton_step']
      cap = 0.75*engine.y0
      if abs(step) > cap:
        step *= cap/abs(step)
      obs['accepted_step'] = step
      z += step
      r, v, centre = state['r'], state['v'], obs['mean_eigenvalue']
    # Reevaluate, then query independent ARPACK modes of the full operator.
    final, state = engine.evaluate(z, state['r'], state['v'],
                                   obs['mean_eigenvalue'], candidates=True)
    checks = root_checks(final, engine.g)
    out.update(z_ep=z, eigenvalue_ep=final['mean_eigenvalue'], final=final,
               convergence_checks=checks, checks_passed=all(checks.values()),
               linear_prediction_error=z-engine.z0,
               error_scale_note='Near EP, a small backward error yields a square-root '
                 'eigenvalue splitting. The raw gap is reported, not used alone to '
                 'declare defectivity. H, nonzero K, K^2, and subspace residuals are retained.')
    np.savez_compressed(OUT/f'{name}_vectors.npz',
                        z=z, R=state['r'], V=state['v'], W=state['w'],
                        H=state['h'], B=np.array(final['B']),
                        candidate_values=state['candidate_values'],
                        candidate_right=state['candidate_right'],
                        candidate_left_values=state['candidate_left_values'],
                        candidate_left=state['candidate_left'],
                        selected=state['selected'])
    print(f'N={n}{tag} checks_passed={out["checks_passed"]} checks={checks}', flush=True)
    out['seconds'] = time.monotonic()-start
    save(f'{name}.json', out)
    return out, engine, state
  except Exception as exc:
    out.update(error=repr(exc), traceback=traceback.format_exc(),
               seconds=time.monotonic()-start)
    save(f'{name}.json', out)
    print(out['traceback'], flush=True)
    return out, None, None


def ordered_assignment(values, reference):
  distances = abs(np.array(reference)[:, None] - values[None, :])
  costs = sorted((float(distances[0, a]+distances[1, b]), (a, b))
                 for a, b in itertools.permutations(range(len(values)), 2))
  return np.array(costs[0][1]), dict(reference=reference,
    distance_matrix=distances, best_cost=costs[0][0], second_best_cost=costs[1][0],
    second_best_pair=costs[1][1], assignment_cost_margin=costs[1][0]-costs[0][0])


def continue_loop(root, engine, state, steps, tag):
  zep = root['z_ep']
  radius = 0.30*abs(zep.imag)
  out = dict(N=engine.n, seed=engine.seed, shift_offset=engine.shift_offset,
             centre=zep, radius=radius, steps_per_turn=steps, turns=2,
             conjugate_EP_distance=2*abs(zep.imag),
             radius_over_conjugate_EP_distance=0.15,
             radius_rationale='Radius 0.30 Im(z_EP), much smaller than the distance '
               '2 Im(z_EP) to its conjugate. Local neighbor separation, assignment '
               'margins, and the winding of the updated discriminant are checked '
               'at every sample. This is not a global exclusion theorem for other EPs.',
             provenance=metadata(), points=[], checks_passed=False)
  name = f'ep_N{engine.n}_loop_{tag}.json'
  r, v, centre = state['r'], state['v'], root['eigenvalue_ep']
  previous_values, previous_vectors, previous2, previous_z = None, None, None, None
  try:
    for j in range(2*steps+1):
      theta = 2*np.pi*j/steps
      z = zep+radius*np.exp(1j*theta)
      obs, state = engine.evaluate(z, r, v, centre, candidates=True)
      values = state['candidate_values']
      if previous_values is None:
        ref = la.eigvals(state['h'])
        ref = ref[np.lexsort((ref.imag, ref.real))]
        predictor = 'Projected eigenvalues at first point, ordered by real part'
      elif previous2 is None:
        ref = previous_values
        predictor = 'Previous full-operator eigenvalues'
      else:
        ref = previous_values + (previous_values-previous2)*(z-previous_z)/(previous_z-previous2_z)
        predictor = 'Complex secant predictor from previous two full-operator eigenvalues'
      selected, matching = ordered_assignment(values, ref)
      selected_values = values[selected]
      matching.update(method=predictor, selected_candidate_indices=selected,
                      agrees_with_projected_pair=set(selected)==set(state['selected']))
      if previous_vectors is not None:
        matching['right_overlap_matrix'] = abs(mm(previous_vectors.conj().T,
                                                   state['candidate_right']))
        matching['max_step_over_previous_pair_gap'] = float(
          max(abs(selected_values-previous_values))/abs(previous_values[1]-previous_values[0]))
      obs.update(theta=theta, full_operator_branches=selected_values, matching=matching)
      out['points'].append(obs)
      previous2, previous_values = previous_values, selected_values
      previous2_z, previous_z = previous_z, z
      previous_vectors = state['candidate_right'][:, selected]
      r, v, centre = state['r'], state['v'], obs['mean_eigenvalue']
      if j % steps == 0:
        print(f'N={engine.n} {tag}: {j/steps:g} turns, values={selected_values}', flush=True)
    branches = np.array([p['full_operator_branches'] for p in out['points']])
    discriminants = np.array([p['discriminant'] for p in out['points']])
    phases = np.unwrap(np.angle(discriminants))
    initial_gap = abs(branches[0, 1]-branches[0, 0])
    swap_error = float(max(abs(branches[steps]-branches[0, ::-1])))
    reset_error = float(max(abs(branches[-1]-branches[0])))
    unswapped_error = float(max(abs(branches[steps]-branches[0])))
    min_sep = min(p['spectrum']['nearest_sampled_external_mode_distance'] for p in out['points'])
    max_resid = max(max(p['spectrum'][key]) for p in out['points']
                    for key in ('right_residuals', 'transpose_left_residuals'))
    max_step = max(p['matching'].get('max_step_over_previous_pair_gap', 0)
                   for p in out['points'])
    margins = [p['matching']['assignment_cost_margin'] for p in out['points']]
    winding1 = float((phases[steps]-phases[0])/(2*np.pi))
    winding2 = float((phases[-1]-phases[0])/(2*np.pi))
    checks = dict(one_turn_swap=swap_error < LIMITS['closure_error'],
                  one_turn_is_not_reset=unswapped_error > 0.9*initial_gap,
                  two_turn_reset=reset_error < LIMITS['closure_error'],
                  discriminant_winding=abs(winding1-1)<1e-7 and abs(winding2-2)<1e-7,
                  resolved_discriminant_phase=float(max(abs(np.diff(phases))))<np.pi/2,
                  full_operator_residuals=max_resid < LIMITS['mode_residual'],
                  subspaces_converged=all(p['inverse_iteration_converged'] for p in out['points']),
                  isolated_sampled_pair=min_sep > LIMITS['external_separation'],
                  continuous_small_steps=max_step < LIMITS['circle_step_over_gap'],
                  unambiguous_assignment=min(margins) > 1e-6,
                  subspace_pair_agreement=all(p['matching']['agrees_with_projected_pair']
                                               for p in out['points']))
    out.update(summary=dict(one_turn_swap_error=swap_error, two_turn_reset_error=reset_error,
                            one_turn_identity_error=unswapped_error, initial_gap=initial_gap,
                            one_turn_discriminant_winding=winding1,
                            two_turn_discriminant_winding=winding2,
                            min_external_separation=min_sep,
                            max_mode_residual=max_resid, max_step_over_gap=max_step,
                            min_assignment_cost_margin=min(margins)),
               checks=checks, checks_passed=all(checks.values()))
  except Exception as exc:
    out.update(error=repr(exc), traceback=traceback.format_exc())
    print(out['traceback'], flush=True)
  save(name, out)
  print(f'Loop {tag} checks_passed={out["checks_passed"]}', flush=True)
  return out


def derivative_check(root, engine, state):
  z = root['z_ep']
  dz = 0.002*abs(z.imag)*(1+0.5j)
  values = []
  for sign in (-1, 1):
    obs, _ = engine.evaluate(z+sign*dz, state['r'], state['v'], root['eigenvalue_ep'])
    values.append(obs)
  finite = (values[1]['discriminant']-values[0]['discriminant'])/(2*dz)
  analytic = root['final']['discriminant_derivative']
  error = abs(finite-analytic)/abs(analytic)
  return dict(complex_step=dz, finite_difference=finite, analytic=analytic,
              relative_error=error, samples=values, checks_passed=error < 1e-5)


def main():
  p = argparse.ArgumentParser(description=__doc__)
  p.add_argument('--sizes', type=int, nargs='+', choices=(24, 32, 40), default=[24, 32, 40])
  p.add_argument('--loop-size', type=int, choices=(24, 32, 40), default=32)
  p.add_argument('--loop-steps', type=int, default=24)
  p.add_argument('--stability', action='store_true')
  p.add_argument('--skip-loop', action='store_true')
  args = p.parse_args()
  OUT.mkdir(parents=True, exist_ok=True)
  started = time.monotonic()
  summary = dict(provenance=metadata(), arguments=vars(args), roots=[], loops=[])
  loop_item = None
  with warnings.catch_warnings(record=True) as captured:
    warnings.simplefilter('always')
    for n in args.sizes:
      root, engine, state = solve_root(n)
      summary['roots'].append({k:root.get(k) for k in
        ('N', 'z_ep', 'eigenvalue_ep', 'checks_passed', 'error', 'seconds')})
      if n == args.loop_size and root['checks_passed']:
        loop_item = root, engine, state
    if loop_item and not args.skip_loop:
      base_loop = continue_loop(*loop_item, args.loop_steps, f's{17}_steps{args.loop_steps}')
      summary['loops'].append(dict(tag='baseline', checks_passed=base_loop['checks_passed'],
                                   summary=base_loop.get('summary')))
    if args.stability and loop_item:
      base_root, base_engine, base_state = loop_item
      other_root, other_engine, other_state = solve_root(
        args.loop_size, seed=83, shift_offset=-0.035+0.021j,
        perturb=0.13+0.09j, tag='_stability')
      stab = dict(seed=83, shift_offset=-0.035+0.021j,
                  initial_perturbation_in_units_y0=0.13+0.09j,
                  root_checks_passed=other_root['checks_passed'])
      if other_root['checks_passed']:
        difference = abs(other_root['z_ep']-base_root['z_ep'])
        stab.update(root_difference=difference,
                    root_repeat_agrees=difference < LIMITS['root_repeat_difference'],
                    derivative_check=derivative_check(base_root, base_engine, base_state))
        if not args.skip_loop:
          other_loop = continue_loop(other_root, other_engine, other_state,
                                    2*args.loop_steps, f's83_steps{2*args.loop_steps}')
          a = np.array([p['full_operator_branches'] for p in base_loop['points']])
          b = np.array([p['full_operator_branches'] for p in other_loop['points'][::2]])
          common_error = min(float(np.max(abs(a-b))), float(np.max(abs(a-b[:, ::-1]))))
          stab.update(common_node_branch_difference=common_error,
                      common_nodes_agree=common_error < 2e-8)
          summary['loops'].append(dict(tag='stability', checks_passed=other_loop['checks_passed'],
                                       summary=other_loop.get('summary')))
      stab['checks_passed'] = (stab.get('root_checks_passed', False)
        and stab.get('root_repeat_agrees', False)
        and stab.get('derivative_check', {}).get('checks_passed', False)
        and (args.skip_loop or stab.get('common_nodes_agree', False)))
      summary['stability'] = stab
    summary['warnings'] = sorted(set(str(w.message) for w in captured))
  summary['seconds'] = time.monotonic()-started
  summary['checks_passed'] = (all(r['checks_passed'] for r in summary['roots'])
    and (args.skip_loop or bool(summary['loops']))
    and all(r['checks_passed'] for r in summary['loops'])
    and (not args.stability or summary.get('stability', {}).get('checks_passed', False)))
  save('ep_summary.json', summary)
  print(json.dumps(serial(summary), indent=2, allow_nan=False), flush=True)
  return 0 if summary['checks_passed'] else 1


if __name__ == '__main__':
  sys.exit(main())
