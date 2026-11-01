"""Independent controls for the fixed-preparation physical-effect audit.

Explicit full spins, alternative time integration, phase quadrature refinement,
and archived-response reproduction. Writes new evidence, not original samples.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import expm_multiply
from biorthogonal_audit import mm
from biorthogonal_generator import generators
from biorthogonal_validate import kron_all
from physical_effects import (ROOT, DATA, TIMES, preparation_parameters,
  preparations, classical_signals, raw_dynamics)
from pra_modes import NF, V


def full_spin(N, detuning, interaction, parameters):
  one = np.eye(2)
  n = np.diag([1., 0.])
  sx = np.array([[0., 1.], [1., 0.]])
  sm = np.array([[0., 0.], [1., 0.]])
  def local(op, i):
    return kron_all([op if j == i else one for j in range(N)])
  Ne = sum(local(n, i) for i in range(N))
  eye = np.eye(2**N)
  H = -detuning*Ne+sum(local(sx, i) for i in range(N))
  H += interaction/(2*N)*(mm(Ne, Ne)-Ne)
  L = -1j*(np.kron(eye, H)-np.kron(H.T, eye))
  for i in range(N):
    jump = local(sm, i)
    square = mm(jump.conj().T, jump)
    L += np.kron(jump.conj(), jump)-.5*(np.kron(eye, square)+np.kron(square.T, eye))
  def product(state):
    pop, u, v = state
    return kron_all([np.array([[pop, u+1j*v], [u-1j*v, 1-pop]])]*N)
  folded = product(parameters['fold'])
  numbers = np.diag(Ne)
  dephased = folded*(numbers[:, None] == numbers[None, :])
  plus, minus = product(parameters['high_plus']), product(parameters['high_minus'])
  sources = np.column_stack([(folded-dephased).ravel(order='F'), (plus-minus).ravel(order='F')])
  distances = [float(np.sum(abs(np.linalg.eigvalsh(x)))/2) for x in (folded-dephased, plus-minus)]
  selected_times = np.array([0., .1, .5, 1.1, 2., 6., 12.])
  def rhs(t, y):
    return mm(L, y.reshape(-1, 2)).ravel()
  ode = solve_ivp(rhs, (0., 12.), sources.ravel(), method='DOP853', t_eval=selected_times,
    rtol=3e-12, atol=3e-14)
  if not ode.success:
    raise RuntimeError(ode.message)
  obs = Ne.T.ravel(order='F')/N
  full = np.einsum('i,ijt->tj', obs, ode.y.reshape(len(L), 2, -1), optimize=False)
  L0, L1, space = generators(N, interaction=interaction)
  piobs = space['trace']*space['ne'][space['row']]/N
  piinput, pidistances = preparations(parameters, space)
  piL = L0+detuning*L1
  projected = np.array([np.dot(piobs, expm_multiply(t*piL, piinput,
    traceA=t*piL.diagonal().sum())) for t in selected_times])
  error = float(np.max(abs(projected-full)))
  distance_error = float(np.max(abs(np.array(distances)-pidistances)))
  if error > 1e-9 or distance_error > 1e-11:
    raise RuntimeError(f'Full-spin control failed: {N}, {interaction}, {error}, {distance_error}')
  return dict(N=N, detuning=detuning, interaction=interaction,
    maximum_density_error=error, maximum_trace_distance_error=distance_error,
    explicit_trace_distances=distances, evaluations=ode.nfev)


def ode_check(scan):
  central = next(r for r in scan['records'] if r['offset'] == 0.)
  N, detuning = scan['N'], central['detuning']
  L0, L1, space = generators(N)
  L = L0+detuning*L1
  source, _ = preparations(scan['preparations'], space)
  def rhs(t, y):
    return (L@y.reshape(-1, 2)).ravel()
  ode = solve_ivp(rhs, (0., 12.), source.ravel(), method='DOP853', t_eval=TIMES,
    rtol=1e-10, atol=1e-13)
  if not ode.success:
    raise RuntimeError(ode.message)
  obs = space['trace']*space['ne'][space['row']]/N
  signal = np.einsum('i,ijt->tj', obs, ode.y.reshape(len(obs), 2, -1), optimize=False)
  reference = np.asarray(central['full'])
  error = float(np.max(abs(signal-reference)))
  relative = np.linalg.norm(signal-reference, axis=0)/np.linalg.norm(reference, axis=0)
  if error > 1e-8 or max(relative) > 1e-7:
    raise RuntimeError(f'Independent ODE failed: {error}, {relative}')
  return dict(N=N, evaluations=ode.nfev, maximum_density_error=error,
    relative_l2_errors=relative.tolist())


def phase_check(scan):
  checks = []
  for record in scan['records']:
    refined, _ = classical_signals(scan['preparations'], record['detuning'],
      V*(scan['N']-1)/scan['N'], phases=128)
    error = float(np.max(abs(refined-np.asarray(record['meanfield_phase_ensemble']))))
    check = dict(offset=record['offset'], max_density_difference=error,
      coarse_passes_2e_minus7=bool(error <= 2e-7))
    if error > 2e-7:
      finer, _ = classical_signals(scan['preparations'], record['detuning'],
        V*(scan['N']-1)/scan['N'], phases=256)
      fine_error = float(np.max(abs(finer-refined)))
      if fine_error > 2e-9:
        raise RuntimeError(f'Refined phase quadrature failed: {fine_error}')
      check.update(accepted_phases=256, refinement_128_to_256_error=fine_error,
        accepted_meanfield_phase_ensemble=finer.tolist(),
        interpretation='The archived 64-phase record is retained; use this refined control at this point.')
    else:
      check['accepted_phases'] = 128
    checks.append(check)
  return dict(N=scan['N'], coarse_phases=64, refined_phases=128, records=checks)


def baseline_check(parameters, detuning):
  L0, L1, space = generators(4, interaction=0.)
  source, _ = preparations(parameters, space)
  obs = space['trace']*space['ne'][space['row']]/4
  pi, _ = raw_dynamics(L0+detuning*L1, source, obs, space)
  single, _ = classical_signals(parameters, detuning, 0., phases=4)
  error = float(np.max(abs(pi-single)))
  if error > 2e-10:
    raise RuntimeError(f'Noninteracting baseline failed: {error}')
  return dict(N=4, interaction=0., maximum_density_error=error,
    explanation='Global dephasing retains correlations, but the independent density observable is one-body.')


def main(args):
  path, output = DATA/args.input, DATA/args.output
  if output.exists():
    raise FileExistsError(f'Choose a new output: {output}')
  data = json.loads(path.read_text())
  mismatches = [p for p, sha in data['source_sha256'].items()
    if hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != sha]
  if mismatches:
    raise RuntimeError(f'Input sources changed: {mismatches}')
  scans = data['scans']
  control = data['fixed_detuning_control']
  detuning = control['detuning']
  parameters = control['preparations']
  spin = [full_spin(N, detuning, interaction, parameters) for N, interaction in ((1, 10.), (4, 10.), (4, 0.))]
  baseline = baseline_check(parameters, detuning)
  phases = [phase_check(s) for s in scans]
  center32 = next(s for s in scans if s['N'] == 32)
  ode = ode_check(center32)
  old_path = ROOT/'data/coherence_response/response_2026-10-30.json'
  old = json.loads(old_path.read_text())
  old32 = next(r for r in old['records'] if r['N'] == 32)
  center = next(r for r in center32['records'] if r['offset'] == 0.)
  old_signal = np.asarray(old32['dynamics']['full'])
  signal_error = float(np.max(abs(np.asarray(center['full'])[:, 0]-old_signal)))
  old_residues = np.asarray(old32['reference_residues'])[:, 2]
  residue_error = float(np.max(abs(np.asarray(center['modes']['residues'])[:, 0]-old_residues)))
  if signal_error > 1e-10 or residue_error > 1e-9:
    raise RuntimeError('Previous response evidence did not reproduce')
  slopes = [abs(s['initial_slope'][0]-NF) for s in scans]
  if max(slopes) > 1e-10:
    raise RuntimeError('Exact coherence-population initial slope failed')
  result = dict(created='2026-11-01', input=str(path.relative_to(ROOT)),
    input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    validation_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    frozen_sources_match=True, explicit_spin=spin, independent_baseline=baseline,
    phase_quadrature=phases, independent_ode=ode,
    previous_response=dict(path=str(old_path.relative_to(ROOT)),
      sha256=hashlib.sha256(old_path.read_bytes()).hexdigest(),
      maximum_signal_error=signal_error, maximum_residue_error=residue_error),
    exact_initial_slope_maximum_error=max(slopes))
  output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
  print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--input', default='effects_2026-11-01.json')
  parser.add_argument('--output', default='validation_2026-11-01.json')
  main(parser.parse_args())
