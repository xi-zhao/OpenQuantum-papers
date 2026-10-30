"""Independent checks of the coherence-response preparation and dynamics.

Uses explicit tensor-product spins for global dephasing at N=1--4, and an
adaptive ODE integrator for the N=32 signal previously computed by expm_multiply.
Reads existing samples without replacing them. Writes a new validation record.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import expm_multiply
from biorthogonal_generator import generators
from biorthogonal_validate import kron_all
from biorthogonal_audit import mm
from pra_modes import product_state, fixed, DF, NF

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/coherence_response'


def full_spin_control():
  qfold = fixed(DF, NF)
  n, u, v = qfold
  rho1 = np.array([[n, u+1j*v], [u-1j*v, 1-n]])
  rows = []
  for N in (1, 2, 3, 4):
    detuning = 3.2
    eye = np.eye(2**N)
    def local(op, site):
      return kron_all([op if k == site else np.eye(2) for k in range(N)])
    Ne = sum(local(np.diag([1., 0.]), k) for k in range(N))
    H = -detuning*Ne+sum(local(np.array([[0., 1.], [1., 0.]]), k) for k in range(N))
    H += 5/N*(mm(Ne, Ne)-Ne)
    fullL = -1j*(np.kron(eye, H)-np.kron(H.T, eye))
    for k in range(N):
      jump = local(np.array([[0., 0.], [1., 0.]]), k)
      jj = mm(jump.conj().T, jump)
      fullL += np.kron(jump.conj(), jump)-.5*(np.kron(eye, jj)+np.kron(jj.T, eye))
    rho = kron_all([rho1]*N)
    numbers = np.diag(Ne)
    dephased = rho*(numbers[:, None] == numbers[None, :])
    minimum_eigenvalue = float(np.linalg.eigvalsh(dephased).min())
    population_error = float(np.max(abs(np.diag(rho)-np.diag(dephased))))
    if minimum_eigenvalue < -1e-13 or population_error > 1e-13:
      raise RuntimeError('Global-dephasing preparation is not physical')
    delta = (rho-dephased).ravel(order='F')
    times = np.array([0., .1, .5, 1.1, 2.])
    ode = solve_ivp(lambda t, y: np.einsum('ij,j->i', fullL, y, optimize=False),
      (0., times[-1]), delta, method='DOP853', t_eval=times, rtol=3e-13, atol=1e-15)
    if not ode.success:
      raise RuntimeError(ode.message)
    observable = Ne.T.ravel(order='F')/N
    reference = np.dot(observable, ode.y)
    L0, L1, s = generators(N)
    piL = L0+detuning*L1
    folded = product_state(qfold, s)
    initial = np.where(s['trace'].astype(bool), 0., folded)
    piobs = s['trace']*s['ne'][s['row']]/N
    signal = [np.dot(piobs, expm_multiply(t*piL, initial, traceA=t*piL.diagonal().sum()))
      for t in times]
    error = float(np.max(abs(np.asarray(signal)-reference)))
    slope = np.dot(observable, np.einsum('ij,j->i', fullL, delta, optimize=False))
    if error > 2e-11 or abs(slope-NF) > 1e-12:
      raise RuntimeError(f'Explicit-spin global-dephasing check failed: {N}, {error}, {slope}')
    rows.append(dict(N=N, detuning=detuning, times=times.tolist(),
      maximum_density_error=error, initial_slope=float(slope.real),
      initial_population_error=population_error,
      dephased_minimum_eigenvalue=minimum_eigenvalue))
  return rows


def ode_control(record):
  N, detuning = record['N'], record['detuning']
  L0, L1, s = generators(N)
  L = L0+detuning*L1
  folded = product_state(fixed(DF, NF), s)
  initial = np.where(s['trace'].astype(bool), 0., folded)
  times = np.array(record['dynamics']['times'])
  ode = solve_ivp(lambda t, y: L@y, (0., times[-1]), initial, method='DOP853',
    t_eval=times, rtol=1e-10, atol=1e-13)
  if not ode.success:
    raise RuntimeError(ode.message)
  observable = s['trace']*s['ne'][s['row']]/N
  signal = np.dot(observable, ode.y)
  reference = np.array(record['dynamics']['full'])
  maximum_error = float(np.max(abs(signal-reference)))
  relative_error = float(np.linalg.norm(signal-reference)/np.linalg.norm(reference))
  if maximum_error > 1e-8 or relative_error > 1e-7:
    raise RuntimeError(f'Independent time integrator failed: {maximum_error}, {relative_error}')
  return dict(N=N, evaluations=ode.nfev, rtol=1e-10, atol=1e-13,
    maximum_density_error=maximum_error, relative_l2_error=relative_error)


def main(args):
  source = DATA/args.input
  output = DATA/args.output
  if output.exists():
    raise FileExistsError(f'Preserve previous evidence; choose a new --output: {output}')
  data = json.loads(source.read_text())
  mismatches = [name for name, sha in data['source_sha256'].items()
    if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != sha]
  if mismatches:
    raise RuntimeError(f'Source snapshot changed: {mismatches}')
  full_spin = full_spin_control()
  large = next(r for r in data['records'] if 'dynamics' in r)
  result = dict(created='2026-10-30', input=str(source.relative_to(ROOT)),
    input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    validation_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    frozen_source_hashes_match=True, explicit_spin_dephasing=full_spin,
    independent_time_integrator=ode_control(large))
  output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
  print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--input', default='response_2026-10-30.json')
  parser.add_argument('--output', default='validation_2026-10-30.json')
  main(parser.parse_args())
