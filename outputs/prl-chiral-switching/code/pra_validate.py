"""Validate PRA preparations, stored observables, and repeatability.

The small-N reference uses explicit tensor-product spins, independent of the
Dicke product-state formula. Stored samples are read without being rewritten.
Only --repeat-modes recomputes the N=32 step/seed check; its new record is kept
separately from the original mode_stability.json exploration.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import platform
import numpy as np
import scipy
import qutip
from scipy.linalg import eig
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import expm_multiply
from biorthogonal_generator import generators
from biorthogonal_audit import project, mm
from biorthogonal_validate import kron_all
from pra_modes import product_state

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data/pra_spectral'


def load(path):
  return json.loads((ROOT/path).read_text())


def full_spin_preparations():
  state = np.array([.29, .13, -.18])
  n, u, v = state
  rho1 = np.array([[n, u+1j*v], [u-1j*v, 1-n]])
  records = []
  for N in (1, 2, 3, 4):
    one = np.eye(2)
    sx = np.array([[0, 1], [1, 0]])
    ne1 = np.diag([1., 0.])
    lower = np.array([[0., 0.], [1., 0.]])
    def local(op, site):
      return kron_all([op if i == site else one for i in range(N)])
    Ne = sum(local(ne1, i) for i in range(N))
    H = -3.25*Ne+sum(local(sx, i) for i in range(N))+5/N*(mm(Ne, Ne)-Ne)
    eye = np.eye(2**N)
    L = -1j*(np.kron(eye, H)-np.kron(H.T, eye))
    for site in range(N):
      jump = local(lower, site)
      jj = mm(jump.conj().T, jump)
      L += np.kron(jump.conj(), jump)-.5*(np.kron(eye, jj)+np.kron(jj.T, eye))
    rho = kron_all([rho1]*N).ravel(order='F')
    L0, L1, s = generators(N)
    pi = (L0+3.25*L1).tocsc()
    initial = product_state(state, s)
    obs = s['trace']*s['ne'][s['row']]/N
    times = (0., .25, .75, 1.5)
    reference = solve_ivp(lambda t, y: np.einsum('ij,j->i', L, y, optimize=False),
      (0., times[-1]), rho, method='DOP853', t_eval=times, rtol=2.3e-14, atol=1e-15)
    assert reference.success
    samples = []
    for index, time in enumerate(times):
      full = reference.y[:, index].reshape((2**N, 2**N), order='F')
      stored = expm_multiply(time*pi, initial, traceA=time*pi.diagonal().sum())
      matrix = np.zeros((s['d'], s['d']), complex)
      matrix[s['row'], s['col']] = stored
      purity = 0.
      min_eigenvalue = 1.
      for tj in range(N, -1, -2):
        at = np.flatnonzero(s['twoj'] == tj)
        block = matrix[np.ix_(at, at)]
        k = (N-tj)//2
        multiplicity = math.comb(N, k)-(math.comb(N, k-1) if k else 0)
        purity += np.trace(mm(block, block)).real/multiplicity
        min_eigenvalue = min(min_eigenvalue, float(np.linalg.eigvalsh(block).min()/multiplicity))
      density_error = float(abs(np.trace(mm(Ne, full))/N-np.dot(obs, stored)))
      purity_error = float(abs(np.trace(mm(full, full)).real-purity))
      trace_error = float(abs(np.dot(s['trace'], stored)-1))
      assert max(density_error, purity_error, trace_error) < 2e-12
      assert min_eigenvalue > -2e-12
      samples.append(dict(time=time, density_error=density_error, purity_error=purity_error,
        trace_error=trace_error, minimum_physical_block_eigenvalue=min_eigenvalue))
    records.append(dict(N=N, samples=samples))
  return records


def stored_checks():
  records = []
  for N, suffix in ((32, ''), (40, '_quick')):
    source = load(f'data/pra_spectral/modes_N{N}{suffix}.json')
    rows = source['records']
    assert len(rows) == (15 if N == 32 else 3)
    for row in rows:
      lam, a = np.array(row['eigenvalues']), np.array(row['residues'])
      qhi = np.array(row['high_meanfield_state'])
      direction = np.array(row['high_perturbation_direction'])
      step = row['high_perturbation_step']
      radii = [float(np.linalg.norm([2*q[1], -2*q[2], 2*q[0]-1]))
               for q in (qhi+step*direction, qhi-step*direction, source['fold']['state'])]
      assert max(radii) < 1
      expected_initial = [2*step*direction[0], source['fold']['state'][0]-row['stationary_density']]
      assert np.max(abs(np.array(expected_initial)-row['initial_signal'])) < 1e-12
      assert lam[0] < lam[1] < 0
      assert row['steady_state_residual'] < 1e-9
      fractions = []
      for profile in row['eigenmode_population_profiles']:
        p = abs(np.array(profile))
        fractions.append(float(p[np.arange(N+1)/N < .22].sum()/p.sum()))
      stored_fraction = [p['low_region_absolute_fraction'] for p in row['eigenmode_population_location']]
      assert np.max(abs(np.array(fractions)-stored_fraction)) < 2e-14
      check = dict(N=N, offset=row['offset'], preparation_bloch_radii=radii,
        invariant_right_error=row['diagnostics']['invariant_right_error'],
        invariant_left_error=row['diagnostics']['invariant_left_error'])
      if 'dynamics' in row:
        d = row['dynamics']
        t = np.array(d['times'])
        exact = np.array(d['exact'])
        pair = np.einsum('ti,ij->tj', np.exp(t[:, None]*lam), a, optimize=False)
        ol = np.array(d['oscillatory_eigenvalues'])
        ol = ol[:, 0]+1j*ol[:, 1]
        oa = np.array(d['oscillatory_residues'])
        oa = oa[:, :, 0]+1j*oa[:, :, 1]
        osc = np.einsum('ti,ij->tj', np.exp(t[:, None]*ol), oa, optimize=False)
        assert np.max(abs(osc.imag)) < 1e-9
        four = pair+osc.real
        assert np.max(abs(pair-np.array(d['two_modes']))) < 1e-12
        assert np.max(abs(four-np.array(d['four_modes']))) < 1e-12
        assert np.max(abs(exact[0]-row['initial_signal'])) < 1e-12
        windows = {}
        for lower in (0., 2., 4., 6.):
          mask = t >= lower
          errors = {}
          for name, signal in [('relative_l2_error', pair), ('four_mode_relative_l2_error', four)]:
            error = np.linalg.norm(exact[mask]-signal[mask], axis=0)/np.linalg.norm(exact[mask], axis=0)
            assert np.max(abs(error-d['windows'][str(lower)][name])) < 1e-12
            errors[name] = error.tolist()
          windows[str(lower)] = errors
        check['verified_windows'] = windows
        assert d['trace_error'] < 2e-12
      records.append(check)
  ep = load('data/pra_spectral/ep_summary.json')
  assert ep['checks_passed'] and all(r['checks_passed'] for r in ep['roots']+ep['loops'])
  assert ep['stability']['checks_passed']
  return records


def repeat_modes():
  original = load('data/pra_spectral/modes_N32.json')
  old = load('data/pra_spectral/mode_stability.json')
  L0, L1, s = generators(32)
  obs = s['trace']*s['ne'][s['row']]/32
  records = []
  for offset in (-.06, 0., .06):
    row = min(original['records'], key=lambda r: abs(r['offset']-offset))
    d = row['detuning']
    audit, Rd, Wd, _ = project(L0, L1, s, d, 1122, -.8)
    lam, transform = eig(np.array(audit['H']))
    order = np.argsort(lam.real)
    lam, transform = lam[order], transform[:, order]
    R = mm(Rd, transform)
    W = mm(Wd, np.linalg.inv(transform).T)
    values = []
    for step in (.005, .0025):
      state = np.array(row['high_meanfield_state'])
      direction = np.array(row['high_perturbation_direction'])
      delta = product_state(state+step*direction, s)-product_state(state-step*direction, s)
      residue = np.dot(obs, R)*np.dot(W.T, delta)/(2*step)
      assert np.max(abs(residue.imag)) < 1e-10
      values.append(residue.real)
    relative = float(np.linalg.norm(values[1]-values[0])/np.linalg.norm(values[0]))
    eig_error = float(np.max(abs(lam-row['eigenvalues'])))
    previous = min(old, key=lambda r: abs(r['offset']-offset))
    assert relative < 8e-6 and eig_error < 2e-10
    assert abs(relative-previous['relative_step_change']) < 2e-9
    assert np.max(abs(np.array(values)-previous['residues_per_unit_step'])) < 2e-9
    records.append(dict(offset=offset, steps=[.005, .0025], residues_per_unit_step=np.array(values).tolist(),
      relative_step_change=relative, seed_shift_eigenvalue_error=eig_error))
  return records


def main(repeat=False):
  inputs = [Path(__file__), ROOT/'code/pra_modes.py', ROOT/'code/biorthogonal_generator.py',
    ROOT/'code/biorthogonal_audit.py', DATA/'modes_N32.json', DATA/'modes_N40_quick.json',
    DATA/'mode_stability.json', DATA/'ep_summary.json']
  result = dict(versions=dict(python=platform.python_version(), numpy=np.__version__,
    scipy=scipy.__version__, qutip=qutip.__version__),
    source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    full_spin_preparations=full_spin_preparations(), stored_records=stored_checks(),
    mode_repeat=repeat_modes() if repeat else None,
    scope='Numerical consistency and independent small-N reference; not an independent manuscript review.')
  result['checks_passed'] = True
  path = DATA/('validation_repeat.json' if repeat else 'validation.json')
  path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
  print(json.dumps(dict(checks_passed=True, output=str(path), mode_repeat=result['mode_repeat'])))


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--repeat-modes', action='store_true')
  main(parser.parse_args().repeat_modes)
