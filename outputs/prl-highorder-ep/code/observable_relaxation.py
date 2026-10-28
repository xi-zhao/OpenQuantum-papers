"""Conditional one-time density observables, with direct dynamics verification.

No fitted coefficients enter theoretical curves. Population evolution retains two
marked modes and the occupation count of the remaining modes. Every initial state
is a diagonal positive density matrix; no odd-parity state is prepared.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply

from spectral_analysis import boolean_generator, chain, full_liouvillian, perron_log

ROOT = Path(__file__).resolve().parents[1]


def population_chain(n, d, z):
  matrix = np.zeros((n + 1, n + 1))
  for k in range(n + 1):
    if k:
      matrix[k - 1, k] = k
    if k + d <= n:
      matrix[k + d, k] = z * math.comb(n - k, d)
  return matrix


def marked_generator(n, d, z):
  """Aggregate probabilities p[a,b,k], not per-configuration coefficients."""
  states = list(itertools.product([0, 1], [0, 1], range(n - 1)))
  index = {state: i for i, state in enumerate(states)}
  matrix = np.zeros((len(states), len(states)))
  for source, (a, b, k) in enumerate(states):
    if a:
      matrix[index[0, b, k], source] += 1
    if b:
      matrix[index[a, 0, k], source] += 1
    if k:
      matrix[index[a, b, k - 1], source] += k
    for u in range(2 - a):
      for v in range(2 - b):
        r = d - u - v
        if 0 <= r <= n - 2 - k:
          matrix[index[a + u, b + v, k + r], source] += z * math.comb(n - 2 - k, r)
  return states, matrix


def initial_probabilities(states, n, filling=.75, bias=.6):
  # rho = rho_filling [1 + bias (n_1-n_2)], normalized because Tr(rho_filling D)=0.
  result = np.array([math.comb(n - 2, k) * filling ** (a + b + k)
                     * (1 - filling) ** (n - a - b - k) * (1 + bias * (a - b))
                     for a, b, k in states])
  assert np.all(result > 0) and abs(result.sum() - 1) < 1e-12
  return result


def rates(n, d, z):
  values = sorted(la.eigvals(population_chain(n, d, z)), key=lambda x: -x.real)
  r = float(values[0].real)
  r1 = float(max(la.eigvals(population_chain(n - 1, d, z)).real))
  r2 = float(max(la.eigvals(population_chain(n - 2, d, z)).real))
  q = z ** (1 / (d + 1))
  certified = q * perron_log(n, d)["value"]
  assert abs(r - certified) / r < 1e-10
  internal = r - values[1].real
  return dict(perron=r, lambda_star=r - n, density=float(internal),
              imbalance=r - r2, odd_candidate=r - r1,
              gap_all=float(min(internal, r - r1)),
              gap_even=float(min(internal, r - r2)), next_density_imag=float(values[1].imag))


def perron_probability(matrix):
  values, vectors = la.eig(matrix)
  top = int(np.argmax(values.real))
  vector = vectors[:, top].real
  vector *= np.sign(vector.sum())
  vector /= vector.sum()
  assert vector.min() > 0
  residual = np.max(abs(np.einsum("ij,j->i", matrix, vector, optimize=False) - values[top].real * vector))
  assert residual < 1e-10
  return values[top].real, vector


def fit_tail(times, signal, start=5., end=8.):
  mask = (times >= start) & (times <= end) & (np.abs(signal) > 2e-12)
  assert mask.sum() >= 5
  coefficient = np.polyfit(times[mask], np.log(np.abs(signal[mask])), 1)
  residual = np.log(np.abs(signal[mask])) - np.polyval(coefficient, times[mask])
  return dict(rate=float(-coefficient[0]), amplitude=float(np.exp(coefficient[1])),
              time_window=[float(times[mask][0]), float(times[mask][-1])],
              max_log_residual=float(np.max(abs(residual))))


def dynamics(n=12, d=3, zeta=1., filling=.75, bias=.6):
  z = zeta * n ** (1 - d)
  assert z <= n / math.comb(n, d)
  theory = rates(n, d, z)
  states, matrix = marked_generator(n, d, z)
  initial = initial_probabilities(states, n, filling, bias)
  times = np.linspace(0, 12, 241)
  shifted = sp.csr_matrix(matrix - theory["perron"] * np.eye(len(states)))
  raw = expm_multiply(shifted, initial, start=times[0], stop=times[-1], num=len(times))
  trace = raw.sum(axis=1)
  probabilities = raw / trace[:, None]
  assert probabilities.min() > -1e-13
  occupations = np.array([a + b + k for a, b, k in states])
  imbalance_weight = np.array([a - b for a, b, k in states])
  density = np.einsum("tk,k->t", probabilities, occupations, optimize=False) / n
  imbalance = np.einsum("tk,k->t", probabilities, imbalance_weight, optimize=False)
  _, stationary = perron_probability(population_chain(n, d, z))
  filling_star = float(np.dot(np.arange(n + 1), stationary) / n)
  density_deviation = density - filling_star
  density_fit = fit_tail(times, density_deviation)
  imbalance_fit = fit_tail(times, imbalance)
  for name, fit in [("density", density_fit), ("imbalance", imbalance_fit)]:
    fit["theory"] = theory[name]
    fit["relative_error"] = abs(fit["rate"] - theory[name]) / theory[name]
    assert fit["relative_error"] < .01, (n, d, name, fit)
  # Independent unmarked propagation checks cancellation of the entire bias.
  unmarked = population_chain(n, d, z)
  total_initial = np.array([math.comb(n, k) * filling ** k * (1 - filling) ** (n - k)
                            for k in range(n + 1)])
  total_raw = expm_multiply(sp.csr_matrix(unmarked - theory["perron"] * np.eye(n + 1)),
                           total_initial, start=times[0], stop=times[-1], num=len(times))
  total_density = np.einsum("tk,k->t", total_raw, np.arange(n + 1), optimize=False) / total_raw.sum(axis=1) / n
  symmetry_error = float(np.max(abs(total_density - density)))
  assert symmetry_error < 2e-12
  traces_error = float(np.max(abs(total_raw.sum(axis=1) - trace)))
  assert traces_error < 2e-11
  imbalance_initial = np.array([2 * bias * math.comb(n - 2, k) * filling ** (k + 1)
                                * (1 - filling) ** (n - k - 1) for k in range(n - 1)])
  imbalance_raw = expm_multiply(
    sp.csr_matrix(population_chain(n - 2, d, z) - theory['perron'] * np.eye(n - 1)),
    imbalance_initial, start=times[0], stop=times[-1], num=len(times))
  imbalance_sector_error = float(np.max(abs(imbalance_raw.sum(axis=1) / trace - imbalance)))
  assert imbalance_sector_error < 2e-12
  survival = np.exp(theory["lambda_star"] * times) * trace
  assert np.all(np.diff(survival) <= 1e-12) and abs(survival[0] - 1) < 1e-12
  # A small readout imbalance does not change the generator or the initial state.
  readouts = {}
  # Derivatives come from the generator, not finite differences of small tails.
  derivative = np.einsum("ij,tj->ti", matrix, probabilities, optimize=False)
  derivative -= probabilities * derivative.sum(axis=1)[:, None]
  density_derivative = np.einsum("tk,k->t", derivative, occupations, optimize=False) / n
  imbalance_derivative = np.einsum("tk,k->t", derivative, imbalance_weight, optimize=False)
  for eta in [0., .01, .1, 1.]:
    deviation = density_deviation + eta * imbalance / n
    deriv = density_derivative + eta * imbalance_derivative / n
    effective = -deriv / deviation
    readouts[str(eta)] = dict(deviation=deviation.tolist(), effective_rate=effective.tolist())
  return dict(n=n, d=d, z=z, zeta=zeta, filling=filling, bias=bias, theory=theory,
              time=times.tolist(), density=density.tolist(), stationary_filling=filling_star,
              density_deviation=density_deviation.tolist(), imbalance=imbalance.tolist(),
              survival=survival.tolist(), readouts=readouts,
              density_tail=density_fit, imbalance_tail=imbalance_fit,
              checks=dict(initial_min_aggregate_probability=float(initial.min()),
                          min_evolved_probability=float(probabilities.min()),
                          density_symmetry_error=symmetry_error, survival_symmetry_error=traces_error,
                          imbalance_sector_error=imbalance_sector_error))


def selection_identity_checks():
  """Integer Boolean intertwiners certify selection even at defective roots."""
  rows = []
  for n in range(2, 8):
    subsets = [frozenset(s) for k in range(n + 1)
               for s in itertools.combinations(range(n), k)]
    total = np.array([len(s) for s in subsets])
    detector = np.array([int(0 in s) - int(1 in s) for s in subsets])
    symmetric = (total[:, None] == np.arange(n + 1)).astype(np.int64)
    antisymmetric = (detector[:, None] *
                     (total[:, None] == np.arange(1, n))).astype(np.int64)
    assert np.array_equal(antisymmetric.sum(axis=1), detector)
    assert np.all(antisymmetric.sum(axis=0) == 0)
    assert np.all(np.einsum('i,ij->j', total, antisymmetric) == 0)
    assert np.array_equal(np.einsum('i,ij->j', detector, antisymmetric),
                          [2 * math.comb(n - 2, k) for k in range(n - 1)])
    for d in range(1, n + 1):
      G = boolean_generator(n, d, 1).astype(np.int64)
      residuals = []
      for T, m in [(symmetric, n), (antisymmetric, n - 2)]:
        B = chain(m, d, 1).astype(np.int64)
        residual = np.einsum('ij,jk->ik', G, T) - np.einsum('ij,jk->ik', T, B)
        residuals.append(int(np.max(abs(residual))))
      assert max(residuals) == 0
      rows.append(dict(n=n, d=d, integer_residuals=residuals))
  return rows


def direct_full_check(n, d, zeta=1.):
  """Independently evolve the full CAR superoperator and a positive density matrix."""
  z = zeta * n ** (1 - d)
  energy = np.linspace(-.6, .8, n)
  full = sp.csr_matrix(full_liouvillian(n, d, z, np.diag(energy)))
  states, small = marked_generator(n, d, z)
  initial = initial_probabilities(states, n)
  basis = list(itertools.product([0, 1], repeat=n))
  probabilities = np.array([.75 ** sum(b) * .25 ** (n - sum(b)) * (1 + .6 * (b[0] - b[1]))
                            for b in basis])
  dimension = 2 ** n
  rho = np.diag(probabilities).reshape(-1).astype(complex)
  times = np.linspace(0, 4, 41)
  r = rates(n, d, z)["perron"]
  raw_full = expm_multiply(full - (r - n) * sp.identity(4 ** n), rho,
                          start=0, stop=times[-1], num=len(times))
  raw_small = expm_multiply(sp.csr_matrix(small - r * np.eye(len(states))), initial,
                           start=0, stop=times[-1], num=len(times))
  max_error = 0.
  min_eigenvalue = 1.
  max_offdiagonal = 0.
  for f, s in zip(raw_full, raw_small):
    fm = f.reshape(dimension, dimension)
    max_offdiagonal = max(max_offdiagonal, float(np.max(abs(fm - np.diag(fm.diagonal())))))
    aggregate = np.zeros(len(states))
    index = {state: i for i, state in enumerate(states)}
    for b, value in zip(basis, fm.diagonal()):
      aggregate[index[b[0], b[1], sum(b[2:])]] += value.real
    max_error = max(max_error, float(np.max(abs(aggregate - s))))
    min_eigenvalue = min(min_eigenvalue, float(np.min(fm.diagonal().real) / fm.trace().real))
  assert max_error < 1e-11 and max_offdiagonal < 1e-13 and min_eigenvalue >= -1e-13
  parity = np.array([(-1) ** sum(b) for b in basis])
  operator_parity = (parity[:, None] * parity[None, :]).reshape(-1)
  even = np.flatnonzero(operator_parity == 1)
  odd = np.flatnonzero(operator_parity == -1)
  assert full[even][:, odd].nnz == 0 and full[odd][:, even].nnz == 0
  even_values = la.eigvals(full[even][:, even].toarray())
  odd_values = la.eigvals(full[odd][:, odd].toarray())
  all_real = np.sort(np.r_[even_values.real, odd_values.real])[::-1]
  even_real = np.sort(even_values.real)[::-1]
  gap_all = float(all_real[0] - all_real[1])
  gap_even = float(even_real[0] - even_real[1])
  prediction = rates(n, d, z)
  assert abs(gap_all - prediction['gap_all']) < 1e-8
  # At n=4,d=3 the even edge includes a size-three nilpotent block.
  # Dense eig splits it by ~5e-6: preserve this diagnostic, not a false pass.
  defective_even_edge = n - 2 < d and abs(prediction['gap_even'] - r) < 1e-10
  if not defective_even_edge:
    assert abs(gap_even - prediction['gap_even']) < 1e-8
  assert np.max(abs(rho[odd])) == 0
  return dict(n=n, d=d, z=z, max_full_vs_reduced_error=max_error,
              max_offdiagonal=max_offdiagonal, min_density_eigenvalue=min_eigenvalue,
              initial_even_parity=True, h_diagonal=True,
              odd_initial_overlap=0., parity_mixing_matrix_elements=0,
              full_gap_direct=gap_all, even_gap_direct=gap_even,
              full_gap_predicted=prediction['gap_all'], even_gap_predicted=prediction['gap_even'],
              even_gap_direct_error=abs(gap_even - prediction['gap_even']),
              even_gap_direct_check=('not used: defective zero-block splitting; matrix identities checked separately'
                                     if defective_even_edge else 'passed'))


def visible_modes(n, d, zeta=1., filling=.75, bias=.6):
  """Exact spectral amplitudes for the two observable numerators."""
  z = zeta * n ** (1 - d)
  B = population_chain(n, d, z)
  values, left, right = la.eig(B, left=True, right=True)
  order = np.argsort(-values.real)
  values, left, right = values[order], left[:, order], right[:, order]
  initial = np.array([math.comb(n,k) * filling**k * (1-filling)**(n-k) for k in range(n+1)])
  # Residual zero Jordan blocks make the complete eigenvector matrix singular.
  # Only simple nonzero edge roots are projected with their own left vectors.
  coordinates = np.zeros(n + 1, dtype=complex)
  for j in [0, 1]:
    coordinates[j] = np.vdot(left[:, j], initial) / np.vdot(left[:, j], right[:, j])
  stationary = right[:, 0].real / right[:, 0].real.sum()
  nu = np.dot(np.arange(n+1), stationary) / n
  detector = np.arange(n+1)/n - nu
  density_coefficients = np.einsum("k,kj->j", detector, right, optimize=False) * coordinates
  trace_coefficient = right[:, 0].sum() * coordinates[0]
  secondary = 1
  density_amplitude = density_coefficients[secondary] / trace_coefficient
  # Antisymmetric two-mode pair times symmetric spectators is the spin-(n-2)/2 irrep.
  # Its probability basis differs from per-configuration coefficients only by binomial factors.
  B2 = population_chain(n-2, d, z)
  vals2, left2, vec2 = la.eig(B2, left=True, right=True)
  top = np.argmax(vals2.real)
  init2 = np.array([2*bias*math.comb(n-2,k)*filling**(k+1)*(1-filling)**(n-k-1)
                    for k in range(n-1)])
  leading2 = np.vdot(left2[:, top], init2) / np.vdot(left2[:, top], vec2[:, top])
  imbalance_amplitude = vec2[:,top].sum()*leading2/trace_coefficient
  assert abs(density_amplitude.imag)<1e-9 and abs(imbalance_amplitude.imag)<1e-9
  return dict(n=n,d=d,density_asymptotic_amplitude=float(density_amplitude.real),
              imbalance_asymptotic_amplitude=float(imbalance_amplitude.real),
              leading_density_overlap_abs=float(abs(density_coefficients[0])),
              second_density_overlap_abs=float(abs(density_coefficients[1])),
              trace_leading_coefficient=float(trace_coefficient.real))


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--output", default=str(ROOT / "data/verified_observables.json"))
  args = parser.parse_args()
  curves=[]
  for n,d in [(8,2),(12,3),(16,4)]:
    row=dynamics(n,d)
    row['modal_amplitudes']=visible_modes(n,d)
    amplitudes = row['modal_amplitudes']
    AN = amplitudes['density_asymptotic_amplitude']
    AD = amplitudes['imbalance_asymptotic_amplitude']
    times = np.array(row['time'])
    density_model = AN * np.exp(-row['theory']['density'] * times)
    imbalance_model = AD * np.exp(-row['theory']['imbalance'] * times)
    mask = (times >= 6) & (times <= 9)
    row['checks']['asymptotic_density_relative_error'] = float(np.max(abs(density_model[mask] / np.array(row['density_deviation'])[mask] - 1)))
    row['checks']['asymptotic_imbalance_relative_error'] = float(np.max(abs(imbalance_model[mask] / np.array(row['imbalance'])[mask] - 1)))
    assert row['checks']['asymptotic_density_relative_error'] < .003
    assert row['checks']['asymptotic_imbalance_relative_error'] < .003
    for eta in [.01, .1, 1.]:
      crossing = math.log(n * abs(AN) / (eta * abs(AD))) / (row['theory']['density'] - row['theory']['imbalance'])
      row['readouts'][str(eta)]['asymptotic_crossover'] = crossing
    curves.append(row)
    print(n,d,row['theory'],row['density_tail'],row['imbalance_tail'],flush=True)
  checks=[direct_full_check(n,d) for n,d in [(4,2),(4,3),(5,2)]]
  result=dict(curves=curves,full_checks=checks,selection_identities=selection_identity_checks(),
              units="gamma=1; all quoted rates are asymptotic conditional exponents",
              initial_state="rho_p [1+epsilon(n_1-n_2)], p=.75, epsilon=.6",
              readout="N/n + eta (n_1-n_2)/n; no change to dynamics",
              parity_scope="one-time even density observables; odd two-time correlators not excluded")
  Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
  print(args.output)


if __name__ == '__main__':
  main()
