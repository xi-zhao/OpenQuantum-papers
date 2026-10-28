"""Reproducible analysis for the PRA draft; no external services required.

Conventions: N acts on coefficients of unnormalised particle-number projectors.
All eigenvalues below are shifted by +n*gamma and divided by gamma.
Exploratory t*.py scripts are deliberately not imported.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import platform
import time
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse as sp
from scipy.special import gammaln

ROOT = Path(__file__).resolve().parents[1]


def chain(n, d, z=1.0):
  if n < 0 or d < 1 or z < 0:
    raise ValueError("Require n >= 0, d >= 1, z >= 0")
  out = np.zeros((n + 1, n + 1))
  for k in range(1, n + 1):
    out[k - 1, k] = n - k + 1
  for k in range(n + 1 - d):
    out[k + d, k] = z * math.comb(k + d, d)
  return out


def annihilators(n):
  """Jordan-Wigner operators, ordered Fock basis, mode 0 the highest bit."""
  identity = sp.identity(2, format="csr")
  parity = sp.csr_matrix(np.diag([1., -1.]))
  lowering = sp.csr_matrix([[0., 1.], [0., 0.]])
  operators = []
  for mode in range(n):
    factors = [parity] * mode + [lowering] + [identity] * (n - mode - 1)
    operator = factors[0]
    for factor in factors[1:]:
      operator = sp.kron(operator, factor, format="csr")
    operators.append(operator)
  return operators


def full_liouvillian(n, d, z, h):
  """Independent full superoperator; gamma=1, row-major vectorization."""
  if h.shape != (n, n) or not np.allclose(h, h.conj().T):
    raise ValueError("h must be an n-by-n Hermitian single-particle matrix")
  c = annihilators(n)
  identity = sp.identity(2 ** n, format="csr")
  h_many = sum(h[i, j] * (c[i].getH() @ c[j]) for i in range(n) for j in range(n))
  result = -1j * (sp.kron(h_many, identity) - sp.kron(identity, h_many.T))
  result -= n * sp.identity(4 ** n)
  for operator in c:
    result += sp.kron(operator, operator.conj())
  for subset in itertools.combinations(range(n), d):
    operator = identity
    for mode in subset:
      operator = operator @ c[mode].getH()
    result += z * sp.kron(operator, operator.conj())
  return result.toarray()


def boolean_generator(m, d, z):
  """Boolean lattice with one removal or d simultaneous additions."""
  subsets = [frozenset(s) for k in range(m + 1)
             for s in itertools.combinations(range(m), k)]
  index = {subset: i for i, subset in enumerate(subsets)}
  result = np.zeros((2 ** m, 2 ** m))
  universe = frozenset(range(m))
  for subset in subsets:
    for mode in subset:
      result[index[subset - {mode}], index[subset]] += 1
    for added in itertools.combinations(sorted(universe - subset), d):
      result[index[subset | set(added)], index[subset]] += z
  return result


def constants(d):
  alpha = 2 * d / (d + 1)
  a = (2.0 ** (-(d - 1)) / math.factorial(d - 1)) ** (1 / (d + 1))
  omega = math.sqrt((d + 1) / 2)
  return dict(alpha=alpha, beta=alpha - 1, C=a / alpha, a=a,
              omega=omega, first_correction=a * (1 - omega),
              internal_gap=2 * a * omega)


def polynomial_coefficients(n, d):
  """Integer W_j in det(mu I-N)=sum_j (-z)^j W_j mu^(n+1-(d+1)j).

  The leading-principal-minor recurrence adds one isolated vertex or a
  length-(d+1) directed cycle. Coefficients are exact Python integers.
  """
  polys = [[1]]
  for length in range(1, n + 2):
    current = polys[-1].copy()
    if length >= d + 1:
      start = length - d - 1
      weight = math.comb(start + d, d) * math.prod(
        n - start - r + 1 for r in range(1, d + 1))
      old = polys[length - d - 1]
      current += [0] * max(0, len(old) + 1 - len(current))
      for j, value in enumerate(old):
        current[j + 1] += weight * value
    polys.append(current)
  return polys[-1]


def perron_log(n, d, rtol=2e-12, maxiter=250000):
  """Positive log-space iteration, with Collatz bounds in floating arithmetic.

  Calculate at z=n^(-(d-1)), then use exact homogeneity to return z=1.
  Every component remains represented even if the vector spans >300 orders
  of magnitude. Bounds are numerical certificates, not interval arithmetic.
  """
  if n < d:
    return dict(n=n, d=d, value=0., lower=0., upper=0., iterations=0,
                converged=True, relative_width=0., log_vector_span=0.)
  z = n ** (-(d - 1.0))
  q = z ** (1 / (d + 1))
  shift = float(n + 1)
  log_shift = math.log(shift)
  loss = np.log(np.arange(n, 0, -1, dtype=float))
  k = np.arange(d, n + 1, dtype=float)
  gain = math.log(z) + gammaln(k + 1) - gammaln(k - d + 1) - gammaln(d + 1)
  lv = np.zeros(n + 1)
  lower = upper = float("nan")
  for iteration in range(1, maxiter + 1):
    lw = lv + log_shift
    lw[:-1] = np.logaddexp(lw[:-1], loss + lv[1:])
    lw[d:] = np.logaddexp(lw[d:], gain + lv[:-d])
    if iteration % 10 == 0:
      ratios = lw - lv
      lower = (math.exp(float(ratios.min())) - shift) / q
      upper = (math.exp(float(ratios.max())) - shift) / q
      width = (upper - lower) / max(abs(lower), 1e-100)
      if lower > 0 and width < rtol:
        break
    lv = lw - lw.max()
  else:
    raise RuntimeError(f"Perron iteration did not converge: n={n}, d={d}, [{lower}, {upper}]")
  return dict(n=n, d=d, value=(upper + lower) / 2, lower=lower, upper=upper,
              iterations=iteration, converged=True, relative_width=width,
              log_vector_span=float(lv.max() - lv.min()))


def mp_spectrum(n, d, digits=60):
  """Independent arbitrary-precision dense eigenproblem (z=1)."""
  with mp.workdps(digits):
    mat = mp.matrix(n + 1)
    for k in range(1, n + 1):
      mat[k - 1, k] = n - k + 1
    for k in range(n + 1 - d):
      mat[k + d, k] = math.comb(k + d, d)
    values = sorted(mp.eig(mat, left=False, right=False), key=lambda x: -mp.re(x))
    return dict(n=n, d=d, digits=digits, perron=float(mp.re(values[0])),
                gap=float(mp.re(values[0] - values[1])),
                next_imag=float(mp.im(values[1])))


def check_algebra():
  import sympy as sy
  mu, z = sy.symbols("mu z")
  cases = []
  for n, d in [(3, 1), (4, 2), (7, 2), (8, 3), (9, 3), (7, 4)]:
    matrix = sy.zeros(n + 1)
    for k in range(1, n + 1):
      matrix[k - 1, k] = n - k + 1
    for k in range(n + 1 - d):
      matrix[k + d, k] = z * math.comb(k + d, d)
    coeffs = polynomial_coefficients(n, d)
    expected = sum((-z) ** j * value * mu ** (n + 1 - (d + 1) * j)
                   for j, value in enumerate(coeffs))
    actual = matrix.charpoly(mu).as_expr()
    assert sy.expand(actual - expected) == 0, (n, d)
    remainder = (n + 1) % (d + 1)
    nullity = n + 1 - matrix.subs(z, 1).rank()
    assert nullity == int(remainder > 0)
    cases.append(dict(n=n, d=d, coefficients=coeffs, remainder=remainder,
                      nullity=nullity, exact_polynomial=True))
  homogeneity = []
  for n, d, zvalue in [(7, 2, 1e-6), (9, 3, 1e-4), (12, 4, .03)]:
    q = zvalue ** (1 / (d + 1))
    scale = q ** np.arange(n + 1)
    actual = chain(n, d, zvalue) * scale[None, :] / scale[:, None]
    relative = float(np.max(abs(actual - q * chain(n, d))) / np.max(abs(actual)))
    assert relative < 2e-14
    homogeneity.append(dict(n=n, d=d, z=zvalue, matrix_relative_error=relative))
  d1 = []
  for n in [4, 12, 30]:
    values = np.sort(la.eigvals(chain(n, 1)).real)
    error = float(np.max(abs(values - np.arange(-n, n + 1, 2))))
    assert error < 1e-8
    d1.append(dict(n=n, exact_spectrum_error=error))
  trace_domain = []
  for n, d in [(3, 2), (4, 3), (4, 4)]:
    upper = n / math.comb(n, d)
    for fraction in [.2, 1., 1.01]:
      strength = fraction * upper
      operator = full_liouvillian(n, d, strength, np.diag(np.arange(n, dtype=float)))
      dimension = 2 ** n
      # Avoid backend floating-status warnings in BLAS matrix multiplication.
      actual = np.einsum("i,ij->j", np.eye(dimension).reshape(-1),
                         operator, optimize=False).reshape((dimension, dimension))
      occupations = [i.bit_count() for i in range(dimension)]
      expected = np.array([-n + k + strength * (math.comb(n - k, d) if n - k >= d else 0)
                           for k in occupations])
      error = float(np.max(abs(actual - np.diag(expected))))
      assert np.isfinite(actual).all() and error < 1e-12
      assert (expected.max() <= 1e-12) == (fraction <= 1.)
      trace_domain.append(dict(n=n, d=d, fraction_of_boundary=fraction,
                               matrix_error=error, max_trace_growth=float(expected.max())))
  return dict(polynomial=cases, homogeneity=homogeneity, d1=d1,
              trace_domain=trace_domain)


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("mode", choices=["algebra", "perron", "precision"])
  args = parser.parse_args()
  start = time.time()
  if args.mode == "algebra":
    result = check_algebra()
  elif args.mode == "perron":
    ns = [16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]
    result = dict(constants={d: constants(d) for d in [1, 2, 3, 4]}, points=[])
    for d in [2, 3, 4]:
      for n in ns:
        row = perron_log(n, d)
        previous = perron_log(n - 1, d)
        before = perron_log(n - 2, d)
        row.update(gap_one=row["value"] - previous["value"],
                   gap_two=row["value"] - before["value"],
                   gap_one_lower=row["lower"] - previous["upper"],
                   gap_one_upper=row["upper"] - previous["lower"],
                   gap_two_lower=row["lower"] - before["upper"],
                   gap_two_upper=row["upper"] - before["lower"],
                   adjacent_certificates=[previous, before])
        result["points"].append(row)
        print(f"d={d} n={n}: mu={row['value']:.12g}; width={row['relative_width']:.2g}; iterations={row['iterations']}", flush=True)
  else:
    result = dict(points=[])
    for d in [2, 3, 4]:
      for n in [16, 32, 48, 64]:
        low = mp_spectrum(n, d, digits=50)
        row = mp_spectrum(n, d, digits=70)
        row["precision_comparison"] = dict(perron_difference=abs(row["perron"] - low["perron"]),
                                           gap_difference=abs(row["gap"] - low["gap"]))
        cert = perron_log(n, d)
        row["positive_iteration_error"] = abs(row["perron"] - cert["value"]) / row["perron"]
        assert row["positive_iteration_error"] < 2e-11
        assert row["precision_comparison"]["gap_difference"] < 1e-10
        result["points"].append(row)
        print(row, flush=True)
  result["provenance"] = dict(python=platform.python_version(), numpy=np.__version__,
                              scipy=scipy.__version__, mpmath=mp.__version__,
                              elapsed_seconds=time.time() - start)
  path = ROOT / "data" / f"verified_{args.mode}.json"
  path.write_text(json.dumps(result, indent=2) + "\n")
  print(path, flush=True)


if __name__ == "__main__":
  main()
