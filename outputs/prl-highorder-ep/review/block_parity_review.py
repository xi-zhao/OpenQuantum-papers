#!/usr/bin/env python3
r"""Independent exact block/parity audit; no existing file is modified.

Run (only pre-existing NumPy/SciPy are needed):
  python3 -I outputs/prl-highorder-ep/review/block_parity_review.py

ASSUMPTIONS AND PROOF
=====================
There are n CAR modes, gamma>0, z>0, integer 1<=d<=n, and an arbitrary
Hermitian one-particle h.  Define c_I^dagger=c_i1^dagger...c_id^dagger for
increasing I, and c_I as its adjoint.  All d-subsets have the same weight:

 L(rho)=-i[sum h_ij c_i^dagger c_j,rho]
        +gamma sum c_i rho c_i^dagger-n gamma rho
        +z gamma sum_{|I|=d} c_I^dagger rho c_I.

1. Covariance under a number-conserving mode rotation.
Let h=U diag(epsilon) U^dagger and c_i=sum_a U_ia f_a.  The loss sum is
invariant by U^dagger U=1.  The pump creation products transform by
 c_I^dagger=sum_J det(U_{I,J})^* f_J^dagger.
The matrix of minors is the unitary exterior power wedge^d U, so summing
c_I^dagger rho c_I eliminates cross terms and gives exactly the same pump
sum in f modes.  Fock-space rotation has matrix V_{X,Y}=det(U_{X,Y}) for
|X|=|Y| (vacuum minor=1); row-major operator rotation is W=V kron V^*.
This statement requires the complete, equally weighted d-subset family;
a restricted/geometric pump set is not generally covariant.  It does not
assert covariance under a Bogoliubov transformation.

2. Frozen sets and the complete fermion-sign gauge.
Use ordered Fock states |X>=prod_{i in X increasing} f_i^dagger|vac>.
Every matrix unit has the unique form |A union C><B union C|, with A,B,C
pairwise disjoint.  Fix A,B and let F=A union B, R=[n]\F, m=|R|.
Loss removes a common occupied i from C; pump adds a d-subset I of R\C.
Both preserve A,B.  The Hamiltonian contributes -i(E_A-E_B) throughout.
Put s_F(i)=(-1)^{|{f in F:f<i}|} and g_F(C)=prod_{i in C}s_F(i).
For loss, the product of ket/bra CAR signs is s_F(i).  For pump by I,
apply the ordered creation operators successively: common occupations
and occupations added in earlier steps occur twice in the ket/bra sign,
leaving prod_{i in I}s_F(i).  Consequently the gauged basis
 e_C=g_F(C)|A union C><B union C|
turns every allowed loss amplitude into +1 and every pump amplitude into
+z.  Explicitly, the transformed sign is
 g_F(C)/g_F(C\{i}) * s_F(i)=1,
 g_F(C)/g_F(C union I) * prod_{i in I}s_F(i)=1.
The block is exactly (-n gamma-i(E_A-E_B))1+gamma G_m, with
 G_m=J_-+z (J_+)^d/d!,  J_+=sum_{i in R} sigma_i^+,
 J_-=J_+^dagger.  Distinct raising operators commute and square to zero,
so division by d! counts each d-subset once.  The diagonal gauge is
unitary on Hilbert-Schmidt operator space but need not preserve positivity.

3. Exact SU(2) decomposition and the pblock.py convention.
The Boolean tensor product is (C^2)^{tensor m}.  For j=0,...,floor(m/2),
spin s=(m-2j)/2 occurs M(m,j)=binom(m,j)-binom(m,j-1) times (binom(m,-1)=0).
This follows by subtracting dimensions of successive J_z weight spaces.
Put q=m-2j.  On one irrep choose a lowest-weight vector v_0 annihilated
by J_-, and v_k=(J_+)^k v_0/k!, k=0,...,q.  Then
 J_- v_k=(q-k+1)v_{k-1},  J_+ v_k=(k+1)v_{k+1}.
Thus G_m acts on this chain by N_{q,d}(z), with column-source entries
 N[k-1,k]=q-k+1, N[k+d,k]=z binom(k+d,d), all others zero.
These are precisely analysis.chain(q,d,z), including q<d where N is nilpotent.
This is an exact similarity decomposition, not merely equality of spectra.

4. Spectrum, multiplicities and operator parity.
For every disjoint A,B, every j=0,...,floor(m/2), and every eigenvalue mu
of N_{m-2j,d}(z), L has eigenvalue
 -n gamma-i(E_A-E_B)+gamma mu
with multiplicity M(m,j) times its algebraic multiplicity in that N block.
Accidental coincidences between frozen sets/spins must be added.
The dimension identity sum_{m=0}^n binom(n,m)2^{n-m}2^m=4^n checks
completeness.  Operator parity P rho P equals (-1)^{|A|+|B|} on the block.
The physical even operator space therefore includes exactly |A|+|B| even,
equivalently m congruent to n mod 2.  It is not just A=B=empty, and it is
not the restriction of rho to even many-body states: odd-to-odd matrix
units are also physical even operators.  Its complex dimension is 4^n/2
for n>=1.  SU(2) descendants have q congruent to m mod 2.

5. Perron monotonicity and rigorous gap formulas.
For m>=d, the directed Boolean graph is strongly connected.  From any
subset, losses reach empty.  From empty, add a d-set; if size k>=d is not
full, retain d-1 occupied sites, remove them, then pump those d-1 sites
and one new empty site.  This raises k by one.  Reach full and then remove
unwanted sites.  G_m is irreducible and nonnegative.  Its positive Perron
vector is permutation-invariant by uniqueness, hence lies in maximal
spin: r_m=r(G_m)=r(N_{m,d}).  For m<d the graph only lowers, so r_m=0.
The subspace of subsets omitting one new site gives G_m as a principal
submatrix of G_{m+1}.  Therefore r_m is nondecreasing, and is strictly
increasing once m+1>=d (strict Perron comparison with an irreducible
nonnegative matrix).  The leading eigenvalue of L is unique, real and
lambda_0=-n gamma+gamma r_n; every Hamiltonian frozen-set shift is purely
imaginary.  Let b_n=max Re(spec N_{n,d}\{r_n}), excluding the simple Perron
root.  The exact spectral separations in real part are
 Delta_full=gamma [r_n-max(b_n,r_{n-1})],
 Delta_even=gamma [r_n-max(b_n,r_{n-2})]  (n>=2).
The n=1 even gap is simply gamma(r_1-b_1).  All other blocks are bounded
by these candidates using monotonicity; the r_{n-1} and r_{n-2} candidates
actually occur.  No sorting of internal roots of N_n has been assumed.
Let t=z^{1/(d+1)} and D[k,k]=t^k.  Then D^{-1}N(z)D=t N(1).
With kappa_m=r(N_{m,d}(1)), beta_n=max Re(spec N_{n,d}(1)\{kappa_n}),
 delta_n=kappa_n-beta_n, the same formulas are
 Delta_full=gamma t min(delta_n,kappa_n-kappa_{n-1}),
 Delta_even=gamma t min(delta_n,kappa_n-kappa_{n-2}).
These hold for all 1<=d<=n, not only large n.  If n<d, pumping vanishes,
all eigenvalues have real part -n gamma and a unique leading relaxation
gap is unavailable.  No conclusion about a globally second root can be
obtained from a local Holstein-Primakoff ladder alone.

A concrete warning against dropping the internal-root candidate is
n=d=4,z=1: N_4 is a weighted five-cycle, N_4^5=4! 1, so its roots are
(4!)^{1/5} exp(2 pi i k/5).  The largest non-Perron real part is
(4!)^{1/5} cos(2 pi/5)>0, whereas kappa_3=kappa_2=0.  The internal roots
therefore control both full and even gaps in this example.

6. d=1 exact check and physical meaning of the gap.
G_m=sum_i(sigma_i^-+z sigma_i^+).  The local 2x2 eigenvalues are +/-sqrt(z),
so its spectrum is (m-2l)sqrt(z), l=0,...,m, multiplicity binom(m,l).
On one spin-q/2 irrep N_{q,1}, the same ladder is simple, l=0,...,q.
Consequently Delta_full=gamma sqrt(z), Delta_even=2 gamma sqrt(z),
including n=1 for the even formula.  The positive leading operator is
R_0 proportional to z^{Nhat/2}; its normalized filling is
sqrt(z)/(1+sqrt(z)).  At z=1 this is 1/2.
The map as written is generally not trace preserving:
 Tr L(rho)=gamma Tr[(Nhat-n+z binom(n-Nhat,d))rho].
For n>=1, gamma>0, 1<=d<=n it is trace preserving iff d=1,z=1.  Thus the
formulas refer to separation relative to lambda_0, not automatically to
the decay gap of a trace-preserving Lindbladian with stationary eigenvalue
zero.  Positive normalized conditional dynamics has R_0 as its leading
state; the even gap controls its slowest allowed asymptotic exponent,
subject to actual overlap and possible subleading Jordan polynomial factors.

CODE-SPECIFIC CHECK
==================
The script deliberately avoids eigvalue matching of high-order defective
blocks.  It verifies matrix covariance, explicit sign-gauge block equality,
SU(2) intertwining, principal submatrices, and exact integer CAR signs.
Only the explicitly diagonalizable d=1 ladder is eigvalue-checked.
t11.py's Hm.T.conj().T equals Hm.T for Hermitian Hm, so it is correct under
the stated input contract despite the opaque notation.  pblock.py's d=1
semiclassical gap equals the internal/even gap, not the full operator gap.
t19.py's filling root equation becomes identically zero for d=1,zeta=1;
brentq returns the endpoint 0 and does not verify the correct filling 1/2.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path
import sys
import warnings

import numpy as np
from scipy import sparse
from scipy.linalg import null_space
from scipy.optimize import brentq

# Import only the explicitly audited source modules, without writing pycache.
sys.dont_write_bytecode = True
CODE_DIR = Path(__file__).resolve().parent.parent / "code"


def load_module(name: str, filename: str):
  spec = importlib.util.spec_from_file_location(name, CODE_DIR / filename)
  module = importlib.util.module_from_spec(spec)
  sys.modules[name] = module
  spec.loader.exec_module(module)
  return module


analysis = load_module("spectral_analysis", "spectral_analysis.py")


def subsets(items):
  items = tuple(items)
  return [tuple(s) for k in range(len(items) + 1)
          for s in itertools.combinations(items, k)]


def fock_index(occupied, n):
  return sum(1 << (n - 1 - i) for i in occupied)


def dense_product(a, b):
  # Avoid backend-dependent BLAS floating-status warnings seen on this host.
  return np.einsum("ij,jk->ik", a, b, optimize=False)


def conjugate_sparse(matrix, basis):
  result = basis.conj().T @ (sparse.csr_matrix(matrix) @ basis)
  return result.toarray()


def fock_rotation(unitary):
  n = len(unitary)
  dimension = 2 ** n
  result = np.zeros((dimension, dimension), dtype=complex)
  all_sets = subsets(range(n))
  for occupied in all_sets:
    for modes in all_sets:
      if len(occupied) == len(modes):
        result[fock_index(occupied, n), fock_index(modes, n)] = (
          np.linalg.det(unitary[np.ix_(occupied, modes)]) if occupied else 1
        )
  return result


def frozen_block_matrix(n, d, z, epsilon):
  dimension = 2 ** n
  columns = []
  rows = []
  signs = []
  blocks = []
  parities = []
  for frozen in subsets(range(n)):
    remaining = tuple(i for i in range(n) if i not in frozen)
    m = len(remaining)
    for a_set in subsets(frozen):
      b_set = tuple(i for i in frozen if i not in a_set)
      energy = sum(epsilon[list(a_set)]) - sum(epsilon[list(b_set)])
      blocks.append((-n - 1j * energy) * np.eye(2 ** m) + analysis.boolean_generator(m, d, z))
      for common in subsets(remaining):
        ket = set(a_set) | set(common)
        bra = set(b_set) | set(common)
        rows.append(fock_index(ket, n) * dimension + fock_index(bra, n))
        columns.append(len(columns))
        signs.append((-1) ** sum(f < i for f in frozen for i in common))
        parities.append(len(frozen) % 2)
  basis = sparse.csr_matrix((signs, (rows, columns)),
                            shape=(dimension * dimension,) * 2)
  return basis, sparse.block_diag(blocks, format="csr").toarray(), parities


def residual(actual, expected, label, tolerance=2e-11):
  error = float(np.max(np.abs(actual - expected)))
  if not np.isfinite(error) or error > tolerance:
    raise AssertionError(f"{label}: residual {error} exceeds {tolerance}")
  return error


def check_random_complex_h(max_n, seed):
  rng = np.random.default_rng(seed)
  reports = []
  for n in range(1, max_n + 1):
    for d in range(1, n + 1):
      raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
      hamiltonian = (raw + raw.conj().T) / 2
      epsilon, unitary = np.linalg.eigh(hamiltonian)
      v = fock_rotation(unitary)
      w = np.kron(v, v.conj())
      w_sparse = sparse.csr_matrix(w)
      z = 0.37
      exact = analysis.full_liouvillian(n, d, z, hamiltonian)
      rotated = conjugate_sparse(exact, w_sparse)
      diagonal = analysis.full_liouvillian(n, d, z, np.diag(epsilon))
      gauge, blocks, parities = frozen_block_matrix(n, d, z, epsilon)
      gauged = conjugate_sparse(rotated, gauge)
      even_dimension = sum(p == 0 for p in parities)
      assert even_dimension == 2 ** (2 * n - 1)
      reports.append({
        "n": n, "d": d,
        "fock_unitarity_residual": residual(
          dense_product(v.conj().T, v), np.eye(2 ** n), "Fock rotation"),
        "covariance_residual": residual(rotated, diagonal, "mode covariance"),
        "sign_gauge_block_residual": residual(gauged, blocks, "frozen blocks"),
        "even_operator_dimension": even_dimension,
      })
  return reports


def car_create(occupied, indices):
  occupied = set(occupied)
  sign = 1
  for i in reversed(indices):
    if i in occupied:
      return None, 0
    sign *= (-1) ** sum(j < i for j in occupied)
    occupied.add(i)
  return occupied, sign


def car_annihilate(occupied, i):
  occupied = set(occupied)
  if i not in occupied:
    return None, 0
  sign = (-1) ** sum(j < i for j in occupied)
  occupied.remove(i)
  return occupied, sign


def check_integer_signs(max_n=5):
  count = 0
  for n in range(1, max_n + 1):
    for ket in subsets(range(n)):
      for bra in subsets(range(n)):
        frozen = set(ket) ^ set(bra)
        common = set(ket) & set(bra)
        def gauge(occupied):
          return (-1) ** sum(f < i for f in frozen for i in occupied)
        for i in common:
          new_ket, left_sign = car_annihilate(ket, i)
          new_bra, right_sign = car_annihilate(bra, i)
          assert gauge(common) * left_sign * right_sign == gauge(new_ket & new_bra)
          count += 1
        empty = tuple(i for i in range(n) if i not in set(ket) | set(bra))
        for d in range(1, n + 1):
          for pump in itertools.combinations(empty, d):
            new_ket, left_sign = car_create(ket, pump)
            new_bra, right_sign = car_create(bra, pump)
            assert gauge(common) * left_sign * right_sign == gauge(new_ket & new_bra)
            count += 1
  return count


def lowering_matrix(m):
  all_sets = subsets(range(m))
  indices = {frozenset(s): k for k, s in enumerate(all_sets)}
  lower = np.zeros((2 ** m,) * 2)
  for col, occupied in enumerate(all_sets):
    for i in occupied:
      lower[indices[frozenset(occupied) - {i}], col] = 1
  return all_sets, lower


def check_spin_intertwiners(max_m=7):
  reports = []
  for m in range(max_m + 1):
    all_sets, lower = lowering_matrix(m)
    raise_op = lower.T
    chains = []
    q_values = []
    for j in range(m // 2 + 1):
      q = m - 2 * j
      sector = np.array([i for i, occupied in enumerate(all_sets) if len(occupied) == j])
      if j == 0:
        kernel = np.ones((1, 1))
      else:
        previous = np.array([i for i, occupied in enumerate(all_sets)
                             if len(occupied) == j - 1])
        kernel = null_space(lower[np.ix_(previous, sector)])
      multiplicity = math.comb(m, j) - (math.comb(m, j - 1) if j else 0)
      assert kernel.shape[1] == multiplicity
      for copy in range(multiplicity):
        vector = np.zeros(2 ** m)
        vector[sector] = kernel[:, copy]
        chain = [vector]
        for k in range(1, q + 1):
          chain.append(np.einsum("ij,j->i", raise_op, chain[-1], optimize=False) / k)
        chains.append(np.column_stack(chain))
        q_values.append(q)
    basis = np.column_stack(chains)
    assert basis.shape == (2 ** m, 2 ** m)
    norms = np.linalg.norm(basis, axis=0)
    orthonormal = basis / norms
    orthogonality = residual(dense_product(orthonormal.T, orthonormal),
                             np.eye(2 ** m), "spin-chain orthogonality")
    errors = []
    for d in range(1, max(1, m) + 1):
      z = .37
      predicted = sparse.block_diag(
        [analysis.chain(q, d, z) for q in q_values], format="csr")
      actual = sparse.csr_matrix(analysis.boolean_generator(m, d, z)) @ basis
      expected = (predicted.T @ basis.T).T
      errors.append(residual(actual, expected, "SU(2) intertwiner"))
    reports.append({"m": m, "irreducible_copies": len(chains),
                    "orthogonality_residual": orthogonality,
                    "max_intertwiner_residual": max(errors)})
  return reports


def check_principal_submatrices(max_m=6):
  count = 0
  for m in range(max_m + 1):
    larger_sets = subsets(range(m + 1))
    retained = [i for i, occupied in enumerate(larger_sets) if m not in occupied]
    for d in range(1, max_m + 2):
      larger = analysis.boolean_generator(m + 1, d, .37)
      smaller = analysis.boolean_generator(m, d, .37)
      assert np.array_equal(larger[np.ix_(retained, retained)], smaller)
      count += 1
  return count


def check_d1(max_n=10):
  z = .37
  errors = []
  for n in range(1, max_n + 1):
    eigenvalues = np.linalg.eigvals(analysis.chain(n, 1, z))
    expected = (n - 2 * np.arange(n + 1)) * np.sqrt(z)
    errors.append(residual(np.sort(eigenvalues.real), np.sort(expected),
                           "d=1 exact ladder", tolerance=1e-10))
    assert np.max(np.abs(eigenvalues.imag)) < 1e-10
  n = 3
  rng = np.random.default_rng(328)
  raw = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
  hamiltonian = (raw + raw.conj().T) / 2
  liouvillian = analysis.full_liouvillian(n, 1, z, hamiltonian)
  full_eigenvalues = np.linalg.eigvals(liouvillian)
  even = [fock_index(x, n) * 2 ** n + fock_index(y, n)
          for x in subsets(range(n)) for y in subsets(range(n))
          if (len(x) + len(y)) % 2 == 0]
  even_eigenvalues = np.linalg.eigvals(liouvillian[np.ix_(even, even)])
  leading = -n + n * np.sqrt(z)
  full_gap = leading - np.sort(full_eigenvalues.real)[-2]
  even_gap = leading - np.sort(even_eigenvalues.real)[-2]
  assert abs(full_gap - np.sqrt(z)) < 1e-10
  assert abs(even_gap - 2 * np.sqrt(z)) < 1e-10
  probabilities = np.array([z ** (i.bit_count() / 2) for i in range(2 ** n)])
  probabilities /= probabilities.sum()
  right_operator = np.diag(probabilities).reshape(-1)
  right_residual = residual(np.einsum("ij,j->i", liouvillian, right_operator,
                                     optimize=False),
                            leading * right_operator, "d=1 positive eigenoperator")
  endpoint = brentq(lambda x: x + (1 - x) - 1, 0, 1)
  return {"max_ladder_error": max(errors), "full_gap": float(full_gap),
          "even_gap": float(even_gap), "expected_full_gap": float(np.sqrt(z)),
          "expected_even_gap": float(2 * np.sqrt(z)),
          "right_operator_residual": right_residual,
          "t19_d1_zeta1_brentq_endpoint": endpoint,
          "exact_d1_zeta1_filling": .5}


def check_counterexample():
  matrix = analysis.chain(4, 4, 1)
  power = np.eye(5)
  for _ in range(5):
    power = dense_product(power, matrix)
  assert np.array_equal(power, 24 * np.eye(5))
  radius = 24 ** (1 / 5)
  internal = radius * math.cos(2 * math.pi / 5)
  assert internal > 0
  return {"n": 4, "d": 4, "z": 1, "kappa_n": radius,
          "max_internal_nonperron_real_part": internal,
          "kappa_n_minus_1": 0, "kappa_n_minus_2": 0,
          "full_and_even_gap_for_gamma1": radius - internal}


def main():
  parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
  parser.add_argument("--max-n", type=int, default=4,
                      help="largest full-matrix n, 1..4 (default 4)")
  parser.add_argument("--seed", type=int, default=20261028)
  args = parser.parse_args()
  if args.max_n not in range(1, 5):
    parser.error("--max-n must be 1..4; the audit intentionally avoids large dense files")
  warnings.filterwarnings("error", category=RuntimeWarning)
  report = {
    "random_complex_h_matrix_identities": check_random_complex_h(args.max_n, args.seed),
    "exact_integer_sign_checks": check_integer_signs(),
    "spin_intertwiners": check_spin_intertwiners(),
    "principal_submatrix_checks": check_principal_submatrices(),
    "d1_exact_checks": check_d1(),
    "internal_root_counterexample": check_counterexample(),
  }
  print(json.dumps(report, indent=2))


if __name__ == "__main__":
  main()
