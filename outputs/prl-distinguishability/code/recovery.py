#!/usr/bin/env python3
"""Exact Clifford recovery, small-system SDPs, and information-spectrum bounds.

No state-vector surrogate or independent-error-image approximation is used.
Run --simulate once; the default regenerates figures from saved encoder samples.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.stats import binom, norm, unitary_group

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261028
SIZES = ((4, 1), (8, 2), (12, 3), (16, 4))
SAMPLES = 64
R_VALUES = np.array([0.0, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.34, 0.42, 0.55, 1.0])
BLUE, ORANGE, INK = "#2a78d6", "#eb6834", "#171717"
PAULIS = (
    np.eye(2, dtype=complex),
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.diag([1, -1]).astype(complex),
)


def mm(a, b):
    return np.einsum("ij,jk->ik", a, b, optimize=False)


def entropy_sigma(r):
    p = np.array([1 - 3 * r / 4, r / 4, r / 4, r / 4])
    p = p[p > 0]
    s = -np.log2(p)
    h = float(p @ s)
    return h, math.sqrt(float(p @ (s - h) ** 2))


def spectrum_cdf(n, r, threshold):
    """P[-log2 p_P < threshold], retaining the strict inequality."""
    if r == 0:
        return float(threshold > 0)
    if r == 1:
        return float(threshold > 2 * n)
    q = 3 * r / 4
    p0, p1 = 1 - q, r / 4
    base = -n * math.log2(p0)
    step = math.log2(p0 / p1)
    cutoff = min(n, max(-1, math.ceil((threshold - base) / step) - 1))

    def below(weight):
        surprisal = -(n - weight) * math.log2(p0) - weight * math.log2(p1)
        guard = 16 * np.finfo(float).eps * max(1.0, abs(surprisal), abs(threshold))
        if abs(surprisal - threshold) > guard:
            return surprisal < threshold
        # Resolve lattice boundaries without rounding the requested threshold.
        from decimal import Decimal, localcontext
        with localcontext() as context:
            context.prec = 80
            log_two = Decimal(2).ln()

            def log2_exact_power(value):
                numerator, denominator = value.as_integer_ratio()
                if numerator & (numerator - 1) == 0:
                    return Decimal(numerator.bit_length() - denominator.bit_length())
                return Decimal.from_float(value).ln() / log_two

            precise = -(n - weight) * log2_exact_power(p0) - weight * log2_exact_power(p1)
            return precise < Decimal.from_float(float(threshold))

    while cutoff >= 0 and not below(cutoff):
        cutoff -= 1
    while cutoff < n and below(cutoff + 1):
        cutoff += 1
    return float(binom.cdf(cutoff, n, q))


def spectrum_bracket(n, k, r, a=None):
    if a is None:
        a = 2 * math.log2(n)
    lower = spectrum_cdf(n, r, n - k - a) - 8 * 2.0 ** (-a)
    residual = math.exp2(a - k - 2) if a - k - 2 > -1074 else 0.0
    upper = spectrum_cdf(n, r, n - k + 2 * a) + 1.25 * 2.0 ** (-a) + residual
    return max(0.0, lower), min(1.0, upper)


def pauli_matrices(n):
    matrices = []
    words = list(itertools.product(range(4), repeat=n))
    for word in words:
        matrix = np.ones((1, 1), complex)
        for digit in word:
            matrix = np.kron(matrix, PAULIS[digit])
        matrices.append(matrix)
    return np.array(matrices), words


def fwht(values):
    """Unnormalized Walsh-Hadamard transform, in place."""
    result = np.asarray(values, dtype=float).copy()
    half = 1
    while half < len(result):
        blocks = result.reshape(-1, 2 * half)
        left = blocks[:, :half].copy()
        right = blocks[:, half:].copy()
        blocks[:, :half] = left + right
        blocks[:, half:] = left - right
        half *= 2
    return result


def retained_fourier_weights(inverse_symplectic, n, k):
    """Weights of dual characters for y=(X_aux,X_logical,Z_logical).

    Qiskit uses row Pauli bits (x_0,...,x_n-1,z_0,...,z_n-1).
    Multiplication by U^dagger's symplectic matrix implements U^dagger P U.
    Low bits of each retained integer are auxiliary-X syndrome bits.
    """
    retained = list(range(k, n)) + list(range(k)) + list(range(n, n + k))
    matrix = np.asarray(inverse_symplectic, dtype=np.uint8)[:, retained]
    masks = []
    for column in matrix.T:
        masks.append(sum(int(bit) << j for j, bit in enumerate(column)))
    images = np.zeros(1 << (n + k), dtype=np.uint64)
    length = 1
    for mask in masks:
        images[length:2 * length] = images[:length] ^ np.uint64(mask)
        length *= 2
    physical_support = (images & np.uint64((1 << n) - 1)) | (images >> np.uint64(n))
    return np.bitwise_count(physical_support).astype(np.uint8)


def clifford_distribution(weights, n, k, r):
    """Exact syndrome/logical Pauli probabilities by a classical Fourier sum."""
    characteristic = (1.0 - r) ** weights
    probability = fwht(characteristic) / len(characteristic)
    if float(probability.min()) < -2e-12:
        raise ArithmeticError("Fourier probability has a substantial negative entry")
    probability = np.maximum(probability, 0)
    probability /= probability.sum()
    return probability.reshape(1 << k, 1 << k, 1 << (n - k))


def optimal_clifford_scores(probability):
    """Axes are Z_logical, X_logical, syndrome; degeneracy is already summed."""
    quantum = float(probability.max(axis=(0, 1)).sum())
    classical = float(probability.sum(axis=0).max(axis=0).sum())
    assert -1e-12 <= quantum <= classical + 1e-12 <= 1 + 2e-12
    return classical, quantum


def direct_clifford_distribution(clifford, p1, k):
    """Independent small-n enumeration using Qiskit's Pauli conjugation."""
    from qiskit.quantum_info import Pauli
    n = clifford.num_qubits
    out = np.zeros((1 << k, 1 << k, 1 << (n - k)))
    for digits in itertools.product(range(4), repeat=n):
        label = "".join("IXYZ"[digit] for digit in digits)
        transformed = Pauli(label).evolve(clifford, frame="h")
        x = sum(int(v) << i for i, v in enumerate(transformed.x))
        z = sum(int(v) << i for i, v in enumerate(transformed.z))
        weight = math.prod(p1[digit] for digit in digits)
        out[z & ((1 << k) - 1), x & ((1 << k) - 1), x >> k] += weight
    return out


def output_states(V, matrices, probabilities):
    blocks = [mm(P, V) for P in matrices]
    states = [sum(p * np.outer(B[:, j], B[:, j].conj())
                  for p, B in zip(probabilities, blocks)) for j in range(V.shape[1])]
    return blocks, states


def solve_optimal_scores(V, matrices, probabilities):
    """Independent primal and dual SDPs; includes feasibility/gap checks."""
    import cvxpy as cp
    d, M = V.shape
    blocks, states = output_states(V, matrices, probabilities)
    effects = [cp.Variable((d, d), hermitian=True) for _ in range(M)]
    classical = cp.Problem(
        cp.Maximize(cp.real(sum(cp.trace(E @ rho) for E, rho in zip(effects, states))) / M),
        [E >> 0 for E in effects] + [sum(effects) == np.eye(d)],
    )
    Y = cp.Variable((d, d), hermitian=True)
    classical_dual = cp.Problem(cp.Minimize(cp.real(cp.trace(Y))), [Y - rho / M >> 0 for rho in states])
    vectors = [B.T.conj().ravel() for B in blocks]
    C = sum(p * np.outer(v, v.conj()) for p, v in zip(probabilities, vectors)) / M ** 2
    J = cp.Variable((M * d, M * d), hermitian=True)
    partial = sum(J[j * d:(j + 1) * d, j * d:(j + 1) * d] for j in range(M))
    quantum = cp.Problem(cp.Maximize(cp.real(cp.trace(C @ J))), [J >> 0, partial == np.eye(d)])
    Z = cp.Variable((d, d), hermitian=True)
    quantum_dual = cp.Problem(cp.Minimize(cp.real(cp.trace(Z))), [cp.kron(np.eye(M), Z) - C >> 0])
    problems = (classical, classical_dual, quantum, quantum_dual)
    attempts = []
    for name, problem in zip(("classical_primal", "classical_dual", "quantum_primal", "quantum_dual"), problems):
        problem.solve(solver="CLARABEL", tol_gap_abs=1e-8, tol_feas=1e-8, tol_gap_rel=1e-8)
        attempts.append({"problem": name, "solver": "CLARABEL", "status": problem.status,
                         "value": float(problem.value)})
        if problem.status == "optimal_inaccurate":
            problem.solve(solver="SCS", eps=1e-9, max_iters=100000)
            attempts.append({"problem": name, "solver": "SCS", "status": problem.status,
                             "value": float(problem.value)})
        if problem.status != "optimal":
            raise ArithmeticError(f"Uncertified SDP status for {name}, d={d}, M={M}: {problem.status}")
    gaps = [classical_dual.value - classical.value, quantum_dual.value - quantum.value]
    residuals = {
        "classical_completeness": float(np.max(np.abs(sum(E.value for E in effects) - np.eye(d)))),
        "quantum_trace_preservation": float(np.max(np.abs(partial.value - np.eye(d)))),
        "classical_min_eigenvalue": float(min(np.linalg.eigvalsh(E.value).min() for E in effects)),
        "quantum_min_eigenvalue": float(np.linalg.eigvalsh(J.value).min()),
        "classical_dual_min_eigenvalue": float(min(np.linalg.eigvalsh(Y.value - rho / M).min() for rho in states)),
        "quantum_dual_min_eigenvalue": float(np.linalg.eigvalsh(np.kron(np.eye(M), Z.value) - C).min()),
    }
    assert max(abs(x) for x in gaps) < 2e-7
    assert residuals["classical_completeness"] < 2e-7
    assert residuals["quantum_trace_preservation"] < 2e-7
    assert min(value for name, value in residuals.items() if "eigenvalue" in name) > -2e-7
    assert quantum.value <= classical.value + 2e-7
    return {"classical": float(classical.value), "quantum": float(quantum.value),
            "dual_gaps": [float(x) for x in gaps], "residuals": residuals,
            "statuses": [problem.status for problem in problems], "solver_attempts": attempts}


def recovery_check(V, matrices, probabilities, head):
    d, M = V.shape
    blocks, _ = output_states(V, matrices, probabilities)
    S = sum(mm(blocks[q], blocks[q].conj().T) for q in head)
    eigenvalues, eigenvectors = np.linalg.eigh(S)
    factors = np.zeros(d)
    support = eigenvalues > 1e-10
    factors[support] = 1 / np.sqrt(eigenvalues[support])
    inverse = np.einsum("ik,k,jk->ij", eigenvectors, factors, eigenvectors.conj(), optimize=False)
    recoveries = [mm(blocks[q].conj().T, inverse) for q in head]
    # Include a CPTP completion on the kernel; it can only increase fidelity.
    for vector in eigenvectors[:, ~support].T:
        R = np.zeros((M, d), complex)
        R[0] = vector.conj()
        recoveries.append(R)
    completeness = sum(mm(R.conj().T, R) for R in recoveries)
    assert np.max(np.abs(completeness - np.eye(d))) < 2e-8
    fidelity = sum(p * sum(abs(np.trace(mm(R, B))) ** 2 for R in recoveries)
                   for p, B in zip(probabilities, blocks)) / M ** 2
    vectors = np.column_stack([blocks[q][:, j] for q in head for j in range(M)])
    gram = mm(vectors.conj().T, vectors)
    cross = np.abs(gram) ** 2
    np.fill_diagonal(cross, 0)
    loss = sum(probabilities[q] * cross[i * M:(i + 1) * M].sum() / M for i, q in enumerate(head))
    lower = probabilities[head].sum() - 8 * loss
    assert fidelity + 2e-8 >= lower
    return {"head_size": len(head), "fidelity": float(fidelity), "pointwise_lower": float(lower),
            "trace_preservation_residual": float(np.max(np.abs(completeness - np.eye(d))))}


def validate():
    import cvxpy
    import clarabel
    import scs
    import qiskit
    from qiskit.quantum_info import Clifford, random_clifford
    rng = np.random.default_rng(SEED)
    cases = []
    largest_fourier_error = 0.0
    for n in (1, 2, 3):
        for k in range(1, n + 1):
            c = random_clifford(n, seed=int(rng.integers(2 ** 31)))
            weights = retained_fourier_weights(c.adjoint().symplectic_matrix, n, k)
            for r in (0.0, 0.13, 0.4, 1.0):
                exact = direct_clifford_distribution(c, [1 - 3 * r / 4, r / 4, r / 4, r / 4], k)
                fast = clifford_distribution(weights, n, k, r)
                error = float(np.max(np.abs(exact - fast)))
                largest_fourier_error = max(error, largest_fourier_error)
                assert error < 2e-12
            identity = Clifford(np.eye(2 * n, dtype=bool))
            iw = retained_fourier_weights(identity.adjoint().symplectic_matrix, n, k)
            pc, fq = optimal_clifford_scores(clifford_distribution(iw, n, k, 0.2))
            assert abs(pc - 0.9 ** k) < 1e-12
            assert abs(fq - 0.85 ** k) < 1e-12
    for n, k, ensemble in ((2, 1, "clifford"), (2, 2, "clifford"), (2, 1, "haar"), (3, 1, "haar")):
        matrices, words = pauli_matrices(n)
        for trial, r in enumerate((0.0, 0.2, 0.6, 1.0)):
            seed = int(rng.integers(2 ** 31))
            if ensemble == "clifford":
                c = random_clifford(n, seed=seed)
                V = c.to_matrix()[:, :1 << k]
            else:
                V = unitary_group.rvs(1 << n, random_state=np.random.default_rng(seed))[:, :1 << k]
            p1 = [1 - 3 * r / 4, r / 4, r / 4, r / 4]
            probabilities = np.array([math.prod(p1[x] for x in word) for word in words])
            result = solve_optimal_scores(V, matrices, probabilities)
            if ensemble == "clifford":
                dist = direct_clifford_distribution(c, p1, k)
                pc, fq = optimal_clifford_scores(dist)
                assert abs(pc - result["classical"]) < 2e-7
                assert abs(fq - result["quantum"]) < 2e-7
                result["exact_classical"], result["exact_quantum"] = pc, fq
            if r == 0:
                assert abs(result["classical"] - 1) < 2e-7
                assert abs(result["quantum"] - 1) < 2e-7
            if r == 1:
                assert abs(result["classical"] - 2.0 ** (-k)) < 2e-7
                assert abs(result["quantum"] - 2.0 ** (-2 * k)) < 2e-7
            checks = []
            for L in (1, 2, 4):
                head = np.argsort(-probabilities, kind="stable")[:L]
                check = recovery_check(V, matrices, probabilities, head)
                assert check["fidelity"] <= result["quantum"] + 2e-7
                checks.append(check)
            cases.append({"n": n, "k": k, "ensemble": ensemble, "seed": seed, "r": r,
                          "sdp": result, "explicit_recovery": checks})
    cdf_error = 0.0
    for n in (1, 2, 4):
        _, words = pauli_matrices(n)
        for r in (0.07, 0.2, 0.6):
            p1 = [1 - 3 * r / 4, r / 4, r / 4, r / 4]
            p = np.array([math.prod(p1[x] for x in word) for word in words])
            s = -np.log2(p)
            for threshold in np.linspace(0.123, 4 * n + 0.123, 31):
                error = abs(float(p[s < threshold].sum()) - spectrum_cdf(n, r, threshold))
                cdf_error = max(cdf_error, error)
                assert error < 1e-12
    lattice_checks = []
    for n in (1, 2, 5, 12):
        # At r=1/2 every nonidentity symbol has exact surprisal 3.
        boundary = float(3 * n)
        for threshold, expected in (
            (math.nextafter(boundary, -math.inf), 1 - (3 / 8) ** n),
            (boundary, 1 - (3 / 8) ** n),
            (math.nextafter(boundary, math.inf), 1.0),
        ):
            actual = spectrum_cdf(n, 0.5, threshold)
            assert abs(actual - expected) < 1e-12
            lattice_checks.append({"n": n, "r": 0.5, "threshold": threshold,
                                   "expected": expected, "actual": actual})
    return {
        "seed": SEED, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                     "qiskit": qiskit.__version__, "cvxpy": cvxpy.__version__,
                     "clarabel": clarabel.__version__, "scs": scs.__version__},
        "fourier_enumeration_max_error": largest_fourier_error, "cdf_max_error": cdf_error,
        "cdf_lattice_boundary_checks": lattice_checks,
        "sdp_cases": cases, "sdp_case_count": len(cases),
        "scope": "Small-system consistency tests; analytic bounds at large n, not a large-n quantum simulation.",
    }


def simulate(samples=SAMPLES):
    from qiskit.quantum_info import random_clifford
    rng = np.random.default_rng(SEED + 1)
    rows, encoders = [], []
    for n, k in SIZES:
        for sample in range(samples):
            seed = int(rng.integers(2 ** 31))
            c = random_clifford(n, seed=seed)
            inverse = c.adjoint().symplectic_matrix.astype(np.uint8)
            weights = retained_fourier_weights(inverse, n, k)
            encoders.append({"n": n, "k": k, "sample": sample, "seed": seed,
                             "inverse_symplectic": inverse.tolist()})
            for r in R_VALUES:
                pc, fq = optimal_clifford_scores(clifford_distribution(weights, n, k, float(r)))
                rows.append({"n": n, "k": k, "sample": sample, "seed": seed, "r": float(r),
                             "classical": pc, "quantum": fq})
            if (sample + 1) % 8 == 0:
                print(f"n={n}, k={k}: {sample + 1}/{samples} exact Clifford encoders", flush=True)
    write_csv(ROOT / "data" / "recovery_clifford_raw.csv", rows)
    (ROOT / "data" / "recovery_encoders.json").write_text(json.dumps(encoders, separators=(",", ":")) + "\n")
    return rows


def write_csv(path, rows):
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_raw():
    with (ROOT / "data" / "recovery_clifford_raw.csv").open() as file:
        return [{key: int(value) if key in ("n", "k", "sample", "seed") else float(value)
                 for key, value in row.items()} for row in csv.DictReader(file)]


def summarize(rows):
    result = []
    for n, k in SIZES:
        for r in R_VALUES:
            selected = [row for row in rows if row["n"] == n and row["k"] == k and abs(row["r"] - r) < 1e-10]
            row = {"n": n, "k": k, "r": float(r), "samples": len(selected)}
            for task in ("classical", "quantum"):
                values = np.array([item[task] for item in selected])
                row[task] = float(values.mean())
                row[task + "_se"] = float(values.std(ddof=1) / math.sqrt(len(values)))
            gaps = np.array([item["classical"] - item["quantum"] for item in selected])
            row["gap"] = float(gaps.mean())
            row["gap_se"] = float(gaps.std(ddof=1) / math.sqrt(len(gaps)))
            result.append(row)
    write_csv(ROOT / "data" / "recovery_clifford_summary.csv", result)
    return result


def critical_data():
    rows = []
    r = 0.1
    h, sigma = entropy_sigma(r)
    for n in (4096, 65536, 1048576):
        for z in np.linspace(-2.5, 2.5, 101):
            k = round(n * (1 - h) + math.sqrt(n) * sigma * z)
            actual_z = (k - n * (1 - h)) / (math.sqrt(n) * sigma)
            lower, upper = spectrum_bracket(n, k, r)
            rows.append({"n": n, "k": k, "r": r, "z": actual_z, "lower": lower,
                         "upper": upper, "normal": float(norm.cdf(-actual_z))})
    write_csv(ROOT / "data" / "recovery_critical_brackets.csv", rows)
    return rows


def plot(summary, critical):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "text.color": INK,
                         "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})
    fig, axes = plt.subplots(1, 3, figsize=(7.05, 2.55), layout="constrained")
    threshold = brentq(lambda r: entropy_sigma(r)[0] - 0.75, 0.01, 0.4)
    ax = axes[0]
    for n, marker, alpha in ((4, "o", 0.45), (16, "s", 1.0)):
        data = [row for row in summary if row["n"] == n and row["r"] <= 0.55]
        x = [row["r"] for row in data]
        for task, color, style, label in (("classical", BLUE, "-", "Classical"), ("quantum", ORANGE, "--", "Quantum")):
            ax.errorbar(x, [row[task] for row in data], yerr=[row[task + "_se"] for row in data],
                        color=color, linestyle=style, marker=marker, markersize=3.5, linewidth=1.1,
                        markeredgecolor="white", markeredgewidth=.4, alpha=alpha,
                        label=f"{label}, $n={n}$", capsize=1.5)
    ax.axvline(threshold, color="#888888", linestyle=":", linewidth=.8)
    ax.set(xlabel="Depolarizing strength $r$", ylabel="Optimal performance", ylim=(-.02, 1.04), xlim=(0, .56))
    ax.legend(fontsize=6.3, loc="upper right", frameon=False, handlelength=1.5)
    ax.text(.04, .06, "$k/n=1/4$", transform=ax.transAxes, fontsize=7)
    ax.set_title("(a) Exact finite-size decoding", loc="left", fontsize=9)
    ax = axes[1]
    for n, style in ((4096, "--"), (1048576, "-")):
        data = [row for row in critical if row["n"] == n]
        x = np.array([row["z"] for row in data])
        lo, hi = [row["lower"] for row in data], [row["upper"] for row in data]
        ax.fill_between(x, lo, hi, color=BLUE, alpha=.07 if n == 4096 else .18)
        ax.plot(x, lo, color=BLUE, linestyle=style, linewidth=.85, label=f"Bounds, $n=2^{{{int(math.log2(n))}}}$")
        ax.plot(x, hi, color=BLUE, linestyle=style, linewidth=.85)
    x = np.linspace(-2.5, 2.5, 200)
    ax.plot(x, norm.cdf(-x), color=INK, linewidth=1.1, label="$\\Phi(-z)$")
    ax.set(xlabel="Scaled message number $z$", ylabel="Common performance bracket", ylim=(-.02, 1.04), xlim=(-2.5, 2.5))
    ax.legend(fontsize=6.3, loc="lower left", frameon=False)
    ax.set_title("(b) Rigorous critical bracket", loc="left", fontsize=9)
    ax = axes[2]
    rgrid = np.linspace(.0001, .6, 300)
    random_rate = np.array([max(1 - entropy_sigma(float(r))[0], 0) for r in rgrid])
    binary = lambda x: -x * math.log2(x) - (1 - x) * math.log2(1 - x)
    capacity = np.array([1 - binary(float(r) / 2) for r in rgrid])
    ax.plot(rgrid, capacity, color=BLUE, linewidth=1.4, label="Optimized classical code")
    ax.plot(rgrid, random_rate, color=ORANGE, linestyle="--", linewidth=1.4, label="Random code: both tasks")
    ax.set(xlabel="Depolarizing strength $r$", ylabel="Asymptotic rate (bits/qubit)", xlim=(0, .6), ylim=(-.02, 1.02))
    ax.legend(fontsize=6.2, loc="upper right", frameon=False)
    ax.set_title("(c) Cost of random encoding", loc="left", fontsize=9)
    fig.savefig(ROOT / "figures" / "recovery.pdf")
    fig.savefig(ROOT / "figures" / "recovery.png", dpi=200)
    plt.close(fig)
    rate_rows = [{"r": float(r), "classical_capacity": float(c), "random_rate": float(q)}
                 for r, c, q in zip(rgrid, capacity, random_rate)]
    write_csv(ROOT / "data" / "recovery_rates.csv", rate_rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulate", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--samples", type=int, default=SAMPLES)
    args = parser.parse_args()
    if args.check_only:
        print(json.dumps(validate(), indent=2))
        return
    rows = simulate(args.samples) if args.simulate else read_raw()
    summary = summarize(rows)
    critical = critical_data()
    plot(summary, critical)
    report = validate()
    report["samples_per_size"] = args.samples if args.simulate else len({row["sample"] for row in rows if row["n"] == SIZES[0][0]})
    report["raw_rows"] = len(rows)
    report["data_sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in sorted((ROOT / "data").glob("recovery_*"))
                             if p.suffix in (".csv", ".json") and p.name != "recovery_verification.json"}
    assert set(report["data_sha256"]) == {
        "recovery_clifford_raw.csv", "recovery_clifford_summary.csv", "recovery_encoders.json",
        "recovery_critical_brackets.csv", "recovery_rates.csv",
    }
    (ROOT / "data" / "recovery_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"raw_rows": len(rows), "sdp_cases": report["sdp_case_count"],
                      "fourier_error": report["fourier_enumeration_max_error"]}, indent=2))


if __name__ == "__main__":
    main()
