#!/usr/bin/env python3
"""Reproduce the revised manuscript's finite-size evidence; no surrogate model."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import binom, norm

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261027
PAIRS = 64
N_SIM = 8
R_VALUES = np.round(np.linspace(0.05, 0.75, 15), 10)
BLUE = "#2a78d6"
ORANGE = "#eb6834"
INK = "#0b0b0b"
PAULIS = (
    np.eye(2, dtype=complex),
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.diag([1, -1]).astype(complex),
)


def entropy(probabilities):
    p = np.asarray(probabilities, dtype=float)
    p = p[p > 0]
    return float(-np.dot(p, np.log2(p)))


def depolarizing_distribution(r):
    return np.array([1 - 3 * r / 4, r / 4, r / 4, r / 4])


def entropy_sigma(r):
    p = depolarizing_distribution(r)
    p = p[p > 0]
    surprisal = -np.log2(p)
    h = float(p @ surprisal)
    variance = float(p @ (surprisal - h) ** 2)
    return h, math.sqrt(variance)


def critical_r():
    return brentq(lambda r: entropy_sigma(r)[0] - 1, 0.05, 0.5, xtol=1e-14)


def haar_state(d, rng):
    vector = rng.normal(size=d) + 1j * rng.normal(size=d)
    return vector / np.linalg.norm(vector)


def depolarize(matrix, n, r):
    """Apply independent replacement noise, linearly to any matrix."""
    if matrix.shape != (2**n, 2**n) or not 0 <= r <= 1:
        raise ValueError("Inconsistent dimension or noise strength")
    result = matrix.copy()
    for k in range(n):
        tensor = result.reshape(2**k, 2, 2 ** (n-k-1), 2**k, 2, 2 ** (n-k-1))
        partial_trace = tensor[:, 0, :, :, 0, :] + tensor[:, 1, :, :, 1, :]
        output = (1-r) * tensor
        output[:, 0, :, :, 0, :] += r * partial_trace / 2
        output[:, 1, :, :, 1, :] += r * partial_trace / 2
        result = output.reshape(2**n, 2**n)
    return result


def trace_distance(first, second):
    return float(np.abs(np.linalg.eigvalsh(first-second)).sum() / 2)


def cutoff_bounds(n, r):
    """Projector lower bound and exact-second-moment two-class upper bound.

    Optimizes over all head cutoffs k<=K, including the empty head.
    Intended for small n; large n uses the binomial information-spectrum bound.
    """
    if n > 32:
        raise ValueError("Use spectrum_bounds for large systems")
    d = 2.0**n
    q = 3*r/4
    count = np.array([math.comb(n, k)*3**k for k in range(n+1)], dtype=float)
    weights = np.array([(1-q)**(n-k)*(q/3)**k for k in range(n+1)])
    mass = count*weights
    collision = count*weights**2
    head_mass = np.r_[0.0, np.cumsum(mass)]
    head_count = np.r_[0.0, np.cumsum(count)]
    tail_mass = np.r_[np.cumsum(mass[::-1])[::-1], 0.0]
    tail_collision = np.r_[np.cumsum(collision[::-1])[::-1], 0.0]
    lower = max(0.0, float(np.max(head_mass-head_count/d)))
    square = (d*d*tail_collision-tail_mass**2)/(2*(d+1))
    upper = min(1.0, float(np.min(head_mass+np.sqrt(np.maximum(0.0, square)))))
    return lower, upper


def surprisal_cdf(n, r, threshold):
    """P[s<threshold], with strict inequality also at lattice points."""
    q = 3*r/4
    if r == 0:
        return float(0 < threshold)
    if r == 1:
        return float(2*n < threshold)
    s0 = -math.log2(1-q)
    step = math.log2(3*(1-q)/q)
    k = int(math.ceil((threshold-n*s0)/step)-1)
    k = max(-1, min(n, k))
    while k >= 0 and n*s0+k*step >= threshold:
        k -= 1
    while k < n and n*s0+(k+1)*step < threshold:
        k += 1
    return float(binom.cdf(k, n, q))


def spectrum_bounds(n, r, a):
    lower = max(0.0, surprisal_cdf(n, r, n-a)-2.0**(-a))
    upper = min(1.0, surprisal_cdf(n, r, n+a)+2.0**(-a/2)/math.sqrt(2))
    return lower, upper


def r_at_scaled_coordinate(n, x):
    return brentq(
        lambda r: math.sqrt(n)*(1-entropy_sigma(r)[0])/entropy_sigma(r)[1]-x,
        0.02, 0.8, xtol=1e-13,
    )


def mub_distance(n):
    d = 2**n
    return (d-1)*(2*d+1)/(d*(d+1)**2)


def haar_dephasing_distance(n):
    d = 2**n
    return (d-1)/(2*d-1)


def pauli_strings(n):
    matrices = []
    for word in itertools.product(range(4), repeat=n):
        matrix = np.ones((1, 1), dtype=complex)
        for symbol in word:
            matrix = np.kron(matrix, PAULIS[symbol])
        matrices.append(matrix)
    return matrices


def complete_mub_states(n):
    """Explicit complete MUB designs for n=1,2, for independent self-checks."""
    if n == 1:
        return [v for matrix in PAULIS[1:] for v in np.linalg.eigh(matrix)[1].T]
    if n == 2:
        # Each pair commutes and generates one of five disjoint Pauli classes.
        classes = [((3, 0), (0, 3)), ((1, 0), (0, 1)), ((2, 0), (0, 2)),
                   ((1, 2), (2, 3)), ((1, 3), (2, 1))]
        states = []
        for a, b in classes:
            observable = np.kron(PAULIS[a[0]], PAULIS[a[1]])
            observable += 2*np.kron(PAULIS[b[0]], PAULIS[b[1]])
            states.extend(np.linalg.eigh(observable)[1].T)
        return states
    raise ValueError("The explicit check is only implemented for n=1,2")


def mean_pair_distance(outputs):
    return float(np.mean([trace_distance(a, b) for a in outputs for b in outputs]))


def verification():
    """Meaningful identities and independent implementations, not fit checks."""
    rng = np.random.default_rng(SEED+1)
    report = {"channel_errors": {}, "mub_checks": {}, "correlated_pauli": {}}
    for n in (1, 2, 3):
        d = 2**n
        a, b = haar_state(d, rng), haar_state(d, rng)
        delta = np.outer(a, a.conj())-np.outer(b, b.conj())
        matrices = pauli_strings(n)
        for r in (0.0, 0.37, 1.0):
            p = depolarizing_distribution(r)
            probabilities = [math.prod(p[i] for i in word)
                             for word in itertools.product(range(4), repeat=n)]
            direct = sum(w*(P@delta@P.conj().T) for w, P in zip(probabilities, matrices))
            efficient = depolarize(delta, n, r)
            error = float(np.max(np.abs(direct-efficient)))
            assert error < 2e-13
            report["channel_errors"][f"n{n}_r{r}"] = error
            if r == 1:
                assert np.linalg.norm(efficient) < 2e-13
    for n in (1, 2):
        d = 2**n
        states = complete_mub_states(n)
        densities = [np.outer(v, v.conj()) for v in states]
        moment = sum(np.kron(rho, rho) for rho in densities)/len(densities)
        swap = np.zeros((d*d, d*d))
        for i in range(d):
            for j in range(d):
                swap[i*d+j, j*d+i] = 1
        target = (np.eye(d*d)+swap)/(d*(d+1))
        moment_error = float(np.max(np.abs(moment-target)))
        assert moment_error < 2e-13
        outputs = [np.diag(np.diag(rho)) for rho in densities]
        exact = mean_pair_distance(outputs)
        assert abs(exact-mub_distance(n)) < 2e-13
        report["mub_checks"][str(n)] = {"design_error": moment_error, "distance": exact}
    # A genuinely correlated two-qubit Pauli distribution, not a product law.
    matrices = pauli_strings(2)
    p = 0.45*rng.dirichlet(np.ones(16))
    p[0] += 0.55
    densities = [np.outer(v, v.conj()) for v in complete_mub_states(2)]
    outputs = [sum(w*(P@rho@P.conj().T) for w, P in zip(p, matrices)) for rho in densities]
    exact_distance = mean_pair_distance(outputs)
    order = np.argsort(-p)
    subsets = [order[:4], order[4:10], order[10:]]
    upper = 0.0
    for subset in subsets:
        m, v = float(p[subset].sum()), float((p[subset]**2).sum())
        partial = [sum(p[j]*(matrices[j]@rho@matrices[j].conj().T) for j in subset)
                   for rho in densities]
        moment = float(np.mean([np.trace((a-b)@(a-b)).real for a in partial for b in partial]))
        predicted = 8*v/5-2*m*m/20
        assert abs(moment-predicted) < 2e-13
        upper += min(m, math.sqrt(max(0.0, (16*v-m*m)/10)))
    lower = max([0.0]+[float(p[order[:k]].sum())-k/4 for k in range(1, 17)])
    assert lower > 0.25
    assert lower <= exact_distance+2e-13 <= upper+4e-13
    for a in np.linspace(0, 8, 33):
        s = -np.log2(p)
        lo = p[s < 2-a].sum()-2**(-a)
        hi = p[s < 2+a].sum()+2**(-a/2)/math.sqrt(2)
        assert lo <= exact_distance+2e-13 and exact_distance <= hi+2e-13
    report["correlated_pauli"] = {"distance": exact_distance, "lower": lower, "upper": upper}
    # Check binomial CDF evaluation independently by enumerating error weights.
    for n in (1, 3, 8):
        for r in (0.1, 0.4, 0.8):
            q = 3*r/4
            for a in (0., 1.5, 3.):
                threshold = n+a
                direct = sum(math.comb(n,k)*q**k*(1-q)**(n-k)
                             for k in range(n+1)
                             if -(n-k)*math.log2(1-q)-k*math.log2(q/3) < threshold)
                assert abs(direct-surprisal_cdf(n,r,threshold)) < 2e-13
    report["critical_r"] = critical_r()
    rate_stab = 1-math.log2(1.3)
    rate_bound = (1-math.log2(1.03))/2
    assert rate_stab > rate_bound
    report["stabilizer_rate_example"] = {"lower_bound": rate_stab, "privacy_amplification_bound": rate_bound}
    # Independent one-dimensional minimization of the depolarizing type formula.
    qcrit = 3*critical_r()/4
    errors = []
    for r in (0.3, 0.5, 0.7, 0.9):
        p = depolarizing_distribution(r)
        def renyi_exponent(s):
            if s == 0:
                return 0.0
            h = -math.log2(float(np.sum(p**(1+s))))/s
            return s/(1+s)*(h-1)
        optimum = minimize_scalar(lambda s:-renyi_exponent(s), bounds=(1e-8,1), method="bounded")
        value = max(0.0, renyi_exponent(1), -optimum.fun)
        q = 3*r/4
        escort = (q*q/3)/((1-q)**2+q*q/3)
        expected = (qcrit*math.log2(qcrit/q)+(1-qcrit)*math.log2((1-qcrit)/(1-q))
                    if escort < qcrit else (-math.log2((1-q)**2+q*q/3)-1)/2)
        errors.append(abs(value-expected))
    assert max(errors) < 1e-8
    report["renyi_identity_max_error"] = max(errors)
    # The projector exponent and the achievable exponent agree to leading order
    # as r decreases to r_c from above; neither proves a general rate limit.
    kappa = 9/(32*math.log(2)*qcrit*(1-qcrit))
    near = []
    for dr in (0.01, 0.002, 0.0004):
        q = 3*(critical_r()+dr)/4
        s0 = -math.log2(1-q)
        step = math.log2(3*(1-q)/q)
        f0 = (1-s0)/step
        def binary_relative_entropy(f):
            return f*math.log2(f/q)+(1-f)*math.log2((1-f)/(1-q))
        achievable = binary_relative_entropy(qcrit)
        projector = binary_relative_entropy(f0)
        assert 0 < achievable <= projector
        near.append({"delta_r":dr,"achievable_over_delta2":achievable/dr**2,
                     "projector_over_delta2":projector/dr**2})
    assert abs(near[-1]["achievable_over_delta2"]-kappa) < 0.01
    assert abs(near[-1]["projector_over_delta2"]-kappa) < 0.01
    report["near_threshold"] = {"kappa":kappa,"values":near}
    report["checks_passed"] = True
    return report


def simulate():
    rng = np.random.default_rng(SEED)
    rows = []
    for pair in range(PAIRS):
        a, b = haar_state(2**N_SIM, rng), haar_state(2**N_SIM, rng)
        delta = np.outer(a, a.conj())-np.outer(b, b.conj())
        for r in R_VALUES:
            output = depolarize(delta, N_SIM, float(r))
            value = float(np.abs(np.linalg.eigvalsh(output)).sum()/2)
            rows.append((pair, r, value))
        if (pair+1) % 16 == 0:
            print(f"Simulated {pair+1}/{PAIRS} Haar pairs", flush=True)
    return np.array(rows)


def figure_data(raw):
    small = []
    for r in R_VALUES:
        values = raw[np.isclose(raw[:,1],r,rtol=0,atol=1e-12),2]
        if len(values) != PAIRS:
            raise ValueError("Raw sample count does not match the fixed design")
        lo,hi = cutoff_bounds(N_SIM,float(r))
        small.append((r,float(values.mean()),float(values.std(ddof=1)/math.sqrt(PAIRS)),lo,hi))
    small = np.array(small)
    large = []
    for n in (1024,65536):
        for x in np.linspace(-3,3,121):
            r = r_at_scaled_coordinate(n,float(x))
            lo,hi = spectrum_bounds(n,r,2*math.log2(n))
            assert 0 <= lo <= hi <= 1
            large.append((n,x,r,lo,hi,float(norm.cdf(x))))
    degenerate = np.array([(n,mub_distance(n),haar_dephasing_distance(n)) for n in range(1,17)])
    return small,np.array(large),degenerate


def plot_figures(small,large,degenerate):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    plt.rcParams.update({
        "font.family":"DejaVu Sans", "font.size":9, "axes.labelsize":10,
        "axes.titlesize":10, "legend.fontsize":8, "lines.linewidth":1.5,
        "pdf.fonttype":42, "ps.fonttype":42, "text.color":INK,
        "axes.labelcolor":INK, "xtick.color":"#52514e", "ytick.color":"#52514e",
        "axes.spines.top":False, "axes.spines.right":False,
        "axes.edgecolor":"#c3c2b7", "axes.linewidth":0.6,
        "grid.color":"#e1e0d9", "grid.linewidth":0.5,
        "savefig.facecolor":"white", "figure.facecolor":"white",
    })
    fig,axes=plt.subplots(1,2,figsize=(7.1,2.75),layout="constrained")
    ax=axes[0]
    rs=np.linspace(.05,.75,241)
    bounds=np.array([cutoff_bounds(N_SIM,float(r)) for r in rs])
    ax.fill_between(rs,bounds[:,0],bounds[:,1],color=BLUE,alpha=.12)
    ax.plot(rs,bounds[:,0],color=BLUE)
    ax.plot(rs,bounds[:,1],color=BLUE,ls="--")
    ax.errorbar(small[:,0],small[:,1],small[:,2],fmt="o",markersize=6,
                color=INK,mec="white",mew=.7,capsize=2,elinewidth=.9,zorder=4)
    rc=critical_r()
    ax.axvline(rc,color="#898781",lw=.9,ls=":")
    ax.text(rc+.009,.055,r"$r_c$",fontsize=9)
    ax.annotate("lower",xy=(.40,cutoff_bounds(N_SIM,.40)[0]),xytext=(.41,.045),
                fontsize=8,arrowprops={"arrowstyle":"-","color":"#898781","lw":.7})
    ax.annotate("upper",xy=(.60,cutoff_bounds(N_SIM,.60)[1]),xytext=(.62,.36),
                fontsize=8,arrowprops={"arrowstyle":"-","color":"#898781","lw":.7})
    ax.legend(handles=[Patch(facecolor=BLUE,alpha=.15,label="rigorous bounds"),
                       Line2D([0],[0],marker="o",color=INK,ls="none",label="Haar sample mean")],
              frameon=False,loc="upper right")
    ax.set(xlim=(.04,.76),ylim=(0,1.03),xlabel="Depolarizing strength $r$",ylabel=r"$\overline{T}$")
    ax.set_title("(a) Finite system: $n=8$",loc="left",pad=9)
    ax=axes[1]
    for n,color,ls in [(1024,BLUE,"--"),(65536,ORANGE,"-")]:
        table=large[large[:,0]==n]
        ax.fill_between(table[:,1],table[:,3],table[:,4],color=color,alpha=.1)
        ax.plot(table[:,1],table[:,3],color=color,ls=ls)
        ax.plot(table[:,1],table[:,4],color=color,ls=ls)
    xx=np.linspace(-3,3,301)
    ax.plot(xx,norm.cdf(xx),color=INK,lw=1.4)
    ax.legend(handles=[Line2D([0],[0],color=BLUE,ls="--",label="$n=1{,}024$ bounds"),
                       Line2D([0],[0],color=ORANGE,label="$n=65{,}536$ bounds"),
                       Line2D([0],[0],color=INK,label=r"$\Phi(x)$")],
              frameon=False,loc="upper left")
    ax.set(xlim=(-3,3),ylim=(0,1.03),xlabel=r"$x=\sqrt{n}\,[1-H]/\sigma$",ylabel=r"bounds on $\overline{T}$")
    ax.set_title("(b) Regular critical window",loc="left",pad=9)
    for ax in axes:
        ax.yaxis.grid(True)
        ax.set_axisbelow(True)
    out=ROOT/"figures"; out.mkdir(exist_ok=True)
    fig.savefig(out/"bounds.pdf",metadata={"Title":"Finite-size bounds on noisy-state distinguishability","CreationDate":None})
    fig.savefig(out/"bounds.png",dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(3.4,2.4),layout="constrained")
    ax.plot(degenerate[:,0],degenerate[:,1],"o-",color=BLUE,ms=6,mec="white",mew=.7,label="complete MUB design")
    ax.plot(degenerate[:,0],degenerate[:,2],"s--",color=ORANGE,ms=6,mec="white",mew=.7,label="Haar ensemble")
    ax.axhline(.5,color="#898781",ls=":",lw=.8)
    ax.text(10.5,.055,"MUB",fontsize=8)
    ax.text(10.5,.43,"Haar",fontsize=8)
    ax.set(xlim=(.5,16.5),ylim=(-.015,.53),xticks=[1,4,8,12,16],xlabel="Qubits $n$",ylabel=r"$\overline{T}$ after dephasing")
    ax.legend(frameon=False,loc="center right")
    ax.yaxis.grid(True); ax.set_axisbelow(True)
    fig.savefig(out/"dephasing.pdf",metadata={"Title":"Design dependence at complete dephasing","CreationDate":None})
    fig.savefig(out/"dephasing.png",dpi=180)
    plt.close(fig)
    return matplotlib.__version__


def save_table(path,table,header):
    np.savetxt(path,table,delimiter=",",header=header,comments="",fmt="%.17g")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resimulate",action="store_true",help="Regenerate raw Haar samples with the fixed seed")
    parser.add_argument("--check-only",action="store_true",help="Run identity checks without writing files")
    args=parser.parse_args()
    checks=verification()
    print("Independent channel, 2-design, moment, bound and entropy checks passed.",flush=True)
    if args.check_only:
        print(json.dumps(checks,indent=2))
        return
    data=ROOT/"data"; data.mkdir(exist_ok=True)
    raw_path=data/"haar_n8_raw.csv"
    if args.resimulate or not raw_path.exists():
        raw=simulate()
        save_table(raw_path,raw,"pair,r,trace_distance")
    else:
        raw=np.loadtxt(raw_path,delimiter=",",skiprows=1)
        if raw.shape != (PAIRS*len(R_VALUES),3) or not np.all(np.isfinite(raw)):
            raise ValueError("Unexpected raw sample format")
    small,large,degenerate=figure_data(raw)
    # This is a statistical consistency check, not a bound on every realization.
    checks["haar_means_within_bounds_plus_5se"] = bool(np.all(small[:,1]+5*small[:,2]>=small[:,3]) and np.all(small[:,1]-5*small[:,2]<=small[:,4]))
    assert checks["haar_means_within_bounds_plus_5se"]
    save_table(data/"haar_n8_summary.csv",small,"r,mean,standard_error,lower,upper")
    save_table(data/"critical_brackets.csv",large,"n,x,r,lower,upper,normal_cdf")
    save_table(data/"dephasing_exact.csv",degenerate,"n,mub_distance,haar_distance")
    mpl_version=plot_figures(small,large,degenerate)
    checks["environment"]={"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"matplotlib":mpl_version}
    checks["sampling"]={"seed":SEED,"pairs":PAIRS,"qubits":N_SIM,"r_values":R_VALUES.tolist(),"same_pairs_across_r":True,"standard_error":"sample standard deviation (ddof=1) / sqrt(pairs)"}
    checks["evidence_scope"]="Exact finite-dimensional channel evaluation on sampled Haar inputs; large-n curves are rigorous bounds, not simulated distances."
    checks["raw_sha256"]=hashlib.sha256(raw_path.read_bytes()).hexdigest()
    checks["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (data/"verification.json").write_text(json.dumps(checks,indent=2,allow_nan=False)+"\n")
    print("Wrote two figures, four numeric tables, and data/verification.json.",flush=True)


if __name__ == "__main__":
    main()
