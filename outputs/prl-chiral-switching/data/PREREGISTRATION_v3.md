# Pre-registration v3: coupling exponent = marginal quasipotential at the fold population

Written 2026-10-28, before any computation at Ω = 2.5 or Ω = 3.0. The hypothesis comes from exploratory data at Ω = 2 (`code/xel.py`, `prof2.py`, `prof3.py`; `data/phiEP.json`), which **does not count as evidence** for this test.

## Hypothesis (H3)
At the avoided crossing Δ_min(N) between the high-branch intrabasin mode and the ghost-escape mode, the half gap g_N satisfies

  ln g_N = c − N · Φ_N(n_f; Δ_min(N)) + o(N),

where
- n_f(Ω) is the MF population at the left fold (the double root of the cubic), fixed by Ω with no freedom;
- Φ_N(n; Δ) = −(1/N)[ln P_ss(⌊nN⌉) − max_k ln P_ss(k)] is the marginal quasipotential from the exact steady-state excitation-number distribution at the same N and Δ, linearly interpolated in n.

Mechanism: the right eigenvector of the intrabasin mode carries the high-basin WKB weight exp(−NΦ); the ghost mode is O(1) at n_f; the coupling is their overlap at the ghost core.

## Inputs at the new parameters
- Ω ∈ {2.5, 3.0}, V = 10, γ = 1; N ∈ {24, 32, 40, 48}.
- Δ_min(N) and g_N = (minimum real gap)/2 between the two real eigenvalues adjacent to the MF high-branch Jacobian real eigenvalue at the fold. These come from a single job; the g_N values are written to file and **not inspected** until the Φ values are computed.

## Test (per Ω)
1. Slope of the linear fit of ln g_N against N·Φ_N over the 4 values of N: PASS if within 1 ± 0.15.
2. Local exponent at the middle N (32 → 40): −Δ ln g/ΔN vs Δ(NΦ_N)/ΔN, PASS if they agree within 15%.
H3 passes only if both criteria pass at both Ω values.

## Not allowed
Changing n_f, the definition of Φ, the N set, or the fit after seeing g_N. If the avoided crossing cannot be located for some N, that N is dropped, and this is reported.
- Setup note: at Ω=2.5 the fold-finding bracket [3,5] crossed both folds, so the Ω=2.5 runs crashed. The bracket upper end is changed to 4.4 for Ω=2.5. Nothing else changes.

## Result (g_N unsealed after `v3_prediction_frozen.json` was written)
| Ω | slope of ln g vs NΦ (1±0.15) | local exponent 32→40: observed / predicted (±15%) |
|---|---|---|
| 2.5 | 0.957 PASS | 0.1408 / 0.1492 = 0.944 PASS |
| 3.0 | 0.878 PASS | 0.1187 / 0.1352 = 0.878 PASS |
**H3: PASS** under the pre-registered criteria. Robustness checks (null n-choice and grid integrity) follow below; they are not part of the registered test.

## Robustness and caveats (post-test; not part of the registered criteria)
- **Sensitivity to the evaluation point.** Scanning n in place of n_f: the fitted slope changes by about 0.08 per 0.01 in n (steep). The slope would equal 1 at n ≈ n_f + 0.005 (Ω=2.5) and n_f + 0.014 (Ω=3.0), by linear interpolation of the scan. The pre-registered n_f lands within ~0.02 of the best point at both Ω, with no tuning. Points 0.04 away give slopes 0.65–1.29 and fail.
- **Systematic.** At both Ω the prediction slightly overshoots the local exponent (ratios 0.94, 0.88). Combined with the near-unit slopes (0.96, 0.88), this suggests a sub-leading prefactor or an N-dependent shift of the effective ghost point toward larger n, consistent with the n_f + 0.005–0.014 offset.
- **What is and is not shown.** Shown: the exponential N-dependence of the EP/avoided-crossing coupling is quantitatively set by the exact steady-state marginal quasipotential at the MF fold population, at two parameter values not used to form the hypothesis. Not shown: an a-priori computation of Φ(n_f) from the WKB Hamiltonian. Φ_N here is still read from the exact steady state. The validated Lin solver could supply it in principle; this has not been done.
