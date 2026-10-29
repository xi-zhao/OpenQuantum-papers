# Pre-registration v4: a-priori WKB prediction of the coupling exponent

Written 2026-10-28, before any computation at Ω = 2.25 or Ω = 2.75 (spectra, steady states, or WKB).

## Derived formula (fixed now)
In the 4-mode Liouville-boson representation, P(k) = Tr(Π_k ρ_N) is the Fock coefficient at occupation
x = (x_gg, x_ee, x_ge, x_eg) = (1−n, n, 0, 0) times the multinomial factor, with n = k/N. In Doi / coherent-state variables
the occupations are x = X̂∘ρ̂ (component-wise product of conjugate pairs). The WKB rate function for P(k) is therefore the
Cole–Hopf action to that occupation point:

  Φ_WKB(n) = S_Doi + (1−n) ln X̂_gg(end) + n ln X̂_ee(end),

where S_Doi = ∫ p·dq is taken along the zero-energy heteroclinic of the reduced Hamiltonian (`code/wkbr.py`). The heteroclinic leaves the high fixed point in its 3D unstable manifold and ends on {x(q,p) = (1−n, n, 0, 0)}. That endpoint set has two branches:
- (A) p_u = p_v = 0, v = −n_ρ/Ω, with n_ρ(1+(1−n_ρ)p_n) = n. Here the endpoint term reduces to KL(n‖n_ρ).
- (B) u = v = 0, h = 0.

Φ_WKB is the smallest positive action over all converged solutions on either branch. The exact steady state is not used anywhere.

## Gate A′ (validation at Ω = 2; must pass before Gate B is interpreted)
Points: (Δ, n) ∈ {(3.25, n_f = 0.15693), (3.7, 0.25), (3.7, 0.30), (3.7, 0.35)}.
Reference: exact marginal quasipotential, extrapolated as a + b/N from N = 30, 40, 50.
PASS if |Φ_WKB/Φ_ref − 1| ≤ 10% at at least 3 of the 4 points, with none worse than 20%.

## Gate B (prediction at unseen Ω ∈ {2.25, 2.75}, V = 10, N ∈ {24, 32, 40, 48})
- The only exact input is the crossing location Δ_min(N), found as in v3 (a spectral quantity, not the steady state). g_N, and the exact Φ_N for diagnostics, are written to sealed files that are not opened until `v4_prediction_frozen.json` exists.
- n_f(Ω) = MF left-fold population.
- Prediction: ln g_N = c − N·Φ_WKB(n_f; Δ_min(N)).
- PASS per Ω: the slope of ln g_N vs N·Φ_WKB is within 1 ± 0.15, **and** the local exponent (32→40) is within ±15% of [40Φ(Δ_40) − 32Φ(Δ_32)]/8.
- v4 passes only if both criteria pass at both Ω values.

## Known, declared systematics
Φ_WKB is the N→∞ rate function. Prefactor corrections of order (ln N)/N are not included.

## Not allowed after unsealing
Changing n_f, the endpoint branches, the selection rule (minimum positive action), the N set, or the fit.
- Solver note (before any Gate A′ number): integrations are now stopped when n leaves [−0.02, 1.02], i.e. when a ray becomes unphysical, and max_step is set to 0.5. One stiff backward ray had taken more than 40 s. The quantity, endpoints, selection rule and criteria are unchanged.
- Solver note: the search is reduced to 2 sections (0.5, 0.3), 3 seeds per branch, max_nfev=150, for compute reasons. This makes the minimum-over-solutions *less* exhaustive, and that is declared here before any number.
- Solver change (before any Gate A′ number): Lin section-matching is replaced by backward-to-A shooting from the endpoint family, minimising the closest approach to the high fixed point over the 2 endpoint parameters (codim-2 condition, well-posed). Converged means closest approach < 1e-5. The quantity, branches, selection rule and criteria are unchanged.

## Gate A′ result (recorded before any Gate B WKB run; v4_SEALED.jsonl unopened)
| (Δ, n) at Ω=2 | Φ_WKB (a priori) | Φ_exact (N→∞) | ratio |
|---|---|---|---|
| (3.25, 0.15693) | 0.13690 | 0.13717 | 0.998 |
| (3.7, 0.25) | 0.06224 | 0.06605 | 0.942 |
| (3.7, 0.30) | 0.02918 | 0.03103 | 0.940 |
| (3.7, 0.35) | 0.00882 | 0.00924 | 0.955 |
**PASS** (4/4 within 10%, worst 6%). Branch A gave every converged minimum; branch B produced no valid rays.
- Freeze note: the prediction was frozen from branch-A minima (all 8 converged to closest approach ≈1.3e-6) while the jobs were still scanning branch B. At Ω=2, branch B produced no valid rays. If a branch-B solution with a smaller positive Φ turns up later, it will be reported separately and will **not** replace the frozen numbers.

## Gate B result (unsealed after v4_prediction_frozen.json, 09:47)
| Ω | Φ_WKB(n_f; Δ_min(N)), N=24..48 | slope of ln g vs NΦ_WKB (1±0.15) | local 32→40: obs/pred (±15%) |
|---|---|---|---|
| 2.25 | 0.1336 0.1372 0.1395 0.1411 | 0.958 PASS | 0.1433/0.1487 = 0.963 PASS |
| 2.75 | 0.1310 0.1336 0.1353 0.1365 | 0.929 PASS | 0.1327/0.1420 = 0.934 PASS |
**v4: PASS.** Using only Δ_min(N) from the exact spectrum and no steady-state input, the a-priori WKB rate function at the MF fold population predicts the exponential N-dependence of the EP coupling at two unseen Ω. Diagnostic: Φ_WKB exceeds the finite-N exact Φ_N by 1–3%, which is the expected finite-N approach to the rate function. Post-freeze: branch B produced no converged solution at any of the 8 points (jobs stopped after this check), so the frozen numbers stand.
