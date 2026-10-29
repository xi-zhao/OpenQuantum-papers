# Pre-registration v2: a-priori prediction of the EP coupling exponent

Written 2026-10-27, **before** any exact-Hamiltonian prediction run. Earlier attempts (v1 and the post-hoc free-endpoint gMAM) are recorded in memory and in the data files `gmam_runs.jsonl` and `gmam_free.jsonl`. They do not enter this test.

## Gate A: solver validation (must pass first)
- Solver: Lin's method on the reduced exact WKB Hamiltonian (`code/wkbr.py`, `code/lin.py`).
- Quantity: the action S_hl(Δ) of the heteroclinic from the high fixed point to the saddle at Ω=2, V=10.
- Reference: the exact switching exponent s_hl(Δ) = −d ln(λ_sw(1−q))/dN, fitted over N=40–60 (`data/g2_rates.json`).
- Pass criterion: |S_hl/s_hl − 1| ≤ 5% at every Δ in {3.7, 3.8, 3.9, 4.0} where the solver converges (residual ≤ 1e−8), with at least 3 such Δ.
- If Gate A fails, no prediction is reported as a test.

## Gate B: the prediction
- Mechanism (fixed now): the EP coupling g_N between the high-branch intrabasin mode and the ghost-escape mode is set by the large-deviation cost of reaching the ghost region from the high branch. The ghost region is the bottleneck of the MF slow manifold of the vanished low branch.
- Operational definition (fixed now; no free choice afterwards): at each N, take Δ = Re Δ_EP(N) and compute the minimum action from the high fixed point to the hyperplane n = n_b(Δ). Here n_b(Δ) is the bottleneck, i.e. the minimiser of |dn/dt| along the MF slow manifold g(n) = r(n)(1−2n) − n restricted to n < n_hi. The endpoint has transversality p_u = p_v = 0 and is free in (u, v); this is the standard condition for a minimum action to a hyperplane.
- Prediction: σ_pred = slope of N·S_b(N) vs N. For the comparison, take the linear fit of ln g_N vs N·S_b(N) over N ∈ {24, 32, 40, 48}.
- Pass criterion: the fitted slope lies within 1 ± 0.10, and the local exponent d(N·S_b)/dN agrees with −d ln g_N/dN within 10% at N=40.
- Anything else counts as a failure and will be reported as such.

## What is NOT allowed after seeing results
- Changing the endpoint (fold state, hyperplane value, saddle).
- Changing the N range or the fit.
- Using LNA/gMAM values in place of exact-H values.

## Gate A result (recorded before any Gate B run)
Lin's method, exact reduced H. S_hl / exact s_hl: Δ=3.8: 0.9994, Δ=3.9: 1.0111, Δ=4.0: 1.0204, Δ=4.2: 1.0282. Matching residual is 1e-8 to 1e-7. At Δ=3.9 the action is independent of the matching section (n_sec = 0.315, 0.346, 0.377 all give 0.15358). **PASS** (Δ=4.2 is outside the declared set but is shown anyway; 4 converged points, all within 5%; Δ=3.7 pending).
Bug fixed before this result: an incorrectly constructed unstable-subspace basis (complex pairs doubled), and the sign of the backward action piece.

## Gate B solver note (before any valid Gate B number)
Plain forward shooting to the hyperplane did not converge: N=24 and 32 had no hits, and N=40 had transversality residual 0.47. That run gives no number and is discarded. The solver is switched to Lin matching, which is also used in Gate A: forward rays from W^u(hi), backward rays from the endpoint manifold {n=n_b, p_u=p_v=0, h=0, non-trivial p_n root}, matched at the section n_sec=(n_b+n_hi)/2. Quantity, endpoint, N range and criteria are unchanged.

## Gate B result: definition is ill-posed (recorded as a FAILURE of the pre-registered prediction)
Restricted to p_u = p_v = 0, the Hamiltonian is h = −p_n (Ω v + n) (symbolic, exact). It is linear in p_n, so the zero-energy endpoint condition forces p_n = 0, or else the endpoint lies on the line v = −n/Ω. With p_n = 0 the endpoint is a point of zero momentum: the transversality condition p_u = p_v = 0 together with h = 0 leaves no non-trivial instanton ending on the hyperplane n = n_b. Forward shooting found no admissible hit with p_u = p_v = 0; residual ≥ 0.47 at N = 40 and no hits at N = 24, 32. The pre-registered quantity S_b(N) **does not exist** as defined. This is reported as a failed prediction. No replacement endpoint is adopted in this registration.

## CORRECTION to the entry above (same session, before any valid Gate B number)
The "ill-posed" conclusion was wrong and is withdrawn. h = −p_n(Ωv + n) vanishes for p_n ≠ 0 on the line v = −n/Ω, which is where the MF n-velocity is zero. The endpoint manifold {n = n_b, v = −n_b/Ω, u free, p_n free, p_u = p_v = 0} is a valid 2-parameter family. The earlier code assumed h was quadratic in p_n and discarded every endpoint, which is a solver bug. The endpoint definition, N range and criteria are unchanged. Gate B is rerun with backward rays started from this family.

## Gate B final result: **FAIL**
Converged matching residuals are about 2e-11 for all N. S_b(N) = {24: −1.220, 32: 0.475, 40: 0.445, 48: 0.425}.
- At N=24 the action is negative: the solver found an unphysical branch of the endpoint family (p_n = 13.2).
- At N=32–48 the action is S_b ≈ 0.43–0.47, about three times the coupling exponent 0.14 (same order as the failed v1 fold-state prediction).
- Local exponent at N≈40: predicted d(N·S_b)/dN = 0.33, observed −d ln g/dN = 0.14 (factor 2.3). The registered fit of ln g vs N·S_b over all 4 N gives slope −0.05 (criterion 1±0.10).
The pre-registered mechanism (coupling = large-deviation cost of reaching the MF ghost bottleneck from the high branch) is **refuted** at this level. No post-hoc variant is adopted.
