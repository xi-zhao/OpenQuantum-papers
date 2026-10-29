# Observable selection: independent derivation and verification record

2026-10-28. The coordinating agent assembled this record from an independent mathematical review and the saved checks in `data/verified_observables.json`. This is internal derivation review, not external human peer review. The independent reviewer did not modify the manuscript or computation files.

## Exact result

Let `m=n−2`, `D=n1−n2`, and `Q_k=D P_(k+1)`. In diagonal Boolean space, the first two spins form the singlet `|10>_B−|01>_B`, while the remaining m spins are symmetric. Hence

- `J_- Q_k=(m−k+1)Q_(k−1)`;
- `J_+ Q_k=(k+1)Q_(k+1)`;
- `Tr Q_k=Tr(N Q_k)=0`;
- `Tr(D Q_k)=2 binom(m,k)`.

These identities assign D to a fixed spin-(n/2−1) copy. Trace and total number belong only to the maximal-spin block. Invariance under both J_+ and J_- ensures that the orthogonal complement does not leak into the measured copy. Non-diagonal energy-basis coherences stay in separate frozen-set blocks and do not contribute.

For arbitrary initial density matrix, aggregate projections are `u_k=Tr(P_k rho0)` and `v_k=Tr(D P_(k+1) rho0)`. With `E_M=exp(gamma t B_(M,d)(z))`, where B acts on aggregate populations,

`<N>=k^T E_n u/(1^T E_n u)`,

`<D>=1^T E_(n−2) v/(1^T E_n u)`,

`P_surv=exp(−n gamma t) 1^T E_n u`.

For the chosen full-rank initial state `rho_p(1+epsilon D)`,

`u_k=binom(n,k)p^k(1−p)^(n−k)`,

`v_k=2 epsilon binom(n−2,k)p^(k+1)(1−p)^(n−k−1)`.

Changing epsilon leaves the entire number and survival curves unchanged; D is proportional to epsilon. At p=.75, epsilon=.6, the initial imbalance is .225. No odd-parity state is prepared.

## Conditions that the manuscript must retain

- A physical nonzero positive initial state has a nonzero Perron denominator overlap. No generic-overlap assumption is needed for this denominator.
- When n−2>=d, the selected biased positive state guarantees a nonzero D Perron overlap. Its exact exponent is `gamma q(kappa_n−kappa_(n−2))`.
- The number signal uses the normalized detector `N/n−nu_* I`. Its overlap with a simple nonleading mode j is `(k/n−nu_*1)^T r_j l_j^dagger u`. This can vanish even for a positive initial state. The internal gap is the number exponent only if the slowest relevant overlap is nonzero.
- Complex roots produce oscillations; zero Jordan blocks produce polynomial prefactors. If n−2<d, the chosen state gives `D(t)~A t^(n−2) exp(−gamma q kappa_n t)`.
- Gamma_D is not universally the even gap. At n=d=4,z=1, the physical even gap is `gamma 24^(1/5)[1−cos(2pi/5)]`, smaller than the imbalance exponent `gamma 24^(1/5)`; the latter also has a t² prefactor.
- Odd sectors are excluded only from the stated one-time even-state/even-observable protocol. An insertion in a fermionic two-time Green function can be odd. No universal unobservability claim is warranted.
- Readout mixing is the ordinary competition of two exponentials. Its application here is quantitative; it is not a new general dynamical principle.

## Reproducible checks

The final `code/observable_relaxation.py` includes 27 exact integer matrix-intertwining cases, all n=2,...,7 and d=1,...,n. Both maximal-spin and selected singlet-sector residuals are exactly zero, including nilpotent small blocks. Detector weights and cancellation against trace/number are checked as integer identities.

Full CAR-superoperator propagation at (n,d)=(4,2),(4,3),(5,2) agrees with marked-mode reduction with maximum absolute error 1.3322676295501878e−15. All three preserve diagonal positivity; initial odd components and even–odd mixing elements are exactly zero. Unrestricted gaps agree with their predictions. Direct even gaps agree within 3.2e−14 in the two nondefective comparison cases.

The n=4,d=3 even-gap comparison initially failed: double-precision diagonalization spuriously split a size-three zero Jordan block by about 5.53e−6. This was not resolved by widening the acceptance tolerance. The data retain the discrepancy and explicitly mark that direct eigenvalue comparison as unsuitable. Selection identities, the prior full-block matrix-identity checks, and direct time propagation provide the evidence in this case.

For (n,d)=(8,2),(12,3),(16,4), the full marked-population reduction is compared independently with the two observable-channel propagators. The largest number-curve error is 3.4e−16 and the largest D-curve error is 1.2e−16. Spectral amplitudes, obtained from left/right projections of simple roots, agree with the propagated long-time signals within .21% on gamma t in [6,9]. Separate tail fits on [5,8] agree with predicted exponents within .23%; these fitted amplitudes and rates are not used in the theoretical figure curves.

For n=12,d=3,z=1/144:

| Quantity | Value in units of gamma |
|---|---:|
| Unrestricted gap | 0.5719628040560201 |
| Imbalance exponent / even gap | 1.1176163958545908 |
| Number exponent | 1.678630296546845 |
| Leading trace exponent | −7.49749344707678 |

At gamma t=5, the survival probability is 1.1829e−15; at gamma t=8, it is 2.0184e−25. The revised main figure includes this cost. The asymptotic exponent separation is verified, but no practical long-time precision advantage or efficient post-selection protocol is claimed.

## Prior-art and proof disposition

van Caspel and Gritsev, PRA 97, 052106 (2018), is cited explicitly for observable-sector selection. The present contribution is the fermionic model's exact two-channel conditional dynamics, its relation to the complete parity-resolved spectrum, and explicit rate/readout consequences. The semiclassical hierarchy remains conditional on radial-edge localization and has no proved uniform global spectral remainder. The original polynomial clarification is retained in the supplement, with source tests and known mathematical precedents respected.
