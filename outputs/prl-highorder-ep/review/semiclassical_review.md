# Independent semiclassical review record

2026-10-28. This record was assembled by the coordinating agent from a separate mathematical review of `paper/manuscript.tex`, `code/spectral_analysis.py`, `data/verified_precision.json`, and `data/verified_perron.json`. It is an internal independent derivation check, not external human peer review. The reviewer did not change the manuscript, code, or data.

## Coefficients checked

In the orthonormal spin representation, set

`T_s = J_- + J_+^d/(d s^(d-1))`, `Q = b†b = (X²+P²−1)/2`.

Keeping noncommutative ordering in the expansion about +x gives

`T_s = s(1+1/d) + H_0 + s^(-1/2)V_3 + O(s^(-1))`,

`H_0 = 1 − P² − (d+1)X²/2`,

`V_3 = −i(d−1)[{X,Q}/2 + (d−2)X³/6]`.

The constant +1 comes from `−2Q = 1−X²−P²`; it is not a missing normal-order correction. `V_3` is odd under oscillator parity and has vanishing diagonal expectation in every eigenstate of `H_0`. Its second-order effect enters at the same order as the next even operator correction.

Using `e^η = z^(1/(d+1)) a_d n^β` and oscillator eigenvalues `1−(2l+1)Ω_d` reproduces the manuscript's leading constant, first finite-size correction, and local level-spacing coefficient. No specific coefficient error was found.

## Scope of the asymptotic statement

The operator expansion is in a fixed low-excitation sector; it is not a norm estimate on the whole spin space. The review did not prove a global uniform spectral remainder bound. The manuscript explicitly labels the formula as a formal fixed-d semiclassical expansion and identifies the association of local l=0,1 levels with the largest two positive radial roots as the edge-localization assumption supported by the numerical tests.

McMillen's exact double-band theorem excludes off-axis eigenvalues close to the Perron edge on the local spacing scale. It does not by itself prove the radial-edge identification or the stated asymptotic remainder. The exact min formulas for unrestricted and parity-even spectral gaps do not require that assumption; their large-n simplifications do.

## Numerical interpretation

The saved arbitrary-precision computations cover d=2,3,4 and n=16,32,48,64. At n=64, the actual real-part gap divided by the predicted local spacing is approximately 0.9977918, 0.9987756, and 1.0024913, respectively. Perron data through n=512 support both the leading term and its first correction.

The 50- and 70-digit eigenproblem outputs are converted to Python floats before their saved differences are compared. A recorded zero difference means equality after double-precision rounding; it does not certify 50 accurate digits. The manuscript claims only that the reported spacings agree within 1e-10, which is consistent with the saved data. Floating-point Collatz bounds likewise provide convergence diagnostics, not interval-arithmetic proofs.

## Physical interpretation checked

For hole number h, trace nonincrease requires `z binom(h,d) <= h`. Since `h/binom(h,d)=d/binom(h−1,d−1)` is nonincreasing, the most stringent inequality is at h=n, yielding exactly `z<=n/binom(n,d)`.

The separation relative to the leading eigenvalue controls normalized conditional dynamics, whereas the survival weight has exponent equal to that leading eigenvalue. These are distinct quantities. The manuscript's separation of the unrestricted operator gap, physical parity-even gap, and symmetric-sector particle-number relaxation is consistent with this interpretation.

## Disposition

The observable-focused revision retains the explicit fixed-low-excitation qualification and radial-edge localization limitation in its main text; the full derivation and original spectral tests are now in `paper/supplement.tex`. The coefficients themselves are unchanged. No further mandatory coefficient or physical correction was identified. The formal asymptotic evidence has not been promoted to a rigorous global spectral theorem.
