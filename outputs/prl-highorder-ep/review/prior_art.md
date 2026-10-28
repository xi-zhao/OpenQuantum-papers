# Prior-art source record

Checked on 2026-10-28. This record lists sources and directly verified evidence for the PRA manuscript. Bibliographic metadata were checked against arXiv and/or Crossref; mathematical statements were checked against the downloaded full texts. It is not a claim of exhaustive literature coverage. No unpublished manuscript was uploaded to external services.

## van Caspel and Gritsev: observable-sector selection is established prior art

Moos van Caspel and Vladimir Gritsev, “Symmetry-protected coherent relaxation of open quantum systems,” *Physical Review A* **97**, 052106 (2018). DOI: [10.1103/PhysRevA.97.052106](https://doi.org/10.1103/PhysRevA.97.052106).

- Full text: https://arxiv.org/abs/1802.05293 ; checked v3 dated 2018-05-07.
- Eq. (4) resolves an expectation value into the product of an initial-state overlap, an observable overlap, and the mode exponential. The text explicitly states that different observables can have different decay rates because they see different symmetry sectors.
- Section III develops weak-symmetry sector selection; Section V, including Eq. (14), applies it to magnetization in a spin chain.
- This directly precedes the present manuscript's central organizing idea. The present work cannot claim that observable-selected relaxation or exclusion of a slow mode by symmetry is a new general principle.
- The checked example concerns weak discrete spin-chain symmetries, not fermionic parity superselection. It does not give the present uniform multiparticle-pump Boolean-spin reduction, exact N/D conditional channels, or the stated fixed-d rate hierarchy.
- Local full-text evidence: `review/symmetry-source-fn5l8aoc/paper.txt`, Eq. (4) and Sections III/V. The downloaded full text is excluded from the source bundle and Git submission.

The revised manuscript's narrower contribution is a quantitative realization in the Xu–Yi model: exact assignment of N and D to two spin blocks, physically admissible simultaneous preparation, explicit normalized dynamics and readout mixing, and a semiclassical rate hierarchy. The readout crossover of two exponentials is elementary and is not claimed as an independent general result. These distinctions strengthen the physical question without making a claim of PRL-level novelty or guaranteeing PRA acceptance.

## Xu and Yi: source model and characteristic-polynomial approximation

Mingtao Xu and Wei Yi, “Higher-order Liouvillian exceptional points in the dissipative dynamics of quadratic fermions,” *Physical Review A* **113**(6), 062455 (2026). DOI: [10.1103/5sk1-qjm1](https://doi.org/10.1103/5sk1-qjm1).

- Crossref records online publication on 2026-06-23. Metadata: https://api.crossref.org/works/10.1103/5sk1-qjm1
- arXiv: https://arxiv.org/abs/2602.00486 ; v1 submitted 2026-01-31, v2 submitted 2026-07-18.
- v1 title: “Exactly solvable higher-order Liouvillian exceptional points in dissipative fermionic systems.”
- The relevant approximation is Eq. (49) in v1 and Eq. (54) in v2. The v2 equation was visually checked on PDF page 7 and uses an approximate-equality sign. Its stated consequence is the fractional splitting exponent, not an explicit claim that all residual roots remain at zero.
- Eq. (54) is a first-order truncation in the perturbation parameter. On the splitting scale, the omitted higher powers can have the same asymptotic order as the retained terms. Criticism should be framed as nonuniform truncation, while retaining the splitting exponent.
- Higher powers first occur when n+1 >= 2(d+1). The displayed n=4 cases with d=2,3,4 do not meet this condition. Appendix C separately solves d=1 and states the full Liouvillian gap gamma*sqrt(z).
- Formal publication metadata were verified, but access to the APS publication PDF failed. Equation-by-equation checks use arXiv v2, not a verified comparison against the typeset journal PDF.
- Local evidence: `/tmp/nov_KnVd/2602.00486v2.txt`, lines 1208–1275 and 1686–1719; layout extraction `/tmp/wy2602.xVtB/p.txt`, lines 538–589; page image `/tmp/xuyi54page.2nJbEb/page-07.jpg`.

## Ma and Edelman: full leading polynomial and Hessenberg rings

Yanyuan Ma and Alan Edelman, “Nongeneric eigenvalue perturbations of Jordan blocks,” *Linear Algebra and its Applications* **273**(1–3), 45–63 (1998). DOI: [10.1016/S0024-3795(97)00342-X](https://doi.org/10.1016/S0024-3795(97)00342-X).

- Author-hosted full text: https://people.stat.sc.edu/yanyuanma/papers/me.pdf
- Theorem 1 and Remark 1, pp. 52–53, give the generic k-ring and remainder-ring result with explicit nonvanishing conditions.
- The proof and Eqs. (1)–(3), pp. 53–54, retain all terms epsilon^q*lambda^(m-qk) contributing to the leading splitting polynomial. Coefficients are sums of products over disjoint sections. This is directly relevant to the nonuniform truncation in Xu–Yi.
- The generic remainder-ring statement must not be substituted for the exact zero-root statement for a pure positive double-band matrix. Use McMillen for the latter.
- Local evidence: `/tmp/maedelman.sc.8hNI0e/paper.txt`, lines 372–457; visually checked pages `/tmp/ma-thm.HgKK7N/page-09.jpg` and `page-10.jpg`.

## Graefe et al.: higher-order spin perturbations and eigenvalue rings

E. M. Graefe, U. Günther, H. J. Korsch, and A. E. Niederle, “A non-Hermitian PT-symmetric Bose-Hubbard model: eigenvalue rings from unfolding higher-order exceptional points,” *Journal of Physics A: Mathematical and Theoretical* **41**(25), 255206 (2008). DOI: [10.1088/1751-8113/41/25/255206](https://doi.org/10.1088/1751-8113/41/25/255206).

- arXiv: https://arxiv.org/abs/0802.3164 ; checked v2.
- Section 5.1, Eqs. (39)–(41), treats triplet unfolding and explicitly invokes Ma–Edelman.
- Eqs. (69)–(70) and the following paragraph generalize to higher powers of the spin operator, (k+1)-rings, and residual eigenvalues.
- Ring splitting and its connection to perturbation bandwidth are prior results. The paper does not provide the present general fixed-d spectral-edge constant or local spacing coefficient; its outlook identifies the large-N analysis as future work.
- Local evidence: `/tmp/nov_KnVd/0802.3164.txt`, lines 1340–1478, 2380–2496, and 2592–2611.

## Shapiro and Tater: the d=2 leading spectral edge is already known

Boris Shapiro and Miloš Tater, “On spectral asymptotic of quasi-exactly solvable quartic potential,” *Analysis and Mathematical Physics* **12**, article 2 (2022). DOI: [10.1007/s13324-021-00612-2](https://doi.org/10.1007/s13324-021-00612-2).

- Online publication: 2021-11-01; journal issue year: 2022.
- arXiv: https://arxiv.org/abs/1412.3026 ; checked v2 dated 2021-10-11.
- Eq. (1.3) defines T_n = partial_x^2 - (x^2-a) partial_x + nx. At a=0, its matrix in Eq. (1.4) is exactly N_{n,2}(2)^T in the manuscript's normalization.
- Theorem 1(i) gives spectral radius (3/4)*n^(4/3)*(1+o(1)). Exact scaling in z converts this to C_2 = 3/2^(7/3), the same leading constant as in the manuscript.
- Theorem 1(ii) gives the limiting three-ray spectral distribution; Proposition 3 supplies the modulo-three polynomial structure. The leading d=2 edge must be attributed as a recovered known result.
- The checked paper does not state the manuscript's general fixed-d first finite-size correction or its local edge-spacing formula.
- Local evidence: `/tmp/qesquartic.9vEMOg/paper.txt`, lines 58–125 and 330–345; visually checked theorem page `/tmp/qes-theorem.etINAf/page-02.jpg`.

## McMillen: exact rays, simple nonzero roots, and zero multiplicity

Tyler McMillen, “On the eigenvalues of double band matrices,” *Linear Algebra and its Applications* **431**(10), 1890–1897 (2009). DOI: [10.1016/j.laa.2009.06.026](https://doi.org/10.1016/j.laa.2009.06.026).

- Preprint: https://arxiv.org/abs/0811.0023 ; preprint title “A note on the eigenvalues of double band matrices.”
- Theorem 1.1 applies to two positive bands with offsets b and -k. For coprime b,k, dimension M=qp+r, and p=b+k, the zero eigenvalue has algebraic multiplicity r; the remaining eigenvalues are rho_s*exp(2*pi*i*j/p), with distinct positive rho_s.
- Taking M=n+1, b=d, k=1 applies directly to N_{n,d}(z) for z>0 and n>=d. Thus every nonzero eigenvalue is simple and lies on one of the d+1 exact rays. Zero algebraic multiplicity does not assert semisimplicity.
- This theorem controls off-axis competitors in a gap argument. It does not establish the local Holstein–Primakoff expansion, a uniform remainder estimate, or which spin sector determines the complete Liouvillian gap.
- Local evidence: `/tmp/mcmillen2009.fYZohO/paper.txt`, lines 23–99, especially Theorem 1.1 at lines 67–72.

## Prosen: pre-existing tensor-Jordan and quadratic-fermion decomposition

Tomaž Prosen, “Spectral theorem for the Lindblad equation for quadratic open fermionic systems,” *Journal of Statistical Mechanics: Theory and Experiment* **2010**(07), P07020 (2010). DOI: [10.1088/1742-5468/2010/07/P07020](https://doi.org/10.1088/1742-5468/2010/07/P07020).

- arXiv: https://arxiv.org/abs/1005.0763
- Author-hosted journal PDF: https://chaos.fmf.uni-lj.si/wp-content/uploads/2017/06/ProsenJSTAT2010.pdf
- Section 3.1, Lemma 3.1, Eq. (43), gives the Jordan decomposition of a tensor sum of two Jordan blocks.
- Section 4.1, Theorem 4.1, gives the complete spectrum, invariant subspaces, and largest Jordan blocks for quadratic open fermionic systems. The subsequent discussion distinguishes this from writing every Jordan block in the most general case.
- These establish prior art for the tensor-Jordan machinery. They do not by themselves supply the present uniformly weighted d-particle creation-channel reduction for d>1.
- Local evidence: `/tmp/prosen2010.8HaUHq/paper.txt`, lines 329–350 and 634–743.

## Shiralieva et al.: general no-jump multi-block exceptional points

Aysel Shiralieva, Grigory A. Starkov, and Björn Trauzettel, “Multi-block exceptional points in open quantum systems,” *Physical Review A* **113**, 052202 (2026). DOI: [10.1103/gtw6-yn95](https://doi.org/10.1103/gtw6-yn95).

- arXiv: https://arxiv.org/abs/2509.11856 ; checked v3.
- Eq. (18) explicitly gives the 2n-1, 2n-3, ..., 1 block structure of the effective no-jump Liouvillian associated with an n-dimensional Hamiltonian Jordan block.
- This is relevant prior art for multi-block terminology; the checked models and result do not supply the present uniform d-fermion creation-channel decomposition.
- Local evidence: `/tmp/nov_KnVd/2509.11856.txt`, lines 421–481.

## Direct asymptotic neighbors checked for applicability

- Alexandre Eremenko and Andrei Gabrielov, “Quasi-exactly solvable quartic: elementary integrals and asymptotics,” *Journal of Physics A: Mathematical and Theoretical* **44**, 312001 (2011). DOI: [10.1088/1751-8113/44/31/312001](https://doi.org/10.1088/1751-8113/44/31/312001). arXiv: https://arxiv.org/abs/1104.2305 . Sections 8–9 take a large potential parameter with fixed polynomial degree, not the manuscript's fixed-d, large-n limit. Evidence: `/tmp/eremenko2011.hMliVC/paper.txt`, lines 497–605.
- Boris Shapiro, “Van Vleck spectra of high-order Heun operators: finite-band universality and exterior asymptotics,” arXiv:2607.02700 (2026), https://arxiv.org/abs/2607.02700 . The operator class requires deg Q_k=k+1 for the highest derivative coefficient. For d>1, the manuscript's transposed differential operator (z/d!)*partial_x^d - x^2*partial_x + nx has a constant highest derivative coefficient and falls outside that hypothesis. Evidence: `/tmp/heun2026.mN8VBU/paper.txt`, lines 70–103.

All temporary full-text paths above are local verification pointers; the DOI, arXiv, and author-hosted links are the durable source references.
