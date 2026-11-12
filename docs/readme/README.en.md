# OpenQuantum Papers

<p align="center"><a href="../../README.md">简体中文</a> · <strong>English</strong> · <a href="README.ja.md">日本語</a> · <a href="README.ko.md">한국어</a> · <a href="README.es.md">Español</a> · <a href="README.fr.md">Français</a> · <a href="README.de.md">Deutsch</a> · <a href="README.pt.md">Português</a> · <a href="README.ru.md">Русский</a> · <a href="README.ar.md">العربية</a></p>

OpenQuantum’s six research manuscripts range from a shared recovery law for random quantum codes to a quantum circuit simulator compiled into Qwen’s weights. Read the results, derivations and validation records.

Use OpenQuantum to recompute a two-qubit ground state locally and check the energy against an independent exact solution. The example needs no model key.

**[Read the six papers](#papers)** · **[Recompute a two-qubit ground state][openquantum-first-run]** · [Install OpenQuantum][openquantum-start]

<a id="papers"></a>

## Six papers

### QwenQPU · A quantum circuit simulator inside Qwen

Four standard Qwen layers with fixed weights run a classical state-vector simulation of six-qubit circuits, updating complex amplitudes and reading out measurement probabilities. Enabling the executor in the same language model raised the pass count on 48 fixed tasks **from 4/48 to 48/48**. The maximum absolute probability error against the Qiskit reference was **1.42 × 10⁻⁵**.

[Paper and appendices][qwenqpu] · [Implementation and evaluation data][qwenqpu-code]

![The original QwenQPU banner shows circuit instructions, 64 logical amplitude addresses and probabilities from a classical six-qubit simulation.](../../assets/qwenqpu-banner.gif)

*The amplitude and final probability readouts on the right come from a run of the standard Qwen layers.*

### Random quantum codes · A shared critical curve for classical and quantum recovery

For random codes drawn from an exact unitary 2-design, the paper gives common finite-blocklength upper and lower bounds on the ensemble-average optima for whole-message decoding and entanglement recovery under Pauli noise. For i.i.d. noise satisfying the stated regularity conditions, both optima approach **the same Gaussian curve** in the positive-rate critical window.

[Paper][random-codes] · [Supplement][random-codes-si]

![A blue band bounds the ensemble-average optima for both tasks. As blocklength grows, it narrows toward the Gaussian limit shown as a dashed white line.](../../assets/random-code-banner.gif)

*In Fig. 1b, the bound width near the critical center narrows from about 31.8 to 3.4 percentage points across the three blocklengths.*

### Dissipative Majorana systems · Selecting energy relaxation with a four-time correlation

For fixed-order interactions and uniform dissipation, the four fermion insertions form two short-time pairs, each summed over all components. When both pairs have duration ε, the energy signal is of order ε² and the other decay modes contribute at order O(ε⁴). The paper bounds the error over a finite time window.

[Paper][majorana] · [Supplement][majorana-si]

![Summing both the source and readout fermion pairs over all components selects the energy mode. The right panel shows maximum relative contamination as the pair duration decreases.](../../assets/majorana-banner.gif)

*In the 18 finite systems of Fig. 1c at εJ = 0.16 and gJ = 0.7, the largest contribution from non-energy modes relative to the energy signal was 0.0983%.*

### Fermionic Liouvillian exceptional points · Two decay rates from the same initial state

Total particle number and mode occupation imbalance select different decay channels from the same initial state. When the readout returns the retained system or an erasure flag, short vacuum trials attain the optimal scaling of the number of preparations needed to discriminate pump strengths, counting every preparation.

[Paper][fermion-ep] · [Supplement][fermion-ep-si]

![Particle number and mode imbalance give different decay rates from one initial state. The right panel shows the conditional signal and survival fraction.](../../assets/fermionic-ep-banner.gif)

*In the 12-mode example in Fig. 1, the relaxation rates for total particle number and mode occupation imbalance are about 1.679γ and 1.118γ, respectively. At γt = 5, the survival fraction is 1.18 × 10⁻¹⁵.*

### E₇ spectroscopy in Rydberg chains · Four candidate lines from two probes

Uniform and staggered probes put the even- and odd-parity channels on the same energy scale. In one parameter window of the strict blockade model, four low-energy candidate peaks have mass ratios **within 3%** of E₇ predictions. For the shallow A₄ candidate, the binding energy at the 36-atom working point is about 62% above the ideal value, quantifying the finite-chain deviation.

[Paper][e7] · [Supplement][e7-si]

![The staggered and uniform probes each reveal two E₇ candidate peaks. The right panel shows both spectra and the shallow binding energy for 36 atoms.](../../assets/rydberg-e7-banner.gif)

*Figures 2a and 3a show the 36-atom calculation. The A₄ binding energy is 4.92% of the lowest gap, compared with 3.04% in the E₇ prediction.*

### Open Rydberg ensembles · From decay-mode weight exchange to the density signal

Near a mean-field fold, two decay modes exchange weights across an avoided crossing. Eliminating coherences recovers the spectral roots and observable weights. In the 48-atom example with the preparation held fixed, a single mode’s weight changes by about **8.5 × 10⁻⁴**, while the largest change in the full density signal is about **1.9 × 10⁻⁶**.

[Paper][rydberg-decay] · [Supplement][rydberg-decay-si]

![Two decay modes exchange spectral weight at fixed preparation. A common log scale compares one mode’s weight change and the largest difference in the full density signal at the same two detunings.](../../assets/rydberg-decay-banner.gif)

*The right panel compares the mode weight change in Fig. 5c and the full density difference in Fig. 5d in the same density units. Their magnitudes differ by about a factor of 450 between the marked detunings.*

## Run your own calculations with OpenQuantum

OpenQuantum is an open-source AI platform for quantum research available in Codex and as a standalone desktop app. Start with a quantum circuit, simulate it in Qiskit, verify the result independently, and continue working with the same code, parameters and outputs.

**[Install OpenQuantum][openquantum-start]** · [Browse runnable examples][openquantum-examples] · [OpenQuantum on GitHub][openquantum]

## Research records

[Manuscript versions and validation][manuscript-history] · [Research archive][research-archive] · [QwenQPU manuscript source and figures][qwenqpu-tex]

[openquantum]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.en.md
[openquantum-start]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.en.md#quick-start
[openquantum-examples]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.en.md#algorithm-learning-and-method-comparison
[openquantum-first-run]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.en.md#no-model-key-needed-recompute-a-fixed-case
[qwenqpu]: ../../outputs/qwenqpu/paper/main.pdf
[qwenqpu-code]: https://github.com/xi-zhao/QwenQPU/tree/caae45327b46458163d0840a4c1466e5575a7b4f
[qwenqpu-tex]: https://github.com/xi-zhao/OpenQuantum-papers/tree/1749862c959b44dcc2810319ffeed74a725b08f9/outputs/qwenqpu
[random-codes]: ../../outputs/prl-distinguishability/paper.pdf
[random-codes-si]: ../../outputs/prl-distinguishability/supplement.pdf
[majorana]: ../../outputs/prl-dissipative-syk/paper/main.pdf
[majorana-si]: ../../outputs/prl-dissipative-syk/paper/supplement.pdf
[fermion-ep]: ../../outputs/prl-highorder-ep/paper/manuscript.pdf
[fermion-ep-si]: ../../outputs/prl-highorder-ep/paper/supplement.pdf
[e7]: ../../outputs/e7-rydberg-prl/paper/e7_rydberg.pdf
[e7-si]: ../../outputs/e7-rydberg-prl/paper/e7_rydberg_supplement.pdf
[rydberg-decay]: ../../outputs/prl-chiral-switching/paper/main.pdf
[rydberg-decay-si]: ../../outputs/prl-chiral-switching/paper/supplement.pdf
[manuscript-history]: https://github.com/xi-zhao/OpenQuantum-papers/blob/e9acbf31ee1391e5eb3890450744147d917461c0/README.md
[research-archive]: https://github.com/xi-zhao/OpenQuantum-papers/tree/c/import-research-papers
