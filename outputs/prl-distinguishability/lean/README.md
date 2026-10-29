# Lean formalization of the earlier two-state manuscript

**Scope (2026-10-28):** This directory checks results from the earlier two-state
distinguishability manuscript. All theorem and equation numbers below refer to
that version. It does **not** formalize the current classical/quantum recovery
bracket, its probability bound, or the common Gaussian coding curve in
`../paper.tex` and `../supplement.tex`. No new Lean build was needed for the
current revision; the build statement below records the earlier validation.

Lean 4 (v4.36.0-rc1) + Mathlib (tag `v4.36.0-rc1`). `lake build` succeeds.
No `sorry`, `admit`, or new `axiom`. `#print axioms` in `Audit.lean`
shows that every result depends only on Lean's standard axioms
(`propext`, `Classical.choice`, `Quot.sound`).

## What is proved, and from what

The proofs fall into two groups. **Fully proved** results use Mathlib only.
**Conditional** results take the three quantum facts that Mathlib cannot yet
express as named hypotheses, and Lean checks every step from those hypotheses
to the paper's statement.

| Paper | Lean name | File | Status |
|---|---|---|---|
| Thm. 1, Eq. (4) | `theorem1_spectrum` | `Theorem.lean` | conditional on `QuantumInputs` |
| Cor. 1, Eq. (9) | `corollary_gaussian` | `Theorem.lean` | conditional on a Berry–Esseen hypothesis |
| `QuantumInputs` is consistent | `quantumInputs_satisfiable` | `Theorem.lean` | fully proved |
| \|G\| < d·2^{-a} for G={s<n−a} | `card_lt_of_mass`, `card_le_of_mass` | `Combinatorics.lean` | fully proved |
| tail v ≤ 2^{-n-a} | `sum_sq_le_of_le` | `Combinatorics.lean` | fully proved |
| tail term ≤ 2^{-a/2}/√2 | `tail_term_le` | `Combinatorics.lean` | fully proved |
| Eq. (7), variance of the difference | `moment_difference` | `Combinatorics.lean` | fully proved (algebra) |
| Eq. (7) → square root in Eq. (6) | `jensen_rhs`, `jensen_rhs_le` | `Combinatorics.lean` | fully proved (algebra) |
| ∑p² for depolarizing noise | `depol_sum_sq` | `Examples.lean` | fully proved |
| purity threshold (1−r)² < 1/3 | `purity_decay_iff` | `Examples.lean` | fully proved |
| 0.422 < 1 − 1/√3 < 0.423 | `purity_threshold_bounds` | `Examples.lean` | fully proved |
| H(3/16) < 1 bit < H(1/5) | `depolEntropy_below`, `depolEntropy_above` | `Examples.lean` | fully proved |
| hashing threshold exists in (3/16, 1/5) | `hashing_threshold_exists` | `Examples.lean` | fully proved |
| global depolarization scales differences by θ | `global_depol_difference`, `global_depol_weights_sum` | `Examples.lean` | fully proved |
| MUB average, Eq. (12), and → 0 | `mub_average`, `mub_average_le` | `Examples.lean` | fully proved |
| stabilizer Pauli-weight sum, Eq. (13) | `stabilizer_weight_sum` | `Examples.lean` | fully proved |
| stabilizer rate 0.621 > 0.479 | `stabilizer_rate_gt` | `Examples.lean` | fully proved |

## The assumed inputs

`QuantumInputs p d Tbar` in `Theorem.lean` packages two statements, each the
direct consequence of quantum facts that Mathlib does not yet contain:

- `lower`: `Tbar ≥ p(G) − |G|/d` for every set `G` of Pauli strings (paper
  Eq. (5)). It rests on the variational formula `½‖A−B‖₁ ≥ tr Π(A−B)` and the
  first moment `E ψ₂ = I/d`.
- `upper`: `Tbar ≤ p(S) + √((d²v − m²)/(2(d+1)))` for every head `S`, where
  `m, v` are the first and second moments of `p` on the complement (paper
  Eq. (6) with two classes). It rests on Pauli orthogonality, the 2-design
  second moments, Jensen's inequality, and `‖A‖₁ ≤ √d‖A‖₂`. The algebra from
  the moment identities to this square root is proved (`moment_difference`,
  `jensen_rhs`).

`corollary_gaussian` takes the Berry–Esseen bounds on the two CDF values and a
Lipschitz bound for Φ as hypotheses.

Not formalized: the type-class decay exponent, Eq. (10), and the
near-threshold rate expansion, Eq. (11). Both need large-deviation theory that
Mathlib does not yet provide.

## Building

The toolchain and Mathlib cache (about 11 GB) live on the external disk at
`/Volumes/Elements/CodexArtifacts/Projects/openQuantum/2026-10-27--lean-formalization`.
`.lake` in this directory is a symbolic link to `lake-packages` there.

```bash
export ELAN_HOME=/Volumes/Elements/CodexArtifacts/Projects/openQuantum/2026-10-27--lean-formalization/elan
```

```bash
export PATH="$ELAN_HOME/bin:$PATH"
```

```bash
lake build
```

Imports are restricted to the Mathlib modules actually used; importing all of
Mathlib from the USB disk took about 13 minutes per file, against about one
minute now.
