import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Distinguishability.Combinatorics

/-!
# Theorem 1 and Corollary 1, from explicitly stated assumptions

Mathlib does not yet provide Schatten norms of matrices, state designs over the
Haar measure, or the Berry–Esseen theorem. The paper uses exactly three such
facts. They enter here as hypotheses of the structure `QuantumInputs` and the
Berry–Esseen hypothesis of `corollary_gaussian`. Every other step is proved.

The abstract setting: a finite set `ι` of Pauli strings with a probability
distribution `p`, a dimension `d > 0`, and the average distance `Tbar`.
-/

open Finset BigOperators

namespace Distinguishability

/-- The three quantum facts used in the proof of Theorem 1, stated at the level
of their consequences for `Tbar`.

* `lower` — the variational formula for the trace distance plus the first moment
  of the second state (`E ψ₂ = I/d`) give `Tbar ≥ p(G) - |G|/d` for every set `G`
  (paper Eq. (5)).
* `upper` — Pauli orthogonality and the second moments of a 2-design (Eq. (7)),
  together with Jensen's inequality and `‖A‖₁ ≤ √d ‖A‖₂`, give for every split
  of the Pauli strings into a head `S` and its complement the bound
  `Tbar ≤ p(S) + √((d² v - m²)/(2(d+1)))`, where `m, v` are the first and second
  moments of `p` on the complement (paper Eq. (6) with two classes).

The moment algebra that converts Eq. (7) into the square-root expression is
proved in `Combinatorics.lean` (`moment_difference`, `jensen_rhs`). -/
structure QuantumInputs {ι : Type*} [Fintype ι] [DecidableEq ι]
    (p : ι → ℝ) (d : ℝ) (Tbar : ℝ) : Prop where
  lower : ∀ G : Finset ι, (∑ i ∈ G, p i) - (G.card : ℝ) / d ≤ Tbar
  upper : ∀ S : Finset ι,
    Tbar ≤ (∑ i ∈ S, p i) +
      Real.sqrt ((d ^ 2 * (∑ i ∈ Sᶜ, p i ^ 2) - (∑ i ∈ Sᶜ, p i) ^ 2) / (2 * (d + 1)))

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- **Theorem 1, Eq. (4).** With `d = 2^n` and surprisal `s(P) = -log₂ p_P`, for
every real `a` (the paper states `a ≥ 0`; the proof does not need it):
`P[s < n - a] - 2^{-a} ≤ Tbar ≤ P[s < n + a] + 2^{-a/2}/√2`.

The event `s < n ∓ a` is written as `p_P > 2^{-(n ∓ a)}`. -/
theorem theorem1_spectrum (p : ι → ℝ) (n a : ℝ) (Tbar : ℝ)
    (hp0 : ∀ i, 0 ≤ p i) (hp1 : ∑ i, p i = 1)
    (hQ : QuantumInputs p ((2:ℝ) ^ n) Tbar) :
    (∑ i ∈ Finset.univ.filter (fun i => (2:ℝ) ^ (-(n - a)) < p i), p i) - (2:ℝ) ^ (-a)
        ≤ Tbar ∧
    Tbar ≤ (∑ i ∈ Finset.univ.filter (fun i => (2:ℝ) ^ (-(n + a)) < p i), p i)
        + (2:ℝ) ^ (-a / 2) / Real.sqrt 2 := by
  have h2 : (0:ℝ) < 2 := by norm_num
  have hd : (0:ℝ) < (2:ℝ) ^ n := by positivity
  constructor
  · -- Lower side: take G = {p > 2^{-(n-a)}}; then |G| < 2^{n-a}, so |G|/d ≤ 2^{-a}.
    set G := Finset.univ.filter (fun i => (2:ℝ) ^ (-(n - a)) < p i)
    have hsumG : ∑ i ∈ G, p i ≤ 1 := by
      rw [← hp1]
      exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
        (fun i _ _ => hp0 i)
    have hcard := card_le_of_mass G p ((2:ℝ) ^ (-(n - a))) (by positivity)
      (fun i hi => (Finset.mem_filter.mp hi).2) hsumG
    have hratio : (G.card : ℝ) / (2:ℝ) ^ n ≤ (2:ℝ) ^ (-a) := by
      rw [div_le_iff₀ hd]
      calc (G.card : ℝ) ≤ 1 / (2:ℝ) ^ (-(n - a)) := hcard
        _ = (2:ℝ) ^ (-a) * (2:ℝ) ^ n := by
          rw [one_div, ← Real.rpow_neg h2.le, ← Real.rpow_add h2]; ring_nf
    have := hQ.lower G
    linarith
  · -- Upper side: split at p > 2^{-(n+a)}; the tail has v ≤ 2^{-n-a}.
    set S := Finset.univ.filter (fun i => (2:ℝ) ^ (-(n + a)) < p i)
    have hup := hQ.upper S
    have htail_le : ∀ i ∈ Sᶜ, p i ≤ (2:ℝ) ^ (-(n + a)) := by
      intro i hi
      simp only [S, Finset.mem_compl, Finset.mem_filter, Finset.mem_univ, true_and,
        not_lt] at hi
      exact hi
    have hsumT : ∑ i ∈ Sᶜ, p i ≤ 1 := by
      rw [← hp1]
      exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
        (fun i _ _ => hp0 i)
    have hv := sum_sq_le_of_le Sᶜ p ((2:ℝ) ^ (-(n + a))) (fun i _ => hp0 i) htail_le
      hsumT (by positivity)
    have hv0 : 0 ≤ ∑ i ∈ Sᶜ, p i ^ 2 := Finset.sum_nonneg (fun i _ => sq_nonneg _)
    -- the square-root term is at most √(d v / 2) ≤ 2^{-a/2}/√2
    have hrad := jensen_rhs_le ((2:ℝ) ^ n) (∑ i ∈ Sᶜ, p i ^ 2) (∑ i ∈ Sᶜ, p i) hd hv0
    have hv' : ∑ i ∈ Sᶜ, p i ^ 2 ≤ (2:ℝ) ^ (-n - a) := by
      rw [show -n - a = -(n + a) by ring]; exact hv
    have hsq2 := tail_term_le n a (∑ i ∈ Sᶜ, p i ^ 2) hv0 hv'
    have hroot : Real.sqrt ((((2:ℝ) ^ n) ^ 2 * (∑ i ∈ Sᶜ, p i ^ 2) - (∑ i ∈ Sᶜ, p i) ^ 2)
        / (2 * ((2:ℝ) ^ n + 1))) ≤ (2:ℝ) ^ (-a / 2) / Real.sqrt 2 :=
      (Real.sqrt_le_sqrt hrad).trans hsq2
    linarith [hup, hroot]

/-! ## Corollary 1 (Gaussian crossover), from a Berry–Esseen hypothesis

Let `F(u)` be the probability that the surprisal is below `u`. Berry–Esseen
(assumed) says `|F(n ± a) - Φ(x ∓ a/(√n σ))| ≤ b`. Combined with Theorem 1 and
the Lipschitz bound `|Φ(x) - Φ(y)| ≤ L|x - y|` (with `L = 1/√(2π)`), this gives
`|Tbar - Φ(x)| ≤ L·a/(√n σ) + b + max(2^{-a}, 2^{-a/2}/√2)`. Choosing
`a = 2 log₂ n` makes every term `O(log n / √n)`. -/

theorem corollary_gaussian (Φ : ℝ → ℝ) (L Tbar x δ b Flo Fhi elo ehi : ℝ)
    (hδ : 0 ≤ δ)
    (hLip : ∀ u v, |Φ u - Φ v| ≤ L * |u - v|)
    -- Theorem 1, in the form `Flo - elo ≤ Tbar ≤ Fhi + ehi`
    (hlow : Flo - elo ≤ Tbar) (hhigh : Tbar ≤ Fhi + ehi)
    -- Berry–Esseen (assumed): the two CDF values are within `b` of Φ at `x ∓ δ`
    (hBElo : |Flo - Φ (x - δ)| ≤ b) (hBEhi : |Fhi - Φ (x + δ)| ≤ b) :
    |Tbar - Φ x| ≤ L * δ + b + max elo ehi := by
  have h1 := hLip x (x - δ)
  have h2 := hLip (x + δ) x
  rw [show x - (x - δ) = δ by ring, abs_of_nonneg hδ] at h1
  rw [show x + δ - x = δ by ring, abs_of_nonneg hδ] at h2
  have hm1 : elo ≤ max elo ehi := le_max_left _ _
  have hm2 : ehi ≤ max elo ehi := le_max_right _ _
  rw [abs_le] at hBElo hBEhi h1 h2 ⊢
  constructor <;> linarith [hBElo.1, hBElo.2, hBEhi.1, hBEhi.2, h1.1, h1.2, h2.1, h2.2]

/-! ## The assumptions are satisfiable

`QuantumInputs` is not contradictory: complete depolarization (uniform `p` over
`d²` Pauli strings, every output equal to `I/d`, so `Tbar = 0`) satisfies it.
Theorem 1 is therefore not vacuously true. -/

theorem quantumInputs_satisfiable (k : ℕ) (hk : 1 ≤ k) :
    QuantumInputs (fun _ : Fin (k * k) => (1:ℝ) / (k * k)) (k:ℝ) 0 := by
  have hk0 : (0:ℝ) < k := by exact_mod_cast hk
  refine ⟨fun G => ?_, fun S => ?_⟩
  · simp only [Finset.sum_const, nsmul_eq_mul]
    rw [sub_nonpos, mul_one_div, div_le_div_iff₀ (by positivity) hk0]
    have hc : (0:ℝ) ≤ G.card := Nat.cast_nonneg _
    have hk1 : (1:ℝ) ≤ k := by exact_mod_cast hk
    nlinarith [mul_nonneg hc hk0.le, mul_le_mul_of_nonneg_left hk1 (mul_nonneg hc hk0.le)]
  · have h1 : (0:ℝ) ≤ ∑ i ∈ S, (fun _ : Fin (k * k) => (1:ℝ) / (k * k)) i :=
      Finset.sum_nonneg (fun _ _ => by positivity)
    linarith [Real.sqrt_nonneg ((((k:ℝ)) ^ 2 * (∑ i ∈ Sᶜ, ((1:ℝ) / (k * k)) ^ 2)
      - (∑ i ∈ Sᶜ, (1:ℝ) / (k * k)) ^ 2) / (2 * ((k:ℝ) + 1)))]

end Distinguishability
