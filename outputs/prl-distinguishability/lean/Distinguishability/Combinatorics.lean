import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Algebra.Order.BigOperators.Group.Finset

/-!
# Elementary facts used in the paper

These statements are proved in full, from Mathlib, without extra assumptions.
Section references are to `paper.tex`.
-/

open Finset BigOperators

namespace Distinguishability

/-! ## Counting the likely errors (Theorem 1, lower side)

If every element of a finite set `G` has probability strictly larger than `c > 0`
and the total probability is at most one, then `G` has fewer than `1 / c`
elements. With `c = 2^(a-n)` this is `|G| < d 2^{-a}`. -/

theorem card_lt_of_mass {ι : Type*} (G : Finset ι) (p : ι → ℝ) (c : ℝ) (hc : 0 < c)
    (hp : ∀ i ∈ G, c < p i) (hsum : ∑ i ∈ G, p i ≤ 1) (hne : G.Nonempty) :
    (G.card : ℝ) < 1 / c := by
  have hlt : (G.card : ℝ) * c < ∑ i ∈ G, p i := by
    have := Finset.sum_lt_sum_of_nonempty hne hp
    simpa [Finset.sum_const, nsmul_eq_mul] using this
  rw [lt_div_iff₀ hc]
  linarith

/-- `card_lt_of_mass` also gives the non-strict bound when `G` may be empty. -/
theorem card_le_of_mass {ι : Type*} (G : Finset ι) (p : ι → ℝ) (c : ℝ) (hc : 0 < c)
    (hp : ∀ i ∈ G, c < p i) (hsum : ∑ i ∈ G, p i ≤ 1) :
    (G.card : ℝ) ≤ 1 / c := by
  rcases G.eq_empty_or_nonempty with h | h
  · subst h; simp; positivity
  · exact (card_lt_of_mass G p c hc hp hsum h).le

/-! ## The tail second moment (Theorem 1, upper side)

If every element of `C` has probability at most `c ≥ 0` and the total mass of `C`
is at most one, then `∑_{C} p_P^2 ≤ c`. With `c = 2^{-n-a}` this is
`v ≤ 2^{-n-a}`. -/

theorem sum_sq_le_of_le {ι : Type*} (C : Finset ι) (p : ι → ℝ) (c : ℝ)
    (hp0 : ∀ i ∈ C, 0 ≤ p i) (hpc : ∀ i ∈ C, p i ≤ c) (hsum : ∑ i ∈ C, p i ≤ 1)
    (hc : 0 ≤ c) : ∑ i ∈ C, p i ^ 2 ≤ c := by
  calc ∑ i ∈ C, p i ^ 2 ≤ ∑ i ∈ C, c * p i := by
        apply Finset.sum_le_sum
        intro i hi
        have := hp0 i hi
        nlinarith [hpc i hi]
    _ = c * ∑ i ∈ C, p i := by rw [Finset.mul_sum]
    _ ≤ c * 1 := by exact mul_le_mul_of_nonneg_left hsum hc
    _ = c := mul_one c

/-- The square-root tail term of Eq. (6): with `v ≤ 2^{-n-a}` and `d = 2^n`,
`√(d v / 2) ≤ 2^{-a/2} / √2`. -/
theorem tail_term_le (n a v : ℝ) (_hv0 : 0 ≤ v) (hv : v ≤ (2:ℝ) ^ (-n - a)) :
    Real.sqrt ((2:ℝ) ^ n * v / 2) ≤ (2:ℝ) ^ (-a / 2) / Real.sqrt 2 := by
  have h2 : (0:ℝ) < 2 := by norm_num
  have hdv : (2:ℝ) ^ n * v ≤ (2:ℝ) ^ (-a) := by
    calc (2:ℝ) ^ n * v ≤ (2:ℝ) ^ n * (2:ℝ) ^ (-n - a) :=
          mul_le_mul_of_nonneg_left hv (by positivity)
      _ = (2:ℝ) ^ (-a) := by
          rw [← Real.rpow_add h2]; ring_nf
  have hsq : (2:ℝ) ^ (-a) = ((2:ℝ) ^ (-a / 2)) ^ 2 := by
    rw [← Real.rpow_natCast, ← Real.rpow_mul h2.le]; norm_num
  rw [div_eq_mul_inv ((2:ℝ) ^ (-a / 2)), ← Real.sqrt_inv]
  rw [show (2:ℝ) ^ (-a / 2) = Real.sqrt (((2:ℝ) ^ (-a / 2)) ^ 2) by
        rw [Real.sqrt_sq (by positivity)]]
  rw [← Real.sqrt_mul (by positivity)]
  apply Real.sqrt_le_sqrt
  rw [← hsq, div_eq_mul_inv]
  exact mul_le_mul_of_nonneg_right hdv (by positivity)

/-! ## The second-moment algebra (Eq. (7))

Given the three moment identities, the variance of the difference follows by
algebra. The identities themselves are the assumed 2-design input
(see `Assumptions.lean`). -/

theorem moment_difference (d v m A B : ℝ) (hd : 0 < d)
    (hA : A = (d * v + m ^ 2) / (d + 1)) (hB : B = m ^ 2 / d) :
    2 * A - 2 * B = 2 * d / (d + 1) * v - 2 * m ^ 2 / (d * (d + 1)) := by
  subst hA hB
  have hd1 : d + 1 ≠ 0 := by linarith
  have hd0 : d ≠ 0 := hd.ne'
  field_simp
  ring

/-- The square-root expression in Eq. (6) is what Jensen produces from Eq. (7):
`(d/4)·(2d v/(d+1) - 2m²/(d(d+1))) = (d²v - m²)/(2(d+1))`. -/
theorem jensen_rhs (d v m : ℝ) (hd : 0 < d) :
    d / 4 * (2 * d / (d + 1) * v - 2 * m ^ 2 / (d * (d + 1)))
      = (d ^ 2 * v - m ^ 2) / (2 * (d + 1)) := by
  have hd1 : d + 1 ≠ 0 := by linarith
  have hd0 : d ≠ 0 := hd.ne'
  field_simp
  ring

/-- The square-root term never exceeds `√(d v / 2)`, the form used for the tail. -/
theorem jensen_rhs_le (d v m : ℝ) (hd : 0 < d) (hv : 0 ≤ v) :
    (d ^ 2 * v - m ^ 2) / (2 * (d + 1)) ≤ d * v / 2 := by
  rw [div_le_div_iff₀ (by positivity) (by norm_num)]
  nlinarith [sq_nonneg m, mul_nonneg hd.le hv]

end Distinguishability
