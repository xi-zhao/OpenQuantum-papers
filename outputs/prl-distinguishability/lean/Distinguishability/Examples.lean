import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Data.Nat.Choose.Sum

/-!
# The examples and numbers in the paper

All statements here are proved in full from Mathlib, without extra assumptions.
-/

open Real

namespace Distinguishability

/-! ## Local depolarizing noise

For `N_r(ρ) = (1-r)ρ + r I/2` the single-qubit Pauli distribution is
`p = (1 - q, q/3, q/3, q/3)` with `q = 3r/4`. -/

/-- Collision probability of single-qubit depolarizing noise, used for the purity
bound: `∑ p_i² = (1 + 3(1-r)²)/4`. -/
theorem depol_sum_sq (r : ℝ) :
    (1 - 3 * r / 4) ^ 2 + 3 * (r / 4) ^ 2 = (1 + 3 * (1 - r) ^ 2) / 4 := by
  ring

/-- The purity bound `d ∑ p²` decays with `n` exactly when `(1-r)² < 1/3`. -/
theorem purity_decay_iff (r : ℝ) :
    (1 + 3 * (1 - r) ^ 2) / 2 < 1 ↔ (1 - r) ^ 2 < 1 / 3 := by
  constructor <;> intro h <;> linarith

/-- The purity threshold `r = 1 - 1/√3` lies in `(0.422, 0.423)`. -/
theorem purity_threshold_bounds :
    (0.422 : ℝ) < 1 - 1 / Real.sqrt 3 ∧ 1 - 1 / Real.sqrt 3 < 0.423 := by
  have h3 : (0:ℝ) < Real.sqrt 3 := by positivity
  have hs : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  have hlo : (1.7320 : ℝ) < Real.sqrt 3 := by nlinarith [Real.sqrt_nonneg 3]
  have hhi : Real.sqrt 3 < 1.7321 := by nlinarith [Real.sqrt_nonneg 3]
  have hinv_lo : 1 / Real.sqrt 3 < 1 / 1.7320 := one_div_lt_one_div_of_lt (by norm_num) hlo
  have hinv_hi : 1 / 1.7321 < 1 / Real.sqrt 3 := one_div_lt_one_div_of_lt h3 hhi
  constructor
  · have : (1:ℝ) / 1.7320 < 0.578 := by norm_num
    linarith
  · have : (0.577 : ℝ) < 1 / 1.7321 := by norm_num
    linarith

/-! ## Hashing threshold for depolarizing noise

The single-qubit error entropy is `H(q) = -(1-q) log₂(1-q) - q log₂(q/3)`. The
hashing threshold `H(q_c) = 1` lies between `q = 3/16` and `q = 1/5`, which give
`r_c ∈ (0.25, 0.2667)`. Both brackets reduce to integer inequalities. -/

/-- Single-qubit error entropy of depolarizing noise, in nats. -/
noncomputable def depolEntropyNat (q : ℝ) : ℝ :=
  -(1 - q) * Real.log (1 - q) - q * Real.log (q / 3)

/-- At `q = 3/16` the error entropy is below one bit. -/
theorem depolEntropy_below : depolEntropyNat (3 / 16) < Real.log 2 := by
  unfold depolEntropyNat
  have e1 : (1:ℝ) - 3 / 16 = 13 / 16 := by norm_num
  have e2 : (3 / 16 : ℝ) / 3 = 1 / 16 := by norm_num
  rw [e1, e2]
  -- 13·log(16/13) < 4·log 2  ⇔  16^13 < 2^4 · 13^13
  have hint : ((16:ℝ) / 13) ^ 13 < 2 ^ 4 := by norm_num
  have hlog : (13:ℝ) * Real.log (16 / 13) < 4 * Real.log 2 := by
    have h := Real.log_lt_log (by positivity) hint
    rw [Real.log_pow, Real.log_pow] at h
    push_cast at h
    linarith
  have hl13 : Real.log (13 / 16 : ℝ) = -Real.log (16 / 13) := by
    rw [← Real.log_inv]; norm_num
  have hl16 : Real.log (1 / 16 : ℝ) = -(4 * Real.log 2) := by
    rw [show (1 / 16 : ℝ) = (2 ^ 4)⁻¹ by norm_num, Real.log_inv, Real.log_pow]
    push_cast; ring
  rw [hl13, hl16]
  nlinarith [hlog]

/-- At `q = 1/5` the error entropy is above one bit. -/
theorem depolEntropy_above : Real.log 2 < depolEntropyNat (1 / 5) := by
  unfold depolEntropyNat
  have e1 : (1:ℝ) - 1 / 5 = 4 / 5 := by norm_num
  have e2 : (1 / 5 : ℝ) / 3 = 1 / 15 := by norm_num
  rw [e1, e2]
  -- 4·log(5/4) + log 15 > 5·log 2  ⇔  (5/4)^4 · 15 > 2^5
  have hint : (2:ℝ) ^ 5 < ((5:ℝ) / 4) ^ 4 * 15 := by norm_num
  have hlog : 5 * Real.log 2 < 4 * Real.log (5 / 4) + Real.log 15 := by
    have h := Real.log_lt_log (by positivity) hint
    rw [Real.log_mul (by positivity) (by norm_num), Real.log_pow, Real.log_pow] at h
    push_cast at h
    linarith
  have hl45 : Real.log (4 / 5 : ℝ) = -Real.log (5 / 4) := by
    rw [← Real.log_inv]; norm_num
  have hl15 : Real.log (1 / 15 : ℝ) = -Real.log 15 := by
    rw [show (1 / 15 : ℝ) = (15:ℝ)⁻¹ by norm_num, Real.log_inv]
  rw [hl45, hl15]
  nlinarith [hlog]

/-- The hashing threshold exists in `(3/16, 1/5)`: the error entropy crosses one
bit there, by continuity (intermediate value theorem). -/
theorem hashing_threshold_exists :
    ∃ q ∈ Set.Ioo (3 / 16 : ℝ) (1 / 5), depolEntropyNat q = Real.log 2 := by
  have hcont : ContinuousOn depolEntropyNat (Set.Icc (3 / 16) (1 / 5)) := by
    unfold depolEntropyNat
    apply ContinuousOn.sub
    · apply ContinuousOn.mul
      · exact (continuousOn_const.sub continuousOn_id).neg
      · apply ContinuousOn.log (continuousOn_const.sub continuousOn_id)
        intro x hx
        simp only [Set.mem_Icc] at hx
        show (1:ℝ) - x ≠ 0
        linarith [hx.2]
    · apply ContinuousOn.mul continuousOn_id
      apply ContinuousOn.log (continuousOn_id.div_const _)
      intro x hx
      simp only [Set.mem_Icc] at hx
      show x / 3 ≠ 0
      have : 0 < x := by linarith [hx.1]
      positivity
  obtain ⟨q, hq, hval⟩ := intermediate_value_Icc (by norm_num) hcont
    ⟨depolEntropy_below.le, depolEntropy_above.le⟩
  refine ⟨q, ⟨?_, ?_⟩, hval⟩
  · rcases hq.1.eq_or_lt with h | h
    · subst h; exact absurd hval depolEntropy_below.ne
    · exact h
  · rcases hq.2.eq_or_lt with h | h
    · subst h; exact absurd hval depolEntropy_above.ne'
    · exact h

/-! ## Global depolarization (correlated-noise example)

`N(ρ) = θρ + (1-θ) I/d` acts linearly, so it scales every difference of states
by `θ`, and hence the trace distance by `θ`. -/

theorem global_depol_difference {d : ℕ} (θ : ℂ) (ρ σ : Matrix (Fin d) (Fin d) ℂ)
    (hρ : ρ.trace = 1) (hσ : σ.trace = 1) :
    (θ • ρ + (1 - θ) • (ρ.trace / d) • (1 : Matrix (Fin d) (Fin d) ℂ))
      - (θ • σ + (1 - θ) • (σ.trace / d) • (1 : Matrix (Fin d) (Fin d) ℂ))
      = θ • (ρ - σ) := by
  rw [hρ, hσ, smul_sub]; abel

/-- Pauli weights of global depolarization: `p_I = θ + (1-θ)/d²` and
`p_P = (1-θ)/d²` for the other `d² - 1` strings sum to one. -/
theorem global_depol_weights_sum (θ D : ℝ) (hD : D ≠ 0) :
    (θ + (1 - θ) / D) + (D - 1) * ((1 - θ) / D) = 1 := by
  field_simp; ring

/-! ## Complete dephasing and mutually unbiased bases (Eq. (12))

Of the `d(d+1)` MUB states, the `d` computational-basis states survive dephasing
and the others map to `I/d`. An ordered pair contributes trace distance `1` if
both are distinct basis states, `1 - 1/d` if exactly one is a basis state, and
`0` otherwise. Averaging over `d²(d+1)²` ordered pairs gives Eq. (12). -/

theorem mub_average (d : ℝ) (hd : 1 ≤ d) :
    (d * (d - 1) * 1 + 2 * d * (d * d) * (1 - 1 / d)) / (d * (d + 1)) ^ 2
      = (d - 1) * (2 * d + 1) / (d * (d + 1) ^ 2) := by
  have hd0 : d ≠ 0 := by linarith
  have hd1 : d + 1 ≠ 0 := by linarith
  field_simp
  ring

/-- Eq. (12) tends to zero: it is at most `3/d`. -/
theorem mub_average_le (d : ℝ) (hd : 1 ≤ d) :
    (d - 1) * (2 * d + 1) / (d * (d + 1) ^ 2) ≤ 3 / d := by
  have hd0 : 0 < d := by linarith
  rw [div_le_div_iff₀ (by positivity) hd0]
  nlinarith [sq_nonneg d, mul_pos hd0 hd0]

/-! ## Random stabilizer states (Eq. (13))

`∑_{P ≠ I} (1-r)^{w(P)} = (1 + 3(1-r))^n - 1`, by the binomial theorem over the
Hamming weight `w`: there are `C(n,w) 3^w` Pauli strings of weight `w`. -/

theorem stabilizer_weight_sum (n : ℕ) (x : ℝ) :
    ∑ w ∈ Finset.range (n + 1), (n.choose w : ℝ) * 3 ^ w * x ^ w = (1 + 3 * x) ^ n := by
  rw [add_comm (1:ℝ), add_pow]
  apply Finset.sum_congr rfl
  intro w _
  rw [one_pow, mul_one, mul_pow]
  ring

/-- The stabilizer decay rate at `r = 0.9` exceeds the achievable exponent:
`1 - log₂ 1.3 > (1 - log₂ 1.03)/2`, i.e. `1.3² < 2 · 1.03`. -/
theorem stabilizer_rate_gt :
    (1 - Real.log 1.3 / Real.log 2) > (1 - Real.log 1.03 / Real.log 2) / 2 := by
  have hl2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hint : (1.3 : ℝ) ^ 2 < 2 * 1.03 := by norm_num
  have h := Real.log_lt_log (by positivity) hint
  rw [Real.log_pow, Real.log_mul (by norm_num) (by norm_num)] at h
  push_cast at h
  -- h : 2 log 1.3 < log 2 + log 1.03; divide by log 2 > 0
  have key : 2 * Real.log 1.3 / Real.log 2 < 1 + Real.log 1.03 / Real.log 2 := by
    rw [div_lt_iff₀ hl2, add_mul, one_mul, div_mul_cancel₀ _ hl2.ne']
    linarith
  have e : 2 * Real.log 1.3 / Real.log 2 = 2 * (Real.log 1.3 / Real.log 2) := by ring
  rw [e] at key
  linarith

end Distinguishability
