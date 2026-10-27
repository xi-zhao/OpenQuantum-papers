import Distinguishability.Combinatorics
import Distinguishability.Examples
import Distinguishability.Theorem

/-! Prints the axioms each main result depends on. Only Lean's three standard
axioms (`propext`, `Classical.choice`, `Quot.sound`) should appear. -/

open Distinguishability

#print axioms theorem1_spectrum
#print axioms corollary_gaussian
#print axioms quantumInputs_satisfiable
#print axioms card_le_of_mass
#print axioms sum_sq_le_of_le
#print axioms tail_term_le
#print axioms moment_difference
#print axioms jensen_rhs
#print axioms jensen_rhs_le
#print axioms depol_sum_sq
#print axioms purity_decay_iff
#print axioms purity_threshold_bounds
#print axioms depolEntropy_below
#print axioms depolEntropy_above
#print axioms hashing_threshold_exists
#print axioms global_depol_difference
#print axioms global_depol_weights_sum
#print axioms mub_average
#print axioms mub_average_le
#print axioms stabilizer_weight_sum
#print axioms stabilizer_rate_gt
