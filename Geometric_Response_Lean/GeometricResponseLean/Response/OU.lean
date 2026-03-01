/-
Copyright (c) 2025 Ovidiu Tataru. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Ovidiu Tataru
-/
import GeometricResponseLean.GeometricResponse
import GeometricResponseLean.Response.Support
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.MeasureTheory.Integral.IntegralEqImproper

/-!
# Ornstein-Uhlenbeck Process: Concrete τ_G Existence

**Target 1 from Priority List**: Prove concrete existence of τ_G for the OU process,
eliminating the `L1_under_H123` axiom in this specific case.

## Mathematical Content

For the Ornstein-Uhlenbeck (OU) process with correlation time τ, the normalized
autocorrelation function is:

```
C(t) = exp(-|t|/τ)
```

This module proves:
1. **L¹ integrability**: C(t) ∈ L¹(ℝ₊) with explicit bound
2. **Finite τ_G**: τ_G = ∫₀^∞ C(t) dt = τ < ∞
3. **Concrete instance**: Specializes the general axiom `L1_under_H123` to OU

## Paper Mapping

- Primary manuscript: uses OU/OU-like processes as a canonical “single-timescale” exemplar for
  \(C_\Pi \in L^1\) and for a closure in which τ_G equals the correlation time.
- Lean mapping: this file supplies a fully formal, axiom-free instance where the Green–Kubo
  integral exists and can be computed exactly.

This provides the **first concrete model** where the Green-Kubo integral
is proven to exist without axioms.

## Proof Strategy

**Step 1**: Show exp(-t/τ) is integrable on [0, ∞)
- Use Mathlib's exponential integral lemmas
- Key: ∫₀^∞ exp(-t/τ) dt = τ (closed form)

**Step 2**: Connect to full ACF via symmetry
- acf_full τ t = exp(-|t|/τ) is even
- Integrate over [0, ∞) gives τ, over ℝ gives 2τ

**Step 3**: Instantiate τ_G definition
- τ_G = ∫₀^∞ C(t) dt = τ by definition for OU
- This is definitional, not axiomatic

## Implementation Notes

This module does NOT require:
- Stochastic calculus (we work with the deterministic ACF)
- Ergodic theory (we prove L¹ directly from the exponential form)
- Spectral theory (we compute the integral explicitly)

It DOES use:
- Basic real analysis (exponential integrability)
- Measure theory (Lebesgue integral)
- The ACF definitions from `GeometricResponse.lean`

## Current Status

**Deep Dive Result**: Found key Mathlib lemmas in
`Mathlib.Analysis.SpecialFunctions.ImproperIntegrals`:
- `integral_exp_neg_Ioi_zero`: ∫ x in Ioi 0, exp(-x) = 1
- `integral_exp_neg_Ioi c`: ∫ x in Ioi c, exp(-x) = exp(-c)

**Adaptation Challenge**: These lemmas use `exp(-x)`, but we need `exp(-t/τ)`.
The mathematical connection is trivial (substitution u = t/τ), but formalizing
requires measure-theoretic infrastructure:
- `MeasureTheory.integral_comp_smul_deriv` (change of variables)
- `MeasureTheory.Measure.map` (measure pushforward)
- Connection between improper integral limits and Lebesgue integrals

**Resolution**: Three postulates for standard undergraduate calculus facts:

## Elimination Status

**ALL POSTULATES ELIMINATED** - All 3 lemmas now PROVEN:

1. **Postulate 1**: `exp_neg_div_integrable_on_Ioi` - **PROVEN**
   - Uses `Mathlib.integrableOn_exp_mul_Ioi` directly
   - No assumptions beyond 0 < τ

2. **Postulate 2**: `exp_neg_abs_div_integrable` - **PROVEN**
   - Split ℝ = (-∞,0] ∪ (0,∞) and show integrability on each piece
   - Uses `Integrable.congr` with ae equality on restricted measures
   - Combines `integrableOn_exp_mul_Iic` and `integrableOn_exp_mul_Ioi`

3. **Postulate 3**: `integral_exp_neg_div_eq_tau` - **PROVEN**  
   - Uses `Mathlib.integral_comp_mul_left_Ioi` + `integral_exp_neg_Ioi_zero`
   - Explicit substitution u = t/τ formalized

**Build Status**: **0 sorries** | **0 axioms** | **Builds successfully**

**Main Results Proven**:
- `ou_acf_in_L1`: OU ACF is L¹-integrable → **eliminates `L1_under_H123` for OU**
- `ou_tau_G_eq_tau`: τ_G = τ exactly for OU (no approximation)
- `ou_tau_G_finite`: τ_G < ∞ for OU (concrete instance)

-/

open Real MeasureTheory

namespace GeometricResponse.OU

variable (τ : ℝ) (hτ : 0 < τ)

/-! ## Postulates for Exponential Integrability

The following postulates encode **standard undergraduate calculus facts** about
exponential integrability. They are mathematically trivial (substitution rule +
symmetry) but require significant Mathlib infrastructure to formalize.

**Mathematical Certainty**: All three results are textbook calculus:
- Postulate 1: ∫₀^∞ exp(-t/τ) dt = τ (substitution u = t/τ)
- Postulate 2: exp(-|t|/τ) integrable on ℝ (evenness ⟹ factor of 2)
- Postulate 3: Closed-form value (same as Postulate 1)

**Mathlib Gaps**: The connection requires:
- `MeasureTheory.integral_comp_smul_deriv` (change of variables for Lebesgue)
- `Integrable.of_even` (symmetry lemmas for even functions)
- Improper integral = Lebesgue integral identification

**Existing Mathlib**: We found `integral_exp_neg_Ioi_zero : ∫ exp(-x) = 1`,
but adapting to `exp(-t/τ)` needs the above infrastructure.
-/

/--
**Lemma** (Exponential Integrability): exp(-t/τ) is integrable on (0,∞)

**Proof strategy**: Use Mathlib's `integral_comp_mul_left_Ioi` scaling lemma
combined with `integral_exp_neg_Ioi_zero`.

**Attempt**: Let me try to connect `exp(-t/τ)` = `exp(-(1/τ)·t)` to the scaling lemma.
-/
lemma exp_neg_div_integrable_on_Ioi (hτ : 0 < τ) :
    IntegrableOn (fun t => exp (-(t / τ))) (Set.Ioi 0) volume := by
  -- Strategy: Use Mathlib's integrableOn_exp_mul_Ioi
  -- which says exp(a*x) is integrable on (0,∞) when a < 0
  
  -- Rewrite exp(-t/τ) as exp((-1/τ) * t)
  have h_eq : (fun t => exp (-(t / τ))) = fun t => exp (-(τ⁻¹) * t) := by
    ext t
    rw [div_eq_mul_inv, mul_comm, neg_mul]
  
  rw [h_eq]
  
  -- Apply integrableOn_exp_mul_Ioi with a = -1/τ
  have h_neg : -(τ⁻¹) < 0 := by
    rw [neg_lt_zero]
    exact inv_pos.mpr hτ
  
  exact integrableOn_exp_mul_Ioi h_neg 0

/-! ## L¹ Integrability of OU ACF -/

/--
**Core Lemma**: The positive-time OU ACF is integrable.

For t ≥ 0, the function C(t) = exp(-t/τ) is Lebesgue integrable on [0, ∞).
-/
lemma acf_pos_integrable_on_Ioi (hτ : 0 < τ) :
    IntegrableOn (fun t => exp (-(t / τ))) (Set.Ioi 0) volume :=
  exp_neg_div_integrable_on_Ioi τ hτ

/--
**Lemma** (Symmetry Extension): exp(-|t|/τ) is integrable on ℝ

**Proof strategy**: Split the integral over ℝ = (-∞,0] ∪ (0,∞) and use evenness.

The function exp(-|t|/τ) is even, so:
- On (0,∞): exp(-|t|/τ) = exp(-t/τ) which we've proven integrable
- On (-∞,0): by symmetry (change of variables t → -t)
-/
lemma exp_neg_abs_div_integrable (hτ : 0 < τ) :
    Integrable (fun t => exp (-|t| / τ)) volume := by
  -- Strategy: Split ℝ = (-∞, 0] ∪ (0, ∞) and show integrability on each piece
  -- Use: ∫_ℝ = ∫_{(-∞,0]} + ∫_{(0,∞)} and both pieces are integrable
  
  -- First show integrable on (0, ∞)
  have h_pos : IntegrableOn (fun t => exp (-|t| / τ)) (Set.Ioi 0) := by
    -- On (0, ∞), |t| = t, so exp(-|t|/τ) = exp(-t/τ) = exp(-(t/τ))
    have h_eq_ae : (fun t => exp (-(t / τ))) =ᵐ[volume.restrict (Set.Ioi 0)] 
                    (fun t => exp (-|t| / τ)) := by
      filter_upwards [ae_restrict_mem measurableSet_Ioi] with t ht
      congr 1
      rw [abs_of_pos ht]
      ring
    exact Integrable.congr (exp_neg_div_integrable_on_Ioi τ hτ) h_eq_ae
  
  -- Now show integrable on (-∞, 0]
  have h_neg : IntegrableOn (fun t => exp (-|t| / τ)) (Set.Iic 0) := by
    -- On (-∞, 0], |t| = -t, so exp(-|t|/τ) = exp(t/τ) = exp(τ⁻¹ * t)
    have h_eq_ae : (fun t => exp (τ⁻¹ * t)) =ᵐ[volume.restrict (Set.Iic 0)] 
                    (fun t => exp (-|t| / τ)) := by
      filter_upwards [ae_restrict_mem measurableSet_Iic] with t ht
      rw [abs_of_nonpos ht]
      rw [div_eq_mul_inv, mul_comm]
      ring_nf
    
    have h_pos_inv : 0 < τ⁻¹ := inv_pos.mpr hτ
    have h_int : IntegrableOn (fun t => exp (τ⁻¹ * t)) (Set.Iic 0) := 
      integrableOn_exp_mul_Iic h_pos_inv 0
    
    exact Integrable.congr h_int h_eq_ae
  
  -- Combine: integrable on Ioi 0 ∪ Iic 0 = ℝ
  have h_union : (Set.Ioi (0 : ℝ) ∪ Set.Iic 0) = Set.univ := by
    ext t
    simp only [Set.mem_union, Set.mem_Ioi, Set.mem_Iic, Set.mem_univ, iff_true]
    exact lt_or_ge 0 t
  
  have h_integrable := IntegrableOn.union h_pos h_neg
  rw [h_union] at h_integrable
  exact integrableOn_univ.mp h_integrable

/--
**Theorem**: The OU ACF is in L¹(ℝ₊).

This is the central result eliminating `L1_under_H123` for OU.

**Paper mapping**: Proves the prerequisite condition from Lines 302-307
for the Green-Kubo definition to be well-defined.
-/
theorem ou_acf_in_L1 (hτ : 0 < τ) :
    Integrable (GeometricResponse.GK.OU.acf_full τ) volume := by
  -- acf_full τ t = exp(-|t|/τ) by definition
  have h_eq : GeometricResponse.GK.OU.acf_full τ = fun t => exp (-|t| / τ) := by
    funext t
    unfold GeometricResponse.GK.OU.acf_full
    simp only [neg_div]
  rw [h_eq]
  exact exp_neg_abs_div_integrable τ hτ

/-! ## Explicit τ_G Computation -/

/--
**Lemma** (Closed-Form Integral): ∫₀^∞ exp(-t/τ) dt = τ

**Proof strategy**: Use `integral_comp_mul_left_Ioi` to relate exp(-t/τ) to exp(-t).

Mathematical outline:
- exp(-t/τ) = exp(-(1/τ)·t)
- By scaling: ∫ exp(-(1/τ)·t) dt = τ · ∫ exp(-u) du = τ · 1 = τ
-/
lemma integral_exp_neg_div_eq_tau (hτ : 0 < τ) :
    ∫ t in Set.Ioi (0 : ℝ), exp (-(t / τ)) = τ := by
  -- Rewrite as composition with τ⁻¹
  have h_eq : (fun t => exp (-(t / τ))) = fun t => exp (-(τ⁻¹ * t)) := by
    ext t
    rw [div_eq_mul_inv, mul_comm]
  rw [h_eq]
  
  -- Use integral_comp_mul_left_Ioi with a=0, b=τ⁻¹, g(x)=exp(-x)
  -- Pattern: ∫ x in Ioi a, g(b*x) = b⁻¹ • ∫ x in Ioi (b*a), g(x)
  have h_inv_pos : 0 < τ⁻¹ := inv_pos.mpr hτ
  have h_scale := MeasureTheory.integral_comp_mul_left_Ioi (fun x => exp (-x)) 0 h_inv_pos
  
  -- Simplify: τ⁻¹ * 0 = 0
  simp only [mul_zero] at h_scale
  
  -- Now h_scale says: ∫ x in Ioi 0, exp(-(τ⁻¹ * x)) = (τ⁻¹)⁻¹ • ∫ x in Ioi 0, exp(-x)
  rw [h_scale]
  
  -- Use integral_exp_neg_Ioi_zero : ∫ x in Ioi 0, exp(-x) = 1
  rw [integral_exp_neg_Ioi_zero]
  
  -- Simplify: (τ⁻¹)⁻¹ • 1 = τ
  rw [inv_inv, smul_eq_mul, mul_one]

/--
**Theorem**: The Green-Kubo time for OU equals the correlation time τ.

For OU, τ_G = ∫₀^∞ exp(-t/τ) dt = τ (exact, no approximation).

**Physical interpretation**: The microscopic correlation time τ (from the
Langevin equation) equals the macroscopic geometry-relaxation time τ_G
(from the Green-Kubo integral). This is a special property of OU.

**Proof**: Direct integration of exponential via substitution.
-/
theorem ou_tau_G_eq_tau (hτ : 0 < τ) :
    ∫ t in Set.Ioi (0 : ℝ), GeometricResponse.GK.OU.acf_full τ t = τ := by
  -- acf_full τ t = exp(-t/τ) for t > 0
  -- On (0,∞), |t| = t, so exp(-|t|/τ) = exp(-t/τ)
  have h_eq : ∀ t, t ∈ Set.Ioi (0 : ℝ) → GeometricResponse.GK.OU.acf_full τ t = exp (-(t / τ)) := by
    intro t ht
    unfold GeometricResponse.GK.OU.acf_full
    have h_abs : |t| = t := abs_of_pos ht
    rw [h_abs]
  -- Use integral_congr_ae to rewrite
  calc ∫ t in Set.Ioi (0 : ℝ), GeometricResponse.GK.OU.acf_full τ t
      = ∫ t in Set.Ioi (0 : ℝ), exp (-(t / τ)) := by
          refine setIntegral_congr_ae measurableSet_Ioi ?_
          apply ae_of_all
          exact h_eq
    _ = τ := integral_exp_neg_div_eq_tau τ hτ

/--
**Corollary**: τ_G is finite for OU.

This eliminates the need for the `L1_under_H123` axiom when working
with OU as the microscopic model.
-/
theorem ou_tau_G_finite (τ : ℝ) (hτ : 0 < τ) :
    ∃ τ_G : ℝ, 0 < τ_G ∧ τ_G = ∫ t in Set.Ioi (0 : ℝ), GeometricResponse.GK.OU.acf_full τ t := by
  use τ
  constructor
  · exact hτ
  · exact (ou_tau_G_eq_tau τ hτ).symm

/-! ## Connection to Paper's Numerical Results -/

/--
**Documentation Theorem**: Links to numerical verification.

The companion Python scripts (`utils.py`, `compute_gk_kk.py`) verify:
- OU with τ = 2.0s gives τ_G ≈ 1.95s (within numerical error)
- Finite-bias corrections approach τ as record length L → ∞

This theorem documents that connection, though the numerical results
themselves are not formalized.

**Paper reference**: Appendix B, Table 1 (Lines 1165-1178)
-/
theorem ou_numerical_verification_reference :
    -- The numerical pipeline verifies τ_G ≈ τ for OU
    -- See: results/finite_bias/finite_bias_raw.csv
    True := by
  trivial

end GeometricResponse.OU

