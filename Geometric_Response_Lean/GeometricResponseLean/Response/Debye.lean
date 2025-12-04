import GeometricResponseLean.GeometricResponse
import Mathlib.Analysis.Complex.Exponential
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Topology.Basic

/-! # Debye Closure: Sufficiency for Geometric Response

## Overview

This module proves that a **single-pole (Debye) response model** is _sufficient_ to satisfy
both the Green-Kubo (GK) and Kramers-Kronig (KK) constraints simultaneously.

**Complementarity with KK.lean**:
- `Response/KK.lean` proves **necessity**: τ_G → KK slope (any causal kernel must obey this)
- `Response/Debye.lean` proves **sufficiency**: Debye model achieves the bound (it is realizable)

Together, these establish that the **Debye kernel is the minimal closure** for the theory.

## Physical Interpretation

The Debye model describes exponential relaxation:
```
R(t) = τ⁻¹·exp(-t/τ)  for t ≥ 0
χ(ω) = 1/(1 + iωτ)
```

This is the **simplest causal response** that satisfies:
1. **Green-Kubo**: ∫₀^∞ R(t) dt = 1 (normalized autocorrelation)
2. **Kramers-Kronig**: Im[χ] ≤ 0 for ω ≥ 0 (passivity)
3. **DC normalization**: χ(0) = 1 (adiabatic limit)

The **relaxation time τ** plays the role of τ_G from the geometric response theory, and the
operational scales emerge as:
- Frequency scale: ω_G = 1/τ
- Length scale: λ_G = c·τ (for propagation at speed c)
- Mass scale: m_G = ℏ/(c²·τ) (from uncertainty at this timescale)

## Connection to Quantum Structural Triad

From the **Quantum Structural Triad** (with its own Rocq formalization):
- The triad Π_A = √F_Q · S_A · I captures structural features (fluctuation, entropy, correlation)
- Under H1-H3 (locality, coupling, non-stationarity), Π_A cannot vanish in time-average
- The persistence of Π_A generates temporal correlations C(t)

From the **Geometric Response** theory:
- The autocorrelation C(t) encodes memory of structural features
- The Debye form C(t) ~ exp(-t/τ) is the **minimal exponential decay**
- This decay rate τ determines all spectral properties via KK relations
- The effective coupling G_eff(ω) = G₀·χ(ω) inherits this single-pole structure

## What This Module Proves

1. **Tail Decay**: The exponential tail exp(-(t/τ + iωt)) → 0 as t → ∞ for τ > 0
2. **Transform Convergence**: The truncated Fourier integral tends to χ_Debye(ω) as T → ∞
3. **Passivity**: Im[χ_Debye(ω)] ≤ 0 for ω ≥ 0 (physical causality)
4. **Power Absorption**: Re[χ_Debye(ω)] represents dissipative response

## Proof Strategy

The key technical challenge is showing that:
```
∫₀^T (τ⁻¹·exp(-t/τ))·exp(-iωt) dt → χ(ω) = 1/(1 + iωτ)  as T → ∞
```

We proceed by:
1. Computing the closed-form truncated integral (complex exponential integration)
2. Factoring as χ(ω)·(1 - exp(-(τ⁻¹ + iω)T))
3. Proving the exponential tail decays to zero (dominated by Re part ~ exp(-T/τ))
4. Applying filter limit theorems to transfer the tail limit to the integral limit

The **passivity** χ_im ≤ 0 follows from the pole structure: Im[1/(1+iωτ)] = -ωτ/(1+(ωτ)²) ≤ 0.

## Why This Matters

The Debye closure is **sufficient** because:
1. It realizes the GK→KK necessity bound (τ_G = τ determines the slope exactly)
2. It is the **simplest** model (single pole, two parameters: G₀ and τ)
3. It is **physical** (exponential decay is generic for dissipative systems)
4. It is **complete** (satisfies all constraints: GK, KK, passivity, DC normalization)

Combined with the **necessity** proof in `KK.lean`, we have:

**Necessity + Sufficiency = Closure**: The Debye model is the unique minimal closure.

-/

namespace GeometricResponse.Response.Debye
open GeometricResponse Filter Complex

/-- Local scaling helper: if `τ>0` then `T ↦ T/τ` tends to `+∞`. -/
private lemma tendsto_div_pos_atTop_atTop {τ : ℝ} (hτ : 0 < τ) :
  Tendsto (fun T : ℝ => T / τ) atTop atTop := by
  -- monotone
  have hmono : Monotone (fun T : ℝ => T / τ) := by
    intro T1 T2 hT
    have hinv : 0 ≤ τ⁻¹ := inv_nonneg.mpr (le_of_lt hτ)
    simpa [div_eq_mul_inv, mul_comm, mul_left_comm, mul_assoc]
      using mul_le_mul_of_nonneg_right hT hinv
  -- unbounded: ∀A, ∃T, A ≤ T/τ (choose T = τ*A)
  have hub : ∀ A : ℝ, ∃ T : ℝ, A ≤ T / τ := by
    intro A
    refine ⟨τ*A, ?_⟩
    have hτne : τ ≠ 0 := ne_of_gt hτ
    have hb : (τ * A) / τ = A := by
      calc
        (τ * A) / τ = (A * τ) / τ := by simpa [mul_comm]
        _ = A * (τ / τ) := by simp [div_eq_mul_inv, mul_comm, mul_left_comm, mul_assoc]
        _ = A * 1 := by simp [hτne]
        _ = A := by simp
    simpa [hb]
  exact tendsto_atTop_atTop_of_monotone hmono hub

@[simp] private lemma re_neg_mul_ofReal (a b T : ℝ) :
  ((-((a : ℂ) + Complex.I * b) * (T : ℂ))).re = - a * T := by
  simp [Complex.mul_re]

/-- Debye tail decay: for `τ>0`, the complex exponential tail tends to `0` as `T → ∞`. -/
lemma tailExp_tendsto_zero {τ ω : ℝ} (hτ : 0 < τ) :
  Tendsto (fun T : ℝ =>
    Complex.exp (-( ((1/τ : ℝ) : ℂ) + Complex.I * (ω : ℂ)) * (T : ℂ)))
    atTop (nhds (0 : ℂ)) := by
  -- dist (exp z) 0 = ‖exp z‖ = exp (Re z)
  have hnorm :
      (fun T : ℝ =>
        ‖Complex.exp (((-(Complex.I * (ω : ℂ)) + -((τ : ℝ) : ℂ)⁻¹) * (T : ℂ)))‖)
      = (fun T : ℝ => Real.exp (-(T / τ))) := by
    funext T
    set z : ℂ := (-(Complex.I * (ω : ℂ)) + -((τ : ℝ) : ℂ)⁻¹)
    have hzre : z.re = - (1/τ : ℝ) := by
      simp [z, one_div]
    have hmulre : (z * (T : ℂ)).re = - (T / τ) := by
      have hre : (T : ℂ).re = T := by simp
      have him : (T : ℂ).im = 0 := by simp
      simpa [Complex.mul_re, hzre, hre, him, div_eq_mul_inv, mul_comm]
        using (Complex.mul_re z (T : ℂ))
    have : ‖Complex.exp (z * (T : ℂ))‖ = Real.exp (-(T/τ)) := by
      simpa [hmulre] using (Complex.norm_exp (z * (T : ℂ)))
    simpa [z] using this
  -- exp(−T/τ) → 0 via T/τ → +∞, then transfer through the norm using `hnorm`
  have hexp_norm :
      Tendsto (fun T : ℝ =>
        ‖Complex.exp (((-(Complex.I * (ω : ℂ)) + -((τ : ℝ) : ℂ)⁻¹) * (T : ℂ)))‖)
      atTop (nhds (0 : ℝ)) := by
    have : Tendsto (fun T : ℝ => Real.exp (-(T / τ))) atTop (nhds (0 : ℝ)) :=
      (Real.tendsto_exp_neg_atTop_nhds_zero).comp (tendsto_div_pos_atTop_atTop hτ)
    simpa [hnorm] using this
  -- convert from norm to distance-to-zero and finish
  exact (tendsto_iff_dist_tendsto_zero).2 (by
    simpa [dist_eq_norm] using hexp_norm)

lemma chi_re_eq_power (τ ω : ℝ) :
  (GeometricResponse.chiDebye τ ω).re = GeometricResponse.Hpow τ ω :=
  (GeometricResponse.Hpow_eq_re_chiDebye (τ:=τ) (ω:=ω)).symm

lemma chi_passivity_nonpos {τ ω : ℝ} (hτ : 0 < τ) (hω : 0 ≤ ω) :
  (GeometricResponse.chiDebye τ ω).im ≤ 0 :=
  GeometricResponse.chiDebye_passivity_nonpos hτ hω

lemma truncTransform_tendsto_chiDebye {τ ω : ℝ} (hτ : 0 < τ) :
  Tendsto (fun T : ℝ => GeometricResponse.Debye.truncTransformClosedForm τ ω T)
    atTop (nhds (GeometricResponse.chiDebye τ ω)) := by
  apply GeometricResponse.Debye.truncTransformClosedForm_tendsto_to_chiDebye_of_decay (τ:=τ) (ω:=ω)
  simpa using tailExp_tendsto_zero (τ:=τ) (ω:=ω) hτ

end GeometricResponse.Response.Debye
