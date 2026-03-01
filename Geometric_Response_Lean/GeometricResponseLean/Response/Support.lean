/-
Copyright (c) 2025 Ovidiu Tataru. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Ovidiu Tataru
-/

import Mathlib.MeasureTheory.Integral.IntervalIntegral.Basic
import Mathlib.MeasureTheory.Function.StronglyMeasurable.Basic
import Mathlib.MeasureTheory.Integral.Bochner.Set

/-!
# Supporting Necessity Lemmas for Geometric Response

This file provides the mathematical foundations for the Green-Kubo → Kramers-Kronig
necessity chain formalized in `KK.lean`.

## Structure

**Proven Theorems** (no axioms, no sorries):
1. `paley_zygmund_simple`: Existence form of Paley-Zygmund inequality
2. `integral_pos_of_nonneg_of_pos_measure`: Measure-theoretic positivity
3. `tau_G_positive`: Positivity of Green-Kubo integral (applies previous lemma)

**Physics-Level Axioms** (detailed proofs in companion paper):
4. `Axioms.variance_bound`: Short-time spectral bound τ_G ≥ Var[Π]/⟨Π̇²⟩
5. `Axioms.L1_under_H123`: Ergodic/mixing → L¹ autocorrelation

## Physical Context

Primary manuscript targeted by this Lean development:
- "From Microscopic Structural Interdependence to Causal Dispersive Geometry"
  (GK→KK necessity chain; τ_G from autocorrelation; causality/passivity/regression constraints).

The triad observable Π_A(t) = √F_Q · S_A(t) · I(A:Ā,t) measures structural persistence.
Under hypotheses H1 (locality), H2 (genuine coupling), H3 (ergodicity), the
Green-Kubo time τ_G = ∫₀^∞ C_Π(t) dt is well-defined and positive.

## References

- Paley & Zygmund (1932): On some series of functions
- Kubo (1957): Statistical-mechanical theory of irreversible processes
- Paper mapping: hypotheses H1–H3 and the Green–Kubo/Kramers–Kronig chain (section-level; no line numbers)
-/

noncomputable section

open MeasureTheory Filter Topology ProbabilityTheory
open scoped Topology ENNReal NNReal

namespace GeometricResponse

/-!
## 1. Paley-Zygmund Inequality (Proven - Existence Form)

We prove the **existence form**: for X ≥ 0, X ∈ L², 𝔼[X] > 0,
there exists θ ∈ (0,1) such that ℙ{X ≥ θ·𝔼[X]} > 0.

**Note**: The classical Paley-Zygmund gives an explicit lower bound
ℙ{X ≥ θ𝔼[X]} ≥ (1-θ)²(𝔼[X])²/𝔼[X²]. We only need existence here,
which follows from taking θ = 1/2 and using Markov + Cauchy-Schwarz.

**Application**: Ensures the triad observable Π_A is strictly positive
on a set of positive measure.
-/

variable {Ω : Type*} [MeasureSpace Ω] [IsProbabilityMeasure (volume : Measure Ω)]

/--
Existence form of Paley-Zygmund inequality.

For nonnegative integrable X with positive expectation, there exists θ ∈ (0,1)
such that X exceeds θ·𝔼[X] on a set of positive measure.

**Proof**: Take θ = 1/2. If the set {X ≥ 𝔼[X]/2} had zero measure, then
X < 𝔼[X]/2 almost surely, so ∫X ≤ (1/2)𝔼[X] by `integral_mono_ae`,
contradicting 𝔼[X] > 0.

**Note**: The classical Paley-Zygmund assumes X ∈ L² and gives an explicit bound
ℙ{X ≥ θ𝔼[X]} ≥ (1-θ)²(𝔼[X])²/𝔼[X²]. We only need existence of such θ,
which requires only integrability.
-/
theorem paley_zygmund_simple
    (X : Ω → ℝ)
    (_hX_nonneg : ∀ ω, 0 ≤ X ω)
    (hX_int : Integrable X)
    (hX_pos : 0 < ∫ ω, X ω) :
    ∃ θ : ℝ, 0 < θ ∧ θ < 1 ∧ 0 < volume {ω | θ * (∫ ω', X ω') ≤ X ω} := by
  -- Take θ = 1/2 as witness
  use 1/2
  refine ⟨by norm_num, by norm_num, ?_⟩
  
  -- Proof by contradiction: if {X ≥ (1/2)E[X]} has measure zero,
  -- then X < (1/2)E[X] a.e., so ∫X ≤ (1/2)E[X], contradicting E[X] > 0.
  by_contra h_zero
  push_neg at h_zero
  simp only [le_zero_iff] at h_zero
  
  -- X < (1/2)E[X] almost everywhere
  have h_ae_lt : ∀ᵐ ω, X ω < (1/2) * ∫ ω', X ω' := by
    have : volume {ω | ¬(X ω < (1/2) * ∫ ω', X ω')} = 0 := by
      simp only [not_lt]; exact h_zero
    exact this
  
  -- Therefore ∫X ≤ (1/2)E[X]
  have h_bound : ∫ ω, X ω ≤ (1/2) * ∫ ω', X ω' := by
    have h_ae_le : ∀ᵐ ω, X ω ≤ (1/2) * ∫ ω', X ω' := by
      filter_upwards [h_ae_lt] with ω hω
      exact le_of_lt hω
    set c := (1/2) * ∫ ω', X ω' with hc_def
    calc ∫ ω, X ω 
        ≤ ∫ ω, c := integral_mono_ae hX_int (integrable_const c) h_ae_le
      _ = c := by 
          rw [integral_const]
          have : (volume (Set.univ : Set Ω)) = 1 := measure_univ
          simp [this]
  
  -- But (1/2)E[X] < E[X] when E[X] > 0
  have h_half_lt : (1/2 : ℝ) * ∫ ω', X ω' < ∫ ω, X ω := by
    linarith [hX_pos]
  
  -- Contradiction
  linarith [h_bound, h_half_lt]

/-!
## 2. Measure-Theoretic Positivity (Proven)

Standard integration theory: integral of nonnegative function positive
on a set of positive measure is strictly positive.
-/

/--
Core measure-theoretic lemma: integral of nonnegative integrable function
is strictly positive if the set where it's positive has positive measure.

**Proof**: Use `set_integral_pos_iff_support_of_nonneg_ae` pattern from Mathlib.
-/
lemma integral_pos_of_nonneg_of_pos_measure
    {f : ℝ → ℝ}
    (h_int : Integrable f)
    (h_nonneg : ∀ᵐ t, 0 ≤ f t)
    (h_pos : 0 < volume {t | 0 < f t}) :
    0 < ∫ t, f t := by
  -- Key insight: {t | 0 < f t} ⊆ support f, and if it has positive measure,
  -- then the integral must be positive
  
  -- First, show integral is nonnegative
  have h_ge : 0 ≤ ∫ t, f t := integral_nonneg_of_ae h_nonneg
  
  -- Now show it's not zero
  by_contra h_not_pos
  push_neg at h_not_pos
  
  -- If ∫ f ≤ 0 and ∫ f ≥ 0, then ∫ f = 0
  have h_zero : ∫ t, f t = 0 := le_antisymm h_not_pos h_ge
  
  -- By integral_eq_zero_iff_of_nonneg_ae, f = 0 almost everywhere
  have h_ae_zero : f =ᵐ[volume] 0 := by
    rw [← integral_eq_zero_iff_of_nonneg_ae h_nonneg h_int]
    exact h_zero
  
  -- But then {t | 0 < f t} has measure zero
  have h_meas_zero : volume {t | 0 < f t} = 0 := by
    -- {t | 0 < f t} ⊆ {t | f t ≠ 0}, and the latter has measure zero
    apply le_antisymm _ (zero_le _)
    calc volume {t | 0 < f t}
        ≤ volume {t | f t ≠ 0} := measure_mono (fun t ht => by
          simp only [Set.mem_setOf_eq] at ht ⊢
          linarith)
      _ = 0 := ae_iff.mp h_ae_zero
  
  -- Contradiction with h_pos
  rw [h_meas_zero] at h_pos
  exact absurd rfl (ne_of_gt h_pos)

/--
Positivity of Green-Kubo integral: specialized version of previous lemma.

**Input**: Normalized autocorrelation C_Π
**Output**: τ_G = ∫ C_Π > 0

This is pure measure theory, no physics.
-/
theorem tau_G_positive
    (C_Pi : ℝ → ℝ)
    (h_integrable : Integrable C_Pi)
    (h_nonneg : ∀ᵐ t, 0 ≤ C_Pi t)
    (h_pos_set : 0 < volume {t | 0 < C_Pi t})
    (τ_G : ℝ)
    (hτ_def : τ_G = ∫ t, C_Pi t) :
    0 < τ_G := by
  rw [hτ_def]
  exact integral_pos_of_nonneg_of_pos_measure h_integrable h_nonneg h_pos_set

/-!
## 3. Physics-Level Axioms

The following results require machinery beyond current Mathlib.
They are treated as **explicit axioms** with paper references.
-/

/--
Structure encoding physical hypotheses H1-H3 from the paper.

**NOTE**: This structure is **not used** in any proofs in this file.
It serves as a documentation bridge between the paper's physics hypotheses
and the Lean formalization. Actual proofs use explicit integrability/mixing hypotheses.

**H1 (Locality)**: Bounded generator (finite propagation speed)
**H2 (Genuine coupling)**: Nonzero mutual information on positive measure set
**H3 (Ergodic non-stationarity)**: Time averages converge to ensemble averages
-/
structure TriadHypotheses (Pi_triad : ℝ → ℝ) : Prop where
  /-- H1: Bounded time derivative (locality) - documentation placeholder -/
  bounded_generator : ∃ M, ∀ t, |Pi_triad t| ≤ M
  /-- H2: Nontrivial on a set of positive measure -/
  genuine_coupling : ∃ δ > 0, ∃ S : Set ℝ, MeasurableSet S ∧
                     0 < volume S ∧ ∀ t ∈ S, δ ≤ |Pi_triad t|
  /-- H3: Ergodic limit exists -/
  ergodic : ∃ μ, Tendsto (fun T => (1 / T) * ∫ t in (0)..T, Pi_triad t) atTop (𝓝 μ)

namespace Axioms

/-!
## Axiom Status

**This namespace intentionally contains NO axioms or sorries.**

The axioms `variance_bound` and `L1_under_H123` previously defined here have been
**completely eliminated** from the codebase.

### What happened to them?

1. **`variance_bound`**: τ_G ≥ Var[Π]/⟨Π̇²⟩
   - ✅ Proven for OU case in `Response/OU.lean`
   - ⚠️ General case requires spectral theory (not formalized, not needed for main chain)

2. **`L1_under_H123`**: H1+H2+H3 → C_Π ∈ L¹
   - ✅ Proven for OU case in `Response/OU.lean` (`ou_acf_in_L1`)
   - ⚠️ General case: empirically validated, bypassed in formal proofs

### Current status

All theorems in this project now use **explicit** integrability/mixing hypotheses
or refer to the **concrete OU proofs** (0 axioms, 0 sorries).

If you see an axiom in this namespace, something has gone wrong.
-/

end Axioms

/-!
## 4. Combined Well-Definedness Theorem

Composition: Physics axioms → Measure theory → τ_G > 0
-/

/--
Minimal core theorem: positivity from autocorrelation function assumptions.

Uses only measure theory (no physics narrative).
-/
theorem tau_G_positive_from_ACF
    (C_Pi : ℝ → ℝ)
    (h_int : Integrable C_Pi)
    (h_nonneg : ∀ᵐ t, 0 ≤ C_Pi t)
    (h_pos : 0 < volume {t | 0 < C_Pi t})
    (τ_G : ℝ) (hτ_def : τ_G = ∫ t, C_Pi t) :
    0 < τ_G :=
  tau_G_positive C_Pi h_int h_nonneg h_pos τ_G hτ_def

/-
**REMOVED THEOREM**: `tau_G_well_defined`

This theorem previously used the `L1_under_H123` axiom, which has been eliminated.

For the OU case, the complete proof without axioms is in OU.lean:
- `ou_acf_in_L1`: Complete L¹ proof
- `ou_tau_G_finite`: τ_G < ∞ proven

For general cases requiring this result, use `tau_G_positive_from_ACF` directly
with an explicit integrability hypothesis rather than deriving it from H1-H3.
-/

end GeometricResponse

