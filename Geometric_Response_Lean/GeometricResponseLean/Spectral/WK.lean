import GeometricResponseLean.GeometricResponse

/-! # Wiener-Khinchin Relations: Spectral View of Temporal Correlations

## Overview

This module establishes the connection between **temporal structure** (autocorrelation in time)
and **spectral structure** (power spectral density in frequency).

The Wiener-Khinchin theorem states that the **PSD is the Fourier transform of the autocorrelation**:
```
S(ω) = ∫₋∞^∞ C(t)·exp(-iωt) dt
```

For **real, stationary processes**, the DC (zero-frequency) PSD has a special relationship to
the Green-Kubo integral:
```
S(0) = 2·∫₀^∞ C(t) dt = 2·τ_G
```

This establishes that **spectral weight at DC is determined by temporal memory**.

## Connection to Geometric Response Theory

The Wiener-Khinchin relations provide a **bridge** between:

1. **Green-Kubo (time domain)**: τ_G = ∫₀^∞ C(t) dt
2. **Kramers-Kronig (frequency domain)**: KK slope ~ τ_G via ∫ R(t)·sin(ωt)/ω dt
3. **Spectral density**: S(0) = 2·τ_G

This triad of relations ensures **consistency** across representations:
```
Time Domain         Frequency Domain         DC Limit
   C(t)         ⟷      S(ω)             →     S(0) = 2·τ_G
     ↓                   ↓                          ↓
   τ_G          ⟷   KK slope            →     G_eff(0) = G₀
```

## Connection to Quantum Structural Triad

From the **Quantum Structural Triad** (with its own Rocq formalization):
- The triad Π_A = √F_Q · S_A · I quantifies structural persistence
- Under H1-H3 (locality, coupling, non-stationarity), Π_A cannot vanish in time-average
- Temporal persistence of Π_A manifests as autocorrelation C(t)

From the **Geometric Response** theory:
- C(t) captures memory of structural features over time
- The integral τ_G = ∫ C(t) dt quantifies total memory
- Wiener-Khinchin translates this into spectral language: S(ω)
- The DC value S(0) = 2·τ_G provides a **spectral signature** of structural memory

## Physical Interpretation

The **DC power spectral density** S(0) represents the "low-frequency spectral weight" — how
much of the system's fluctuation energy is concentrated at long timescales.

For systems with **strong memory** (large τ_G):
- S(0) is large (spectral weight at DC)
- χ(ω) varies slowly with ω (broad response)
- Operational scales {ω_G, λ_G, m_G} are pushed to low energy

For systems with **weak memory** (small τ_G):
- S(0) is small (little DC spectral weight)
- χ(ω) varies rapidly (narrow response)
- Operational scales are at high energy

This connects **temporal persistence of quantum structure** (from the Triad) to
**observable spectral features** (via Wiener-Khinchin).

## What This Module Proves

1. **DC Wiener-Khinchin**: If the GK integral converges to L, then the DC PSD converges to 2L
2. **Truncation consistency**: Finite-window approximations respect the 2:1 factor
3. **Limit transfer**: Filter.Tendsto properties are preserved under the PSD transformation

These are **wrapper lemmas** that connect the GK definitions in `GeometricResponse.lean` to
spectral interpretations used in computational verification.

## Proof Strategy

The proof is straightforward because the **PSD is defined** as the Fourier transform of C(t):
```
S₀(T) := 2·∫₀^T C(t) dt
```

So the DC limit is:
```
lim_{T→∞} S₀(T) = 2·lim_{T→∞} ∫₀^T C(t) dt = 2·τ_G
```

The lemmas here simply **package this relationship** in a form suitable for downstream use.

## Why This Matters

The Wiener-Khinchin relations provide **computational accessibility**:

1. **Time-domain simulations** compute C(t) from dynamics
2. **FFT methods** transform C(t) → S(ω) efficiently
3. **DC extraction** yields S(0), which must equal 2·τ_G
4. **Consistency checks** verify that GK and spectral methods agree

This bridges **formal theory** (Lean proofs) and **numerical validation** (Python/Julia scripts).

-/

namespace GeometricResponse.Spectral.WK

open GeometricResponse
open GeometricResponse.Spectral
open Filter

/-- DC Wiener–Khinchin wrapper on truncated forms:
    if the Green–Kubo truncation `∫₀ᵀ C → L` as `T→∞`, then the one-sided
    DC PSD truncation `S₀(T) := 2∫₀ᵀ C` tends to `2L`. -/
lemma psd0_trunc_tendsto_of_gk_trunc
    {C : ℝ → ℝ} {L : ℝ}
    (h : Filter.Tendsto (fun T : ℝ => gkTrunc C T) Filter.atTop (nhds L)) :
    Filter.Tendsto (fun T : ℝ => psd0Trunc C T) Filter.atTop (nhds (2 * L)) :=
  psd0_limit_of_gk_limit (C:=C) (L:=L) h

end GeometricResponse.Spectral.WK


