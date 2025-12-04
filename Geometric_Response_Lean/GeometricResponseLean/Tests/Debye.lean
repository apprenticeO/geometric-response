import GeometricResponseLean.GeometricResponse
import GeometricResponseLean.Response.Debye

open GeometricResponse

/- Basic checks for the Debye susceptibility and kernel bounds.
   Paper mapping (Geometric Response and the Frequency–Dependent Gravitational Coupling.tex):
   - Response kernel and passivity: L795–L821 (χ, H=|χ|^2, Imχ≤0)
   - GK definition and PSD link: L301–L314
   - Scales and λ_G, m_G: L360–L372, L591–L595
   - OU exemplar and Debye closure: L269–L277 -/
namespace GeometricResponse.Tests

-- H(ω) ≥ 0 and ≤ 1
example (τ ω : ℝ) : 0 ≤ Hpow τ ω := Hpow_nonneg τ ω
example (τ ω : ℝ) : Hpow τ ω ≤ 1 := Hpow_le_one τ ω

-- DC checks
example (τ : ℝ) : Hpow τ 0 = 1 := Hpow_zero τ
example (G0 τ : ℝ) : Geff G0 τ 0 = G0 := Geff_zero G0 τ

-- Evenness in frequency
example (τ ω : ℝ) : Hpow τ (-ω) = Hpow τ ω := Hpow_even τ ω
example (G0 τ ω : ℝ) : Geff G0 τ (-ω) = Geff G0 τ ω := Geff_even G0 τ ω

-- Antitonicity in ω²
example (τ ω ω0 : ℝ) (h : ω0 ^ 2 ≤ ω ^ 2) : Hpow τ ω ≤ Hpow τ ω0 :=
  Hpow_le_of_sq_le τ ω ω0 h

-- High-frequency tail bound
example (τ ω : ℝ) (h : 0 < |ω * τ|) : Hpow τ ω ≤ 1 / (ω * τ) ^ 2 :=
  Hpow_le_inv_sq τ ω h

-- OU GK exemplar
example (τ : ℝ) : GeometricResponse.GK.OU.tauG τ = τ := rfl

-- OU ACF sanity
example (τ : ℝ) : GeometricResponse.GK.OU.acf τ 0 = 1 :=
  GeometricResponse.GK.OU.acf_at_zero τ

-- OU integral approximation properties
example {τ T : ℝ} (hτpos : 0 < τ) (hT : 0 ≤ T) :
    0 ≤ GeometricResponse.GK.OU.acfIntApprox τ T :=
  GeometricResponse.GK.OU.acfIntApprox_nonneg hτpos hT

example {τ T : ℝ} (hτ : 0 ≤ τ) :
    GeometricResponse.GK.OU.acfIntApprox τ T ≤ τ :=
  GeometricResponse.GK.OU.acfIntApprox_le_tau hτ

example {τ T1 T2 : ℝ} (hτpos : 0 < τ) (hT : T1 ≤ T2) :
    GeometricResponse.GK.OU.acfIntApprox τ T1 ≤ GeometricResponse.GK.OU.acfIntApprox τ T2 :=
  (GeometricResponse.GK.OU.acfIntApprox_mono hτpos) hT

-- lambdaG basic properties
example (c τ : ℝ) (hc : 0 ≤ c) (hτ : 0 ≤ τ) : 0 ≤ lambdaG c τ :=
  GeometricResponse.lambdaG_nonneg hc hτ

example (τ : ℝ) : lambdaG 0 τ = 0 := GeometricResponse.lambdaG_zero_left τ
example (c : ℝ) : lambdaG c 0 = 0 := GeometricResponse.lambdaG_zero_right c

example (a c τ : ℝ) : lambdaG (a * c) τ = a * lambdaG c τ :=
  GeometricResponse.lambdaG_mul_left a c τ

example (c a τ : ℝ) : lambdaG c (a * τ) = a * lambdaG c τ :=
  GeometricResponse.lambdaG_mul_right c a τ

-- Operational scales extras
example {τ : ℝ} (hτ : 0 < τ) : 0 < GeometricResponse.omegaG τ :=
  GeometricResponse.omegaG_pos hτ

example {ħ c τ : ℝ} (hħ : 0 ≤ ħ) (hc : 0 < c) (hτ : 0 < τ) :
    0 ≤ GeometricResponse.mG ħ c τ :=
  GeometricResponse.mG_nonneg hħ hc hτ

-- Debye small-ω slope (limit form)
example (τ : ℝ) :
    Filter.Tendsto (fun ω : ℝ => - τ / (1 + (ω * τ) ^ 2))
      (nhds 0) (nhds (-τ)) :=
  GeometricResponse.kk_slope_limit_Debye_nhds τ

-- Debye imaginary part and passivity
-- Debye imaginary part and passivity lemmas will be added when the API is pinned.
example {τ ω : ℝ} (hτ : 0 < τ) (hω : 0 ≤ ω) :
    GeometricResponse.chiImDebyePow τ ω ≤ 0 :=
  GeometricResponse.chiImDebyePow_nonpos hτ hω

-- Passivity for the complex Debye susceptibility: Im χ_D(ω) ≤ 0 for ω ≥ 0, τ > 0
example {τ ω : ℝ} (hτ : 0 < τ) (hω : 0 ≤ ω) :
    (chiDebye τ ω).im ≤ 0 :=
  GeometricResponse.chiDebye_passivity_nonpos hτ hω

-- Complex link will be enabled once compatible helpers are pinned.
example (τ ω : ℝ) : (chiDebye τ ω).im = GeometricResponse.chiImDebyePow τ ω :=
  GeometricResponse.chiDebye_im_eq_pow τ ω

-- Power-form Imχ properties
example (τ : ℝ) : GeometricResponse.chiImDebyePow τ 0 = 0 :=
  GeometricResponse.chiImDebyePow_zero τ

example (τ ω : ℝ) : GeometricResponse.chiImDebyePow τ (-ω) = - GeometricResponse.chiImDebyePow τ ω :=
  GeometricResponse.chiImDebyePow_odd τ ω

-- Debye tail decay and improper transform
example {τ ω : ℝ} (hτ : 0 < τ) :
    Filter.Tendsto
      (fun T : ℝ =>
        Complex.exp (-( ((1/τ : ℝ) : ℂ) + Complex.I * (ω : ℂ)) * (T : ℂ)))
      Filter.atTop (nhds (0 : ℂ)) :=
  GeometricResponse.Response.Debye.tailExp_tendsto_zero hτ

example {τ ω : ℝ} (hτ : 0 < τ) :
    Filter.Tendsto
      (fun T : ℝ => GeometricResponse.Debye.truncTransformClosedForm τ ω T)
      Filter.atTop (nhds (GeometricResponse.chiDebye τ ω)) :=
  GeometricResponse.Response.Debye.truncTransform_tendsto_chiDebye hτ

end GeometricResponse.Tests


