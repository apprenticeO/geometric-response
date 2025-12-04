import GeometricResponseLean.GeometricResponse
import GeometricResponseLean.Response.KK

/-!
Lightweight tests/examples for KK low-frequency slope wrappers.
-/

namespace GeometricResponse.Tests.KK

open GeometricResponse
open GeometricResponse.Response

/-- Low-ω slope wrapper: `∂_ω Im G_eff(ω)|₀ = - G0 · τ` (limit form). -/
example (G0 τ : ℝ) :
    Filter.Tendsto (fun ω : ℝ => G0 * (- τ / (1 + (ω * τ) ^ 2)))
      (nhds 0) (nhds (- G0 * τ)) :=
  KK.slope_lowOmega_Geff G0 τ

end GeometricResponse.Tests.KK


