Geometric Response scripts
==========================

Planck‑time DC normalization (primary)
--------------------------------------
- Use `calibrate_G0.py` to set the low‑frequency coupling:
  - Planck‑time route (intended): `python3 scripts/calibrate_G0.py --planck_time --P <P>`
  - This identifies τ_G = t_P √P and sets G0 → G_N in the adiabatic limit (ω→0).
- The nominal‑from‑Π² branch is retained only as a diagnostic on datasets with ⟨Π²⟩; it is not a universal calibration.

Microscopic validation (two‑qubit)
----------------------------------
- Use the two‑qubit demo and pipeline to validate the machinery: Π → C_Π → τ_G (GK) → χ(ω) (KK) and passivity.
- Do not interpret microscopic runs as absolute G0 calibrations.

DC scaling cross‑check (experimental)
-------------------------------------
- `scripts/experimental/pi2_scaling_estimator.py` estimates Π² for a matter model and reports:
  - A diagnostic `G0_nominal ≈ c^5/(ħ Π^2)` (P≈1) and, if provided, a triad‑formula `G0` using ε_F and P.
- This is a DC/adiabatic diagnostic only; do not use it to set the universal Newton constant.

# TRIAD Geometry – Scripts

Minimal utilities to reproduce the GK→KK→Geff pipeline described in the paper.

Outputs follow the filenames referenced in the text:
- `acf_diagnostics.csv` – lags (s) and normalized ACF C(τ); l1_abs integrability diagnostic
- `tau_estimates.csv` – τ_G, τ_bound (short-time floor), optional τ_KK and passivity flag
- `g_eff.csv` – ω, Re/Im χ(ω); Re/Im Geff(ω) (scaled by G0 externally)

## Generate OU time series
```bash
python3 scripts/generate_ou.py --tmax 200 --dt 0.01 --tau0 2.0 --sigma 1.0 --out results/ou_pi.csv --plot --png results/ou_pi.png
```

## Compute τ_G (GK), ACF, optional χ(ω) and τ_KK
Input CSV must contain columns: `t,Pi` (and optionally `R,R_info` for cross-spectral χ).
```bash
python3 scripts/compute_gk_kk.py --in results/ou_pi.csv --outdir ./results --plot
```
Writes `results/acf_diagnostics.csv`, `results/tau_estimates.csv`, and if `R,R_info` present also `results/g_eff.csv` (with τ_KK and passivity in `tau_estimates.csv`).

## Fit Debye single-pole near ω≈0
```bash
python3 scripts/fit_debye.py --in results/g_eff.csv --out results/debye_fit.csv --plot --png results/debye_fit.png
```
Fits Im χ(ω) ≈ −(ω τ)/(1+(ω τ)^2) in the band |ω|τ≤0.2 to get τ_fit.

Notes:
- One-sided PSD convention S(0)=2∫_0^∞ C(τ)dτ is used; τ_G=S(0)/(2 Var Π).
- Passivity check: Im χ(ω)≤0 near ω>0.
- No virtualenv assumed; run with `python3`.

Dimensional notes (scripts alignment with the paper):
- Use ℓ_G:=c τ_G when forming dimensionless operator combinations with □ in SI; e.g., G_eff(□)=G_0 (1+ℓ_G^2 □)^(-1).
- The R_{μν} diagnostic (`export_rmunu_diagnostic.py`) now accepts `--kappa_star` [m^-2] so that (ħ/c^5) Π^2 κ_* completes curvature units; default is 1.0 for relative scaling.

## Phase 1 extras (reliability)
- τ_PSD and Re χ flatness (if available):
```bash
python3 tests/compute_flatness_and_tau_psd.py --results_dir ./results --plot
```
- Multi-window τ_G CV/CI:
```bash
python3 tests/multiwindow_tau.py --in results/ou_pi.csv --results_dir ./results --window_s 400 --step_s 200
```
- Multi-realization stats:
```bash
python3 tests/multi_realization_stats.py --outdir results/multi_real --n 4 --tmax 4000 --dt 0.01 --tau0 2.0
```
- Finite-record bias:
```bash
python3 tests/finite_record_bias.py --outdir results/finite_bias --lengths_s 500 1000 2000 4000
```

## Phase 2 (flatness and violations)
- Single-pole dataset with R_info and R (χ with real + imaginary parts):
```bash
python3 scripts/generate_lindblad_qubit.py --tmax 4000 --dt 0.01 --tau 2.0 --out results/lindblad/lindblad.csv
python3 scripts/compute_gk_kk.py --in results/lindblad/lindblad.csv --outdir results/lindblad --plot --nperseg 131072 --nfft_factor 16
python3 tests/compute_flatness_and_tau_psd.py --results_dir results/lindblad --plot
```
- Violation suite (non-passive, quench):
```bash
python3 tests/violation_cases.py --outdir results/violation --tmax 2000 --dt 0.01 --tau 2.0 --tau1 2.0 --tau2 0.2 --t_quench 1000
```
Outputs a `violation_summary.csv` with τ_G, τ_KK, passivity, and stationarity metrics.


