TRIAD Geometry – Scripts
========================

Minimal utilities to reproduce the GK→KK→G\_eff pipeline described in the paper
*From Microscopic Structural Interdependence to Causal Dispersive Geometry*.

Outputs follow the filenames referenced in the text:
- `acf_diagnostics.csv` – lags (s) and normalized ACF \(C_\Pi(\tau)\); L¹ / integrability diagnostics.
- `tau_estimates.csv` – \(\tau_G\), short‑time floor \(\tau_{\mathrm{bound}}\), optional \(\tau_{\mathrm{KK}}\) and passivity flag.
- `g_eff.csv` – \(\omega\), \(\Re/\Im \chi(\omega)\) and (optionally) \(\Re/\Im G_{\mathrm{eff}}(\omega)\).

## Generate OU time series

```bash
python3 scripts/generate_ou.py \
  --tmax 200 --dt 0.01 --tau0 2.0 --sigma 1.0 \
  --out results/ou_pi.csv --plot --png results/ou_pi.png
```

## Compute \(\tau_G\) (GK), ACF, optional \(\chi(\omega)\) and \(\tau_{\mathrm{KK}}\)

Input CSV must contain columns: `t,Pi` (and optionally `R,R_info` for cross‑spectral \(\chi\)).

```bash
python3 scripts/compute_gk_kk.py \
  --in results/ou_pi.csv --outdir ./results --plot
```

This writes:
- `results/acf_diagnostics.csv`
- `results/tau_estimates.csv`
- and, if `R,R_info` are present, `results/g_eff.csv` (with \(\tau_{\mathrm{KK}}\) and passivity flags in `tau_estimates.csv`).

## Fit Debye single-pole near \(\omega \approx 0\)

```bash
python3 scripts/fit_debye.py \
  --in results/g_eff.csv \
  --out results/debye_fit.csv --plot --png results/debye_fit.png
```

Fits \(\Im\chi(\omega) \approx -(\omega \tau)/(1+(\omega \tau)^2)\) in the band \(|\omega|\tau\le 0.2\) to obtain \(\tau_{\mathrm{fit}}\).

Notes:
- One-sided PSD convention \(S(0)=2\int_0^\infty C_\Pi(\tau)\,d\tau\) is used;
  \(\tau_G = S(0)/(2\,\mathrm{Var}\,\Pi)\).
- Passivity check: \(\Im\chi(\omega)\le 0\) near \(\omega>0\).
- No virtualenv assumed; run with `python3`.

Dimensional notes (alignment with the paper):
- Use \(\ell_G := c\,\tau_G\) when forming dimensionless operator combinations with \(\Box\) in SI,
  e.g. \(G_{\mathrm{eff}}(\Box)=G_0 (1+\ell_G^2 \Box)^{-1}\).
- The \(R_{\mu\nu}\) diagnostic (`export_rmunu_diagnostic.py`) accepts `--kappa_star` [m\(^{-2}\)]
  so that \((\hbar/c^5)\,\Pi^2\,\kappa_\ast\) completes curvature units; default is 1.0 for relative scaling.

## Phase 1 extras (reliability)

- \(\tau_{\mathrm{PSD}}\) and \(\Re\chi\) flatness:

```bash
python3 tests/compute_flatness_and_tau_psd.py --results_dir ./results --plot
```

- Multi-window \(\tau_G\) CV/CI:

```bash
python3 tests/multiwindow_tau.py \
  --in results/ou_pi.csv --results_dir ./results \
  --window_s 400 --step_s 200
```

- Multi-realization stats:

```bash
python3 tests/multi_realization_stats.py \
  --outdir results/multi_real --n 4 \
  --tmax 4000 --dt 0.01 --tau0 2.0
```

- Finite-record bias:

```bash
python3 tests/finite_record_bias.py \
  --outdir results/finite_bias --lengths_s 500 1000 2000 4000
```

## Phase 2 (flatness and violations)

- Single-pole dataset with \(R_{\mathrm{info}}\) and \(R\) (\(\chi\) with real + imaginary parts):

```bash
python3 scripts/generate_lindblad_qubit.py \
  --tmax 4000 --dt 0.01 --tau 2.0 \
  --out results/lindblad/lindblad.csv

python3 scripts/compute_gk_kk.py \
  --in results/lindblad/lindblad.csv \
  --outdir results/lindblad --plot \
  --nperseg 131072 --nfft_factor 16

python3 tests/compute_flatness_and_tau_psd.py \
  --results_dir results/lindblad --plot
```

- Violation suite (non-passive, quench):

```bash
python3 tests/violation_cases.py \
  --outdir results/violation \
  --tmax 2000 --dt 0.01 \
  --tau 2.0 --tau1 2.0 --tau2 0.2 --t_quench 1000
```

Outputs a `violation_summary.csv` with \(\tau_G\), \(\tau_{\mathrm{KK}}\), passivity, and stationarity metrics.


