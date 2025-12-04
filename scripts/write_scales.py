#!/usr/bin/env python3
import argparse
import math
import os
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Compute operational scales (omega_G, lambda_G, m_G) and a pipeline summary from results CSVs")
    ap.add_argument("--results_dir", type=str, default="./results", help="Directory containing tau_estimates.csv, debye_fit.csv")
    ap.add_argument("--scales_out", type=str, default=None, help="Output CSV for scales (default: results/scales.csv)")
    ap.add_argument("--summary_out", type=str, default=None, help="Output CSV for pipeline summary (default: results/pipeline_check.csv)")
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    tau_path = os.path.join(results_dir, "tau_estimates.csv")
    debye_path = os.path.join(results_dir, "debye_fit.csv")
    os.makedirs(results_dir, exist_ok=True)

    df_tau = pd.read_csv(tau_path)
    # Expect single-row file
    row = df_tau.iloc[0]
    tau_g = float(row.get("tau_G_trunc", float("nan")))
    tau_kk = float(row.get("tau_KK", float("nan")))
    passivity_ok = bool(row.get("passivity_ok", False))

    # Physical constants
    c = 299_792_458.0  # m/s
    hbar = 1.054_571_817e-34  # J*s

    # Scales from Section 2
    omega_g = float("nan")
    lambda_g = float("nan")
    m_g = float("nan")
    if math.isfinite(tau_g) and tau_g > 0:
        omega_g = 1.0 / tau_g
        lambda_g = 2.0 * math.pi * c * tau_g
        m_g = hbar / (c * c * tau_g)
        ell_g = c * tau_g
    else:
        ell_g = float("nan")

    scales_out = args.scales_out or os.path.join(results_dir, "scales.csv")
    pd.DataFrame([
        {
            "tau_G_trunc_s": tau_g,
            "omega_G_rad_per_s": omega_g,
            "lambda_G_m": lambda_g,
            "m_G_kg": m_g,
            "ell_G_m": ell_g,
        }
    ]).to_csv(scales_out, index=False)

    # Pipeline summary: τ_G, τ_KK, τ_fit, passivity, L1
    tau_fit = float("nan")
    if os.path.exists(debye_path):
        try:
            df_fit = pd.read_csv(debye_path)
            tau_fit = float(df_fit.iloc[0].get("tau_fit", float("nan")))
        except Exception:
            pass

    l1_abs = float(row.get("l1_abs", float("nan")))
    summary_out = args.summary_out or os.path.join(results_dir, "pipeline_check.csv")
    pd.DataFrame([
        {
            "tau_G_trunc_s": tau_g,
            "tau_KK_s": tau_kk,
            "tau_fit_s": tau_fit,
            "passivity_ok": passivity_ok,
            "l1_abs": l1_abs,
        }
    ]).to_csv(summary_out, index=False)

    print(f"Wrote {scales_out} and {summary_out}")


if __name__ == "__main__":
    main()


