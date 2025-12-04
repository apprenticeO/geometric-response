#!/usr/bin/env python3
"""
Role: DC/adiabatic G0 calibration (primary: Planck-time route).

Use this script to set the low-frequency (ω→0) gravitational coupling G0 in the
geometric-response framework. The intended, theory-native normalization is the
Planck-time route:
  - τ_G = t_P * sqrt(P), with structural factor P (dimensionless)
  - G0 → G_Newton in the adiabatic limit (ω→0)

An optional nominal-from-Π² path is retained only as a diagnostic cross-check
on datasets where ⟨Π²⟩ can be estimated; it is not a universal calibration.
"""
import argparse
import os
import math
import numpy as np
import pandas as pd


def robust_mean_square(x: np.ndarray, trim_frac: float = 0.0) -> float:
    """Return mean(x^2) with optional symmetric trimming by fraction."""
    x = np.asarray(x, dtype=float)
    if trim_frac > 0:
        q_lo = np.quantile(x, trim_frac / 2.0)
        q_hi = np.quantile(x, 1.0 - trim_frac / 2.0)
        m = (x >= q_lo) & (x <= q_hi)
        x = x[m]
    return float(np.mean(x * x))


def main():
    ap = argparse.ArgumentParser(
        description="G0 calibration (DC/adiabatic). Primary: Planck-time (τ_G=t_P√P ⇒ G0=G_N). Optional: nominal-from-Π² diagnostic."
    )
    ap.add_argument("--results_dir", type=str, default="./results/two_qubit_pipeline", help="Directory with tau_estimates.csv and (optionally) two_qubit.csv")
    ap.add_argument("--pi_csv", type=str, default=None, help="CSV with columns t,Pi to compute ⟨Π^2⟩ empirically (default: use two_qubit.csv in results_dir)")
    ap.add_argument("--trim_frac", type=float, default=0.0, help="Symmetric trimming fraction for robust ⟨Π^2⟩ (e.g., 0.02)")
    ap.add_argument("--write_notice", action="store_true", help="Write a NOTICE explaining DC/adiabatic meaning and toy-model caveat")
    ap.add_argument("--planck_time", action="store_true", help="Planck-time route: set τ_G = t_P * sqrt(P) and G0 -> G_N (DC limit)")
    ap.add_argument("--P", type=float, default=1.0, help="Structural factor P=δ ε_S ε_I (dimensionless) for Planck-time route")
    ap.add_argument("--out", type=str, default=None, help="Output CSV path (default: results_dir/calibrate_g0.csv)")
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    os.makedirs(results_dir, exist_ok=True)

    # Physical constants (SI)
    c = 299_792_458.0
    hbar = 1.054_571_817e-34
    G_Newton = 6.67430e-11
    t_P = math.sqrt(hbar * G_Newton / (c ** 5))

    # Read τ metrics (optional)
    tau_path = os.path.join(results_dir, "tau_estimates.csv")
    tau_G = float("nan")
    if os.path.exists(tau_path):
        try:
            td = pd.read_csv(tau_path)
            tau_G = float(td.iloc[0].get("tau_G_trunc", float("nan")))
        except Exception:
            pass

    # Planck-time normalization route: DC limit G0 -> G_N, τ_G = t_P * sqrt(P)
    if args.planck_time:
        P_val = max(0.0, float(args.P))
        tauG_planck = t_P * math.sqrt(P_val)
        out_path = args.out or os.path.join(results_dir, "calibrate_g0.csv")
        pd.DataFrame([{
            "route": "planck_time",
            "P_dimless": P_val,
            "tau_G_trunc_s": tau_G,
            "tau_G_planck_s": tauG_planck,
            "G0_SI": G_Newton,
            "note": "Planck-time route (DC/adiabatic): τ_G=t_P*sqrt(P), G0→G_N."
        }]).to_csv(out_path, index=False)
        print(f"Wrote {out_path}")
        if args.write_notice:
            try:
                notice = (
                    "NOTICE (Planck-time route):\n"
                    "- Adiabatic/DC identification: τ_G = t_P * sqrt(P) and G0 → G_N.\n"
                    "- Use only in the ω→0 limit; microscopic demos are not absolute G0 calibrations.\n"
                )
                with open(os.path.join(results_dir, "NOTICE_G0_calibration.txt"), "w") as f:
                    f.write(notice)
            except Exception:
                pass
        return

    # Otherwise: nominal DC estimate from ⟨Π^2⟩ (P≈1 proxy)
    pi_csv = args.pi_csv or os.path.join(results_dir, "two_qubit.csv")
    if not os.path.exists(pi_csv):
        raise FileNotFoundError(f"Cannot find Pi CSV at {pi_csv}")
    df = pd.read_csv(pi_csv)
    if "Pi" not in df.columns:
        raise ValueError("Input CSV must contain column 'Pi'")
    Pi = df["Pi"].to_numpy(dtype=float)

    # Compute ⟨Π^2⟩ with optional trimming
    mean_Pi2 = robust_mean_square(Pi, trim_frac=max(0.0, float(args.trim_frac)))

    # Nominal calibration using R ~ (ħ/c^5) Π^2 P with P≈1 proxy:
    #    G0^{-1} ~ (ħ/c^5) ⟨Π^2⟩  ⇒  G0_nominal ~ c^5 / (ħ ⟨Π^2⟩)
    if mean_Pi2 <= 0 or not math.isfinite(mean_Pi2):
        raise ValueError("Non-positive or invalid ⟨Π^2⟩; cannot calibrate G0.")
    G0_nominal = (c ** 5) / (hbar * mean_Pi2)
    ratio_to_Newton = G0_nominal / G_Newton

    out_path = args.out or os.path.join(results_dir, "calibrate_g0.csv")
    pd.DataFrame([{
        "tau_G_trunc_s": tau_G,
        "mean_Pi2_s^-2": mean_Pi2,
        "G0_nominal_SI": G0_nominal,
        "G_Newton_SI": G_Newton,
        "ratio_G0_over_G_Newton": ratio_to_Newton,
        "note": "Adiabatic/DC quantity; nominal P≈1. Use Planck-time or realistic matter DC limit; toy models are not absolute calibrations."
    }]).to_csv(out_path, index=False)
    print(f"Wrote {out_path}")

    if args.write_notice:
        try:
            notice = (
                "NOTICE (G0 calibration):\n"
                "- G0 is an adiabatic/DC stiffness (ω→0). This script computes a nominal value from ⟨Π^2⟩ using G0≈c^5/(ħ⟨Π^2⟩) with P≈1.\n"
                "- Do NOT interpret toy microscopic datasets (e.g., two-qubit) as absolute G0 calibrations; they are only for τ_G/KK/passivity checks.\n"
                "- To match Newton’s G in the IR, use the Planck-time route (τ_G≈t_P√P) or supply realistic matter inputs and take the DC limit.\n"
            )
            with open(os.path.join(results_dir, "NOTICE_G0_calibration.txt"), "w") as f:
                f.write(notice)
        except Exception:
            pass


if __name__ == "__main__":
    main()


