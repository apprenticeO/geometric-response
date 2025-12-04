#!/usr/bin/env python3
"""
R_{μν} ∝ (ħ/c^5) Π^2 diagnostic exporter (DC scaling aid).

Reads a time series CSV with columns t, Pi and produces:
  - rmunu_diagnostic.csv with columns: t, Pi, Pi2_s^-2, kappa_star_m^-2, rmunu_scalar_SI
  - rmunu_timeseries.png plotting Π^2(t) and the scaled rmunu_scalar(t)

This is a diagnostic for the geometric→gravitational mapping visibility.
It is NOT a universal calibration of G; use the Planck-time route for that.
"""
import argparse
import os
import numpy as np
import pandas as pd

HBAR = 1.054_571_817e-34  # J*s
C = 299_792_458.0         # m/s


def main():
    ap = argparse.ArgumentParser(description="Export R_{μν} ∝ (ħ/c^5) Π^2 diagnostic CSV and plot")
    ap.add_argument("--in", dest="inp", type=str, required=True, help="Input CSV with columns t,Pi")
    ap.add_argument("--outdir", type=str, required=True, help="Output directory")
    ap.add_argument("--plot", action="store_true", help="Write rmunu_timeseries.png")
    ap.add_argument("--kappa_star", type=float, default=1.0, help="Curvature scale κ_* [m^-2] to complete units to R (default 1.0)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.inp)
    if not {"t", "Pi"}.issubset(df.columns):
        raise ValueError("Input must contain columns: t, Pi")
    t = df["t"].to_numpy(dtype=float)
    pi = df["Pi"].to_numpy(dtype=float)
    pi2 = pi * pi
    # Proportional scalar (units combine to SI via ħ/c^5 factor) and curvature scale κ_*
    kappa_star = float(args.kappa_star)
    rmunu_scalar = (HBAR / (C ** 5)) * pi2 * kappa_star

    out_df = pd.DataFrame({
        "t": t,
        "Pi": pi,
        "Pi2_s^-2": pi2,
        "kappa_star_m^-2": kappa_star,
        "rmunu_scalar_SI": rmunu_scalar,
    })
    out_csv = os.path.join(args.outdir, "rmunu_diagnostic.csv")
    out_df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            fig, ax1 = plt.subplots(figsize=(8, 3))
            ax1.plot(t, pi2, lw=0.9, label="Π²(t) [s⁻²]")
            ax1.set_xlabel("t [s]")
            ax1.set_ylabel("Π² [s⁻²]")
            ax2 = ax1.twinx()
            ax2.plot(t, rmunu_scalar, lw=0.9, color="tab:red", alpha=0.8, label="R_scalar ∝ (ħ/c⁵)Π²")
            ax2.set_ylabel("R_scalar [SI]")
            fig.tight_layout()
            out_png = os.path.join(args.outdir, "rmunu_timeseries.png")
            plt.savefig(out_png, dpi=150)
            plt.close(fig)
            print(f"Wrote {out_png}")
        except Exception as e:
            print(f"Plot skipped: {e}")


if __name__ == "__main__":
    main()


