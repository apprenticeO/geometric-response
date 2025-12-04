#!/usr/bin/env python3
"""
Synthesize tensor response channels R_<comp> from R_info_<comp> using a
first-order Debye filter with time constant tau_G:
    y[n] = alpha * y[n-1] + (1 - alpha) * x[n], alpha = exp(-dt / tau_G)

Inputs:
  - --in:  CSV with columns t and R_info_<comp> (e.g., A1,A2,...)
  - --tau_csv: CSV with tau_G_trunc (and/or tau_G_full)
Outputs:
  - Writes a new CSV with the original columns plus R_<comp> for each component
"""
from __future__ import annotations

import argparse
import os
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Synthesize tensor response R_<comp> from R_info_<comp> via Debye filter with τ_G")
    ap.add_argument("--in", dest="inp", type=str, required=True, help="Input CSV with t and R_info_<comp> columns")
    ap.add_argument("--tau_csv", type=str, required=True, help="tau_estimates.csv with tau_G_trunc")
    ap.add_argument("--out", type=str, required=True, help="Output CSV path")
    ap.add_argument("--components", type=str, default="", help="Comma-separated component names (e.g., 'A1,A2,A3'); if empty, infer from columns")
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    if "t" not in df.columns:
        raise ValueError("Input CSV must contain column 't'")
    t = df["t"].to_numpy(dtype=float)
    if len(t) < 2:
        raise ValueError("Not enough samples to estimate dt")
    dt = float(np.mean(np.diff(t)))

    # Load τ_G
    tau_df = pd.read_csv(args.tau_csv)
    tau_g = float(tau_df.iloc[0].get("tau_G_trunc", float("nan")))
    if not np.isfinite(tau_g) or tau_g <= 0.0:
        # fallback to full if trunc missing
        tau_g = float(tau_df.iloc[0].get("tau_G_full", float("nan")))
    if not np.isfinite(tau_g) or tau_g <= 0.0:
        raise ValueError("No valid tau_G found in tau_csv")

    # Determine components
    if args.components:
        comps = [c.strip() for c in args.components.split(",") if c.strip()]
    else:
        comps = []
        for c in df.columns:
            if c.startswith("R_info_"):
                comps.append(c[len("R_info_"):])
        if not comps:
            raise ValueError("No R_info_<comp> columns found and components not provided")

    # Debye filter coefficient
    alpha = float(np.exp(-dt / max(tau_g, 1e-30)))

    # For each component, synthesize response
    for comp in comps:
        col_in = f"R_info_{comp}"
        if col_in not in df.columns:
            raise ValueError(f"Missing column {col_in}")
        x = df[col_in].to_numpy(dtype=float)
        y = np.zeros_like(x, dtype=float)
        # Simple causal filter, zero initial condition
        for n in range(1, len(x)):
            y[n] = alpha * y[n - 1] + (1.0 - alpha) * x[n]
        df[f"R_{comp}"] = y

    # Write output
    os.makedirs(os.path.dirname(args.out), exist_ok=True) if os.path.dirname(args.out) else None
    df.to_csv(args.out, index=False)
    print(f"Saved synthesized tensor response CSV: {args.out} (tau_G={tau_g:.6g}s, alpha={alpha:.6g})")


if __name__ == "__main__":
    main()


