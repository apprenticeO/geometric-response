#!/usr/bin/env python3
"""
Write an analytic Debye susceptibility CSV:
  chi(ω) = 1 / (1 + i ω τ_G)
given tau_G from tau_estimates.csv. The output has columns:
  omega_rad_s, Re_chi, Im_chi
on a uniform omega grid up to the Nyquist of a supplied dt (or an explicit max_omega).
"""
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
import os


def main():
    ap = argparse.ArgumentParser(description="Write Debye chi(ω)=1/(1+iωτ) from tau_estimates.csv")
    ap.add_argument("--tau_csv", type=str, required=True, help="tau_estimates.csv with tau_G_trunc")
    ap.add_argument("--out_csv", type=str, required=True, help="Output g_eff CSV path")
    ap.add_argument("--dt", type=float, default=None, help="Sampling interval to set Nyquist ω_max=π/dt")
    ap.add_argument("--max_omega", type=float, default=None, help="Override ω_max if provided")
    ap.add_argument("--num_points", type=int, default=4096, help="Number of frequency points")
    args = ap.parse_args()

    tau_df = pd.read_csv(args.tau_csv)
    tau_g = float(tau_df.iloc[0].get("tau_G_trunc", float("nan")))
    if not np.isfinite(tau_g) or tau_g <= 0.0:
        tau_g = float(tau_df.iloc[0].get("tau_G_full", float("nan")))
    if not np.isfinite(tau_g) or tau_g <= 0.0:
        raise ValueError("No valid tau_G found")

    if args.max_omega is not None:
        w_max = float(args.max_omega)
    elif args.dt is not None and args.dt > 0:
        w_max = float(np.pi / args.dt)
    else:
        # Fallback to a conservative band
        w_max = float(10.0 / tau_g)
    w = np.linspace(0.0, w_max, int(max(32, args.num_points)))
    chi = 1.0 / (1.0 + 1j * w * tau_g)
    df = pd.DataFrame({
        "omega_rad_s": w,
        "Re_chi": np.real(chi),
        "Im_chi": np.imag(chi),
    })
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True) if os.path.dirname(args.out_csv) else None
    df.to_csv(args.out_csv, index=False)
    print(f"Saved Debye chi to {args.out_csv} (tau_G={tau_g:.6g}, ω_max={w_max:.6g})")


if __name__ == "__main__":
    main()


