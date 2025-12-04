#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import numpy as np
import pandas as pd


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description="Multi-realization variance of τ_G and τ_KK (OU)")
    ap.add_argument("--outdir", type=str, default="./results/multi_real", help="Output base directory")
    ap.add_argument("--n", type=int, default=4, help="Number of realizations")
    ap.add_argument("--tmax", type=float, default=4000.0, help="Record length [s]")
    ap.add_argument("--dt", type=float, default=0.01, help="Time step [s]")
    ap.add_argument("--tau0", type=float, default=2.0, help="OU τ0 [s]")
    ap.add_argument("--seed0", type=int, default=0, help="Base seed")
    ap.add_argument("--nperseg", type=int, default=131072, help="Welch segment length")
    ap.add_argument("--nfft_factor", type=int, default=16, help="Zero-padding factor")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    project = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    gen = os.path.join(project, "scripts", "generate_ou.py")
    comp = os.path.join(project, "scripts", "compute_gk_kk.py")

    rows = []
    for i in range(args.n):
        seed = args.seed0 + i
        run_dir = os.path.join(args.outdir, f"run_{i:02d}")
        os.makedirs(run_dir, exist_ok=True)
        csv_path = os.path.join(run_dir, "ou_pi.csv")
        png_path = os.path.join(run_dir, "ou_pi.png")

        run([
            sys.executable, gen,
            "--tmax", str(args.tmax),
            "--dt", str(args.dt),
            "--tau0", str(args.tau0),
            "--seed", str(seed),
            "--out", csv_path,
            "--plot", "--png", png_path,
        ])

        run([
            sys.executable, comp,
            "--in", csv_path,
            "--outdir", run_dir,
            "--plot",
            "--nperseg", str(args.nperseg),
            "--nfft_factor", str(args.nfft_factor),
            "--acf_cut", "zero",
        ])

        tau_df = pd.read_csv(os.path.join(run_dir, "tau_estimates.csv"))
        tau_g = float(tau_df["tau_G_trunc"].iloc[0])
        tau_kk = float(tau_df.get("tau_KK", pd.Series([np.nan])).iloc[0])
        rows.append({"run": i, "seed": seed, "tau_G_trunc_s": tau_g, "tau_KK_s": tau_kk})

    agg = pd.DataFrame(rows)
    agg.to_csv(os.path.join(args.outdir, "multi_realization_raw.csv"), index=False)

    def stats(col: str):
        vals = agg[col].to_numpy(dtype=float)
        mean = float(np.nanmean(vals))
        std = float(np.nanstd(vals, ddof=1)) if len(vals) > 1 else np.nan
        cv = float(std / mean) if mean > 0 and np.isfinite(std) else np.nan
        return mean, std, cv

    g_mean, g_std, g_cv = stats("tau_G_trunc_s")
    kk_mean, kk_std, kk_cv = stats("tau_KK_s")
    summ = pd.DataFrame([{ 
        "n": int(args.n),
        "tau_G_mean_s": g_mean, "tau_G_std_s": g_std, "tau_G_cv": g_cv,
        "tau_KK_mean_s": kk_mean, "tau_KK_std_s": kk_std, "tau_KK_cv": kk_cv,
        "tmax_s": float(args.tmax), "dt_s": float(args.dt), "tau0_s": float(args.tau0)
    }])
    summ.to_csv(os.path.join(args.outdir, "multi_realization_summary.csv"), index=False)
    print(f"Wrote multi_realization_summary.csv in {args.outdir}")


if __name__ == "__main__":
    main()


