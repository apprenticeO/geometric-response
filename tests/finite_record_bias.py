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
    ap = argparse.ArgumentParser(description="Finite-record bias: τ_G vs L, fit τ(L)=τ_inf(1-c/L)")
    ap.add_argument("--outdir", type=str, default="./results/finite_bias", help="Output directory")
    ap.add_argument("--lengths_s", type=float, nargs="*", default=[500.0, 1000.0, 2000.0, 4000.0], help="Record lengths [s]")
    ap.add_argument("--dt", type=float, default=0.01, help="Time step [s]")
    ap.add_argument("--tau0", type=float, default=2.0, help="OU τ0 [s]")
    ap.add_argument("--seed", type=int, default=0, help="Seed for the long record")
    ap.add_argument("--nperseg", type=int, default=131072, help="Welch segment length")
    ap.add_argument("--nfft_factor", type=int, default=16, help="Zero-padding factor")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    project = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    gen = os.path.join(project, "scripts", "generate_ou.py")
    comp = os.path.join(project, "scripts", "compute_gk_kk.py")

    Lmax = float(max(args.lengths_s))
    long_csv = os.path.join(args.outdir, "ou_long.csv")
    long_png = os.path.join(args.outdir, "ou_long.png")
    run([
        sys.executable, gen,
        "--tmax", str(Lmax), "--dt", str(args.dt), "--tau0", str(args.tau0),
        "--seed", str(args.seed), "--out", long_csv, "--plot", "--png", long_png,
    ])

    df = pd.read_csv(long_csv)
    t = df["t"].to_numpy(dtype=float)
    dt = float(np.mean(np.diff(t)))

    rows = []
    for L in args.lengths_s:
        n = int(round(L / dt)) + 1
        sub = df.iloc[:n, :].copy()
        sub_csv = os.path.join(args.outdir, f"ou_L{int(L)}.csv")
        sub.to_csv(sub_csv, index=False)
        run([
            sys.executable, comp,
            "--in", sub_csv,
            "--outdir", args.outdir,
            "--nperseg", str(args.nperseg),
            "--nfft_factor", str(args.nfft_factor),
            "--acf_cut", "zero",
        ])
        tau_df = pd.read_csv(os.path.join(args.outdir, "tau_estimates.csv"))
        tau_g = float(tau_df["tau_G_trunc"].iloc[0])
        rows.append({"L_s": float(L), "tau_G_trunc_s": tau_g})

    tab = pd.DataFrame(rows).sort_values("L_s")
    tab.to_csv(os.path.join(args.outdir, "finite_bias_raw.csv"), index=False)

    # Fit τ(L) ≈ τ_inf (1 - c/L) ⇒ τ(L) ≈ a + b*(1/L), with a=τ_inf, b= -τ_inf*c
    L = tab["L_s"].to_numpy(dtype=float)
    invL = 1.0 / L
    tau = tab["tau_G_trunc_s"].to_numpy(dtype=float)
    A = np.vstack([np.ones_like(invL), invL]).T
    coef, _, _, _ = np.linalg.lstsq(A, tau, rcond=None)
    tau_inf = float(coef[0])
    b = float(coef[1])
    c = -b / tau_inf if tau_inf != 0 else np.nan
    pd.DataFrame([{ "tau_inf_s": tau_inf, "c_bias": c }]).to_csv(
        os.path.join(args.outdir, "finite_bias_fit.csv"), index=False
    )
    print(f"Wrote finite_bias_raw.csv and finite_bias_fit.csv in {args.outdir}")


if __name__ == "__main__":
    main()


