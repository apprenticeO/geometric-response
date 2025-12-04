#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import pandas as pd
import subprocess


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def generate_quench_pi(t: np.ndarray, tau1: float, tau2: float, t_quench: float, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    dt = float(np.mean(np.diff(t)))
    pi = np.zeros_like(t)
    a1 = float(np.exp(-dt / tau1))
    a2 = float(np.exp(-dt / tau2))
    n_quench = int(round(t_quench / dt))
    for n in range(1, len(t)):
        a = a1 if n < n_quench else a2
        pi[n] = a * pi[n - 1] + np.sqrt(1 - a * a) * rng.normal()
    return pi


def main():
    ap = argparse.ArgumentParser(description="Generate and analyze violation cases (V-1 non-passive, V-2 quench)")
    ap.add_argument("--outdir", type=str, default="./results/violation", help="Base output directory")
    ap.add_argument("--tmax", type=float, default=2000.0, help="Total time (s)")
    ap.add_argument("--dt", type=float, default=0.01, help="Time step (s)")
    ap.add_argument("--tau", type=float, default=2.0, help="Time constant for non-passive base (s)")
    ap.add_argument("--tau1", type=float, default=2.0, help="Pre-quench tau (s)")
    ap.add_argument("--tau2", type=float, default=0.2, help="Post-quench tau (s)")
    ap.add_argument("--t_quench", type=float, default=1000.0, help="Quench time (s)")
    ap.add_argument("--seed", type=int, default=123, help="Seed")
    ap.add_argument("--nperseg", type=int, default=65536, help="Welch segment length")
    ap.add_argument("--nfft_factor", type=int, default=8, help="Zero-padding factor")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    project = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    gen_lind = os.path.join(project, "scripts", "generate_lindblad_qubit.py")
    comp = os.path.join(project, "scripts", "compute_gk_kk.py")
    stat = os.path.join(project, "tests", "stationarity_drift.py")

    # V-1: Non-passive case (flip sign of response)
    t = np.arange(0.0, args.tmax + 1e-12, args.dt, dtype=float)
    # Base single-pole
    u = np.random.default_rng(args.seed).normal(0.0, 1.0, size=t.shape[0])
    a = float(np.exp(-args.dt / args.tau))
    y = np.zeros_like(u)
    for n in range(1, len(u)):
        y[n] = a * y[n - 1] + (1.0 - a) * u[n]
    # Non-passive: flip the sign of output (gain < 0)
    y_np = -y
    df_np = pd.DataFrame({"t": t, "Pi": y_np, "R_info": u, "R": y_np})
    out_np = os.path.join(args.outdir, "nonpassive")
    os.makedirs(out_np, exist_ok=True)
    csv_np = os.path.join(out_np, "data.csv")
    df_np.to_csv(csv_np, index=False)
    run([sys.executable, comp, "--in", csv_np, "--outdir", out_np, "--nperseg", str(args.nperseg), "--nfft_factor", str(args.nfft_factor), "--acf_cut", "zero", "--plot"])

    # V-2: Quench (non-stationary)
    pi_q = generate_quench_pi(t, args.tau1, args.tau2, args.t_quench, seed=args.seed + 1)
    df_q = pd.DataFrame({"t": t, "Pi": pi_q})
    out_q = os.path.join(args.outdir, "quench")
    os.makedirs(out_q, exist_ok=True)
    csv_q = os.path.join(out_q, "data.csv")
    df_q.to_csv(csv_q, index=False)
    run([sys.executable, comp, "--in", csv_q, "--outdir", out_q, "--nperseg", str(args.nperseg), "--nfft_factor", str(args.nfft_factor), "--acf_cut", "zero", "--plot"])
    run([sys.executable, stat, "--in", csv_q, "--out", os.path.join(out_q, "stationarity.csv"), "--window_s", "200"]) 

    # Summaries
    def read_tau(path: str) -> dict:
        td = pd.read_csv(os.path.join(path, "tau_estimates.csv"))
        return {
            "tau_G_trunc_s": float(td["tau_G_trunc"].iloc[0]),
            "tau_KK_s": float(td.get("tau_KK", pd.Series([np.nan])).iloc[0]) if os.path.exists(os.path.join(path, "g_eff.csv")) else np.nan,
            "passivity_ok": bool(td.get("passivity_ok", pd.Series([np.nan])).iloc[0]) if os.path.exists(os.path.join(path, "g_eff.csv")) else np.nan,
        }

    summ_np = read_tau(out_np)
    summ_q = read_tau(out_q)
    # Stationarity metric for quench
    st_q = pd.read_csv(os.path.join(out_q, "stationarity.csv"))
    var_drift_q = float(st_q["var_drift_pct"].iloc[0])

    rows = [
        {"case": "Non-passive", "tau_G_trunc_s": summ_np["tau_G_trunc_s"], "tau_KK_s": summ_np["tau_KK_s"], "passivity_ok": summ_np["passivity_ok"], "var_drift_pct": np.nan},
        {"case": "Quench", "tau_G_trunc_s": summ_q["tau_G_trunc_s"], "tau_KK_s": summ_q["tau_KK_s"], "passivity_ok": summ_q["passivity_ok"], "var_drift_pct": var_drift_q},
    ]
    pd.DataFrame(rows).to_csv(os.path.join(args.outdir, "violation_summary.csv"), index=False)
    print(f"Wrote violation_summary.csv in {args.outdir}")

    # Assemble a 2x2 summary plot from existing images
    try:
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        fig, axs = plt.subplots(2, 2, figsize=(9, 6))
        for ax in axs.ravel():
            ax.axis('off')
        # non-passive chi_im and acf
        np_chi = os.path.join(out_np, 'chi_im.png')
        np_acf = os.path.join(out_np, 'acf.png')
        # quench chi_im and acf
        q_chi = os.path.join(out_q, 'chi_im.png')
        q_acf = os.path.join(out_q, 'acf.png')
        if os.path.exists(np_chi):
            axs[0, 0].imshow(mpimg.imread(np_chi)); axs[0, 0].set_title('Non-passive: Im χ(ω)')
        if os.path.exists(q_chi):
            axs[0, 1].imshow(mpimg.imread(q_chi)); axs[0, 1].set_title('Quench: Im χ(ω)')
        if os.path.exists(np_acf):
            axs[1, 0].imshow(mpimg.imread(np_acf)); axs[1, 0].set_title('Non-passive: ACF')
        if os.path.exists(q_acf):
            axs[1, 1].imshow(mpimg.imread(q_acf)); axs[1, 1].set_title('Quench: ACF')
        plt.tight_layout()
        summary_png = os.path.join(args.outdir, 'violation_plots.png')
        plt.savefig(summary_png, dpi=150)
        plt.close()
        print(f"Saved {summary_png}")
    except Exception as e:
        print(f"Violation summary plot skipped: {e}")


if __name__ == "__main__":
    main()


