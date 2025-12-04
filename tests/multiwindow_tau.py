#!/usr/bin/env python3
import argparse
import os
import sys
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Compute τ_G across sliding windows (CV and 95% CI)")
    ap.add_argument("--in", dest="inp", type=str, default="./results/ou_pi.csv", help="Input CSV with columns t,Pi")
    ap.add_argument("--results_dir", type=str, default="./results", help="Output directory for CSVs/plots")
    ap.add_argument("--timescale", type=float, default=1.0, help="Seconds per time unit in CSV")
    ap.add_argument("--window_s", type=float, default=400.0, help="Sliding window length [s]")
    ap.add_argument("--step_s", type=float, default=200.0, help="Window step [s]")
    ap.add_argument("--acf_cut", type=str, default="zero", help="'zero' first nonpositive; or numeric epsilon")
    ap.add_argument("--n_boot", type=int, default=1000, help="Bootstrap resamples for CI")
    args = ap.parse_args()

    os.makedirs(args.results_dir, exist_ok=True)

    # import utils from scripts
    sys.path.append(os.path.join(os.path.dirname(__file__), "..", "scripts"))
    from utils import normalized_acf, gk_tau_from_acf  # type: ignore

    df = pd.read_csv(args.inp)
    if not {"t", "Pi"}.issubset(df.columns):
        raise ValueError("Input must contain columns: t, Pi")
    t = df["t"].to_numpy(dtype=float) * float(args.timescale)
    x = df["Pi"].to_numpy(dtype=float)
    dt = float(np.mean(np.diff(t)))

    n = len(t)
    w = int(round(args.window_s / dt))
    s = int(round(args.step_s / dt))
    if w < 10 or w >= n:
        raise ValueError("window_s incompatible with record length")
    starts = np.arange(0, n - w + 1, max(1, s))

    tau_list = []
    rows = []
    for i, start in enumerate(starts):
        end = int(start + w)
        t_win = t[start:end]
        x_win = x[start:end]
        lags, acf = normalized_acf(x_win, dt)
        # truncation rule
        cut_idx = len(lags) - 1
        if args.acf_cut.strip().lower() == "zero":
            nz = np.where(acf <= 0)[0]
            if nz.size > 0:
                cut_idx = int(nz[0])
        else:
            try:
                eps = float(args.acf_cut)
                nz = np.where(np.abs(acf) <= eps)[0]
                if nz.size > 0:
                    cut_idx = int(nz[0])
            except Exception:
                pass
        tau_g = gk_tau_from_acf(lags[: cut_idx + 1], acf[: cut_idx + 1])
        tau_list.append(tau_g)
        rows.append({
            "window_idx": i,
            "t_start_s": float(t_win[0]),
            "t_end_s": float(t_win[-1]),
            "tau_G_trunc_s": float(tau_g)
        })

    perwin_path = os.path.join(args.results_dir, "multiwindow_tau.csv")
    pd.DataFrame(rows).to_csv(perwin_path, index=False)

    taus = np.array(tau_list, dtype=float)
    mean = float(np.nanmean(taus))
    std = float(np.nanstd(taus, ddof=1)) if len(taus) > 1 else np.nan
    cv = float(std / mean) if mean > 0 and np.isfinite(std) else np.nan

    # bootstrap CI
    rng = np.random.default_rng(0)
    boots = []
    for _ in range(int(args.n_boot)):
        idx = rng.integers(0, len(taus), size=len(taus))
        boots.append(np.nanmean(taus[idx]))
    lo, hi = (float(np.nanpercentile(boots, 2.5)), float(np.nanpercentile(boots, 97.5))) if boots else (np.nan, np.nan)

    summary = {
        "num_windows": int(len(taus)),
        "tau_mean_s": mean,
        "tau_std_s": std,
        "tau_cv": cv,
        "tau_ci95_lo_s": lo,
        "tau_ci95_hi_s": hi,
        "window_s": float(args.window_s),
        "step_s": float(args.step_s),
        "acf_cut": args.acf_cut,
    }
    summ_path = os.path.join(args.results_dir, "multiwindow_summary.csv")
    pd.DataFrame([summary]).to_csv(summ_path, index=False)
    print(f"Wrote {perwin_path} and {summ_path}")


if __name__ == "__main__":
    main()


