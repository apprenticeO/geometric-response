#!/usr/bin/env python3
"""
Two-qubit LWSS and passivity checks (script-style test).

- LWSS: rolling-window wide-sense stationarity checks on Π(t):
  * window means should not drift beyond mean_tol_std * global std
  * window variances should stay within var_tol_frac of global variance
- Passivity: require passivity_ok=True from compute_gk_kk outputs.

Usage:
  python3 scripts/test_two_qubit_lwss_passivity.py \
      --results_dir results/two_qubit_pipeline \
      --window_s 20 --step_s 10 \
      --mean_tol_std 0.25 --var_tol_frac 0.5
"""
import argparse
import os
import sys
import math
import numpy as np
import pandas as pd


def rolling_stats(t: np.ndarray, x: np.ndarray, window_s: float, step_s: float) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    dt = float(np.mean(np.diff(t)))
    n = len(x)
    if n < 3:
        return np.array([]), np.array([])
    win = max(2, int(round(window_s / dt)))
    step = max(1, int(round(step_s / dt)))
    means = []
    vars_ = []
    for start in range(0, n - win + 1, step):
        seg = x[start:start + win]
        means.append(float(np.mean(seg)))
        vars_.append(float(np.var(seg)))
    return np.array(means, dtype=float), np.array(vars_, dtype=float)


def main():
    ap = argparse.ArgumentParser(description="LWSS and passivity checks for two-qubit pipeline outputs")
    ap.add_argument("--results_dir", type=str, default="results/two_qubit_pipeline", help="Directory with two_qubit.csv and tau_estimates.csv")
    ap.add_argument("--window_s", type=float, default=20.0, help="Rolling window size [s] for LWSS")
    ap.add_argument("--step_s", type=float, default=None, help="Rolling step [s] (default: window/2)")
    ap.add_argument("--mean_tol_std", type=float, default=0.25, help="Max |window_mean - global_mean| in units of global std")
    ap.add_argument("--var_tol_frac", type=float, default=0.5, help="Max fractional deviation |var_win - var_global|/var_global")
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    pi_csv = os.path.join(results_dir, "two_qubit.csv")
    tau_csv = os.path.join(results_dir, "tau_estimates.csv")

    if not os.path.exists(pi_csv):
        print(f"FAIL: missing {pi_csv}; run the pipeline first.", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(tau_csv):
        print(f"FAIL: missing {tau_csv}; run compute_gk_kk first.", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(pi_csv)
    if not {"t", "Pi"}.issubset(df.columns):
        print("FAIL: two_qubit.csv must contain columns t, Pi", file=sys.stderr)
        sys.exit(1)
    t = df["t"].to_numpy(dtype=float)
    pi = df["Pi"].to_numpy(dtype=float)
    pi = np.asarray(pi, dtype=float)
    if np.allclose(np.var(pi), 0.0):
        print("FAIL: Pi variance ~ 0; cannot test LWSS.", file=sys.stderr)
        sys.exit(1)

    # LWSS checks
    global_mean = float(np.mean(pi))
    global_std = float(np.std(pi))
    global_var = float(np.var(pi))
    step_s = (args.window_s / 2.0) if args.step_s is None else float(args.step_s)
    means, vars_ = rolling_stats(t, pi, window_s=float(args.window_s), step_s=step_s)
    if means.size == 0 or vars_.size == 0:
        print("FAIL: insufficient data for rolling LWSS checks.", file=sys.stderr)
        sys.exit(1)
    mean_dev_std = float(np.max(np.abs(means - global_mean))) / max(global_std, 1e-24)
    var_dev_frac = float(np.max(np.abs(vars_ - global_var))) / max(global_var, 1e-24)
    lwss_ok = (mean_dev_std <= float(args.mean_tol_std)) and (var_dev_frac <= float(args.var_tol_frac))

    # Passivity from KK pipeline
    tau_df = pd.read_csv(tau_csv)
    passivity_ok = bool(tau_df.iloc[0].get("passivity_ok", False))

    # Report
    if lwss_ok and passivity_ok:
        print(f"OK: LWSS and passivity satisfied. mean_dev_std={mean_dev_std:.3g}, var_dev_frac={var_dev_frac:.3g}, passivity_ok={passivity_ok}")
        sys.exit(0)
    else:
        if not lwss_ok:
            print(f"FAIL: LWSS check failed. mean_dev_std={mean_dev_std:.3g} (tol {args.mean_tol_std}), var_dev_frac={var_dev_frac:.3g} (tol {args.var_tol_frac})", file=sys.stderr)
        if not passivity_ok:
            print("FAIL: passivity_ok=False in tau_estimates.csv; adjust drive/response or check passivity_sign.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()


