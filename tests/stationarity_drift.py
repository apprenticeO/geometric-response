#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd


def window_indices(n: int, w: int) -> list[tuple[int, int]]:
    idx = []
    start = 0
    while start + w <= n:
        idx.append((start, start + w))
        start += w
    if not idx and w < n:
        idx.append((0, n))
    return idx


def main():
    ap = argparse.ArgumentParser(description="LWSS drift metrics: rolling mean/variance drift percentages")
    ap.add_argument("--in", dest="inp", type=str, default="./results/ou_pi.csv", help="Input CSV with t,Pi")
    ap.add_argument("--out", type=str, default="./results/stationarity.csv", help="Output CSV path")
    ap.add_argument("--window_s", type=float, default=200.0, help="Window length [s]")
    ap.add_argument("--timescale", type=float, default=1.0, help="Seconds per time unit")
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    if not {"t", "Pi"}.issubset(df.columns):
        raise ValueError("Input must contain columns: t, Pi")
    t = df["t"].to_numpy(dtype=float) * float(args.timescale)
    x = df["Pi"].to_numpy(dtype=float)
    dt = float(np.mean(np.diff(t)))
    n = len(x)
    w = max(2, int(round(args.window_s / dt)))
    idx = window_indices(n, w)
    means, vars_ = [], []
    for a, b in idx:
        seg = x[a:b]
        means.append(float(np.mean(seg)))
        vars_.append(float(np.var(seg)))
    means = np.array(means, dtype=float)
    vars_ = np.array(vars_, dtype=float)

    def drift_percent(arr: np.ndarray) -> float:
        med = float(np.median(arr))
        if med == 0:
            return float(np.nan)
        return float(100.0 * (np.max(arr) - np.min(arr)) / abs(med))

    mean_drift_pct = drift_percent(means)
    var_drift_pct = drift_percent(vars_)

    out = {
        "window_s": float(args.window_s),
        "num_windows": int(len(idx)),
        "mean_drift_pct": mean_drift_pct,
        "var_drift_pct": var_drift_pct,
        "dt_s": dt,
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    pd.DataFrame([out]).to_csv(args.out, index=False)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()


