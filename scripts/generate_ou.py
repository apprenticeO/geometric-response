#!/usr/bin/env python3
import argparse, os
import numpy as np
import pandas as pd


def generate_ou(t: np.ndarray, tau0: float, sigma: float, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    dt = np.mean(np.diff(t))
    alpha = np.exp(-dt / tau0)
    x = np.zeros_like(t)
    # stationary variance of OU: sigma^2 * tau0 / 2 (for dX = -X/tau0 dt + sigma dW)
    var_stat = sigma * sigma * tau0 / 2.0
    x[0] = rng.normal(0.0, np.sqrt(var_stat))
    for n in range(1, len(t)):
        x[n] = alpha * x[n - 1] + np.sqrt(var_stat * (1 - alpha * alpha)) * rng.normal()
    return x


def main():
    ap = argparse.ArgumentParser(description="Generate OU triad time series Pi(t)")
    ap.add_argument("--tmax", type=float, default=100.0, help="Total time (s)")
    ap.add_argument("--dt", type=float, default=0.01, help="Time step (s)")
    ap.add_argument("--tau0", type=float, default=2.0, help="OU correlation time tau0 (s)")
    ap.add_argument("--sigma", type=float, default=1.0, help="OU noise strength (sqrt power)")
    ap.add_argument("--seed", type=int, default=0, help="Random seed")
    ap.add_argument("--out", type=str, default="results/ou_pi.csv", help="Output CSV path")
    ap.add_argument("--plot", action="store_true", help="Save a PNG plot of Pi(t)")
    ap.add_argument("--png", type=str, default="results/ou_pi.png", help="Output PNG path")
    args = ap.parse_args()

    t = np.arange(0.0, args.tmax + 1e-12, args.dt)
    pi = generate_ou(t, args.tau0, args.sigma, args.seed)
    # ensure output directories
    out_dir = os.path.dirname(args.out) or "."
    os.makedirs(out_dir, exist_ok=True)
    if args.plot:
        os.makedirs(os.path.dirname(args.png) or ".", exist_ok=True)

    df = pd.DataFrame({"t": t, "Pi": pi})
    df.to_csv(args.out, index=False)
    print(f"Wrote {args.out} with {len(df)} rows")
    if args.plot:
        try:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(8, 3))
            plt.plot(t, pi, lw=1.0)
            plt.xlabel("t [s]")
            plt.ylabel("Pi(t)")
            plt.title(f"OU Triad Time Series (tau0={args.tau0})")
            plt.tight_layout()
            plt.savefig(args.png, dpi=150)
            plt.close()
            print(f"Saved {args.png}")
        except Exception as e:
            print(f"Plotting skipped: {e}")


if __name__ == "__main__":
    main()


