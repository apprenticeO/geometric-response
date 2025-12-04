#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd


def generate_single_pole(t: np.ndarray, tau: float, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Generate (R_info, R) where R is single-pole response of time-constant tau to input R_info.
    Discrete-time exact decay: y[n] = a*y[n-1] + (1-a)*u[n] with a = exp(-dt/tau).
    """
    rng = np.random.default_rng(seed)
    dt = float(np.mean(np.diff(t)))
    a = float(np.exp(-dt / tau))
    u = rng.normal(0.0, 1.0, size=t.shape[0]).astype(float)
    y = np.zeros_like(u)
    for n in range(1, len(u)):
        y[n] = a * y[n - 1] + (1.0 - a) * u[n]
    return u, y


def main():
    ap = argparse.ArgumentParser(description="Generate Lindblad-like single-pole dataset with R_info and R")
    ap.add_argument("--tmax", type=float, default=4000.0, help="Total time (s)")
    ap.add_argument("--dt", type=float, default=0.01, help="Time step (s)")
    ap.add_argument("--tau", type=float, default=2.0, help="Single-pole time constant τ (s)")
    ap.add_argument("--seed", type=int, default=0, help="Random seed")
    ap.add_argument("--out", type=str, default="results/lindblad.csv", help="Output CSV path")
    args = ap.parse_args()

    t = np.arange(0.0, args.tmax + 1e-12, args.dt, dtype=float)
    rinfo, r = generate_single_pole(t, args.tau, args.seed)
    # Use Pi = R so GK metrics refer to output process
    df = pd.DataFrame({"t": t, "Pi": r, "R_info": rinfo, "R": r})
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"Wrote {args.out} with {len(df)} rows (tau={args.tau}s)")


if __name__ == "__main__":
    main()


