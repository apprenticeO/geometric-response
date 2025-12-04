#!/usr/bin/env python3
"""
Paley–Zygmund (PZ) empirical validation for nonnegative time series X_t and triad products.

Features:
- Tests A–F (bounded i.i.d., correlated legs via copula for X=XYZ, plateau+burst, near-zero mean, bounded heavy tail, physics-driven triad proxy)
- Per-cut and systemic averaging modes for CSV input
- Exceedance rate vs PZ lower bound across theta grid with Wilson CIs and bootstrap CIs
- Blocking analysis for plateau/non-ergodic sectors
- Optional hardware-cap bound M and energy-form Y checks using ΔE ≥ (ħ/2)√F
- Plots: survival curve and bar comparison

Usage examples:
  python3 -u "Hyppopotamus/Scripts Demonstration/15_paley_zygmund_validation.py" --mode A --T 20000
  python3 -u "Hyppopotamus/Scripts Demonstration/15_paley_zygmund_validation.py" --mode F --T 5000 --cap-H-op 1.0 --dA 8 --dBar 8
  python3 -u "Hyppopotamus/Scripts Demonstration/15_paley_zygmund_validation.py" --mode file --input-csv data.csv --theta 0.2 0.4 0.6 0.8
"""

import argparse
import json
import math
import os
import random
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from math import erf, sqrt

try:
    import pandas as pd
except Exception:  # pragma: no cover
    pd = None

try:
    import matplotlib
    matplotlib.use("Agg")  # safe default for headless environments
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None


# -------------------------------
# Utilities
# -------------------------------


def ensure_dir(path: str) -> None:
    if path and not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def normal_cdf(x: np.ndarray) -> np.ndarray:
    if np.isscalar(x):
        return 0.5 * (1.0 + erf(x / sqrt(2.0)))
    x = np.asarray(x, dtype=float)
    return 0.5 * (1.0 + np.vectorize(erf)(x / sqrt(2.0)))


def wilson_ci(num_success: int, num_trials: int, confidence: float = 0.95) -> Tuple[float, float]:
    if num_trials <= 0:
        return (float("nan"), float("nan"))
    # Deterministic z for common confidences; default to 0.95
    z_table = {
        0.90: 1.6448536269514722,
        0.95: 1.959963984540054,
        0.99: 2.5758293035489004,
    }
    z = z_table.get(round(confidence, 2), 1.959963984540054)
    p_hat = num_success / num_trials
    denom = 1.0 + (z * z) / num_trials
    center = (p_hat + (z * z) / (2.0 * num_trials)) / denom
    width = (z * math.sqrt((p_hat * (1 - p_hat)) / num_trials + (z * z) / (4.0 * num_trials * num_trials))) / denom
    lo = max(0.0, center - width)
    hi = min(1.0, center + width)
    return (lo, hi)


def bootstrap_ci(values: np.ndarray, stat_fn, num_bootstrap: int = 500, confidence: float = 0.95, rng: Optional[np.random.Generator] = None) -> Tuple[float, float]:
    if values.size == 0:
        return (float("nan"), float("nan"))
    rng = rng or np.random.default_rng()
    stats = []
    n = values.size
    for _ in range(num_bootstrap):
        idx = rng.integers(0, n, size=n)
        try:
            stats.append(stat_fn(values[idx]))
        except Exception:
            stats.append(float("nan"))
    stats = np.array(stats, dtype=float)
    stats = stats[np.isfinite(stats)]
    if stats.size == 0:
        return (float("nan"), float("nan"))
    alpha = (1 - confidence) / 2
    return (np.quantile(stats, alpha), np.quantile(stats, 1 - alpha))


def block_bootstrap_ci(values: np.ndarray, stat_fn, block_len: int, num_bootstrap: int = 500, confidence: float = 0.95, rng: Optional[np.random.Generator] = None) -> Tuple[float, float]:
    """Block bootstrap CI for dependent series.
    Resamples contiguous blocks of length block_len with wrap-around concatenation.
    """
    values = np.asarray(values, dtype=float)
    if values.size == 0 or block_len <= 0:
        return (float("nan"), float("nan"))
    rng = rng or np.random.default_rng()
    n = values.size
    b = max(1, int(block_len))
    num_blocks = int(math.ceil(n / b))
    stats = []
    for _ in range(num_bootstrap):
        blocks = []
        for _ in range(num_blocks):
            start = int(rng.integers(0, n))
            end = start + b
            if end <= n:
                blocks.append(values[start:end])
            else:
                wrap = np.concatenate((values[start:n], values[0:(end - n)]))
                blocks.append(wrap)
        sample = np.concatenate(blocks)[:n]
        try:
            stats.append(stat_fn(sample))
        except Exception:
            stats.append(float("nan"))
    stats = np.array(stats, dtype=float)
    stats = stats[np.isfinite(stats)]
    if stats.size == 0:
        return (float("nan"), float("nan"))
    alpha = (1 - confidence) / 2
    return (np.quantile(stats, alpha), np.quantile(stats, 1 - alpha))


def paley_zygmund_bound(x: np.ndarray, theta: float) -> float:
    x = np.asarray(x, dtype=float)
    x = np.where(np.isfinite(x), x, 0.0)
    x = np.clip(x, 0.0, None)
    mu = float(np.mean(x))
    m2 = float(np.mean(x * x))
    if m2 <= 0.0 or mu <= 0.0:
        return 0.0
    theta = float(theta)
    theta = min(max(theta, 0.0), 1.0)
    return (1.0 - theta) ** 2 * (mu * mu) / m2


def exceedance_rate(x: np.ndarray, theta: float) -> Tuple[float, int]:
    mu = float(np.mean(x))
    threshold = theta * mu
    hits = int(np.sum(x >= threshold))
    n = x.size
    return (hits / n if n > 0 else float("nan"), hits)


def survival_curve(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x >= 0]
    if x.size == 0:
        return np.array([]), np.array([])
    xs = np.sort(x)
    n = xs.size
    # right-continuous survival: P[X >= u]
    surv = 1.0 - (np.arange(n) / n)
    return xs, surv


# -------------------------------
# Synthetic Generators
# -------------------------------


def gen_bounded_iid(T: int, alpha: float = 2.0, beta: float = 5.0, M: float = 1.0, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    rng = rng or np.random.default_rng()
    x = rng.beta(alpha, beta, size=T)
    return M * x


def erf_vec(arr: np.ndarray) -> np.ndarray:
    return np.vectorize(math.erf)(arr)


def gen_correlated_product(T: int, rho: float = 0.5, Mx: float = 1.0, My: float = 1.0, Mz: float = 1.0, rng: Optional[np.random.Generator] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate X=XYZ with Gaussian-copula dependence among legs.
    Note: Pearson correlations among the resulting uniforms differ from input rho; the copula structure is Gaussian but marginal correlations are distorted by the CDF transform.
    """
    rng = rng or np.random.default_rng()
    rho = max(min(rho, 0.999), -0.999)
    cov = np.array([[1.0, rho, rho], [rho, 1.0, rho], [rho, rho, 1.0]], dtype=float)
    L = np.linalg.cholesky(cov)
    z = rng.standard_normal(size=(T, 3)) @ L.T
    u = 0.5 * (1.0 + erf_vec(z / math.sqrt(2.0)))
    x = Mx * u[:, 0]
    y = My * u[:, 1]
    z_leg = Mz * u[:, 2]
    prod = x * y * z_leg
    return prod, x, y, z_leg


def gen_plateau_burst(T: int, lam: float = 0.1, x0: float = 0.01, M: float = 1.0, alpha: float = 2.0, beta: float = 5.0, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    rng = rng or np.random.default_rng()
    mask = rng.random(T) < lam
    bursts = x0 + (M - x0) * rng.beta(alpha, beta, size=T)
    x = np.full(T, x0, dtype=float)
    x[mask] = bursts[mask]
    return x


def gen_near_zero_mean(T: int, alpha: float = 0.2, beta: float = 5.0, M: float = 1.0, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    return gen_bounded_iid(T, alpha=alpha, beta=beta, M=M, rng=rng)


def gen_bounded_heavy_tail(T: int, sigma_logn: float = 1.5, M: float = 1.0, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    rng = rng or np.random.default_rng()
    y = rng.lognormal(mean=0.0, sigma=sigma_logn, size=T)
    y = np.minimum(y, np.quantile(y, 0.999))  # clip extreme outliers before scaling
    y = y / (np.max(y) + 1e-12)
    return M * y


# Physics-driven proxy for triad legs
def gen_physics_proxy(
    T: int,
    h_op: float = 1.0,
    d_a: int = 4,
    d_bar: int = 4,
    a0: float = 0.05,
    a1: float = 0.8,
    omega: float = 0.2,
    eps_std: float = 0.02,
    b0: float = 0.01,
    b1: float = 0.6,
    tau_s: float = 50.0,
    eta_std: float = 0.01,
    c0: float = 0.01,
    c1: float = 0.6,
    t0: float = 30.0,
    tau_i: float = 60.0,
    zeta_std: float = 0.01,
    dt: float = 1.0,
    seed: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    hbar = 1.0  # unit convention (time-scaled generators)
    cap_f = 2.0 * h_op / hbar
    cap_s = math.log(max(d_a, 1))
    cap_i = 2.0 * math.log(max(min(d_a, d_bar), 1))
    t = np.arange(T, dtype=float) * dt
    sqrtF = a0 + a1 * np.abs(np.sin(omega * t)) + rng.normal(0.0, eps_std, size=T)
    sqrtF = np.clip(sqrtF, 0.0, cap_f)
    S = b0 + b1 * (1.0 - np.exp(-t / max(tau_s, 1e-6))) + rng.normal(0.0, eta_std, size=T)
    S = np.clip(S, 0.0, cap_s)
    I = c0 + c1 * (1.0 - np.exp(-np.maximum(t - t0, 0.0) / max(tau_i, 1e-6))) + rng.normal(0.0, zeta_std, size=T)
    I = np.clip(I, 0.0, cap_i)
    X = sqrtF * S * I
    # Energy form using ΔE ≥ (ħ/2)√F (set ħ=1)
    Y = 0.5 * sqrtF * S * I
    M_sup = cap_f * cap_s * cap_i
    return {"t": t, "sqrtF": sqrtF, "S": S, "I": I, "X": X, "Y": Y, "M": np.full(T, M_sup)}


# -------------------------------
# Core evaluation
# -------------------------------


@dataclass
class PZResult:
    theta: float
    mu_hat: float
    m2_hat: float
    p_hat: float
    p_ci_lo: float
    p_ci_hi: float
    pz_hat: float
    pz_ci_lo: float
    pz_ci_hi: float
    pass_primary: bool
    pass_strict: bool
    mu_near_zero: bool
    time_avg_floor_X: float
    time_avg_floor_Y: float


def evaluate_series(x: np.ndarray, thetas: Iterable[float], num_bootstrap: int = 500, confidence: float = 0.95, rng: Optional[np.random.Generator] = None, block_len: int = 0, mu_floor: float = 1e-12) -> List[PZResult]:
    x = np.asarray(x, dtype=float)
    x = np.where(np.isfinite(x), x, 0.0)
    x = np.clip(x, 0.0, None)
    rng = rng or np.random.default_rng()
    mu_hat = float(np.mean(x))
    m2_hat = float(np.mean(x * x))
    results: List[PZResult] = []
    n = x.size
    for th in thetas:
        p_hat, hits = exceedance_rate(x, th)
        p_lo, p_hi = wilson_ci(hits, n, confidence=confidence)
        pz = paley_zygmund_bound(x, th)
        # Bootstrap CI for PZ bound
        def _pz_stat(sample):
            mu = float(np.mean(sample))
            m2 = float(np.mean(sample * sample))
            return (1.0 - th) ** 2 * (mu * mu) / m2 if (m2 > 0 and mu > 0) else 0.0

        if block_len and block_len > 0:
            pz_lo, pz_hi = block_bootstrap_ci(x, _pz_stat, block_len=int(block_len), num_bootstrap=num_bootstrap, confidence=confidence, rng=rng)
        else:
            pz_lo, pz_hi = bootstrap_ci(x, _pz_stat, num_bootstrap=num_bootstrap, confidence=confidence, rng=rng)
        pass_primary = (p_hat - pz) >= 0.0 or (p_lo > pz) or (p_hat > pz_hi)
        pass_strict = (p_lo >= pz_hi)
        mu_near_zero = (mu_hat < mu_floor)
        floor_x = float(th * (1.0 - th) ** 2 * (mu_hat ** 3) / (m2_hat + 1e-16)) if mu_hat > 0 and m2_hat > 0 else 0.0
        floor_y = 0.5 * floor_x
        results.append(
            PZResult(theta=float(th), mu_hat=mu_hat, m2_hat=m2_hat, p_hat=p_hat, p_ci_lo=p_lo, p_ci_hi=p_hi, pz_hat=pz, pz_ci_lo=pz_lo, pz_ci_hi=pz_hi, pass_primary=bool(pass_primary), pass_strict=bool(pass_strict), mu_near_zero=bool(mu_near_zero), time_avg_floor_X=floor_x, time_avg_floor_Y=floor_y)
        )
    return results


def energy_form_checks(x: np.ndarray, y: Optional[np.ndarray], thetas: Iterable[float], M_sup: Optional[float]) -> Dict[str, List[Tuple[float, float]]]:
    # Returns dict with keys "ineq_mu_m2" and "ineq_M" each as list of (rhs, lhs) pairs per theta
    x = np.asarray(x, dtype=float)
    x = np.clip(np.where(np.isfinite(x), x, 0.0), 0.0, None)
    y = np.asarray(y, dtype=float) if y is not None else 0.5 * x  # Y ≥ (ħ/2) X with ħ=1
    y = np.clip(np.where(np.isfinite(y), y, 0.0), 0.0, None)
    mu = float(np.mean(x))
    m2 = float(np.mean(x * x))
    y_avg = float(np.mean(y))
    out_mu_m2: List[Tuple[float, float]] = []
    out_M: List[Tuple[float, float]] = []
    for th in thetas:
        rhs1 = 0.5 * th * (1.0 - th) ** 2 * (mu ** 3) / (m2 + 1e-16) if mu > 0 and m2 > 0 else 0.0
        out_mu_m2.append((rhs1, y_avg))
        rhs2 = 0.5 * th * (1.0 - th) ** 2 * (mu ** 2) / (float(M_sup) + 1e-16) if (mu > 0 and M_sup and M_sup > 0) else 0.0
        out_M.append((rhs2, y_avg))
    return {"ineq_mu_m2": out_mu_m2, "ineq_M": out_M}


# -------------------------------
# Plotting
# -------------------------------


def plot_survival_and_thresholds(x: np.ndarray, thetas: Iterable[float], out_path: str) -> None:
    if plt is None:
        return
    xs, surv = survival_curve(x)
    if xs.size == 0:
        return
    mu = float(np.mean(x))
    plt.figure(figsize=(7, 4))
    plt.plot(xs, surv, label="empirical survival")
    for th in thetas:
        thr = th * mu
        plt.axvline(thr, linestyle="--", alpha=0.5, label=f"u = {th:.2f}·mu")
    plt.xlabel("u")
    plt.ylabel("P[X ≥ u]")
    plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def plot_pz_bars(results: List[PZResult], out_path: str) -> None:
    if plt is None or not results:
        return
    thetas = [r.theta for r in results]
    p_hats = [r.p_hat for r in results]
    p_lo = [r.p_ci_lo for r in results]
    p_hi = [r.p_ci_hi for r in results]
    pz = [r.pz_hat for r in results]
    pz_lo = [r.pz_ci_lo for r in results]
    pz_hi = [r.pz_ci_hi for r in results]
    x_idx = np.arange(len(thetas))
    width = 0.35
    plt.figure(figsize=(8, 4))
    plt.bar(x_idx - width / 2, p_hats, width=width, label="empirical p", alpha=0.8, color="#4C78A8")
    plt.bar(x_idx + width / 2, pz, width=width, label="PZ bound", alpha=0.8, color="#F58518")
    # Error bars with nonnegative yerr
    p_hats_arr = np.array(p_hats)
    p_lo_arr = np.array(p_lo)
    p_hi_arr = np.array(p_hi)
    pz_arr = np.array(pz)
    pz_lo_arr = np.array(pz_lo)
    pz_hi_arr = np.array(pz_hi)
    yerr_p_lower = np.maximum(0.0, p_hats_arr - p_lo_arr)
    yerr_p_upper = np.maximum(0.0, p_hi_arr - p_hats_arr)
    yerr_pz_lower = np.maximum(0.0, pz_arr - pz_lo_arr)
    yerr_pz_upper = np.maximum(0.0, pz_hi_arr - pz_arr)
    plt.errorbar(x_idx - width / 2, p_hats_arr, yerr=[yerr_p_lower, yerr_p_upper], fmt="none", ecolor="#1F3A93", capsize=3)
    plt.errorbar(x_idx + width / 2, pz_arr, yerr=[yerr_pz_lower, yerr_pz_upper], fmt="none", ecolor="#A33F03", capsize=3)
    plt.xticks(x_idx, [f"{t:.2f}" for t in thetas])
    plt.xlabel("theta")
    plt.ylabel("rate / bound")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


# -------------------------------
# CSV I/O for per-cut and systemic
# -------------------------------


def read_time_series_from_csv(path: str) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    if pd is None:
        raise RuntimeError("pandas not available to read CSV")
    df = pd.read_csv(path)
    # Modes supported:
    # 1) single series in column 'X' -> return X, None
    # 2) wide per-cut columns prefixed with 'X_' -> return None, matrix [T x num_cuts]
    has_X = "X" in df.columns
    x_cols = [c for c in df.columns if c.startswith("X_")]
    has_Xstar = bool(x_cols)
    if has_X and has_Xstar:
        raise ValueError("CSV contains both 'X' and 'X_*' columns; please provide one mode only.")
    if has_X:
        x = df["X"].astype(float).to_numpy()
        return x, None
    if has_Xstar:
        X = df[x_cols].astype(float).to_numpy()
        return np.array([]), X
    # fallback: use all numeric columns
    num_df = df.select_dtypes(include=[np.number])
    if num_df.shape[1] == 1:
        return num_df.iloc[:, 0].to_numpy(dtype=float), None
    return np.array([]), num_df.to_numpy(dtype=float)


def evaluate_per_cut(X_mat: np.ndarray, thetas: Iterable[float], **kwargs) -> Dict[str, object]:
    results_per_cut: Dict[str, List[PZResult]] = {}
    pass_flags: Dict[str, bool] = {}
    for j in range(X_mat.shape[1]):
        r = evaluate_series(X_mat[:, j], thetas, **kwargs)
        results_per_cut[f"A{j+1}"] = r
        pass_flags[f"A{j+1}"] = all(rr.pass_primary for rr in r)
    # systemic average across cuts
    x_bar = np.mean(X_mat, axis=1)
    sys_results = evaluate_series(x_bar, thetas, **kwargs)
    sys_pass = all(rr.pass_primary for rr in sys_results)
    return {"per_cut": results_per_cut, "per_cut_pass": pass_flags, "systemic": sys_results, "systemic_pass": sys_pass}


# -------------------------------
# Main
# -------------------------------


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Paley–Zygmund empirical validation (A15)")
    parser.add_argument("--mode", type=str, required=True, choices=["A", "B", "C", "D", "E", "F", "file"], help="Test mode or CSV file mode")
    parser.add_argument("--input-csv", type=str, default="", help="CSV path (mode=file). Columns: X or X_* per-cut (not both)")
    parser.add_argument("--theta", type=float, nargs="*", default=[1/3, 0.2, 0.4, 0.5, 0.6, 0.8], help="Theta grid in (0,1); includes 1/3 which maximizes theta(1-theta)^2")
    parser.add_argument("--T", type=int, default=10000, help="Time series length for generators")
    parser.add_argument("--seed", type=int, default=12345, help="PRNG seed")
    parser.add_argument("--bootstrap", type=int, default=500, help="Bootstrap resamples for PZ CI")
    parser.add_argument("--confidence", type=float, default=0.95, help="Confidence level for CIs")
    parser.add_argument("--rho", type=float, nargs="*", default=[0.0, 0.5, 0.9], help="Correlation values for mode B sweep")
    parser.add_argument("--lam", type=float, nargs="*", default=[0.01, 0.1, 0.5], help="Burst fraction for mode C sweep")
    parser.add_argument("--block-size", type=int, default=200, help="Block length for blocking analysis")
    parser.add_argument("--plateau-var-th", type=float, default=1e-9, help="Variance threshold to tag plateau blocks")
    parser.add_argument("--cap-H-op", type=float, default=1.0, help="Operator norm cap ||H_A||_op (used for sqrt(F) cap via 2||H||/ħ in time-scaled units)")
    parser.add_argument("--block-bootstrap", type=int, default=0, help="Block length for block bootstrap CI of PZ bound (0 disables)")
    parser.add_argument("--dA", type=int, default=4, help="d_A for caps")
    parser.add_argument("--dBar", type=int, default=4, help="d_\\bar{A} for caps")
    parser.add_argument("--out-dir", type=str, default="Hyppopotamus/Scripts Demonstration", help="Base output directory")
    args = parser.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    # stash block_len for evaluate_series
    args.block_len = int(max(0, args.block_bootstrap))
    base_out = os.path.join(args.out_dir, "A15")
    plots_dir = os.path.join(base_out, "plots")
    results_dir = os.path.join(base_out, "results")
    ensure_dir(plots_dir)
    ensure_dir(results_dir)

    thetas = [t for t in args.theta if 0.0 < t < 1.0]
    if not thetas:
        print("No valid theta in (0,1)", file=sys.stderr)
        return 2

    summary_rows: List[Dict[str, object]] = []

    if args.mode == "A":
        x = gen_bounded_iid(args.T, rng=rng)
        results = evaluate_series(x, thetas, num_bootstrap=args.bootstrap, confidence=args.confidence, rng=rng, block_len=args.block_len)
        plot_survival_and_thresholds(x, thetas, os.path.join(plots_dir, "A_bounded_iid_survival.png"))
        plot_pz_bars(results, os.path.join(plots_dir, "A_bounded_iid_bars.png"))
        for r in results:
            summary_rows.append({"mode": "A", "theta": r.theta, "mu": r.mu_hat, "m2": r.m2_hat, "p_hat": r.p_hat, "p_lo": r.p_ci_lo, "p_hi": r.p_ci_hi, "pz": r.pz_hat, "pz_lo": r.pz_ci_lo, "pz_hi": r.pz_ci_hi, "pass": r.pass_primary, "pass_strict": r.pass_strict, "mu_near_zero": r.mu_near_zero, "time_avg_floor_X": r.time_avg_floor_X, "time_avg_floor_Y": r.time_avg_floor_Y})

    elif args.mode == "B":
        for rho in args.rho:
            prod, x_leg, y_leg, z_leg = gen_correlated_product(args.T, rho=float(rho), rng=rng)
            results = evaluate_series(prod, thetas, num_bootstrap=args.bootstrap, confidence=args.confidence, rng=rng, block_len=args.block_len)
            plot_survival_and_thresholds(prod, thetas, os.path.join(plots_dir, f"B_copula_rho{rho:.2f}_survival.png"))
            plot_pz_bars(results, os.path.join(plots_dir, f"B_copula_rho{rho:.2f}_bars.png"))
            for r in results:
                summary_rows.append({"mode": "B", "rho": rho, "theta": r.theta, "mu": r.mu_hat, "m2": r.m2_hat, "p_hat": r.p_hat, "p_lo": r.p_ci_lo, "p_hi": r.p_ci_hi, "pz": r.pz_hat, "pz_lo": r.pz_ci_lo, "pz_hi": r.pz_ci_hi, "pass": r.pass_primary, "pass_strict": r.pass_strict, "mu_near_zero": r.mu_near_zero, "time_avg_floor_X": r.time_avg_floor_X, "time_avg_floor_Y": r.time_avg_floor_Y})

    elif args.mode == "C":
        for lam in args.lam:
            x = gen_plateau_burst(args.T, lam=float(lam), rng=rng)
            results = evaluate_series(x, thetas, num_bootstrap=args.bootstrap, confidence=args.confidence, rng=rng, block_len=args.block_len)
            plot_survival_and_thresholds(x, thetas, os.path.join(plots_dir, f"C_plateau_lam{lam:.2f}_survival.png"))
            plot_pz_bars(results, os.path.join(plots_dir, f"C_plateau_lam{lam:.2f}_bars.png"))
            for r in results:
                summary_rows.append({"mode": "C", "lam": lam, "theta": r.theta, "mu": r.mu_hat, "m2": r.m2_hat, "p_hat": r.p_hat, "p_lo": r.p_ci_lo, "p_hi": r.p_ci_hi, "pz": r.pz_hat, "pz_lo": r.pz_ci_lo, "pz_hi": r.pz_ci_hi, "pass": r.pass_primary, "pass_strict": r.pass_strict, "mu_near_zero": r.mu_near_zero, "time_avg_floor_X": r.time_avg_floor_X, "time_avg_floor_Y": r.time_avg_floor_Y})

    elif args.mode == "D":
        x = gen_near_zero_mean(args.T, rng=rng)
        results = evaluate_series(x, thetas, num_bootstrap=args.bootstrap, confidence=args.confidence, rng=rng, block_len=args.block_len)
        plot_survival_and_thresholds(x, thetas, os.path.join(plots_dir, "D_near_zero_mean_survival.png"))
        plot_pz_bars(results, os.path.join(plots_dir, "D_near_zero_mean_bars.png"))
        for r in results:
            summary_rows.append({"mode": "D", "theta": r.theta, "mu": r.mu_hat, "m2": r.m2_hat, "p_hat": r.p_hat, "p_lo": r.p_ci_lo, "p_hi": r.p_ci_hi, "pz": r.pz_hat, "pz_lo": r.pz_ci_lo, "pz_hi": r.pz_ci_hi, "pass": r.pass_primary, "pass_strict": r.pass_strict, "mu_near_zero": r.mu_near_zero, "time_avg_floor_X": r.time_avg_floor_X, "time_avg_floor_Y": r.time_avg_floor_Y})

    elif args.mode == "E":
        x = gen_bounded_heavy_tail(args.T, rng=rng)
        results = evaluate_series(x, thetas, num_bootstrap=args.bootstrap, confidence=args.confidence, rng=rng, block_len=args.block_len)
        plot_survival_and_thresholds(x, thetas, os.path.join(plots_dir, "E_heavy_tail_survival.png"))
        plot_pz_bars(results, os.path.join(plots_dir, "E_heavy_tail_bars.png"))
        for r in results:
            summary_rows.append({"mode": "E", "theta": r.theta, "mu": r.mu_hat, "m2": r.m2_hat, "p_hat": r.p_hat, "p_lo": r.p_ci_lo, "p_hi": r.p_ci_hi, "pz": r.pz_hat, "pz_lo": r.pz_ci_lo, "pz_hi": r.pz_ci_hi, "pass": r.pass_primary, "pass_strict": r.pass_strict, "mu_near_zero": r.mu_near_zero, "time_avg_floor_X": r.time_avg_floor_X, "time_avg_floor_Y": r.time_avg_floor_Y})

    elif args.mode == "F":
        proxy = gen_physics_proxy(args.T, h_op=float(args.cap_H_op), d_a=int(args.dA), d_bar=int(args.dBar), seed=args.seed)
        x = proxy["X"]
        y = proxy["Y"]
        M_sup = float(proxy["M"][0])
        results = evaluate_series(x, thetas, num_bootstrap=args.bootstrap, confidence=args.confidence, rng=rng, block_len=args.block_len)
        energy_checks = energy_form_checks(x, y, thetas, M_sup)
        plot_survival_and_thresholds(x, thetas, os.path.join(plots_dir, "F_physics_proxy_survival.png"))
        plot_pz_bars(results, os.path.join(plots_dir, "F_physics_proxy_bars.png"))
        for idx, r in enumerate(results):
            rhs1, lhs1 = energy_checks["ineq_mu_m2"][idx]
            rhs2, lhs2 = energy_checks["ineq_M"][idx]
            summary_rows.append({
                "mode": "F",
                "theta": r.theta,
                "mu": r.mu_hat,
                "m2": r.m2_hat,
                "p_hat": r.p_hat,
                "p_lo": r.p_ci_lo,
                "p_hi": r.p_ci_hi,
                "pz": r.pz_hat,
                "pz_lo": r.pz_ci_lo,
                "pz_hi": r.pz_ci_hi,
                "pass": r.pass_primary,
                "pass_strict": r.pass_strict,
                "mu_near_zero": r.mu_near_zero,
                "Y_avg": float(np.mean(y)),
                "ineq_mu_m2_rhs": rhs1,
                "ineq_mu_m2_lhs": lhs1,
                "ineq_M_rhs": rhs2,
                "ineq_M_lhs": lhs2,
                "M_sup": M_sup,
                "time_avg_floor_X": r.time_avg_floor_X,
                "time_avg_floor_Y": r.time_avg_floor_Y,
            })

    elif args.mode == "file":
        if not args.input_csv:
            print("--input-csv required for mode=file", file=sys.stderr)
            return 2
        x_single, X_mat = read_time_series_from_csv(args.input_csv)
        if X_mat is not None and X_mat.size > 0:
            evald = evaluate_per_cut(X_mat, thetas, num_bootstrap=args.bootstrap, confidence=args.confidence, rng=rng)
            # per-cut
            for k, rlist in evald["per_cut"].items():
                for r in rlist:
                    summary_rows.append({"mode": "file-per-cut", "cut": k, "theta": r.theta, "mu": r.mu_hat, "m2": r.m2_hat, "p_hat": r.p_hat, "p_lo": r.p_ci_lo, "p_hi": r.p_ci_hi, "pz": r.pz_hat, "pz_lo": r.pz_ci_lo, "pz_hi": r.pz_ci_hi, "pass": r.pass_primary, "pass_strict": r.pass_strict, "mu_near_zero": r.mu_near_zero, "time_avg_floor_X": r.time_avg_floor_X, "time_avg_floor_Y": r.time_avg_floor_Y})
            # systemic
            sys_res = evald["systemic"]
            for r in sys_res:
                summary_rows.append({"mode": "file-systemic", "theta": r.theta, "mu": r.mu_hat, "m2": r.m2_hat, "p_hat": r.p_hat, "p_lo": r.p_ci_lo, "p_hi": r.p_ci_hi, "pz": r.pz_hat, "pz_lo": r.pz_ci_lo, "pz_hi": r.pz_ci_hi, "pass": r.pass_primary, "pass_strict": r.pass_strict, "mu_near_zero": r.mu_near_zero, "time_avg_floor_X": r.time_avg_floor_X, "time_avg_floor_Y": r.time_avg_floor_Y})
            # plots for systemic series
            plot_survival_and_thresholds(np.mean(X_mat, axis=1), thetas, os.path.join(plots_dir, "file_systemic_survival.png"))
            plot_pz_bars(sys_res, os.path.join(plots_dir, "file_systemic_bars.png"))
        else:
            x = x_single
            results = evaluate_series(x, thetas, num_bootstrap=args.bootstrap, confidence=args.confidence, rng=rng, block_len=args.block_len)
            for r in results:
                summary_rows.append({"mode": "file-single", "theta": r.theta, "mu": r.mu_hat, "m2": r.m2_hat, "p_hat": r.p_hat, "p_lo": r.p_ci_lo, "p_hi": r.p_ci_hi, "pz": r.pz_hat, "pz_lo": r.pz_ci_lo, "pz_hi": r.pz_ci_hi, "pass": r.pass_primary, "pass_strict": r.pass_strict, "mu_near_zero": r.mu_near_zero, "time_avg_floor_X": r.time_avg_floor_X, "time_avg_floor_Y": r.time_avg_floor_Y})
            plot_survival_and_thresholds(x, thetas, os.path.join(plots_dir, "file_single_survival.png"))
            plot_pz_bars(results, os.path.join(plots_dir, "file_single_bars.png"))

        # Blocking analysis on whichever series we can form
        if X_mat is not None and X_mat.size > 0:
            x_block_base = np.mean(X_mat, axis=1)
        else:
            x_block_base = x_single
        if x_block_base is not None and x_block_base.size > 0 and args.block_size > 0:
            num_blocks = max(1, x_block_base.size // args.block_size)
            blocks = np.array_split(x_block_base, num_blocks)
            block_rows = []
            for b_idx, xb in enumerate(blocks):
                var_b = float(np.var(xb))
                tag = "plateau" if var_b <= args.plateau_var_th else "active"
                rlist = evaluate_series(xb, thetas, num_bootstrap=max(100, args.bootstrap // 5), confidence=args.confidence, rng=rng)
                for r in rlist:
                    block_rows.append({"block": b_idx, "tag": tag, "theta": r.theta, "mu": r.mu_hat, "m2": r.m2_hat, "p_hat": r.p_hat, "pz": r.pz_hat, "pass": r.pass_primary})
            if pd is not None:
                pd.DataFrame(block_rows).to_csv(os.path.join(results_dir, "file_blocking_summary.csv"), index=False)

    else:
        print(f"Unknown mode: {args.mode}", file=sys.stderr)
        return 2

    # Write summary
    if pd is not None and summary_rows:
        pd.DataFrame(summary_rows).to_csv(os.path.join(results_dir, "A15_pz_summary.csv"), index=False)
    # Write JSON snapshot of args
    with open(os.path.join(results_dir, "A15_args.json"), "w") as f:
        json.dump(vars(args), f, indent=2)

    # Optional caps CSV for physics proxy (mode F)
    if args.mode == "F":
        try:
            proxy_caps = {
                "cap_f": float(2.0 * args.cap_H_op / 1.0),
                "cap_s": float(math.log(max(args.dA, 1))),
                "cap_i": float(2.0 * math.log(max(min(args.dA, args.dBar), 1))),
                "M_sup": float(2.0 * args.cap_H_op / 1.0 * math.log(max(args.dA, 1)) * 2.0 * math.log(max(min(args.dA, args.dBar), 1))),
            }
            if pd is not None:
                pd.DataFrame([proxy_caps]).to_csv(os.path.join(results_dir, "A15_caps_summary.csv"), index=False)
        except Exception:
            pass

    print(f"A15 completed. Results: {results_dir}, Plots: {plots_dir}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())


