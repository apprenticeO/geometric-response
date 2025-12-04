#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd


def small_omega_mask(omega: np.ndarray, tau: float, max_frac: float = 0.2) -> np.ndarray:
    return np.isfinite(omega) & (np.abs(omega) * float(tau) <= max_frac) & (omega != 0)


def main():
    ap = argparse.ArgumentParser(description="Compute τ_PSD and Re χ flatness (if available)")
    ap.add_argument("--results_dir", type=str, default="./results", help="Directory with tau_estimates.csv and g_eff.csv")
    ap.add_argument("--plot", action="store_true", help="Save chi_re.png if Re χ is available")
    args = ap.parse_args()

    tau_csv = os.path.join(args.results_dir, "tau_estimates.csv")
    g_csv = os.path.join(args.results_dir, "g_eff.csv")
    out_csv = os.path.join(args.results_dir, "metrics_extras.csv")
    out_png = os.path.join(args.results_dir, "chi_re.png")

    if not os.path.exists(tau_csv):
        raise FileNotFoundError(tau_csv)
    tau_df = pd.read_csv(tau_csv)
    row = tau_df.iloc[0]

    var_pi = float(row.get("var_pi", np.nan))
    s0_raw = float(row.get("psd_S0_one_sided", np.nan))
    s0_trunc = float(row.get("psd_S0_trunc", np.nan))
    tau_g_trunc = float(row.get("tau_G_trunc", np.nan))
    tau_kk = float(row.get("tau_KK", np.nan)) if "tau_KK" in row else np.nan

    tau_psd_raw = np.nan
    if np.isfinite(s0_raw) and np.isfinite(var_pi) and var_pi > 0:
        tau_psd_raw = s0_raw / (2.0 * var_pi)
    tau_psd_trunc = np.nan
    if np.isfinite(s0_trunc) and np.isfinite(var_pi) and var_pi > 0:
        tau_psd_trunc = s0_trunc / (2.0 * var_pi)

    # Defaults for Re χ flatness
    flatness = np.nan
    band_tau = tau_kk if np.isfinite(tau_kk) and tau_kk > 0 else tau_g_trunc

    if os.path.exists(g_csv):
        gdf = pd.read_csv(g_csv)
        if "omega_rad_s" in gdf.columns and "Re_chi" in gdf.columns:
            omega = gdf["omega_rad_s"].to_numpy(dtype=float)
            re_chi = gdf["Re_chi"].to_numpy(dtype=float)
            m = small_omega_mask(omega, band_tau, max_frac=0.2)
            if np.any(m) and np.isfinite(re_chi[m]).sum() > 2:
                vals = re_chi[m]
                vals = vals[np.isfinite(vals)]
                if vals.size > 2:
                    mu = float(np.mean(vals))
                    sigma = float(np.std(vals, ddof=1))
                    if mu != 0:
                        flatness = abs(sigma / mu)
            if args.plot and np.any(m):
                try:
                    import matplotlib.pyplot as plt
                    plt.figure(figsize=(6, 3))
                    plt.plot(omega, re_chi, lw=1.0, label="Re χ(ω)")
                    plt.scatter(omega[m], re_chi[m], s=8, alpha=0.7, label="|ω|τ≤0.2 band")
                    plt.axhline(0.0, color="k", lw=0.5)
                    plt.xlabel("ω [rad/s]")
                    plt.ylabel("Re χ(ω)")
                    title = f"Re χ flatness σ/mean={flatness:.3g}" if np.isfinite(flatness) else "Re χ flatness (N/A)"
                    plt.title(title)
                    plt.legend()
                    plt.tight_layout()
                    plt.savefig(out_png, dpi=150)
                    plt.close()
                except Exception as e:
                    print(f"Plotting Re χ skipped: {e}")

    out = {
        "tau_psd_raw_s": tau_psd_raw,
        "tau_psd_trunc_s": tau_psd_trunc,
        "band_tau_s": band_tau,
        "re_chi_flatness_sigma_over_mean": flatness,
    }
    pd.DataFrame([out]).to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()


