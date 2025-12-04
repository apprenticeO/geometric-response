#!/usr/bin/env python3
import argparse, os
import numpy as np
import pandas as pd
from utils import (
    normalized_acf, gk_tau_from_acf, short_time_bound, one_sided_psd_zero,
    GKResults, welch_cross_spectrum, tau_kk_from_chi, bootstrap_tau_g
)
import json, hashlib


def main():
    ap = argparse.ArgumentParser(description="Compute GK tau_G, ACF diagnostics, and KK slope from Pi(t); optional chi(ω) from R and R_info")
    ap.add_argument("--in", dest="inp", type=str, required=True, help="Input CSV with columns: t, Pi [optional: R, R_info]")
    ap.add_argument("--outdir", type=str, default="./results", help="Output directory")
    ap.add_argument("--nperseg", type=int, default=512, help="Welch segment length")
    ap.add_argument("--timescale", type=float, default=1.0, help="Physical seconds per time-unit in CSV (default=1)")
    ap.add_argument("--acf_cut", type=str, default="zero", help="Truncation rule for GK integral: 'zero' (first zero-crossing) or numeric epsilon (e.g., 0.01)")
    ap.add_argument("--plot", action="store_true", help="Save diagnostic plots (Pi, ACF, chi)")
    ap.add_argument("--nfft_factor", type=int, default=1, help="Zero-padding factor for Welch (nfft = nperseg * factor)")
    ap.add_argument("--loww_band", type=float, default=None, help="Set small-ω band threshold (e.g., 0.2 means |ω|τ<=0.2)")
    ap.add_argument("--passivity_sign", type=int, default=-1, choices=[-1, 1], help="Passivity check uses sign*Imχ(ω) ≤ 0 for ω>0; use -1 for e^{+iωt}, +1 for e^{-iωt}")
    ap.add_argument("--bootstrap", action="store_true", help="Compute sliding-window bootstrap CI for τ_G")
    ap.add_argument("--bootstrap_window_s", type=float, default=20.0, help="Bootstrap window length (seconds)")
    ap.add_argument("--bootstrap_step_s", type=float, default=None, help="Bootstrap step (seconds); default window/2")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.inp)
    if not {"t", "Pi"}.issubset(df.columns):
        raise ValueError("Input must contain columns: t, Pi")
    # sanitize
    keep_cols = ["t", "Pi"] + [c for c in ["R", "R_info"] if c in df.columns]
    df = df[keep_cols].dropna()
    t = df["t"].to_numpy(dtype=float) * float(args.timescale)
    pi = df["Pi"].to_numpy(dtype=float)
    if np.allclose(np.var(pi), 0.0):
        raise ValueError("Pi series has ~zero variance; cannot normalize ACF or estimate τ_G.")
    dt = float(np.mean(np.diff(t)))
    print(f"dt={dt:.6g} s, len={len(t)}, total_time={t[-1]-t[0]:.6g} s")

    # ACF and GK integral
    lags, acf = normalized_acf(pi, dt)
    # GK integral (full)
    tau_g_full = gk_tau_from_acf(lags, acf)
    # Optional truncation to avoid oscillatory tails
    cut_rule = args.acf_cut.strip().lower()
    cut_idx = len(lags) - 1
    cut_val = None
    if cut_rule == "zero":
        # first index where ACF <= 0
        nz = np.where(acf <= 0)[0]
        if nz.size > 0:
            cut_idx = int(nz[0])
    else:
        try:
            eps = float(cut_rule)
            cut_val = eps
            nz = np.where(np.abs(acf) <= eps)[0]
            if nz.size > 0:
                cut_idx = int(nz[0])
        except Exception:
            pass
    tau_g_trunc = gk_tau_from_acf(lags[: cut_idx + 1], acf[: cut_idx + 1])
    tau_bound = short_time_bound(pi, dt)
    var_pi = float(np.var(pi))
    s0 = one_sided_psd_zero(acf, dt)
    s0_trunc = 2.0 * var_pi * tau_g_trunc
    l1_abs = float(np.trapz(np.abs(acf), lags))

    # Write ACF diagnostics
    pd.DataFrame({
        "lag_s": lags,
        "acf": acf
    }).to_csv(f"{args.outdir}/acf_diagnostics.csv", index=False)

    # tau estimates
    pd.DataFrame([{
        "tau_G_full": tau_g_full,
        "tau_G_trunc": tau_g_trunc,
        "acf_cut": cut_rule,
        "cut_lag_s": lags[cut_idx],
        "tau_bound": tau_bound,
        "psd_S0_one_sided": s0,
        "psd_S0_trunc": s0_trunc,
        "var_pi": var_pi,
        "l1_abs": l1_abs,
        "kk_seed": "tau_G_trunc"
    }]).to_csv(f"{args.outdir}/tau_estimates.csv", index=False)

    # optional bootstrap CI
    if args.bootstrap:
        mean_b, low_b, high_b, n_b = bootstrap_tau_g(
            pi, dt, window_s=float(args.bootstrap_window_s), step_s=(None if args.bootstrap_step_s is None else float(args.bootstrap_step_s))
        )
        tau_df = pd.read_csv(f"{args.outdir}/tau_estimates.csv")
        tau_df["tau_G_boot_mean"] = mean_b
        tau_df["tau_G_boot_low95"] = low_b
        tau_df["tau_G_boot_high95"] = high_b
        tau_df["tau_G_boot_n"] = n_b
        tau_df.to_csv(f"{args.outdir}/tau_estimates.csv", index=False)

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(8, 3))
            plt.plot(t, pi, lw=1.0)
            plt.xlabel("t [s]")
            plt.ylabel("Pi(t)")
            plt.title("Triad Time Series")
            plt.tight_layout()
            plt.savefig(f"{args.outdir}/pi_timeseries.png", dpi=150)
            plt.close()

            plt.figure(figsize=(6, 3))
            plt.plot(lags, acf, lw=1.0)
            plt.axvline(lags[cut_idx], ls="--", lw=0.8)
            plt.xlabel("lag τ [s]")
            plt.ylabel("C(τ)")
            plt.title(f"ACF; τ_G(full)≈{tau_g_full:.3g} s, τ_G(trunc)≈{tau_g_trunc:.3g} s; cut@{lags[cut_idx]:.3g}s")
            plt.tight_layout()
            plt.savefig(f"{args.outdir}/acf.png", dpi=150)
            plt.close()
        except Exception as e:
            print(f"Plotting (Pi/ACF) skipped: {e}")

    # Optional: empirical χ(ω) from R and R_info
    if {"R", "R_info"}.issubset(df.columns):
        r = df["R"].to_numpy(dtype=float)
        rinfo = df["R_info"].to_numpy(dtype=float)
        # mean removal to reduce leakage
        r = r - np.mean(r)
        rinfo = rinfo - np.mean(rinfo)
        fs = 1.0 / dt
        nfft = int(args.nperseg * max(1, args.nfft_factor))
        f, sxy = welch_cross_spectrum(rinfo, r, fs=fs, nperseg=args.nperseg, nfft=nfft)
        # Transfer χ(ω) ≈ S_R,Rinfo / S_Rinfo,Rinfo
        _, sxx = welch_cross_spectrum(rinfo, rinfo, fs=fs, nperseg=args.nperseg, nfft=nfft)
        chi = sxy / (sxx + 1e-12)
        omega = 2 * np.pi * f
        chi_re = np.real(chi)
        chi_im = np.imag(chi)
        tau_seed = float(tau_g_trunc if np.isfinite(tau_g_trunc) else tau_g_full)
        tau_kk = tau_kk_from_chi(omega, chi_im, tau_guess=tau_seed, max_frac_start=(None if args.loww_band is None else float(args.loww_band)))
        band = np.isfinite(tau_kk) & (np.abs(omega) * float(tau_kk) <= 0.2)
        sgn = int(args.passivity_sign)
        passivity_ok = bool(np.all(sgn * chi_im[(omega > 0) & band] <= 0.0))

        df_chi = pd.DataFrame({
            "omega_rad_s": omega,
            "Re_chi": chi_re,
            "Im_chi": chi_im,
            "Re_Geff": chi_re,
            "Im_Geff": chi_im
        })
        df_chi.to_csv(f"{args.outdir}/g_eff.csv", index=False)

        # append tau_KK and passivity
        tau_df = pd.read_csv(f"{args.outdir}/tau_estimates.csv")
        tau_df["tau_KK"] = tau_kk
        tau_df["passivity_ok"] = passivity_ok
        tau_df["passivity_sign"] = sgn
        tau_df.to_csv(f"{args.outdir}/tau_estimates.csv", index=False)

        if args.plot:
            try:
                import matplotlib.pyplot as plt
                plt.figure(figsize=(6, 3))
                plt.plot(omega, chi_im, lw=1.0, label="Im χ(ω)")
                if np.isfinite(tau_kk) and np.any(band):
                    plt.scatter(omega[band], chi_im[band], s=8, alpha=0.7, label="|ω|τ≤0.2 band")
                plt.axhline(0.0, color="k", lw=0.5)
                plt.xlabel("ω [rad/s]")
                plt.ylabel("Im χ(ω)")
                plt.title(f"KK τ≈{tau_kk:.3g} s; passivity={bool(passivity_ok)}")
                plt.legend()
                plt.tight_layout()
                plt.savefig(f"{args.outdir}/chi_im.png", dpi=150)
                plt.close()
            except Exception as e:
                print(f"Plotting (chi) skipped: {e}")
    else:
        # Π-only path via FDT-normalized PSD: χ̃_im(ω) = -ω S_Π(ω) / S_Π(0)
        try:
            from scipy.signal import welch, get_window
            fs = 1.0 / dt
            nfft = int(args.nperseg * max(1, args.nfft_factor))
            f, pxx = welch(
                pi - np.mean(pi),
                fs=fs,
                window=get_window('hann', args.nperseg),
                nperseg=args.nperseg,
                noverlap=args.nperseg // 2,
                nfft=nfft,
                return_onesided=True,
            )
            omega = 2 * np.pi * f
            # One-sided S_Π(0) from GK (more robust than extrapolating): S(0)=2 Var Π τ_G
            tau_seed = float(tau_g_trunc if np.isfinite(tau_g_trunc) else tau_g_full)
            S0_norm = 2.0 * var_pi * tau_seed
            chi_im = -(omega * pxx) / (S0_norm + 1e-24)

            tau_kk = tau_kk_from_chi(omega, chi_im, tau_guess=tau_seed, max_frac_start=(None if args.loww_band is None else float(args.loww_band)))
            band = np.isfinite(tau_kk) & (np.abs(omega) * float(tau_kk) <= 0.2)
            sgn = int(args.passivity_sign)
            passivity_ok = bool(np.all(sgn * chi_im[(omega > 0) & band] <= 0.0))

            df_chi = pd.DataFrame({
                "omega_rad_s": omega,
                "Re_chi": np.full_like(omega, np.nan, dtype=float),
                "Im_chi": chi_im,
                "source": "psd"
            })
            df_chi.to_csv(f"{args.outdir}/g_eff.csv", index=False)

            tau_df = pd.read_csv(f"{args.outdir}/tau_estimates.csv")
            tau_df["tau_KK"] = tau_kk
            tau_df["passivity_ok"] = passivity_ok
            tau_df["passivity_sign"] = sgn
            tau_df.to_csv(f"{args.outdir}/tau_estimates.csv", index=False)

            if args.plot:
                import matplotlib.pyplot as plt
                plt.figure(figsize=(6, 3))
                plt.plot(omega, chi_im, lw=1.0, label="Im χ̃(ω) from PSD")
                if np.isfinite(tau_kk) and np.any(band):
                    plt.scatter(omega[band], chi_im[band], s=8, alpha=0.7, label="|ω|τ≤0.2 band")
                plt.axhline(0.0, color="k", lw=0.5)
                plt.xlabel("ω [rad/s]")
                plt.ylabel("Im χ̃(ω)")
                plt.title(f"KK τ≈{tau_kk:.3g} s (PSD); passivity={bool(passivity_ok)}")
                plt.legend()
                plt.tight_layout()
                plt.savefig(f"{args.outdir}/chi_im.png", dpi=150)
                plt.close()
        except Exception as e:
            print(f"PSD/χ̃ path skipped: {e}")

    # provenance JSON
    try:
        prov = {
            "script": "compute_gk_kk.py",
            "input_csv": os.path.abspath(args.inp),
            "outdir": os.path.abspath(args.outdir),
            "params": {
                "nperseg": int(args.nperseg),
                "nfft_factor": int(args.nfft_factor),
                "acf_cut": args.acf_cut,
                "timescale": float(args.timescale),
                "loww_band": (None if args.loww_band is None else float(args.loww_band)),
                "bootstrap": bool(args.bootstrap),
                "bootstrap_window_s": float(args.bootstrap_window_s),
                "bootstrap_step_s": (None if args.bootstrap_step_s is None else float(args.bootstrap_step_s)),
            },
            "stats": {
                "dt": dt,
                "len": int(len(t)),
                "var_pi": var_pi,
                "tau_G_full": tau_g_full,
                "tau_G_trunc": tau_g_trunc,
                "l1_abs": l1_abs,
            }
        }
        with open(os.path.join(args.outdir, "provenance.json"), "w") as f:
            json.dump(prov, f, indent=2)
    except Exception as e:
        print(f"Provenance write skipped: {e}")

    print(f"tau_G(full)={tau_g_full:.6g}s, tau_G(trunc)={tau_g_trunc:.6g}s, l1_abs={l1_abs:.6g}; outputs in {args.outdir}")


if __name__ == "__main__":
    main()


