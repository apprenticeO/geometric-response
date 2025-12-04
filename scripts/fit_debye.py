#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from utils import small_omega_mask


def main():
    ap = argparse.ArgumentParser(description="Fit Debye single-pole chi(ω)=(1+iωτ)^{-1} near ω=0")
    ap.add_argument("--in", dest="inp", type=str, required=True, help="Input g_eff.csv with omega, Re_chi, Im_chi")
    ap.add_argument("--out", type=str, default="debye_fit.csv", help="Output CSV path")
    ap.add_argument("--plot", action="store_true", help="Save a plot of Im χ(ω) and Debye fit")
    ap.add_argument("--png", type=str, default="debye_fit.png", help="Output PNG path")
    ap.add_argument("--tau_guess", type=float, default=None, help="Initial τ guess (s)")
    ap.add_argument("--max_frac", type=float, default=0.2, help="Use |ω|τ<=max_frac for fit band")
    args = ap.parse_args()

    df = pd.read_csv(args.inp)
    # Accept either 'omega_rad_s' (pipeline default) or 'omega' as the frequency column
    if "omega_rad_s" in df.columns:
        omega = df["omega_rad_s"].to_numpy()
    elif "omega" in df.columns:
        omega = df["omega"].to_numpy()
    else:
        raise ValueError("Input CSV must contain 'omega_rad_s' or 'omega' column")
    chi_re = df["Re_chi"].to_numpy()
    chi_im = df["Im_chi"].to_numpy()

    # If no guess, infer from low-ω slope of Im
    tau_guess = args.tau_guess
    if tau_guess is None:
        # robust slope from first 5% positive ω samples
        pos = omega > 0
        sel = np.where(pos)[0][: max(1, int(0.05 * np.sum(pos)))]
        if sel.size == 0:
            tau_guess = 1.0
        else:
            slope = np.median(-chi_im[sel] / omega[sel])
            tau_guess = float(np.abs(slope)) if np.isfinite(slope) else 1.0

    mask = small_omega_mask(omega, tau_guess, max_frac=args.max_frac) & (omega != 0)
    if not np.any(mask):
        raise RuntimeError("No small-omega samples for Debye fit band")
    w = omega[mask]
    y_im = chi_im[mask]

    # Debye model: Im χ = -(ω τ) / (1 + (ω τ)^2)
    def model_im(om, tau):
        return -(om * tau) / (1.0 + (om * tau) ** 2)

    # Simple 1D search over τ near guess
    taus = tau_guess * np.logspace(-2, 2, 200)
    errs = [np.mean((y_im - model_im(w, tau)) ** 2) for tau in taus]
    best_idx = int(np.argmin(errs))
    tau_fit = float(taus[best_idx])

    pd.DataFrame([{"tau_fit": tau_fit, "tau_guess": tau_guess}]).to_csv(args.out, index=False)
    print(f"Debye τ_fit={tau_fit:.6g}s (guess {tau_guess:.6g}s); wrote {args.out}")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(6, 3))
            plt.plot(omega, chi_im, label='Im χ(ω)')
            plt.plot(omega, model_im(omega, tau_fit), '--', label=f'Debye fit τ={tau_fit:.3g}s')
            plt.axhline(0, color='k', lw=0.5)
            plt.xlabel('ω [rad/s]')
            plt.ylabel('Im χ(ω)')
            plt.legend()
            plt.tight_layout()
            plt.savefig(args.png, dpi=150)
            plt.close()
            print(f"Saved {args.png}")
        except Exception as e:
            print(f"Plotting skipped: {e}")


if __name__ == "__main__":
    main()


