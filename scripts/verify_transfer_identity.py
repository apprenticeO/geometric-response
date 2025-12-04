#!/usr/bin/env python3
import argparse
import os
import json
import numpy as np
import pandas as pd

from typing import Tuple


def welch_spectrum(x: np.ndarray, fs: float, nperseg: int, nfft: int) -> Tuple[np.ndarray, np.ndarray]:
    from scipy.signal import welch, get_window
    f, pxx = welch(
        x - np.mean(x),
        fs=fs,
        window=get_window('hann', nperseg),
        nperseg=nperseg,
        noverlap=nperseg // 2,
        nfft=nfft,
        return_onesided=True,
    )
    omega = 2.0 * np.pi * f
    return omega, pxx


def welch_csd(x: np.ndarray, y: np.ndarray, fs: float, nperseg: int, nfft: int) -> Tuple[np.ndarray, np.ndarray]:
    from scipy.signal import csd, get_window
    f, pxy = csd(
        x - np.mean(x),
        y - np.mean(y),
        fs=fs,
        window=get_window('hann', nperseg),
        nperseg=nperseg,
        noverlap=nperseg // 2,
        nfft=nfft,
        return_onesided=True,
    )
    omega = 2.0 * np.pi * f
    return omega, pxy


def main():
    ap = argparse.ArgumentParser(description="Verify R(ω) ≈ χ(ω) R_info(ω) in the slow band |ω τ_G| ≤ 0.2")
    ap.add_argument("--in", dest="inp", type=str, required=True, help="Input CSV with columns: t, R_info, R (plus Pi optional). For tensor mode, supply R_<comp> and R_info_<comp>.")
    ap.add_argument("--g_csv", type=str, required=True, help="g_eff.csv containing omega_rad_s, Re_chi, Im_chi")
    ap.add_argument("--tau_csv", type=str, required=True, help="tau_estimates.csv containing tau_G_trunc (and optionally tau_KK)")
    ap.add_argument("--outdir", type=str, default="./results", help="Output directory for plots/CSV")
    ap.add_argument("--provenance", type=str, default=None, help="Optional provenance.json with nperseg and nfft_factor")
    ap.add_argument("--nperseg", type=int, default=None, help="Override Welch segment length")
    ap.add_argument("--nfft_factor", type=int, default=None, help="Override nfft factor (nfft = nperseg * factor)")
    # Tensor-mode options
    ap.add_argument("--tensor_components", type=str, default="", help="Comma-separated components like 'xx,yy,zz,xy,yz,xz'. Uses columns R_<c>, R_info_<c> if present.")
    ap.add_argument("--projector", type=str, default="sum", choices=["sum", "trace", "custom"], help="Projector for tensor reduction to scalar channel")
    ap.add_argument("--projector_weights", type=str, default="", help="JSON mapping of component -> weight for projector='custom' (e.g., '{\"xx\":1, \"yy\":1, \"zz\":1}')")
    # Band and comparison options
    ap.add_argument("--loww_band_frac", type=float, default=0.2, help="Slow-band threshold: use |ω| τ_G ≤ frac")
    ap.add_argument("--use_complex", action="store_true", help="Also compare complex transfer Sxy/Sxx vs χ(ω) (reports complex nRMSE and phase MAE)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.inp)
    if "t" not in df.columns:
        raise ValueError("Input must contain column 't'")
    t = df["t"].to_numpy(dtype=float)
    # Tensor vs scalar input handling
    tensor_mode = False
    used_components: list[str] = []
    if args.tensor_components:
        comps = [c.strip() for c in args.tensor_components.split(",") if c.strip()]
        # Verify columns exist; collect present ones
        present = []
        for c in comps:
            r_col = f"R_{c}"
            ri_col = f"R_info_{c}"
            if r_col in df.columns and ri_col in df.columns:
                present.append(c)
        if present:
            tensor_mode = True
            used_components = present
        # If none present, fall back to scalar mode
    if tensor_mode:
        # Determine weights
        weights: dict[str, float] = {}
        if args.projector == "sum":
            weights = {c: 1.0 for c in used_components}
        elif args.projector == "trace":
            # weight 1 for diagonal-like components; 0 for off-diagonals
            weights = {c: (1.0 if c in {"xx", "yy", "zz"} else 0.0) for c in used_components}
        else:
            # custom
            try:
                weights = json.loads(args.projector_weights) if args.projector_weights else {}
            except Exception as e:
                raise ValueError(f"Failed to parse projector_weights as JSON: {e}")
            # Default missing weights to 0
            weights = {c: float(weights.get(c, 0.0)) for c in used_components}
        # Build scalar channel time series by weighted sum
        rinfo_series = np.zeros_like(t, dtype=float)
        r_series = np.zeros_like(t, dtype=float)
        for c in used_components:
            w = float(weights.get(c, 0.0))
            rinfo_series += w * df[f"R_info_{c}"].to_numpy(dtype=float)
            r_series += w * df[f"R_{c}"].to_numpy(dtype=float)
        rinfo = rinfo_series
        r = r_series
    else:
        # Scalar mode expects R_info and R
        if not {"R_info", "R"}.issubset(df.columns):
            raise ValueError("Input must contain columns: R_info, R (or tensor columns with --tensor_components)")
        rinfo = df["R_info"].to_numpy(dtype=float)
        r = df["R"].to_numpy(dtype=float)
    dt = float(np.mean(np.diff(t)))
    fs = 1.0 / dt

    # Load χ(ω)
    gdf = pd.read_csv(args.g_csv)
    if "omega_rad_s" not in gdf.columns or "Im_chi" not in gdf.columns or "Re_chi" not in gdf.columns:
        raise ValueError("g_eff.csv must contain columns: omega_rad_s, Re_chi, Im_chi")
    omega_g = gdf["omega_rad_s"].to_numpy(dtype=float)
    chi = gdf["Re_chi"].to_numpy(dtype=float) + 1j * gdf["Im_chi"].to_numpy(dtype=float)

    # Load τ_G
    tdf = pd.read_csv(args.tau_csv)
    tau_g = float(tdf.iloc[0].get("tau_G_trunc", float("nan")))
    if not np.isfinite(tau_g) or tau_g <= 0:
        raise ValueError("tau_estimates.csv missing valid tau_G_trunc")

    # Welch settings
    nperseg = args.nperseg
    nfft_factor = args.nfft_factor
    if (nperseg is None or nfft_factor is None) and args.provenance and os.path.exists(args.provenance):
        try:
            with open(args.provenance, "r") as f:
                prov = json.load(f)
            params = prov.get("params", {})
            if nperseg is None:
                nperseg = int(params.get("nperseg", 131072))
            if nfft_factor is None:
                nfft_factor = int(params.get("nfft_factor", 16))
        except Exception:
            pass
    if nperseg is None:
        # Heuristic: aim for ~10 segments; power of two not strictly required
        target = max(32, len(t) // 10)
        nperseg = int(min(131072, max(16, target)))
    if nfft_factor is None:
        nfft_factor = 8
    # Ensure nperseg is valid for signal length
    if nperseg > len(t):
        nperseg = max(16, len(t) // 2)
    nfft = int(nperseg * nfft_factor)

    # Compute spectra for R_info and R on same grid used by χ(ω)
    # We will compute CSD x: R_info, y: R to infer the discrete ω grid
    omega_xy, sxy = welch_csd(rinfo, r, fs=fs, nperseg=nperseg, nfft=nfft)
    omega_xx, sxx = welch_csd(rinfo, rinfo, fs=fs, nperseg=nperseg, nfft=nfft)
    # Also compute direct Welch of R_info and R magnitudes for overlays
    omega_rinfo, pxx_rinfo = welch_spectrum(rinfo, fs=fs, nperseg=nperseg, nfft=nfft)
    omega_r, pxx_r = welch_spectrum(r, fs=fs, nperseg=nperseg, nfft=nfft)

    # Interpolate χ(ω) onto the (CSD) ω grid if needed
    if omega_g.shape != omega_xy.shape or not np.allclose(omega_g, omega_xy):
        chi_interp_real = np.interp(omega_xy, omega_g, np.real(chi), left=np.nan, right=np.nan)
        chi_interp_imag = np.interp(omega_xy, omega_g, np.imag(chi), left=np.nan, right=np.nan)
        # Guard NaNs in Re χ (e.g., PSD-only pipelines provide Im χ̃)
        chi_interp_real = np.nan_to_num(chi_interp_real, nan=0.0, posinf=0.0, neginf=0.0)
        chi_use = chi_interp_real + 1j * chi_interp_imag
    else:
        chi_use = np.nan_to_num(np.real(chi), nan=0.0, posinf=0.0, neginf=0.0) + 1j * np.imag(chi)

    # Predicted R spectrum from χ(ω)·R_info(ω)
    # Use transfer on the input spectrum: Yn_pred ≈ χ(ω)·Xn, where Xn comes from R_info
    # For robustness here we use the cross-spectral ratio Sxy/Sxx as χ_hat consistency check
    # but to test R ≈ χ R_info we compare magnitudes from Welch of R vs |χ|·|R_info|
    # Align R_info Welch grid to CSD grid for consistent comparison
    if omega_rinfo.shape != omega_xy.shape or not np.allclose(omega_rinfo, omega_xy):
        # Interpolate magnitudes to ω_xy
        mag_rinfo = np.sqrt(np.interp(omega_xy, omega_rinfo, np.maximum(pxx_rinfo, 0.0), left=np.nan, right=np.nan))
        mag_r = np.sqrt(np.interp(omega_xy, omega_r, np.maximum(pxx_r, 0.0), left=np.nan, right=np.nan))
    else:
        mag_rinfo = np.sqrt(np.maximum(pxx_rinfo, 0.0))
        mag_r = np.sqrt(np.maximum(pxx_r, 0.0))

    mag_pred = np.abs(chi_use) * mag_rinfo
    residual = mag_r - mag_pred

    # Slow band mask
    band = np.isfinite(omega_xy) & np.isfinite(mag_r) & np.isfinite(mag_pred) & (np.abs(omega_xy) * tau_g <= float(args.loww_band_frac))
    # Normalized RMSE in band
    denom = np.sqrt(np.mean(mag_r[band] ** 2)) if np.any(band) else np.nan
    nrmse = (np.sqrt(np.mean(residual[band] ** 2)) / denom) if (np.any(band) and denom and np.isfinite(denom)) else np.nan

    # Optional complex comparison using transfer ratio
    nrmse_complex = None
    phase_mae = None
    if args.use_complex:
        # Empirical transfer
        with np.errstate(divide="ignore", invalid="ignore"):
            transfer_emp = sxy / (sxx + 1e-24)
        # Predicted transfer
        transfer_pred = chi_use
        # Interpolate if grids differ (should match omega_xy already)
        if transfer_pred.shape != transfer_emp.shape or not np.allclose(omega_xy, omega_xy):
            # noop: both are already on omega_xy
            pass
        diff = transfer_emp - transfer_pred
        # Complex nRMSE over band
        denom_c = np.sqrt(np.mean(np.abs(transfer_emp[band]) ** 2)) if np.any(band) else np.nan
        nrmse_complex = (np.sqrt(np.mean(np.abs(diff[band]) ** 2)) / denom_c) if (np.any(band) and denom_c and np.isfinite(denom_c)) else np.nan
        # Phase MAE in band (radians)
        phase_emp = np.angle(transfer_emp)
        phase_pred = np.angle(transfer_pred)
        dphi = np.unwrap(phase_emp) - np.unwrap(phase_pred)
        phase_mae = float(np.mean(np.abs(dphi[band]))) if np.any(band) else None

    # Save CSV
    suffix = ""
    if tensor_mode:
        suffix = f"_tensor_{args.projector}"
    out_csv = os.path.join(args.outdir, f"transfer_identity_check{suffix}.csv")
    pd.DataFrame({
        "omega_rad_s": omega_xy,
        "mag_R": mag_r,
        "mag_Rinfo": mag_rinfo,
        "mag_chi": np.abs(chi_use),
        "mag_pred": mag_pred,
        "residual_mag": residual,
        "band_mask": band.astype(int),
    }).to_csv(out_csv, index=False)

    # Save metrics JSON
    metrics: dict[str, object] = {
        "tau_G_s": float(tau_g),
        "nperseg": int(nperseg),
        "nfft_factor": int(nfft_factor),
        "n_points_band": int(np.count_nonzero(band)),
        "nrmse_mag_band": (None if not np.isfinite(nrmse) else float(nrmse)),
        "loww_band_frac": float(args.loww_band_frac),
    }
    if tensor_mode:
        metrics["tensor_components"] = used_components
        metrics["projector"] = args.projector
    if args.use_complex:
        metrics["nrmse_complex_band"] = (None if (nrmse_complex is None or not np.isfinite(nrmse_complex)) else float(nrmse_complex))
        metrics["phase_mae_band"] = phase_mae
    with open(os.path.join(args.outdir, f"transfer_identity_metrics{suffix}.json"), "w") as f:
        json.dump(metrics, f, indent=2)

    # Plots
    try:
        import matplotlib.pyplot as plt
        # Magnitude overlay
        plt.figure(figsize=(7, 3.2))
        plt.plot(omega_xy, mag_r, label="|R(ω)|", lw=1.0)
        plt.plot(omega_xy, mag_pred, "--", label="|χ(ω)|·|R_info(ω)|", lw=1.0)
        if np.any(band):
            plt.scatter(omega_xy[band], mag_r[band], s=6, alpha=0.6, label="band |ω|τ_G≤0.2")
        plt.xlabel("ω [rad/s]")
        plt.ylabel("magnitude")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, f"transfer_identity_mag{suffix}.png"), dpi=150)
        plt.close()

        # Residual magnitude
        plt.figure(figsize=(7, 3.2))
        plt.plot(omega_xy, np.abs(residual), lw=1.0, label="|R| - |χ|·|R_info|")
        if np.any(band):
            plt.scatter(omega_xy[band], np.abs(residual[band]), s=6, alpha=0.6, label="band |ω|τ_G≤0.2")
        plt.xlabel("ω [rad/s]")
        plt.ylabel("residual magnitude")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, f"transfer_residual_mag{suffix}.png"), dpi=150)
        plt.close()
    except Exception as e:
        print(f"Plotting skipped: {e}")

    print(f"Saved {out_csv} and transfer_identity_metrics{suffix}.json; band points={metrics['n_points_band']}, nRMSE={metrics['nrmse_mag_band']}")


if __name__ == "__main__":
    main()


