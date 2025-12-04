#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys

# Suggested run for sine drive and analysis settings:
# python3 scripts/chain_two_qubit_pipeline.py --drive_type sine --drive_freq 0.02 --drive_eps 0.01 --nperseg 32768 --nfft_factor 8
# Response channel notes:
#   --response ObsZA   -> R(t) = ⟨σz_A⟩ = Tr[ρ(t) (σz ⊗ I)]
#   --response ObsXA   -> R(t) = ⟨σx_A⟩ = Tr[ρ(t) (σx ⊗ I)]
#   --response EnergyA -> R(t) = ⟨σx_A⟩ as a local-energy proxy aligned with h(t)=h+εu(t)
#   --response Pi      -> R(t) = Π(t) (useful for Π-only checks; GK still uses Π)


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description="Chain: two-qubit demo -> GK/KK -> Debye fit -> scales")
    ap.add_argument("--outdir", type=str, default="results/two_qubit_pipeline", help="Base output directory")
    # two-qubit params
    ap.add_argument("--J", type=float, default=0.5, help="Coupling for σz⊗σz")
    ap.add_argument("--h", type=float, default=1.0, help="Transverse field amplitude")
    ap.add_argument("--dt", type=float, default=0.01, help="Time step [s]")
    ap.add_argument("--tmax", type=float, default=400.0, help="Total time [s]")
    ap.add_argument("--drive_eps", type=float, default=0.02, help="Drive amplitude ε for h(t)=h+ε u(t)")
    ap.add_argument("--drive_type", type=str, default="white", choices=["white", "sine"], help="Drive type")
    ap.add_argument("--drive_freq", type=float, default=0.05, help="Sine drive frequency [Hz] if drive_type=sine")
    ap.add_argument("--drive_seed", type=int, default=0, help="Random seed for white drive")
    # analysis params
    ap.add_argument("--nperseg", type=int, default=131072, help="Welch segment length")
    ap.add_argument("--nfft_factor", type=int, default=16, help="Zero-padding factor")
    ap.add_argument("--acf_cut", type=str, default="zero", help="ACF truncation rule: 'zero' or epsilon")
    ap.add_argument("--response", type=str, default="ObsZA", choices=["Pi", "ObsZA", "ObsXA", "EnergyA"], help="Response R(t) channel for cross-spectral χ(ω)")
    # auto-frequency discovery
    ap.add_argument("--auto_freq", action="store_true", help="Estimate dominant response frequency and set sine drive at a fraction")
    ap.add_argument("--auto_frac", type=float, default=0.1, help="Fraction of dominant frequency for sine drive (e.g., 0.1)")
    args = ap.parse_args()

    base = os.path.abspath(args.outdir)
    os.makedirs(base, exist_ok=True)
    csv_path = os.path.join(base, "two_qubit.csv")

    project = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    demo = os.path.join(project, "scripts", "run_two_qubit_demo.py")
    comp = os.path.join(project, "scripts", "compute_gk_kk.py")
    fit = os.path.join(project, "scripts", "fit_debye.py")
    scales = os.path.join(project, "scripts", "write_scales.py")
    rmunu = os.path.join(project, "scripts", "export_rmunu_diagnostic.py")

    # Optional 0) Auto-discover dominant response frequency on short un-driven run
    drive_freq = args.drive_freq
    if args.auto_freq:
        tmp_csv = os.path.join(base, "two_qubit_nodrive.csv")
        run([
            sys.executable, demo,
            "--J", str(args.J),
            "--h", str(args.h),
            "--dt", str(args.dt),
            "--tmax", str(min(args.tmax, 200.0)),
            "--out", tmp_csv,
            "--response", args.response,
            "--drive_eps", "0.0",
        ])
        # Estimate peak frequency (Hz) from response R(t)
        import pandas as pd
        import numpy as np
        from scipy.signal import welch, get_window
        df = pd.read_csv(tmp_csv)
        t = df["t"].to_numpy(dtype=float)
        r = df["R"].to_numpy(dtype=float) if "R" in df.columns else df["Pi"].to_numpy(dtype=float)
        r = r - np.mean(r)
        dt = float(np.mean(np.diff(t)))
        fs = 1.0 / dt
        f, pxx = welch(r, fs=fs, window=get_window("hann", min(8192, len(r))), nperseg=min(8192, len(r)), noverlap=min(8192, len(r))//2, nfft=min(8192, len(r)))
        # pick the largest peak at f>0
        if f.size > 1:
            idx = int(np.argmax(pxx[1:]) + 1)
            f_peak = float(max(f[idx], 1e-6))
            drive_freq = float(max(1e-4, args.auto_frac * f_peak))
        print(f"[auto_freq] peak≈{f_peak:.6g} Hz, drive_freq set to {drive_freq:.6g} Hz (frac={args.auto_frac})")

    # 1) Generate two-qubit dataset with drive and response
    run([
        sys.executable, demo,
        "--J", str(args.J),
        "--h", str(args.h),
        "--dt", str(args.dt),
        "--tmax", str(args.tmax),
        "--out", csv_path,
        "--plot",
        "--png_dir", base,
        "--drive_eps", str(args.drive_eps),
        "--drive_type", ("sine" if args.auto_freq else args.drive_type),
        "--drive_freq", str(drive_freq),
        "--drive_seed", str(args.drive_seed),
        "--response", args.response,
    ])

    # 2) Compute GK/KK and diagnostics
    run([
        sys.executable, comp,
        "--in", csv_path,
        "--outdir", base,
        "--plot",
        "--nperseg", str(args.nperseg),
        "--nfft_factor", str(args.nfft_factor),
        "--acf_cut", args.acf_cut,
        "--loww_band", "0.2",
        "--passivity_sign", "1",
    ])

    # 3) Debye fit (optional, if Re/Im χ present)
    g_csv = os.path.join(base, "g_eff.csv")
    if os.path.exists(g_csv):
        run([
            sys.executable, fit,
            "--in", g_csv,
            "--out", os.path.join(base, "debye_fit.csv"),
            "--plot",
            "--png", os.path.join(base, "debye_fit.png"),
        ])

    # 4) Scales and pipeline summary
    run([
        sys.executable, scales,
        "--results_dir", base,
    ])

    # 4b) Optional diagnostic: R_{μν} ∝ (ħ/c^5) Π²
    try:
        run([
            sys.executable, rmunu,
            "--in", csv_path,
            "--outdir", base,
            "--plot",
        ])
    except Exception:
        pass

    # 5) Write a cautionary notice: two-qubit is for τ_G/KK/passivity, not G0
    try:
        notice = (
            "NOTICE: This two-qubit dataset is intended for τ_G (GK), KK low-ω slope, and passivity checks.\n"
            "Do NOT interpret it as an absolute G0 calibration. In the framework, G0 is an adiabatic/DC stiffness,\n"
            "to be set either by Planck-time normalization or by a realistic matter model in the ω→0 limit.\n"
        )
        with open(os.path.join(base, "NOTICE_toy_model.txt"), "w") as f:
            f.write(notice)
    except Exception:
        pass

    print(f"Pipeline completed. Outputs in {base}")


if __name__ == "__main__":
    main()


