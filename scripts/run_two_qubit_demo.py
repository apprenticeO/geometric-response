#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from typing import Tuple

try:
    from scipy.linalg import expm, eigh
except Exception as e:
    raise SystemExit(f"scipy is required for this demo (expm, eigh): {e}")


def paulis() -> Tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return I, X, Y, Z


def kron(*ops: NDArray[np.complex128]) -> NDArray[np.complex128]:
    out = np.array([[1.0 + 0j]], dtype=complex)
    for a in ops:
        out = np.kron(out, a)
    return out


def projector(vec: NDArray[np.complex128]) -> NDArray[np.complex128]:
    v = vec.reshape(-1, 1).astype(complex)
    return v @ v.conj().T


def partial_trace_two_qubit(rho: NDArray[np.complex128], keep: str) -> NDArray[np.complex128]:
    """Trace over the other qubit. keep ∈ {'A','B'}."""
    rho4 = rho.reshape(2, 2, 2, 2)
    if keep == 'A':
        # trace over B: sum over indices (b = d), keep (a,c)
        return np.einsum('abcb->ac', rho4)
    elif keep == 'B':
        # trace over A: sum over indices (a = c), keep (b,d)
        return np.einsum('abad->bd', rho4)
    else:
        raise ValueError("keep must be 'A' or 'B'")


def von_neumann_entropy(rho: NDArray[np.complex128], eps: float = 1e-10) -> float:
    w = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
    w = np.clip(np.real(w), eps, 1.0)
    w = w / np.sum(w)
    return float(-np.sum(w * np.log(w)))


def qfi_spectral(rho: NDArray[np.complex128], H: NDArray[np.complex128], eps: float = 1e-10) -> float:
    # Braunstein–Caves spectral formula: F_Q = 2 Σ_{i,j} ( (p_i - p_j)^2 / (p_i + p_j) |⟨i|H|j⟩|^2 )
    w, U = eigh((rho + rho.conj().T) / 2.0)
    w = np.clip(np.real(w), eps, 1.0)
    w = w / np.sum(w)
    Hm = U.conj().T @ H @ U
    F = 0.0
    for i in range(len(w)):
        for j in range(len(w)):
            denom = w[i] + w[j]
            if denom > eps:
                num = (w[i] - w[j]) ** 2
                F += 2.0 * num / denom * (np.abs(Hm[i, j]) ** 2)
    return float(np.real(F))


def evolve_two_qubit(J: float, h: float, dt: float, tmax: float) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    I, X, Y, Z = paulis()
    # H = J σz⊗σz + h (σx⊗I + I⊗σx)
    H = J * kron(Z, Z) + h * (kron(X, I) + kron(I, X))
    # initial product state |+0>
    plus = (1 / np.sqrt(2)) * np.array([1, 1], dtype=complex)
    zero = np.array([1, 0], dtype=complex)
    psi0 = np.kron(plus, zero)
    rho0 = projector(psi0)

    times = np.arange(0.0, tmax + 1e-12, dt, dtype=float)
    rhos = np.zeros((len(times), 4, 4), dtype=complex)
    for k, t in enumerate(times):
        U = expm(-1j * H * t)
        rhos[k] = U @ rho0 @ U.conj().T
    return times, rhos


def evolve_two_qubit_driven(
    J: float,
    h0: float,
    dt: float,
    tmax: float,
    u: NDArray[np.float64],
) -> Tuple[NDArray[np.float64], NDArray[np.complex128]]:
    """
    Piecewise-constant driven evolution with h_n = h0 + eps*u[n].
    Uses per-step propagators U_n = exp(-i H(h_n) dt).
    """
    I, X, Y, Z = paulis()
    # initial product state |+0>
    plus = (1 / np.sqrt(2)) * np.array([1, 1], dtype=complex)
    zero = np.array([1, 0], dtype=complex)
    psi0 = np.kron(plus, zero)
    rho = projector(psi0)

    times = np.arange(0.0, tmax + 1e-12, dt, dtype=float)
    if len(u) != len(times):
        raise ValueError("drive u length must match number of time samples")
    rhos = np.zeros((len(times), 4, 4), dtype=complex)
    rhos[0] = rho
    for k in range(1, len(times)):
        h_k = h0 + float(u[k])
        H_k = J * kron(Z, Z) + h_k * (kron(X, I) + kron(I, X))
        U = expm(-1j * H_k * dt)
        rho = U @ rho @ U.conj().T
        rhos[k] = rho
    return times, rhos


def main():
    ap = argparse.ArgumentParser(description="Two-qubit microscopic demo: ρ(t) -> ρ_A(t) -> {S,I,F_Q} -> Π(t) with optional weak drive")
    ap.add_argument("--J", type=float, default=0.5, help="Coupling for σz⊗σz")
    ap.add_argument("--h", type=float, default=1.0, help="Transverse field for σx⊗I + I⊗σx")
    ap.add_argument("--dt", type=float, default=0.01, help="Time step (s)")
    ap.add_argument("--tmax", type=float, default=100.0, help="Total time (s)")
    ap.add_argument("--eps", type=float, default=1e-10, help="Rank-lift epsilon for entropy/QFI")
    ap.add_argument("--out", type=str, default="results/two_qubit.csv", help="Output CSV path")
    ap.add_argument("--plot", action="store_true", help="Save PNG plots")
    ap.add_argument("--png_dir", type=str, default="results", help="Directory for PNG outputs")
    # Drive options
    ap.add_argument("--drive_eps", type=float, default=0.0, help="Drive amplitude ε; uses h(t)=h+ε*u(t)")
    ap.add_argument("--drive_type", type=str, default="white", choices=["white", "sine"], help="Drive type for u(t)")
    ap.add_argument("--drive_freq", type=float, default=0.1, help="Drive frequency [Hz] for sine drive (ignored for white)")
    ap.add_argument("--drive_seed", type=int, default=0, help="Random seed for white drive")
    # Response channel options
    ap.add_argument("--response", type=str, default="Pi", choices=["Pi", "ObsZA", "ObsXA", "EnergyA"], help="Output response R(t): 'Pi', ⟨σz_A⟩, ⟨σx_A⟩, or local A-energy proxy")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    if args.plot:
        os.makedirs(args.png_dir or ".", exist_ok=True)

    # Optional weak drive u(t); when drive_eps=0, fall back to constant-H evolution
    use_drive = float(args.drive_eps) != 0.0
    if use_drive:
        t = np.arange(0.0, args.tmax + 1e-12, args.dt, dtype=float)
        if args.drive_type == "white":
            rng = np.random.default_rng(int(args.drive_seed))
            u_sig = rng.normal(0.0, 1.0, size=t.shape[0]).astype(float)
            # normalize to unit variance and zero mean
            u_sig = (u_sig - np.mean(u_sig))
        else:
            # sine drive with unit amplitude
            u_sig = np.sin(2.0 * np.pi * float(args.drive_freq) * t).astype(float)
        # scale by ε
        u_series = float(args.drive_eps) * u_sig
        t, rhos = evolve_two_qubit_driven(args.J, args.h, args.dt, args.tmax, u_series.astype(float))
    else:
        t, rhos = evolve_two_qubit(args.J, args.h, args.dt, args.tmax)
        u_series = np.zeros_like(t, dtype=float)
    I2, X, Y, Z = paulis()
    H_A = Z  # σz on qubit A

    S_A = []
    S_B = []
    S_AB = []
    I_MI = []
    FQ = []
    Pi = []
    for rho in rhos:
        rhoA = partial_trace_two_qubit(rho, 'A')
        rhoB = partial_trace_two_qubit(rho, 'B')
        sA = von_neumann_entropy(rhoA, eps=args.eps)
        sB = von_neumann_entropy(rhoB, eps=args.eps)
        sAB = von_neumann_entropy(rho, eps=args.eps)
        Ival = sA + sB - sAB
        fq = qfi_spectral(rhoA, H_A, eps=args.eps)
        pi = np.sqrt(max(fq, 0.0)) * sA * max(Ival, 0.0)
        S_A.append(sA)
        S_B.append(sB)
        S_AB.append(sAB)
        I_MI.append(Ival)
        FQ.append(fq)
        Pi.append(pi)

    # Choose response observable R(t)
    # Option ObsZA/ObsXA/EnergyA: expectation values on full state ρ(t)
    I2, X, Y, Z = paulis()
    Z_A = kron(Z, I2)
    X_A = kron(X, I2)
    R_series: NDArray[np.float64]
    if args.response == "ObsZA":
        R_vals = []
        for rho in rhos:
            val = np.real(np.trace(rho @ Z_A))
            R_vals.append(float(val))
        R_series = np.array(R_vals, dtype=float)
    elif args.response == "ObsXA":
        R_vals = []
        for rho in rhos:
            val = np.real(np.trace(rho @ X_A))
            R_vals.append(float(val))
        R_series = np.array(R_vals, dtype=float)
    elif args.response == "EnergyA":
        # local energy proxy on A: use ⟨σx_A⟩ which couples to h modulation
        R_vals = []
        for rho in rhos:
            val = np.real(np.trace(rho @ X_A))
            R_vals.append(float(val))
        R_series = np.array(R_vals, dtype=float)
    else:
        R_series = np.array(Pi, dtype=float)

    df = pd.DataFrame({
        "t": t.astype(float),
        "S_A": np.array(S_A, dtype=float),
        "S_B": np.array(S_B, dtype=float),
        "S_AB": np.array(S_AB, dtype=float),
        "I": np.array(I_MI, dtype=float),
        "FQ": np.array(FQ, dtype=float),
        "Pi": np.array(Pi, dtype=float),
        # For GK→KK cross-spectrum: treat input as R_info and output as R=Pi
        "R_info": u_series.astype(float),
        "R": R_series,
    })
    df.to_csv(args.out, index=False)
    print(f"Wrote {args.out} with {len(df)} rows; mean I={np.mean(df['I']):.4g}, mean S_A={np.mean(df['S_A']):.4g}")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            plt.figure(figsize=(7, 3))
            plt.plot(df["t"], df["Pi"], lw=1.0)
            plt.xlabel("t [s]")
            plt.ylabel("Pi(t)")
            plt.title("Two-qubit triad Π(t)")
            plt.tight_layout()
            png = os.path.join(args.png_dir, "two_qubit_pi.png")
            plt.savefig(png, dpi=150)
            plt.close()
            print(f"Saved {png}")
            if use_drive:
                plt.figure(figsize=(7, 2.5))
                plt.plot(df["t"], df["R_info"], lw=0.9, label="R_info(t)")
                plt.xlabel("t [s]")
                plt.ylabel("drive")
                plt.title("Drive R_info(t)")
                plt.tight_layout()
                png2 = os.path.join(args.png_dir, "two_qubit_rinfo.png")
                plt.savefig(png2, dpi=150)
                plt.close()
                print(f"Saved {png2}")
        except Exception as e:
            print(f"Plotting skipped: {e}")


if __name__ == "__main__":
    main()


