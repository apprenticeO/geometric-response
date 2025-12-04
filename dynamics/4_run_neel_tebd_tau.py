#!/usr/bin/env python3
import math, csv, json, os
from dataclasses import dataclass, asdict
from typing import List, Optional, Tuple

import numpy as np
import quimb as qu
import quimb.tensor as qtn

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None


# -------------------------------
# Data structures
# -------------------------------

@dataclass
class Frame:
    t: float
    S_A_nats: float
    dH_A: float
    Pi_t: float
    good_time: int
    Pi_fluct: Optional[float] = None
    trace_norm_d: Optional[float] = None
    # New fields
    I_AB: Optional[float] = None
    Pi_rate: Optional[float] = None
    Pi_eng: Optional[float] = None


@dataclass
class Summary:
    L: int; ell: int; J: float; h: float; dt: float; steps: int; chi: int; cut: float; hbar: float
    sample_k: int; tau_max_samples: int
    T_total: float; n_frames: int; n_tau_samples: int
    Pi_avg: float; Pi_liminf_proxy: float; frac_good_times: float
    Pi_fluct_avg: Optional[float]; Pi_fluct_liminf_proxy: Optional[float]
    tau_A4_avg: Optional[float]; tau_A4_liminf_proxy: Optional[float]
    # Floors and RHS (speed-based)
    epsF: Optional[float] = None; epsS: Optional[float] = None; epsI: Optional[float] = None; delta: Optional[float] = None
    rhs_speed: Optional[float] = None
    # Hamiltonian-capped components
    triad_hat_mean: Optional[float] = None
    rhs_capped: Optional[float] = None


# -------------------------------
# Helpers reused from existing scripts
# -------------------------------

def entanglement_entropy_nats(psi, ell: int) -> float:
    if not (1 <= ell <= psi.L - 1):
        if ell == psi.L:
            return 0.0
        raise ValueError("ell must be in [1, L-1]")
    # quimb entropy default is base-2; convert to nats
    return float(psi.entropy(ell - 1)) * float(np.log(2.0))


def build_dense_HA(L: int, ell: int, J: float = 1.0, h: float = 0.0) -> np.ndarray:
    sx, sy, sz = qu.pauli('X'), qu.pauli('Y'), qu.pauli('Z')
    id2 = np.eye(2, dtype=complex)
    H = np.zeros((2**L, 2**L), dtype=complex)
    if h != 0.0:
        for i in range(ell):
            ops = [id2] * L; ops[i] = sz
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + h * term
    for i in range(L - 1):
        left_in_A = (i < ell); right_in_A = (i + 1 < ell)
        if left_in_A and right_in_A: w = J
        elif (left_in_A and not right_in_A) or (not left_in_A and right_in_A): w = 0.5 * J
        else: w = 0.0
        if w == 0.0: continue
        for A, B in ((sx, sx), (sy, sy), (sz, sz)):
            ops = [id2] * L; ops[i] = A; ops[i + 1] = B
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + w * term
    return H


def build_dense_HA_on_A(ell: int, J: float = 1.0, h: float = 0.0) -> np.ndarray:
    sx, sy, sz = qu.pauli('X'), qu.pauli('Y'), qu.pauli('Z')
    id2 = np.eye(2, dtype=complex)
    H = np.zeros((2**ell, 2**ell), dtype=complex)
    if h != 0.0:
        for i in range(ell):
            ops = [id2] * ell; ops[i] = sz
            term = ops[0]
            for k in range(1, ell): term = np.kron(term, ops[k])
            H = H + h * term
    for i in range(ell - 1):
        for A, B in ((qu.pauli('X'), qu.pauli('X')), (qu.pauli('Y'), qu.pauli('Y')), (qu.pauli('Z'), qu.pauli('Z'))):
            ops = [id2] * ell; ops[i] = A; ops[i + 1] = B
            term = ops[0]
            for k in range(1, ell): term = np.kron(term, ops[k])
            H = H + J * term
    return H


def _partial_trace_rho(rho: np.ndarray, L: int, keep: List[int]) -> np.ndarray:
    rho_t = rho.reshape([2] * (2 * L)); keep_set = set(keep)
    for k in range(L - 1, -1, -1):
        if k not in keep_set:
            half = rho_t.ndim // 2
            rho_t = np.trace(rho_t, axis1=k, axis2=half + k)
    return rho_t.reshape(2 ** len(keep), 2 ** len(keep))


def variance_from_dense_state(psi_vec: np.ndarray, H_dense: np.ndarray) -> float:
    EH = np.vdot(psi_vec, H_dense @ psi_vec)
    EH2 = np.vdot(psi_vec, H_dense @ (H_dense @ psi_vec))
    return math.sqrt(max(float(np.real(EH2 - EH * EH)), 0.0))


def qfi_sld_unitary(rho: np.ndarray, H: np.ndarray) -> float:
    rho = 0.5 * (rho + rho.conj().T)
    vals, vecs = np.linalg.eigh(rho)
    He = vecs.conj().T @ H @ vecs
    F = 0.0
    for i in range(len(vals)):
        for j in range(len(vals)):
            denom = vals[i] + vals[j]
            if denom > 0.0:
                num = (vals[i] - vals[j]) ** 2
                F += 2.0 * (num / denom) * (abs(He[i, j]) ** 2)
    return float(np.real(F))


# -------------------------------
# Core run
# -------------------------------

def run_neel_tau(L: int, ell: int, J: float, h: float, dt: float, steps: int,
                 chi: int, cut: float, hbar: float, sample_k: int, tau_max_samples: int,
                 out_prefix: str) -> Tuple[List[Frame], Summary]:
    # LocalHam1D accepted by TEBD
    H_local = qtn.ham_1d_heis(L=L, j=J, bz=h, S=0.5, cyclic=False)
    # Dense operators used at sampled frames only
    H_A_dense = build_dense_HA(L, ell, J, h)
    # Operator norm and logs for capped RHS and intensive triad
    try:
        vals = np.linalg.eigvalsh(H_A_dense)
        H_opnorm = float(np.max(np.abs(vals)))
    except Exception:
        H_opnorm = float(np.linalg.norm(H_A_dense, ord=2))
    log_dA = float(ell * np.log(2.0))
    log_dmin = float(min(ell, L - ell) * np.log(2.0)) if L - ell > 0 else log_dA
    H_A_on_A = build_dense_HA_on_A(ell, J, h)

    # Néel initial state
    pattern = ''.join('01'[i % 2] for i in range(L))
    psi = qtn.MPS_computational_state(pattern)

    tebd = qtn.TEBD(psi, H_local, dt=dt, split_opts={'max_bond': chi, 'cutoff': cut}, imag=False)

    frames: List[Frame] = []
    A_sites = list(range(ell))
    Abar_sites = list(range(ell, L))

    t = 0.0
    eps = 1e-12

    d4_samples: List[float] = []
    n_tau = 0
    sqrtF_seq: List[float] = []
    S_seq: List[float] = []
    I_seq: List[float] = []

    triad_hat_accum = 0.0
    for step in range(1, steps + 1):
        tebd.step(); t += dt
        if (step % sample_k) != 0:
            continue

        # Compute entropy on MPS
        try:
            S_A = entanglement_entropy_nats(tebd.pt, ell)
        except ValueError:
            S_A = 0.0

        # Dense state for dH, rho, and QFI
        psi_dense = tebd.pt.to_dense()
        dH = variance_from_dense_state(psi_dense, H_A_dense)
        # Build reduced rho, rho_A, rho_Abar
        rho = np.outer(psi_dense, psi_dense.conj())
        rho_A = _partial_trace_rho(rho, L, A_sites)
        rho_Abar = _partial_trace_rho(rho, L, Abar_sites)

        # Von Neumann entropies and mutual information I(A:Ā)
        def S_vn(r):
            r = 0.5*(r + r.conj().T)
            vals = np.clip(np.linalg.eigvalsh(r), 0.0, 1.0)
            nz = vals[vals > 0.0]
            return float(-np.sum(nz * np.log(nz)))
        S_B  = S_vn(rho_Abar)
        # If exactly pure, S_AB=0; otherwise compute
        S_AB = 0.0
        if abs(np.trace(rho @ rho) - 1.0) > 1e-10:
            S_AB = S_vn(rho)
        I_AB = S_A + S_B - S_AB

        # QFI-based instantaneous functional
        F = qfi_sld_unitary(rho_A, H_A_on_A)
        sqrtF = math.sqrt(max(F, 0.0))
        
        # CORRECTED: Canonical triad formulas (no ℏ factor)
        Pi_rate = sqrtF * S_A * I_AB                    # Main triad: √F_Q · S_A · I(A:Ā)
        Pi_eng  = dH * S_A * I_AB                      # Energy product: ΔE_A · S_A · I(A:Ā)
        Pi_t = Pi_rate                                  # Use main triad for Pi_t
        Pi_fl = sqrtF * S_A * I_AB                     # Same as Pi_rate

        # Intensive triad contribution Ψ̂_A for single cut
        denom_op = max(H_opnorm, 1e-12)
        denom_S = max(log_dA, 1e-12)
        denom_I = max(2.0 * log_dmin, 1e-12)
        triad_hat_accum += (hbar * sqrtF) / (2.0 * denom_op) * (S_A / denom_S) * (I_AB / denom_I)

        # Trace-norm distance d(t) sampled with cap to avoid heavy runs
        trn = None
        if n_tau < tau_max_samples:
            rho_Abar = _partial_trace_rho(rho, L, Abar_sites)
            rho_prod = np.kron(rho_A, rho_Abar)
            diff = rho - rho_prod
            # Nuclear norm via SVD (heavy; keep capped by tau_max_samples)
            s = np.linalg.svd(diff, compute_uv=False)
            trn = float(np.sum(np.abs(s)))
            d4_samples.append(trn ** 4)
            n_tau += 1

        good = 0
        frames.append(Frame(t=t, S_A_nats=S_A, dH_A=dH, Pi_t=Pi_t, good_time=good, Pi_fluct=Pi_fl, trace_norm_d=trn,
                          I_AB=I_AB, Pi_rate=Pi_rate, Pi_eng=Pi_eng))

    # Aggregates
    n_frames = len(frames)
    Pi_vals = [f.Pi_t for f in frames]
    tail = max(1, n_frames // 10)
    rm = [np.mean(Pi_vals[i:]) for i in range(0, n_frames - tail + 1)]
    Pi_liminf_proxy = float(np.min(rm)) if rm else float(np.mean(Pi_vals) if Pi_vals else 0.0)

    Pi_fluct_vals = [f.Pi_fluct for f in frames if f.Pi_fluct is not None]
    Pi_fluct_avg = float(np.mean(Pi_fluct_vals)) if Pi_fluct_vals else None
    Pi_fluct_liminf_proxy = None
    if Pi_fluct_vals:
        tail_fl = max(1, len(Pi_fluct_vals) // 10)
        rm_fl = [np.mean(Pi_fluct_vals[i:]) for i in range(0, len(Pi_fluct_vals) - tail_fl + 1)]
        Pi_fluct_liminf_proxy = float(np.min(rm_fl)) if rm_fl else float(np.mean(Pi_fluct_vals))

    tau_A4_avg = float(np.mean(d4_samples)) if d4_samples else None
    tau_A4_liminf_proxy = None
    if d4_samples:
        tail_d = max(1, len(d4_samples) // 10)
        rm_d4 = [float(np.mean(d4_samples[i:])) for i in range(0, len(d4_samples) - tail_d + 1)]
        tau_A4_liminf_proxy = float(np.min(rm_d4)) if rm_d4 else tau_A4_avg

    frac_good = 0.0
    # NOTE / ATTENTION (floors estimation):
    # The blocks below compute proxy floors for (ε_F, ε_S, ε_I) and δ using
    # simplified series. This keeps runtime light but can over/under-estimate
    # the speed-based RHS because RHS_speed ∝ δ·ε_F·ε_S·ε_I.
    # Recommended for production runs:
    #   1) Accumulate full time-series for √F_Q, S_A, and I(A:Ā).
    #   2) Define the active epoch K* = {t : all three legs exceed small baselines}.
    #   3) Compute positive-support quantiles (e.g., q10 over K*) for ε_F, ε_S, ε_I.
    #   4) Set δ = fraction of times in K* (or a rolling-window min to guard transients).
    # This avoids “no floor” pathologies when inactive frames dominate and stabilizes
    # results without changing the theory.
    # Floors and RHS (speed-based)
    # Use canonical √F_Q without extra rescaling; RHS_speed already carries (ħ/2)
    sqrtF_seq = [ math.sqrt(max(qfi_sld_unitary(_partial_trace_rho(np.outer(tebd.pt.to_dense(), tebd.pt.to_dense().conj()), L, A_sites), H_A_on_A), 0.0)) ] if frames else []
    # Use collected sequences for floors
    # In practice, recompute cheaply from stored frames; here use placeholders from Pi_rate and S_A
    S_seq = [fr.S_A_nats for fr in frames]
    I_seq = [fr.I_AB for fr in frames if fr.I_AB is not None]
    def q10(x):
        arr = np.array(x, dtype=float)
        return float(np.quantile(arr, 0.10)) if arr.size else 0.0
    epsF = q10(sqrtF_seq)
    epsS = q10(S_seq)
    epsI = q10(I_seq)
    good_mask = (np.array(sqrtF_seq) >= epsF) & (np.array(S_seq) >= epsS) & (np.array(I_seq) >= epsI)
    delta = float(np.mean(good_mask)) if len(good_mask) else 0.0
    rhs_speed = (hbar / 2.0) * epsF * epsS * epsI * delta  # CORRECTED: (ℏ/2) · δ · ε_F · ε_S · ε_I

    # Hamiltonian-capped RHS
    triad_hat_mean = (triad_hat_accum / len(frames)) if frames else None
    rhs_capped = None
    if triad_hat_mean is not None:
        rhs_capped = 2.0 * (log_dmin ** 2) * H_opnorm * float(triad_hat_mean)

    summary = Summary(L=L, ell=ell, J=J, h=h, dt=dt, steps=steps, chi=chi, cut=cut, hbar=hbar,
                      sample_k=sample_k, tau_max_samples=tau_max_samples,
                      T_total=(frames[-1].t if frames else 0.0), n_frames=n_frames, n_tau_samples=n_tau,
                      Pi_avg=float(np.mean(Pi_vals)) if Pi_vals else 0.0, Pi_liminf_proxy=Pi_liminf_proxy,
                      frac_good_times=frac_good,
                      Pi_fluct_avg=Pi_fluct_avg, Pi_fluct_liminf_proxy=Pi_fluct_liminf_proxy,
                      tau_A4_avg=tau_A4_avg, tau_A4_liminf_proxy=tau_A4_liminf_proxy,
                      epsF=epsF, epsS=epsS, epsI=epsI, delta=delta, rhs_speed=rhs_speed,
                      triad_hat_mean=float(triad_hat_mean) if triad_hat_mean is not None else None,
                      rhs_capped=float(rhs_capped) if rhs_capped is not None else None)

    return frames, summary


# -------------------------------
# I/O helpers
# -------------------------------

def save_csv(frames: List[Frame], prefix: str, L: int, ell: int, chi: int, dt: float, steps: int, sample_k: int) -> str:
    path = f"{prefix}_L{L}_ell{ell}_chi{chi}_dt{dt}_steps{steps}_k{sample_k}.csv"
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(Frame.__annotations__.keys()))
        w.writeheader(); [w.writerow(asdict(fr)) for fr in frames]
    return path


def save_json(summary: Summary, prefix: str, L: int, ell: int, chi: int, dt: float, steps: int, sample_k: int) -> str:
    path = f"{prefix}_L{L}_ell{ell}_chi{chi}_dt{dt}_steps{steps}_k{sample_k}_summary.json"
    with open(path, 'w') as f: json.dump(asdict(summary), f, indent=2)
    return path


def save_plot(frames: List[Frame], prefix: str, L: int, ell: int, chi: int, dt: float, steps: int, sample_k: int) -> Optional[str]:
    if plt is None or not frames: return None
    t = [fr.t for fr in frames]; SA = [fr.S_A_nats for fr in frames]; dH = [fr.dH_A for fr in frames]; Pi = [fr.Pi_t for fr in frames]
    fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    ax[0].plot(t, SA, 'b-'); ax[0].set_ylabel('S_A (nats)'); ax[0].grid(True, alpha=0.3)
    ax[1].plot(t, dH, 'r-'); ax[1].set_ylabel('ΔH_A'); ax[1].grid(True, alpha=0.3)
    ax[2].plot(t, Pi, 'g-'); ax[2].set_ylabel('Π_proxy(t)'); ax[2].set_xlabel('t'); ax[2].grid(True, alpha=0.3)
    fig.suptitle(f"ESSE TEBD Néel+τ (L={L}, ell={ell}, chi={chi}, dt={dt}, k={sample_k})"); plt.tight_layout()
    out_path = f"{prefix}_L{L}_ell{ell}_chi{chi}_dt{dt}_steps{steps}_k{sample_k}_plot.png"; fig.savefig(out_path, dpi=150); plt.close(fig); return out_path


# -------------------------------
# CLI
# -------------------------------
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='ESSE TEBD Néel state with τ sampling')
    ap.add_argument('--L', type=int, default=10)
    ap.add_argument('--ell', type=int, default=2)
    ap.add_argument('--J', type=float, default=1.0)
    ap.add_argument('--h', type=float, default=0.0)
    ap.add_argument('--dt', type=float, default=0.05)
    ap.add_argument('--steps', type=int, default=100)
    ap.add_argument('--chi', type=int, default=128)
    ap.add_argument('--cut', type=float, default=1e-9)
    ap.add_argument('--hbar', type=float, default=1.0)
    ap.add_argument('--sample-k', type=int, default=5, dest='sample_k')
    ap.add_argument('--tau-max-samples', type=int, default=60, dest='tau_max_samples')
    ap.add_argument('--out-prefix', type=str, default='results/esse_tdmrg_neel_tau')
    args = ap.parse_args()

    print('=== ESSE TEBD Néel + τ run ==='); print(vars(args))
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)

    frames, summary = run_neel_tau(args.L, args.ell, args.J, args.h, args.dt, args.steps,
                                   args.chi, args.cut, args.hbar, args.sample_k, args.tau_max_samples,
                                   args.out_prefix)

    csv_path = save_csv(frames, args.out_prefix, args.L, args.ell, args.chi, args.dt, args.steps, args.sample_k)
    json_path = save_json(summary, args.out_prefix, args.L, args.ell, args.chi, args.dt, args.steps, args.sample_k)
    plot_path = save_plot(frames, args.out_prefix, args.L, args.ell, args.chi, args.dt, args.steps, args.sample_k)

    print(f'Saved frames CSV: {csv_path}')
    print(f'Saved summary JSON: {json_path}')
    if plot_path: print(f'Saved plot PNG: {plot_path}')
    print('Summary:'); print(json.dumps(asdict(summary), indent=2)) 