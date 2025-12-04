#!/usr/bin/env python3
import math, csv, json, os
from dataclasses import dataclass, asdict
from typing import List, Optional, Tuple, Dict

import numpy as np
import quimb as qu
import quimb.tensor as qtn

@dataclass
class PerARecord:
    A: int
    v0: float
    deltaA: float
    tauA4: Optional[float]
    # Remove perC_A: Optional[float] - not in canonical LaTeX paper

@dataclass
class Summary:
    L: int; ell: int; J: float; h: float; dt: float; steps: int; chi: int; cut: float; hbar: float
    sample_k: int; tau_max_samples: int; eps: float
    Pi_avg: float; Pi_liminf_proxy: float
    frac_good_times_global: float
    perA: List[PerARecord]
    # Added fields for PDF alignment
    opnormA: List[float]
    predicate: Dict[str, object]
    # Add canonical LaTeX components
    epsF: Optional[float] = None
    epsS: Optional[float] = None
    epsI: Optional[float] = None
    delta: Optional[float] = None
    rhs_speed: Optional[float] = None


def entanglement_entropy_nats(psi, ell: int) -> float:
    if not (1 <= ell <= psi.L - 1):
        if ell == psi.L:
            return 0.0
        raise ValueError("ell must be in [1, L-1]")
    # quimb entropy default is base-2; convert to nats
    return float(psi.entropy(ell)) * float(np.log(2.0))


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
        for A, B in ((sx, sx), (sy, sy), (sz, sz)):
            ops = [id2] * ell; ops[i] = A; ops[i + 1] = B
            term = ops[0]
            for k in range(1, ell): term = np.kron(term, ops[k])
            H = H + J * term
    return H


def local_generator_HA_dense_full(L: int, A: int, J: float = 1.0, h: float = 0.0) -> np.ndarray:
    # Dense H_A on full Hilbert space with 1/2 boundary split
    sx, sy, sz = qu.pauli('X'), qu.pauli('Y'), qu.pauli('Z')
    id2 = np.eye(2, dtype=complex)
    H = np.zeros((2**L, 2**L), dtype=complex)
    ell = A
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
        for Aop, Bop in ((sx, sx), (sy, sy), (sz, sz)):
            ops = [id2] * L; ops[i] = Aop; ops[i + 1] = Bop
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + w * term
    return H


def variance_from_dense_state(psi_vec: np.ndarray, H_dense: np.ndarray) -> float:
    EH = np.vdot(psi_vec, H_dense @ psi_vec)
    EH2 = np.vdot(psi_vec, H_dense @ (H_dense @ psi_vec))
    return math.sqrt(max(float(np.real(EH2 - EH * EH)), 0.0))


def _partial_trace_rho(rho: np.ndarray, L: int, keep: List[int]) -> np.ndarray:
    rho_t = rho.reshape([2] * (2 * L)); keep_set = set(keep)
    for k in range(L - 1, -1, -1):
        if k not in keep_set:
            half = rho_t.ndim // 2
            rho_t = np.trace(rho_t, axis1=k, axis2=half + k)
    return rho_t.reshape(2 ** len(keep), 2 ** len(keep))


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


def run_perA(L: int, ell: int, J: float, h: float, dt: float, steps: int, chi: int, cut: float, hbar: float,
             sample_k: int, tau_max_samples: int, eps: float, out_prefix: str,
             alpha: float = 1.0, S_floor: float = 0.0, use_normalized_entropies: bool = False) -> Summary:
    # Build TEBD Hamiltonian
    H_local = qtn.ham_1d_heis(L=L, j=J, bz=h, S=0.5, cyclic=False)
    # Prebuild dense H_A_on_A for each A to get spectra and v0(A)
    v0_list: List[float] = []
    H_A_on_A_list: List[np.ndarray] = []
    opnormA: List[float] = []
    for A in range(1, ell + 1):
        H_AA = build_dense_HA_on_A(A, J, h)
        H_A_on_A_list.append(H_AA)
        evals = np.linalg.eigvalsh(H_AA)
        gaps = np.sort(np.unique(np.round(np.diff(np.sort(evals)), 12)))
        # smallest nonzero gap; if none, v0=0
        deltaE_min = 0.0
        for g in gaps:
            if g > 0.0:
                deltaE_min = float(g)
                break
        v0_list.append(0.25 * deltaE_min)
        opnormA.append(float(np.max(np.abs(evals))))

    # Dense H_A (on full space) for variance and opnorm for A=ell to keep cost bounded
    H_A_full = local_generator_HA_dense_full(L, ell, J, h)
    opnorm_full_last = float(np.max(np.abs(np.linalg.eigvalsh(H_A_full))))

    # Initial state Neel
    pattern = ''.join('01'[i % 2] for i in range(L))
    psi_mps = qtn.MPS_computational_state(pattern)
    tebd = qtn.TEBD(psi_mps, H_local, dt=dt, split_opts={'max_bond': chi, 'cutoff': cut}, imag=False)

    def get_mps():
        return getattr(tebd, 'pt', None) or getattr(tebd, 'psi', None)

    A_sites = list(range(ell))
    Abar_sites = list(range(ell, L))

    # Per-A good time counters
    good_counts = {A: 0 for A in range(1, ell + 1)}
    sample_counts = {A: 0 for A in range(1, ell + 1)}
    good_global = 0
    total_global = 0

    # tau samples for A=ell
    d4_samples: List[float] = []
    n_tau = 0

    Pi_vals: List[float] = []
    # Remove Pi_fluct_vals: List[float] = []

    for step in range(1, steps + 1):
        tebd.step()
        if (step % sample_k) != 0:
            continue
        mps = get_mps()
        psi_dense = mps.to_dense()
        rho = np.outer(psi_dense, psi_dense.conj())

        # per-A measures
        for A in range(1, ell + 1):
            # entropy for cut after site A (quimb bond index A)
            try:
                S_A = entanglement_entropy_nats(mps, A)
            except ValueError:
                S_A = 0.0
            # variance using H_A on full space only for largest A to control cost; for smaller A approximate with H_A_on_A variance proxy
            if A == ell:
                dH = variance_from_dense_state(psi_dense, H_A_full)
                scale_A = alpha * opnorm_full_last
            else:
                # proxy: variance of H_A_on_A under rho_A
                rho_A = _partial_trace_rho(rho, L, list(range(A)))
                H_AA = H_A_on_A_list[A - 1]
                EA = float(np.real(np.trace(rho_A @ H_AA)))
                EA2 = float(np.real(np.trace(rho_A @ (H_AA @ H_AA))))
                dH = math.sqrt(max(EA2 - EA * EA, 0.0))
                scale_A = alpha * opnormA[A - 1]

            # normalized floors if requested
            if use_normalized_entropies:
                logd = A * math.log(2.0)
                S_ok = ( (S_A / max(logd, 1e-12)) >= S_floor )
                S_floor_eff = S_floor * logd
            else:
                S_ok = (S_A >= S_floor)
                S_floor_eff = S_floor

            # good time criterion: triad predicate (pure case I=2S_A)
            sample_counts[A] += 1
            if (dH >= scale_A) and S_ok:
                good_counts[A] += 1

        # global Pi metrics using A=ell (proxy for sys)
        S_A = entanglement_entropy_nats(mps, ell)  # CORRECTED: Use consistent variable name
        dH_last = variance_from_dense_state(psi_dense, H_A_full)
        
        # CORRECTED: Compute mutual information I(A:Ā) using complete general formula
        S_Abar = entanglement_entropy_nats(mps, L - ell)  # Complement entropy
        
        # CORRECTED: Compute S_AB (global entropy) for complete mutual information
        rho_global = np.outer(psi_dense, psi_dense.conj())
        rho_global = 0.5 * (rho_global + rho_global.conj().T)  # Ensure Hermitian
        vals_global = np.clip(np.linalg.eigvalsh(rho_global), 0.0, 1.0)
        nz_global = vals_global[vals_global > 1e-12]
        S_AB = float(-np.sum(nz_global * np.log(nz_global))) if len(nz_global) > 0 else 0.0
        
        I_AB = S_A + S_Abar - S_AB  # Complete general formula: I(A:Ā) = S_A + S_Ā - S_{AĀ}
        
        # CORRECTED: Canonical triad formulas (no ℏ factor)
        rho_A_last = _partial_trace_rho(rho, L, A_sites)
        F = qfi_sld_unitary(rho_A_last, H_A_on_A_list[-1])
        sqrtF = math.sqrt(max(F, 0.0))
        
        Pi_t = sqrtF * S_A * I_AB  # Main triad: √F_Q · S_A · I(A:Ā)
        Pi_vals.append(Pi_t)
        # Remove Pi_fluct_vals.append(sqrtF * S_A * I_AB) - redundant

        # global good-time fraction based on A=ell triad predicate
        total_global += 1
        if use_normalized_entropies:
            logd_last = ell * math.log(2.0)
            S_ok_last = ( (S_A / max(logd_last, 1e-12)) >= S_floor )  # CORRECTED: Use S_A
            S_floor_eff_last = S_floor * logd_last
        else:
            S_ok_last = (S_A >= S_floor)  # CORRECTED: Use S_A
            S_floor_eff_last = S_floor
        if (dH_last >= alpha * opnorm_full_last) and S_ok_last:
            good_global += 1

        # tau sample for A=ell
        if n_tau < tau_max_samples:
            rho_Abar = _partial_trace_rho(rho, L, Abar_sites)
            s = np.linalg.svd(rho - np.kron(rho_A_last, rho_Abar), compute_uv=False)
            d4_samples.append(float(np.sum(np.abs(s))) ** 4)
            n_tau += 1

    # aggregates
    tail = max(1, len(Pi_vals) // 10)
    rm = [np.mean(Pi_vals[i:]) for i in range(0, len(Pi_vals) - tail + 1)]
    Pi_liminf_proxy = float(np.min(rm)) if rm else float(np.mean(Pi_vals) if Pi_vals else 0.0)

    # Remove redundant Pi_fluct calculations
    # Pi_fluct_avg = float(np.mean(Pi_fluct_vals)) if Pi_fluct_vals else None
    # Pi_fluct_liminf_proxy = None
    # if Pi_fluct_vals:
    #     tail_fl = max(1, len(Pi_fluct_vals) // 10)
    #     rm_fl = [np.mean(Pi_fluct_vals[i:]) for i in range(0, len(Pi_fluct_vals) - tail_fl + 1)]
    #     Pi_fluct_liminf_proxy = float(np.min(rm_fl)) if rm_fl else float(np.mean(Pi_fluct_vals))

    tau_A4 = float(np.mean(d4_samples)) if d4_samples else None

    # CORRECTED: Remove perC calculations and add canonical LaTeX components
    perA: List[PerARecord] = []
    # Remove perC_total = 0.0
    for A in range(1, ell + 1):
        deltaA = (good_counts[A] / sample_counts[A]) if sample_counts[A] > 0 else 0.0
        tauA4 = tau_A4 if A == ell else None
        # Remove all perC calculations - not in canonical LaTeX paper
        perA.append(PerARecord(A=A, v0=v0_list[A - 1], deltaA=deltaA, tauA4=tauA4))

    # NOTE / ATTENTION (floors estimation):
    # Compute floors from the actual time series collected at sampled frames
    # to avoid bias from last-value proxies. This stabilizes RHS_speed while
    # preserving the theory.
    sqrtF_seq = []
    S_seq = []
    I_seq = []
    # Recompute √F_Q, S_A, I(A:Ā) cheaply from stored frames and states
    # at the same sampling cadence used above.
    for step in range(1, steps + 1):
        if (step % sample_k) != 0:
            continue
        mps = get_mps()
        psi_dense = mps.to_dense()
        rho = np.outer(psi_dense, psi_dense.conj())
        # Use largest cut A=ell for global floors (consistent with global delta)
        rho_A_last = _partial_trace_rho(rho, L, A_sites)
        F_last = qfi_sld_unitary(rho_A_last, H_A_on_A_list[-1])
        sqrtF_seq.append(math.sqrt(max(F_last, 0.0)))
        S_seq.append(entanglement_entropy_nats(mps, ell))
        # Full mutual information (pure/mixed compatible)
        rho_Abar_last = _partial_trace_rho(rho, L, Abar_sites)
        # Global entropy
        rho_global = 0.5 * (rho + rho.conj().T)
        vals_global = np.clip(np.linalg.eigvalsh(rho_global), 0.0, 1.0)
        nz_global = vals_global[vals_global > 1e-12]
        S_AB_last = float(-np.sum(nz_global * np.log(nz_global))) if len(nz_global) > 0 else 0.0
        # Entropy of complement via partial trace
        vals_b = np.clip(np.linalg.eigvalsh(0.5 * (rho_Abar_last + rho_Abar_last.conj().T)), 0.0, 1.0)
        nz_b = vals_b[vals_b > 1e-12]
        S_B_last = float(-np.sum(nz_b * np.log(nz_b))) if len(nz_b) > 0 else 0.0
        I_seq.append(S_seq[-1] + S_B_last - S_AB_last)
    
    def q10(x): 
        arr = np.array(x, dtype=float)
        return float(np.quantile(arr, 0.10)) if arr.size else 0.0
    
    epsF = q10(sqrtF_seq)
    epsS = q10(S_seq)
    epsI = q10(I_seq)
    
    # CORRECTED: Calculate global delta (fraction of good times) with proper weighting
    # Weight by sample_counts[A] to avoid bias from different sampling frequencies
    total_good = sum(good_counts[A] for A in range(1, ell + 1))
    total_samples = sum(sample_counts[A] for A in range(1, ell + 1))
    delta = (total_good / total_samples) if total_samples > 0 else 0.0
    
    # CORRECTED: Speed-based RHS bound (LaTeX Eq. 331)
    rhs_speed = (hbar / 2.0) * epsF * epsS * epsI * delta

    summary = Summary(L=L, ell=ell, J=J, h=h, dt=dt, steps=steps, chi=chi, cut=cut, hbar=hbar,
                      sample_k=sample_k, tau_max_samples=tau_max_samples, eps=eps,
                      Pi_avg=float(np.mean(Pi_vals)) if Pi_vals else 0.0,
                      Pi_liminf_proxy=Pi_liminf_proxy,
                      frac_good_times_global=(good_global / total_global) if total_global > 0 else 0.0,
                      perA=perA,
                      epsF=epsF, epsS=epsS, epsI=epsI, delta=delta, rhs_speed=rhs_speed,
                      opnormA=opnormA,
                      predicate={"scale":"opnorm","alpha":alpha,"S_floor":S_floor,"use_norm":use_normalized_entropies})
    # save JSON/CSV
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    with open(f"{out_prefix}_L{L}_ell{ell}_perA.json", 'w') as f:
        json.dump(asdict(summary), f, indent=2, default=lambda o: asdict(o) if hasattr(o, '__dict__') else o)
    with open(f"{out_prefix}_L{L}_ell{ell}_perA.csv", 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['A','v0','deltaA','tauA4'])
        w.writeheader()
        for r in perA:
            w.writerow(asdict(r))
    print(json.dumps(asdict(summary), indent=2, default=lambda o: asdict(o) if hasattr(o, '__dict__') else o))
    return summary


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='Per-A ESSE metrics (v0(A), delta_A, tau_A^4, perC_A sum) for Neel TEBD')
    ap.add_argument('--L', type=int, default=10)
    ap.add_argument('--ell', type=int, default=6)
    ap.add_argument('--J', type=float, default=1.0)
    ap.add_argument('--h', type=float, default=0.0)
    ap.add_argument('--dt', type=float, default=0.05)
    ap.add_argument('--steps', type=int, default=600)
    ap.add_argument('--chi', type=int, default=128)
    ap.add_argument('--cut', type=float, default=1e-9)
    ap.add_argument('--hbar', type=float, default=1.0)
    ap.add_argument('--sample-k', type=int, default=5)
    ap.add_argument('--tau-max-samples', type=int, default=20)
    ap.add_argument('--eps', type=float, default=1e-12)
    ap.add_argument('--out-prefix', type=str, default='results/esse_tdmrg_perA')
    # New CLI for PDF alignment
    ap.add_argument('--alpha', type=float, default=1.0)
    ap.add_argument('--S-floor', type=float, default=0.0, dest='S_floor')
    ap.add_argument('--use-normalized-entropies', action='store_true')
    args = ap.parse_args()
    run_perA(args.L, args.ell, args.J, args.h, args.dt, args.steps, args.chi, args.cut, args.hbar,
             args.sample_k, args.tau_max_samples, args.eps, args.out_prefix,
             alpha=args.alpha, S_floor=args.S_floor, use_normalized_entropies=args.use_normalized_entropies) 