#!/usr/bin/env python3
"""
ESSE lattice ball-K_A proxy (Heisenberg chain, open)
- Build standard H_A (1/2 boundary split) as dense operator (for ΔH_A standard)
- Build weighted K_A^lat on FULL space: on-site and bonds inside A weighted by conformal-like w(i)
  * w(i) = (R^2 - x_i^2) / (2R), sites i in A, with x_i centered around block midpoint
  * internal bond (i,i+1) weight = 0.5*(w(i)+w(i+1))
  * boundary bond to B (i in A, i+1 in B) weight = 0.5*w(i)
- Run TEBD (quimb) with Néel or zeros initial state (configurable)
- Record Π_std=(4/ħ) ΔH_A_std S_A^2, Π_ball=(4/ħ) ΔK_A_ball S_A^2, and Π_fluct_ball via QFI on A-only weighted op
Outputs: CSV, JSON, PNG
"""

import math
import csv
import json
import os
from dataclasses import dataclass, asdict
from typing import List, Optional, Tuple
from typing import Dict

import numpy as np
import quimb as qu
import quimb.tensor as qtn

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

@dataclass
class Frame:
    t: float
    S_A_nats: float
    dH_std: float
    dK_ball: float
    Pi_std: float
    Pi_ball: float
    Pi_fluct_ball: Optional[float] = None
    sqrtF_sld: Optional[float] = None
    sqrtF_bures: Optional[float] = None
    triad_pure: Optional[float] = None
    hatPsi: Optional[float] = None

@dataclass
class Summary:
    L: int; ell: int; J: float; h: float; dt: float; steps: int; chi: int; cut: float; hbar: float
    init: str
    T_total: float; n_frames: int
    Pi_std_avg: float; Pi_std_liminf: float
    Pi_ball_avg: float; Pi_ball_liminf: float
    Pi_fluct_ball_avg: Optional[float]; Pi_fluct_ball_liminf: Optional[float]
    hatPsi_avg: Optional[float] = None
    LHS_avg: Optional[float] = None
    ham_RHS: Optional[float] = None
    speed_RHS: Optional[float] = None
    H_A_op_norm: Optional[float] = None
    log_dA: Optional[float] = None
    log_dmin: Optional[float] = None
    delta_active: Optional[float] = None
    eps_F: Optional[float] = None
    eps_S: Optional[float] = None
    eps_I: Optional[float] = None

# New config dataclass for quench sweep
@dataclass
class QuenchCfg:
    enable: bool = False
    epsilons: Tuple[float, ...] = (1e-4, 5e-4, 1e-3)
    t_linear_max: float = 0.6
    hx_weight_in_KA: float = 0.0  # add X-weight into K_A to align with X quench
    normalize_KA: bool = True


def entanglement_entropy_nats(psi, ell: int) -> float:
    if not (1 <= ell <= psi.L - 1):
        if ell == psi.L:
            return 0.0
        raise ValueError("ell must be in [1, L-1]")
    # quimb entropy default is base-2; convert to nats
    return float(psi.entropy(ell - 1)) * float(np.log(2.0))


def _dense_paulis():
    sx = qu.pauli('X'); sy = qu.pauli('Y'); sz = qu.pauli('Z'); id2 = np.eye(2, dtype=complex)
    return sx, sy, sz, id2


def build_dense_HA_standard(L: int, ell: int, J: float = 1.0, h: float = 0.0) -> np.ndarray:
    sx, sy, sz, id2 = _dense_paulis()
    H = np.zeros((2**L, 2**L), dtype=complex)
    if h != 0.0:
        for i in range(ell):
            ops = [id2] * L; ops[i] = sz
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + h * term
    for i in range(L - 1):
        left_in_A = (i < ell); right_in_A = (i + 1 < ell)
        if left_in_A and right_in_A:
            w = J
        elif (left_in_A and not right_in_A) or (not left_in_A and right_in_A):
            w = 0.5 * J
        else:
            w = 0.0
        if w == 0.0: continue
        for A, B in ((sx, sx), (sy, sy), (sz, sz)):
            ops = [id2] * L; ops[i] = A; ops[i + 1] = B
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + w * term
    return H


def build_dense_KA_ball(L: int, ell: int, J: float = 1.0, h: float = 0.0) -> np.ndarray:
    """Lattice ball-K_A proxy on FULL space with conformal-like weights over A.
    Weights: w(i) = (R^2 - x_i^2) / (2R), where i=0..ell-1 and x_i centered at block midpoint.
    Apply weights to on-site in A; bonds inside A use average of endpoint weights; boundary bonds 1/2*weight(A-end).
    """
    sx, sy, sz, id2 = _dense_paulis()
    H = np.zeros((2**L, 2**L), dtype=complex)
    # coordinates within A
    R = (ell - 1) / 2.0 if ell > 1 else 0.5
    xs = np.array([i - R for i in range(ell)], dtype=float)
    def w_val(idx_in_A: int) -> float:
        x = xs[idx_in_A]
        return float((R*R - x*x) / (2.0 * R)) if R > 0 else 1.0
    # on-site in A
    if h != 0.0:
        for i in range(ell):
            w = w_val(i)
            if abs(w) < 1e-15: continue
            ops = [id2] * L; ops[i] = sz
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + (w * h) * term
    # bonds
    for i in range(L - 1):
        left_in_A = (i < ell); right_in_A = (i + 1 < ell)
        if left_in_A and right_in_A:
            w = 0.5 * (w_val(i) + w_val(i + 1))
            coeff = J * w
        elif left_in_A and not right_in_A:
            coeff = 0.5 * J * w_val(i)
        elif right_in_A and not left_in_A:
            coeff = 0.5 * J * w_val(i + 1)
        else:
            coeff = 0.0
        if coeff == 0.0: continue
        for A, B in ((sx, sx), (sy, sy), (sz, sz)):
            ops = [id2] * L; ops[i] = A; ops[i + 1] = B
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + coeff * term
    return H

def build_KA_mpo_weighted(L: int, ell: int, J: float, two_pi: bool = True,
                          hx_weight: float = 0.0, normalize: bool = True,
                          include_hz: bool = False, hz: float = 0.0) -> Tuple[object, float, float]:
    """Build K_A as an MPO with conformal-like weights and optional 2π factor.
    Returns (MPO, sum_w, factor_used).
    """
    lh = qtn.LocalHam1D(L=L, cyclic=False)
    R = (ell - 1) / 2.0 if ell > 1 else 0.5
    xs = np.array([i - R for i in range(ell)], dtype=float)
    def w(i: int) -> float:
        x = xs[i]
        return float((R * R - x * x) / (2.0 * R)) if R > 0 else 1.0
    wsum = sum(max(w(i), 0.0) for i in range(ell))
    norm = (wsum if (normalize and wsum > 1e-15) else 1.0)
    factor = (2.0 * math.pi) if two_pi else 1.0
    # on-site inside A
    for i in range(ell):
        wi = (w(i) / norm) * factor
        if include_hz and abs(hz) > 0:
            lh += (hz * wi, ('Z', i))
        if abs(hx_weight) > 0:
            lh += (hx_weight * wi, ('X', i))
    # bonds weighted
    for i in range(L - 1):
        left_in = (i < ell); right_in = (i + 1 < ell)
        coeff = 0.0
        if left_in and right_in:
            coeff = J * (0.5 * (w(i) + w(i + 1)) / norm) * factor
        elif left_in and not right_in:
            coeff = 0.5 * J * (w(i) / norm) * factor
        elif right_in and not left_in:
            coeff = 0.5 * J * (w(i + 1) / norm) * factor
        if coeff != 0.0:
            lh += (coeff, ('X', i), ('X', i + 1))
            lh += (coeff, ('Y', i), ('Y', i + 1))
            lh += (coeff, ('Z', i), ('Z', i + 1))
    return lh.build_mpo(S=0.5), float(wsum), float(factor)

def build_dense_KA_ball_with_x(L: int, ell: int, J: float, hz: float = 0.0,
                               hx_weight: float = 0.0, normalize: bool = True,
                               two_pi: bool = True) -> np.ndarray:
    """Dense K_A with conformal-like weights, optional X-weight, optional 2π factor.
    Includes: weighted bonds (XX,YY,ZZ) and optional on-site Z in A.
    """
    sx, sy, sz, id2 = _dense_paulis()
    H = np.zeros((2**L, 2**L), dtype=complex)
    R = (ell - 1) / 2.0 if ell > 1 else 0.5
    xs = np.array([i - R for i in range(ell)], dtype=float)
    def w_val(i: int) -> float:
        x = xs[i]
        return float((R * R - x * x) / (2.0 * R)) if R > 0 else 1.0
    wsum = sum(max(w_val(i), 0.0) for i in range(ell))
    norm = (wsum if (normalize and wsum > 1e-15) else 1.0)
    factor = (2.0 * math.pi) if two_pi else 1.0
    # on-site in A
    if abs(hz) > 0:
        for i in range(ell):
            w = (w_val(i) / norm) * factor
            ops = [id2] * L; ops[i] = sz
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + (w * hz) * term
    # optional X-weight in A
    if abs(hx_weight) > 0:
        for i in range(ell):
            w = (w_val(i) / norm) * factor
            ops = [id2] * L; ops[i] = sx
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + (w * hx_weight) * term
    # bonds
    for i in range(L - 1):
        left_in = (i < ell); right_in = (i + 1 < ell)
        if left_in and right_in:
            w = 0.5 * (w_val(i) + w_val(i + 1)) / norm
            coeff = J * w * factor
        elif left_in and not right_in:
            coeff = 0.5 * J * (w_val(i) / norm) * factor
        elif right_in and not left_in:
            coeff = 0.5 * J * (w_val(i + 1) / norm) * factor
        else:
            coeff = 0.0
        if coeff == 0.0: continue
        for Aop, Bop in ((sx, sx), (sy, sy), (sz, sz)):
            ops = [id2] * L; ops[i] = Aop; ops[i + 1] = Bop
            term = ops[0]
            for k in range(1, L): term = np.kron(term, ops[k])
            H = H + coeff * term
    return H


def variance_dense(psi_vec: np.ndarray, H_dense: np.ndarray) -> float:
    E1 = np.vdot(psi_vec, H_dense @ psi_vec)
    E2 = np.vdot(psi_vec, H_dense @ (H_dense @ psi_vec))
    var = float(np.real(E2 - E1 * E1))
    return math.sqrt(max(var, 0.0))

def variance_mpo(psi_mps, H_mpo, chi: int, cut: float) -> float:
    EH = float(qtn.expec_TN_1D(psi_mps, H_mpo))
    phi = H_mpo.apply(psi_mps, normalize=False, cutoff=cut, max_bond=chi)
    EH2 = float(phi.H @ phi)
    return math.sqrt(max(EH2 - EH * EH, 0.0))


def _partial_trace_rho(rho: np.ndarray, L: int, keep: List[int]) -> np.ndarray:
    rho_t = rho.reshape([2] * (2 * L))
    keep_set = set(keep)
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

def bures_angle_dense(rho1: np.ndarray, rho2: np.ndarray) -> float:
    """Return Bures angle L_B( rho1, rho2 ) via Uhlmann fidelity."""
    rho1 = 0.5 * (rho1 + rho1.conj().T)
    rho2 = 0.5 * (rho2 + rho2.conj().T)
    v, U = np.linalg.eigh(rho1)
    v = np.clip(np.real(v), 0.0, None)
    sqrt_rho1 = U @ np.diag(np.sqrt(v)) @ U.conj().T
    M = sqrt_rho1 @ rho2 @ sqrt_rho1
    M = 0.5 * (M + M.conj().T)
    mvals = np.linalg.eigvalsh(M)
    mvals = np.clip(np.real(mvals), 0.0, None)
    tr_sqrt = float(np.sum(np.sqrt(mvals)))
    amp = min(1.0, max(0.0, tr_sqrt))
    return float(math.acos(amp))


def run_ball_ka(L: int, ell: int, J: float, h: float, dt: float, steps: int, chi: int, cut: float, hbar: float,
                init: str, out_prefix: str) -> Tuple[List[Frame], Summary]:
    H_local = qtn.ham_1d_heis(L=L, j=J, bz=h, S=0.5, cyclic=False)
    H_std_dense = build_dense_HA_standard(L, ell, J, h)
    # Prefer MPO K_A (CHM weights, 2π factor, no onsite); fallback to dense
    KA_mpo = None
    ka_wsum = None
    ka_factor = 2.0 * math.pi
    ka_backend = 'dense'
    try:
        KA_mpo, ka_wsum, ka_factor = build_KA_mpo_weighted(L, ell, J, two_pi=True, hx_weight=0.0, normalize=True, include_hz=False, hz=0.0)
        ka_backend = 'MPO'
    except Exception:
        KA_mpo = None
    KA_dense_fallback = build_dense_KA_ball_with_x(L, ell, J, hz=0.0, hx_weight=0.0, normalize=True, two_pi=True)

    # Initial state
    if init == 'neel':
        pattern = ''.join('01'[i % 2] for i in range(L))
        psi = qtn.MPS_computational_state(pattern)
    else:
        psi = qtn.MPS_computational_state('0' * L)

    tebd = qtn.TEBD(psi, H_local, dt=dt, split_opts={'max_bond': chi, 'cutoff': cut}, imag=False)

    A_sites = list(range(ell))
    # QFI operator on A-only space, unified with CHM K_A choice (no onsite, 2π, normalized)
    H_ball_A = build_dense_KA_ball_with_x(ell, ell, J, hz=0.0, hx_weight=0.0, normalize=True, two_pi=True)
    # Operator-norm cap and log dimensions for intensive triad
    try:
        evals_A = np.linalg.eigvalsh(0.5 * (H_ball_A + H_ball_A.conj().T))
        H_A_op_norm = float(np.max(np.abs(np.real(evals_A))))
    except Exception:
        H_A_op_norm = float(np.linalg.norm(H_ball_A, ord=2))
    log_dA = float(ell * np.log(2.0))
    log_dmin = float(min(ell, L - ell) * np.log(2.0))

    frames: List[Frame] = []
    t = 0.0
    prev_rho_A = None
    sqrtF_list: List[float] = []
    SA_list: List[float] = []
    I_list: List[float] = []
    LHS_list: List[float] = []
    hatPsi_list: List[float] = []
    for step in range(1, steps + 1):
        tebd.step(); t += dt
        S_A = entanglement_entropy_nats(tebd.pt, ell)
        psi_dense = tebd.pt.to_dense()
        dH_std = variance_dense(psi_dense, H_std_dense)
        # ΔK_A via MPO if available, else dense fallback
        try:
            if KA_mpo is not None:
                dK_ball = variance_mpo(tebd.pt, KA_mpo, chi=chi, cut=cut)
            else:
                dK_ball = variance_dense(psi_dense, KA_dense_fallback)
        except Exception:
            dK_ball = variance_dense(psi_dense, KA_dense_fallback)
        Pi_std = (4.0 / hbar) * dH_std * (S_A ** 2)
        Pi_ball = (4.0 / hbar) * dK_ball * (S_A ** 2)
        rho = np.outer(psi_dense, psi_dense.conj())
        rho_A = _partial_trace_rho(rho, L, A_sites)
        F = qfi_sld_unitary(rho_A, H_ball_A)
        sqrtF = math.sqrt(max(F, 0.0))
        # Bures speed finite-difference estimator (unitary): sqrtF_fd ≈ 2 ΔL_B / Δt
        if prev_rho_A is not None:
            try:
                LB = bures_angle_dense(prev_rho_A, rho_A)
                sqrtF_fd = (2.0 * LB) / max(dt, 1e-30)
            except Exception:
                sqrtF_fd = None
        else:
            sqrtF_fd = None
        prev_rho_A = rho_A
        # Mutual information for pure global state: I = 2 S_A
        I_val = 2.0 * S_A
        Pi_fl = (2.0 / hbar) * sqrtF * (S_A ** 2)
        # Triad (pure-state form): sqrtF * S_A * I
        triad_val = sqrtF * S_A * I_val
        # Intensive triad \hat{Psi}
        denomH = max(H_A_op_norm, 1e-30)
        hatPsi = ((hbar * sqrtF) / (2.0 * denomH)) * (S_A / max(log_dA, 1e-30)) * (I_val / (2.0 * max(log_dmin, 1e-30)))
        frames.append(Frame(t=t, S_A_nats=S_A, dH_std=dH_std, dK_ball=dK_ball, Pi_std=Pi_std, Pi_ball=Pi_ball,
                            Pi_fluct_ball=Pi_fl, sqrtF_sld=sqrtF, sqrtF_bures=sqrtF_fd, triad_pure=triad_val, hatPsi=hatPsi))
        # Accumulate for inequalities
        sqrtF_list.append(sqrtF)
        SA_list.append(S_A)
        I_list.append(I_val)
        LHS_list.append(dK_ball * S_A * I_val)
        hatPsi_list.append(hatPsi)

    n = len(frames)
    def liminf(xs):
        if not xs: return 0.0
        arr = np.array(xs, dtype=float)
        tail = max(1, n // 10)
        rm = [float(np.mean(arr[i:])) for i in range(0, n - tail + 1)]
        return float(np.min(rm)) if rm else float(np.mean(arr))

    Pi_std_vals = [f.Pi_std for f in frames]
    Pi_ball_vals = [f.Pi_ball for f in frames]
    Pi_fl_vals = [f.Pi_fluct_ball for f in frames]

    # Inequality aggregates
    hatPsi_avg = float(np.mean(hatPsi_list)) if hatPsi_list else 0.0
    LHS_avg = float(np.mean(LHS_list)) if LHS_list else 0.0
    ham_RHS = 2.0 * (log_dmin ** 2) * H_A_op_norm * hatPsi_avg
    sqrtF_arr = np.array(sqrtF_list, dtype=float)
    SA_arr = np.array(SA_list, dtype=float)
    I_arr = np.array(I_list, dtype=float)
    # Robust floors via 10th percentile on positive support
    def qfloor(arr):
        arrp = arr[arr > 0.0]
        return float(np.quantile(arrp, 0.1)) if arrp.size else 0.0
    eps_F = qfloor(sqrtF_arr)
    eps_S = qfloor(SA_arr)
    eps_I = qfloor(I_arr)
    active_mask = (sqrtF_arr >= eps_F) & (SA_arr >= eps_S) & (I_arr >= eps_I)
    delta_active = float(np.mean(active_mask)) if active_mask.size else 0.0
    speed_RHS = (0.5 * hbar) * delta_active * eps_F * eps_S * eps_I

    summary = Summary(
        L=L, ell=ell, J=J, h=h, dt=dt, steps=steps, chi=chi, cut=cut, hbar=hbar, init=init,
        T_total=frames[-1].t if frames else 0.0, n_frames=n,
        Pi_std_avg=float(np.mean(Pi_std_vals)), Pi_std_liminf=liminf(Pi_std_vals),
        Pi_ball_avg=float(np.mean(Pi_ball_vals)), Pi_ball_liminf=liminf(Pi_ball_vals),
        Pi_fluct_ball_avg=(float(np.mean(Pi_fl_vals)) if Pi_fl_vals else None),
        Pi_fluct_ball_liminf=(liminf(Pi_fl_vals) if Pi_fl_vals else None),
        hatPsi_avg=hatPsi_avg, LHS_avg=LHS_avg, ham_RHS=ham_RHS, speed_RHS=speed_RHS,
        H_A_op_norm=H_A_op_norm, log_dA=log_dA, log_dmin=log_dmin, delta_active=delta_active,
        eps_F=eps_F, eps_S=eps_S, eps_I=eps_I,
    )
    # Attach logging extras
    summary_dict = asdict(summary)
    summary_dict.update({
        'ka_weight_sum': float(ka_wsum) if ka_wsum is not None else None,
        'ka_two_pi_factor': float(ka_factor),
        'ka_normalized': True,
        'ka_hx_weight': 0.0,
        'ka_backend': ka_backend,
        'ka_onsite_included': False,
        'entropy_units': 'nats'
    })
    # Return summary built back from dict for JSON
    # (we will serialize summary_dict directly in save_json call path below if desired)
    return frames, summary

# New: quench sweep with first-law fit
def quench_sweep_first_law(L: int, ell: int, J: float, h: float, dt: float, steps: int, chi: int, cut: float,
                           hbar: float, init: str, out_prefix: str, qcfg: QuenchCfg) -> Optional[Dict[str, float]]:
    if not qcfg.enable:
        return None
    H_local = qtn.ham_1d_heis(L=L, j=J, bz=h, S=0.5, cyclic=False)
    psi0 = qtn.MPS_computational_state(''.join('01'[i % 2] for i in range(L)) if init=='neel' else '0'*L)
    tebd0 = qtn.TEBD(psi0.copy(), H_local, dt=dt, split_opts={'max_bond': chi, 'cutoff': cut}, imag=False)
    # Build K_A with X-weight/normalization
    KA_dense = build_dense_KA_ball_with_x(L, ell, J, h, qcfg.hx_weight_in_KA, qcfg.normalize_KA)
    # Baseline
    states0 = []
    for _ in range(steps): tebd0.step(); states0.append(tebd0.pt.copy())
    # Observable builders
    def S_A(mps):
        try: return float(mps.entropy(ell - 1))
        except: return 0.0
    def exp_KA(mps):
        vec = mps.to_dense(); return float(np.real(np.vdot(vec, KA_dense @ vec)))
    rows = []
    for eps in qcfg.epsilons:
        # Build H1 explicitly as LocalHam1D to avoid API issues
        H1_lh = qtn.LocalHam1D(L=L, cyclic=False)
        # Heisenberg bonds
        for i in range(L - 1):
            H1_lh += (J, ('X', i), ('X', i + 1))
            H1_lh += (J, ('Y', i), ('Y', i + 1))
            H1_lh += (J, ('Z', i), ('Z', i + 1))
        # Onsite Z field
        if abs(h) > 0:
            for i in range(L):
                H1_lh += (h, ('Z', i))
        # Quench: small X field
        for i in range(L):
            H1_lh += (eps, ('X', i))
        H1 = H1_lh.build_mpo(S=0.5)
        tebd1 = qtn.TEBD(psi0.copy(), H1, dt=dt, split_opts={'max_bond': chi, 'cutoff': cut}, imag=False)
        states1 = []
        for _ in range(steps): tebd1.step(); states1.append(tebd1.pt.copy())
        for t_idx in range(steps):
            tb = (t_idx + 1) * dt
            if tb > qcfg.t_linear_max: break
            m0, m1 = states0[t_idx], states1[t_idx]
            dS = S_A(m1) - S_A(m0)
            dK = exp_KA(m1) - exp_KA(m0)
            rows.append((tb, eps, dS, dK))
    if not rows:
        return None
    arr = np.array(rows, dtype=float)
    x = arr[:, 3]; y = arr[:, 2]
    # Linear fit y ~ s x
    denom = float(np.dot(x, x)) if x.size else 1.0
    slope = float(np.dot(x, y) / max(denom, 1e-30))
    yhat = slope * x
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2)) if y.size > 1 else 0.0
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    # Save CSV and scatter
    csv_path = f"{out_prefix}_first_law_quench.csv"
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['t','eps','dS','dK']); [w.writerow(list(r)) for r in rows]
    try:
        import matplotlib.pyplot as plt
        png_path = f"{out_prefix}_first_law_scatter.png"
        plt.figure(figsize=(6,5)); plt.scatter(x, y, s=10, alpha=0.6)
        if x.size>0:
            xx = np.linspace(x.min(), x.max(), 100); plt.plot(xx, slope*xx, 'r--', label=f'slope={slope:.3f}, R2={r2:.3f}')
        plt.xlabel('Δ⟨K_A⟩'); plt.ylabel('ΔS_A'); plt.legend(); plt.grid(True, alpha=0.3)
        plt.tight_layout(); plt.savefig(png_path, dpi=150); plt.close()
    except Exception:
        pass
    return dict(slope=slope, r2=r2)


def save_csv(frames: List[Frame], prefix: str, L: int, ell: int, chi: int, dt: float, steps: int) -> str:
    path = f"{prefix}_L{L}_ell{ell}_chi{chi}_dt{dt}_steps{steps}.csv"
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(Frame.__annotations__.keys()))
        w.writeheader(); [w.writerow(asdict(fr)) for fr in frames]
    return path


def save_json(summary: Summary, prefix: str, L: int, ell: int, chi: int, dt: float, steps: int) -> str:
    path = f"{prefix}_L{L}_ell{ell}_chi{chi}_dt{dt}_steps{steps}_summary.json"
    with open(path, 'w') as f:
        json.dump(asdict(summary), f, indent=2)
    return path


def save_plot(frames: List[Frame], prefix: str, L: int, ell: int, chi: int, dt: float, steps: int) -> Optional[str]:
    if plt is None or not frames: return None
    t = [fr.t for fr in frames]
    SA = [fr.S_A_nats for fr in frames]
    stdP = [fr.Pi_std for fr in frames]
    ballP = [fr.Pi_ball for fr in frames]
    fig, ax = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    ax[0].plot(t, SA, 'b-'); ax[0].set_ylabel('S_A (nats)'); ax[0].grid(True, alpha=0.3)
    ax[1].plot(t, stdP, 'g-', label='Π_std'); ax[1].plot(t, ballP, 'm--', label='Π_ball'); ax[1].set_ylabel('Π'); ax[1].legend(); ax[1].grid(True, alpha=0.3)
    ax[2].plot(t, [fr.dH_std for fr in frames], 'r-', label='ΔH_std'); ax[2].plot(t, [fr.dK_ball for fr in frames], 'k--', label='ΔK_ball'); ax[2].set_ylabel('Δ'); ax[2].set_xlabel('t'); ax[2].legend(); ax[2].grid(True, alpha=0.3)
    fig.suptitle(f"ESSE ball-K_A proxy (L={L}, ell={ell}, chi={chi}, dt={dt}, init=neel)")
    plt.tight_layout()
    out_path = f"{prefix}_L{L}_ell{ell}_chi{chi}_dt{dt}_steps{steps}_plot.png"; fig.savefig(out_path, dpi=150); plt.close(fig); return out_path


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description='ESSE lattice ball-K_A proxy')
    ap.add_argument('--L', type=int, default=10)
    ap.add_argument('--ell', type=int, default=4)
    ap.add_argument('--J', type=float, default=1.0)
    ap.add_argument('--h', type=float, default=0.0)
    ap.add_argument('--dt', type=float, default=0.05)
    ap.add_argument('--steps', type=int, default=300)
    ap.add_argument('--chi', type=int, default=128)
    ap.add_argument('--cut', type=float, default=1e-9)
    ap.add_argument('--hbar', type=float, default=1.0)
    ap.add_argument('--init', choices=['neel','zeros'], default='neel')
    ap.add_argument('--out-prefix', type=str, default='results/esse_ball_ka')
    # quench/first-law options
    ap.add_argument('--quench', action='store_true')
    ap.add_argument('--eps', type=float, nargs='+', default=[1e-4, 5e-4, 1e-3])
    ap.add_argument('--t-lin', type=float, default=0.6)
    ap.add_argument('--ka-hx', type=float, default=0.0)
    ap.add_argument('--ka-norm', action='store_true')
    # CHM/boost and Planck-scale mapping (pure-state diagnostic)
    ap.add_argument('--kappa', type=float, default=1.0, help='Boost surface gravity scale κ')
    ap.add_argument('--c', type=float, default=299792458.0, help='Speed of light (m/s)')
    ap.add_argument('--G', type=float, default=6.67430e-11, help='Newton constant (m^3/kg/s^2)')
    ap.add_argument('--a-meters', type=float, default=0.0, help='Optional lattice spacing a in meters to set κ=c^2/R with R=(ell/2)*a')
    ap.add_argument('--linear-gate', action='store_true', help='Gate κ-test to early-time window t<=t_lin')
    args = ap.parse_args()
    print("=== ESSE ball-K_A proxy run ==="); print(vars(args))
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    # Optional: tie κ to radius if a provided
    if args.a_meters and args.a_meters > 0.0:
        R = 0.5 * args.ell * args.a_meters
        kappa_auto = (args.c ** 2) / max(R, 1e-30)
        print(f"Using κ from R=(ell/2)*a: R={R:.6e} m → κ={kappa_auto:.6e} s^-1")
        args.kappa = kappa_auto
    frames, summary = run_ball_ka(args.L, args.ell, args.J, args.h, args.dt, args.steps, args.chi, args.cut, args.hbar, args.init, args.out_prefix)
    csv_path = save_csv(frames, args.out_prefix, args.L, args.ell, args.chi, args.dt, args.steps)
    json_path = save_json(summary, args.out_prefix, args.L, args.ell, args.chi, args.dt, args.steps)
    plot_path = save_plot(frames, args.out_prefix, args.L, args.ell, args.chi, args.dt, args.steps)
    print(f"Saved frames CSV: {csv_path}")
    print(f"Saved summary JSON: {json_path}")
    if plot_path: print(f"Saved plot PNG: {plot_path}")
    print("Summary:"); print(json.dumps(asdict(summary), indent=2))

    # Report root inequalities (closed-system, pure-state MI assumption I=2 S_A)
    try:
        print("Root inequalities (closed, pure-state MI):")
        print(f"  Hamiltonian-capped: LHS=⟨ΔK_A·S_A·I⟩ ≈ {summary.LHS_avg:.6e}, RHS=2(\log d_min)^2‖H_A‖·⟨\hatΨ⟩ ≈ {summary.ham_RHS:.6e}")
        print(f"    with ‖H_A‖≈{summary.H_A_op_norm:.6e}, log d_A≈{summary.log_dA:.6e}, log d_min≈{summary.log_dmin:.6e}, ⟨\hatΨ⟩≈{summary.hatPsi_avg:.6e}")
        print(f"  Speed-based: RHS=(ħ/2)·δ·ε_F·ε_S·ε_I ≈ {summary.speed_RHS:.6e}")
        print(f"    with δ≈{summary.delta_active:.3f}, ε_F≈{summary.eps_F:.3e}, ε_S≈{summary.eps_S:.3e}, ε_I≈{summary.eps_I:.3e}")
    except Exception:
        pass

    # CHM/boost mapping: check ΔK_A S_A^2 ≥ (2π/κ) ħ (pure-state diagnostic only)
    try:
        rhs = (2.0 * math.pi / max(args.kappa, 1e-30)) * args.hbar
        dK_list = [fr.dK_ball for fr in frames]
        S2_list = [fr.S_A_nats * fr.S_A_nats for fr in frames]
        times = [fr.t for fr in frames]
        slacks = [dK * S2 - rhs for dK, S2 in zip(dK_list, S2_list)]
        pass_frac_full = float(np.mean([1.0 if s >= 0.0 else 0.0 for s in slacks])) if slacks else 0.0
        min_slack_full = float(np.min(slacks)) if slacks else 0.0
        if args.linear_gate:
            mask = [t <= args.t_lin for t in times]
            slacks_g = [s for s, m in zip(slacks, mask) if m]
            pass_frac = float(np.mean([1.0 if s >= 0.0 else 0.0 for s in slacks_g])) if slacks_g else 0.0
            min_slack = float(np.min(slacks_g)) if slacks_g else 0.0
        else:
            pass_frac = pass_frac_full
            min_slack = min_slack_full
        # Planck area floor estimate from bound ΔK_A (A/4ℓ_P^2)^2 ≥ (2π/κ) ħ
        HBAR_SI = 1.054571817e-34  # J·s (CODATA)
        ell_P2 = args.G * HBAR_SI / (args.c ** 3)
        # Natural-units area floors (in ℓ_P^2) and SI (m^2)
        A_min_lP2 = [4.0 * math.sqrt(max((2.0 * math.pi * args.hbar) / (max(args.kappa, 1e-30) * max(dK, 1e-30)), 0.0)) for dK in dK_list]
        A_min_SI = [ell_P2 * a_lP2 for a_lP2 in A_min_lP2]
        A_from_S_lP2 = [4.0 * fr.S_A_nats for fr in frames]
        A_from_S_SI = [ell_P2 * a_lP2 for a_lP2 in A_from_S_lP2]
        print("CHM/boost mapping (pure-state diagnostic):")
        print(f"  RHS=(2π/κ)ħ = {rhs:.6e} (units of ΔK)")
        if args.linear_gate:
            print(f"  Pass fraction (t<=t_lin={args.t_lin}): {pass_frac:.3f}; min slack: {min_slack:.6e}")
            print(f"  Pass fraction (full window): {pass_frac_full:.3f}; min slack: {min_slack_full:.6e}")
        else:
            print(f"  Pass fraction (ΔK_A S_A^2 ≥ RHS): {pass_frac:.3f}; min slack: {min_slack:.6e}")
        # Dimensionless natural-units floor (no κ): A_nat ≥ 4 ℓ_P^2 / √ΔK_A
        sqrt_dK = [math.sqrt(max(dK, 0.0)) for dK in dK_list]
        A_nat_lP2 = [4.0 / max(sd, 1e-30) for sd in sqrt_dK]
        A_nat_SI = [ell_P2 * a for a in A_nat_lP2]
        if A_nat_SI:
            print(f"  A_nat floor (mean): {float(np.mean(A_nat_SI)):.6e} m^2; (min): {float(np.min(A_nat_SI)):.6e} m^2  [ΔK-only]")
            print(f"  A_nat floor (mean): {float(np.mean(A_nat_lP2)):.6e} ℓ_P^2; (min): {float(np.min(A_nat_lP2)):.6e} ℓ_P^2  [ΔK-only]")
        if A_min_SI:
            print(f"  A_min floor (mean): {float(np.mean(A_min_SI)):.6e} m^2; (min): {float(np.min(A_min_SI)):.6e} m^2")
            print(f"  A_min floor (mean): {float(np.mean(A_min_lP2)):.6e} ℓ_P^2; (min): {float(np.min(A_min_lP2)):.6e} ℓ_P^2")
        if A_from_S_SI:
            print(f"  A from S_A (mean): {float(np.mean(A_from_S_SI)):.6e} m^2; {float(np.mean(A_from_S_lP2)):.6e} ℓ_P^2")
    except Exception as _:
        pass

    if args.quench:
        qres = quench_sweep_first_law(args.L, args.ell, args.J, args.h, args.dt, args.steps, args.chi, args.cut, args.hbar,
                                      args.init, args.out_prefix, QuenchCfg(enable=True, epsilons=tuple(args.eps),
                                                                            t_linear_max=args.t_lin, hx_weight_in_KA=args.ka_hx,
                                                                            normalize_KA=bool(args.ka_norm)))
        if qres:
            print(f"First-law fit: slope={qres['slope']:.4f}, R2={qres['r2']:.4f}")
        else:
            print("First-law fit: no data produced (check parameters)")
    print("Summary:"); print(json.dumps(asdict(summary), indent=2)) 