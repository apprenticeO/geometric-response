
#!/usr/bin/env python3
"""
triad_density_test.py

Purpose
-------
Estimate the empirical lower density of the intersection
K = { t : sqrt(F_Q(ρ_A(t); H_A)) >= eps_F  and  S_A(t) >= eps_S  and  I(A:Ā,t) >= eps_I }
for a small locally-coupled spin-1/2 chain evolved under a local Hamiltonian.

This script implements:
- Exact unitary dynamics for a chain (TFIM + XX coupling)
- Choice of subsystem A (contiguous block) with genuine cross-cut coupling
- Reduced state ρ_A(t), von Neumann entropy S_A(t)
- Mutual information I(A:Ā,t) = S_A + S_Ā - S_AB (for pure global |ψ(t)⟩, I = 2 S_A)
- SLD Quantum Fisher Information F_Q(ρ_A; H_A) via spectral formula
- Empirical "intersection density" over sampled times
- Optional sweep over a small non-commuting perturbation ε to test continuity from commuting → noncommuting

Requirements
------------
pip install numpy scipy quimb matplotlib

References
----------
- SLD-QFI spectral formula: Braunstein & Caves, Phys. Rev. Lett. 72, 3439 (1994)
- Mutual information & entropy basics: Nielsen & Chuang, "Quantum Computation and Quantum Information" (2010)
- Locality/Lieb-Robinson background (not computed here): Nachtergaele & Sims, arXiv:1102.0835 (2011)

Usage
-----
python triad_density_test.py

You can tweak PARAMETERS section below for thresholds, system size, times, and epsilon sweep.
"""

import numpy as np
import quimb as qu
import quimb.tensor as qtn
from numpy.linalg import eigh
from scipy.linalg import expm
import os

try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None

# -----------------------
# PARAMETERS (edit these)
# -----------------------
L = 8                  # chain length (2^L Hilbert dimension)
ELL = 3                # size of contiguous block A (1 <= ELL < L)
Jx = 1.0               # XX coupling strength
hz = 0.7               # transverse field strength (non-commuting with XX)
jx_commuting = 0.0     # optional commuting term scale (set nonzero to bias towards commuting dynamics)
cut_start = 2          # start site of block A (0-indexed), A = [cut_start, ..., cut_start+ELL-1]
dt = 0.2               # time step
steps = 150            # number of time steps
eps_F = 0.08           # threshold for sqrt(F_Q)
eps_S = 0.05           # threshold for von Neumann entropy (nats)
eps_I = 0.10           # threshold for mutual information (nats)
seed = 1               # random seed for a non-eigen initial state
do_eps_sweep = True    # sweep ε to study commuting→noncommuting continuity
eps_sweep_vals = [0.0, 0.01, 0.05, 0.1, 0.2]  # strength for a small cross-cut perturbation
use_percentiles = False  # adapt thresholds by quantiles of each series per run
# quantile_q = 0.3  # DISABLED: using fixed thresholds        # e.g., 70th percentile for each leg
use_commuting_baseline = False  # if True, override to commuting baseline (Jx=0, hz=0, jx_commuting>0, eps_nc forced 0)
commuting_jx = 0.0
commuting_hz = 0.0
commuting_jzz = 1.0
use_open_system = True    # NEW: enable local dephasing to avoid I=2S
gamma = 0.05              # dephasing rate per unit time
# Controls and reporting
decoupled_cut = False     # NEW: zero the cross-cut terms (no A|Ā coupling) for control
# quantiles = {"F": 0.7, "S": 0.7, "I": 0.7}  # DISABLED: using fixed thresholds  # per-leg quantiles
num_windows = 5           # split trajectory into windows for lower-density proxy
save_npz = True           # save raw arrays and indicators for reproducibility
assert_controls = True    # NEW: assert S_A=I=0 in decoupled-closed control
out_dir = os.path.join(os.path.dirname(__file__), 'results')
os.makedirs(out_dir, exist_ok=True)

# -----------------------
# Helper: spin operators
# -----------------------
sx = qu.pauli('X')
sy = qu.pauli('Y')
sz = qu.pauli('Z')
id2 = np.eye(2)

def kron_n(ops):
    """Kronecker product for a list of single-site ops."""
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out

def op_on_sites(op, sites, L):
    """Place single-site operator `op` on each site in `sites` of an L-site chain (sum)."""
    H = np.zeros((2**L, 2**L), dtype=complex)
    for s in sites:
        ops = [id2]*L
        ops[s] = op
        H += kron_n(ops)
    return H

def two_site_bond(opA, i, opB, j, L):
    """Two-site operator opA on i tensor opB on j (tensor identity elsewhere)."""
    ops = [id2]*L
    ops[i] = opA
    ops[j] = opB
    return kron_n(ops)

def xx_coupling(L, Jx):
    """Nearest-neighbor XX coupling H = -Jx * sum (X_i X_{i+1} + Y_i Y_{i+1})."""
    H = np.zeros((2**L, 2**L), dtype=complex)
    for i in range(L-1):
        H += -Jx * (two_site_bond(sx, i, sx, i+1, L) + two_site_bond(sy, i, sy, i+1, L))
    return H

def tfim_field(L, hz):
    """Transverse field term H = -hz * sum Z_i (commutes with bond ZZ but not with XX)."""
    return -hz * op_on_sites(sz, range(L), L)

def commuting_term(L, jx):
    """Add a commuting (diagonal in Z) nearest-neighbor ZZ term Hc = -jx * sum Z_i Z_{i+1}."""
    Hc = np.zeros((2**L, 2**L), dtype=complex)
    for i in range(L-1):
        Hc += -jx * two_site_bond(sz, i, sz, i+1, L)
    return Hc

def build_total_H(L, Jx, hz, jx_commuting, eps_nc=0.0):
    """
    Total Hamiltonian: H = H_xx + H_field + H_ZZ_commuting + eps_nc * V_nc,
    where V_nc adds a small cross-cut noncommuting perturbation across A|Ā boundary.
    """
    H = xx_coupling(L, Jx) + tfim_field(L, hz) + commuting_term(L, jx_commuting)
    if eps_nc > 0.0:
        i_b = CUT_END
        if 0 <= i_b < L-1:
            V = two_site_bond(sx, i_b, sz, i_b+1, L)
            H += eps_nc * V
    # Optionally decouple cross-cut bond by zeroing the XX, YY, and ZZ terms at the boundary
    if decoupled_cut:
        i_b = CUT_END
        if 0 <= i_b < L-1:
            H += Jx * (two_site_bond(sx, i_b, sx, i_b+1, L) + two_site_bond(sy, i_b, sy, i_b+1, L))
            H += jx_commuting * two_site_bond(sz, i_b, sz, i_b+1, L)
    return H

# -----------------------
# Reduced states & entropies
# -----------------------

def reduced_state(psi, A_sites, L):
    """Reduced density matrix ρ_A for pure state vector |ψ⟩ on sites A."""
    # psi is length 2^L vector
    A = list(A_sites)
    B = [i for i in range(L) if i not in A]
    dims = [2]*L
    # reshape into tensor with L physical indices
    psi_t = psi.reshape([2]*L)
    # move A indices to front
    perm = A + B
    psi_perm = np.transpose(psi_t, axes=perm)
    dA = 2**len(A); dB = 2**len(B)
    psi_mat = psi_perm.reshape(dA, dB)
    rhoA = psi_mat @ psi_mat.conj().T
    return rhoA

def von_neumann_entropy(rho):
    """Von Neumann entropy S = -Tr ρ ln ρ (nats)."""
    evals = np.clip(np.real(np.linalg.eigvalsh(rho)), 0.0, 1.0)
    nz = evals[evals > 0]
    return float(-np.sum(nz * np.log(nz)))

def mutual_information_pure(SA):
    """For a pure global state, I(A:Ā) = 2 S_A."""
    return 2.0 * SA

def mutual_information_general(rhoA, rhoB, rhoAB):
    """General MI = S(ρA) + S(ρB) - S(ρAB)."""
    return von_neumann_entropy(rhoA) + von_neumann_entropy(rhoB) - von_neumann_entropy(rhoAB)

def reduced_dm(rho_dm, keep_sites, L):
    """Partial trace over complement of keep_sites; returns density matrix on keep_sites."""
    keep = sorted(keep_sites)
    throw = sorted(set(range(L)) - set(keep))
    rho_t = rho_dm.reshape([2]*L + [2]*L)
    perm = keep + throw + [L + i for i in keep] + [L + i for i in throw]
    rho_perm = np.transpose(rho_t, axes=perm)
    dK = 2**len(keep); dT = 2**len(throw)
    rho_perm = rho_perm.reshape(dK, dT, dK, dT)
    rho_keep = np.einsum('ikjk->ij', rho_perm)
    return rho_keep

def vN_entropy_dm(rho_dm):
    evals = np.clip(np.real(np.linalg.eigvalsh((rho_dm + rho_dm.conj().T)/2)), 0.0, 1.0)
    nz = evals[evals > 0]
    return float(-np.sum(nz * np.log(nz)))

def apply_local_dephasing(rho_dm, sites, L, gamma, dt):
    p = 1.0 - np.exp(-gamma * dt)
    E0 = np.sqrt(max(0.0, 1.0 - p)) * np.eye(2)
    E1 = np.sqrt(max(0.0, p)) * np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    for s in sites:
        if s < 0 or s >= L:
            continue
        ops0 = [np.eye(2)] * L; ops0[s] = E0; K0 = kron_n(ops0)
        ops1 = [np.eye(2)] * L; ops1[s] = E1; K1 = kron_n(ops1)
        rho_dm = K0 @ rho_dm @ K0.conj().T + K1 @ rho_dm @ K1.conj().T
    return rho_dm

# -----------------------
# SLD QFI for mixed state ρ and Hermitian H (on A)
# -----------------------

def sld_qfi(rho, H):
    """
    Compute SLD Quantum Fisher Information F_Q(ρ; H)
    using the spectral decomposition formula (Braunstein & Caves, 1994):
    F_Q = 2 * sum_{i,j} ((λ_i - λ_j)^2 / (λ_i + λ_j)) * |⟨i|H|j⟩|^2,
    with terms where λ_i + λ_j > 0. Returns a real float.
    """
    lam, U = eigh((rho + rho.conj().T)/2)  # ensure Hermiticity & numerical stability
    # Clip small negatives to 0
    lam = np.clip(np.real(lam), 0.0, 1.0)
    # Transform H into eigenbasis of ρ
    H_e = U.conj().T @ H @ U
    # Compute contribution
    F = 0.0
    for i in range(len(lam)):
        for j in range(len(lam)):
            denom = lam[i] + lam[j]
            if denom > 1e-14:
                num = (lam[i] - lam[j])**2
                Hij = H_e[i, j]
                F += 2.0 * (num / denom) * (abs(Hij)**2)
    return float(np.real(F))

def local_H_on_A(L, H_total, A_sites):
    """
    Construct a local generator H_A by restricting to terms fully supported in A.
    Here we build H_A by projecting H_total onto A by tracing out Ā in the operator sense.
    For simplicity, we *approximate* H_A by the sum of 1- and 2-site terms that lie entirely in A
    for the XX+field+ZZ model.
    """
    H_A = np.zeros((2**len(A_sites), 2**len(A_sites)), dtype=complex)
    # one-site field terms in A
    for a in A_sites:
        ops = [np.eye(2)] * len(A_sites)
        ai = A_sites.index(a)
        ops[ai] = -hz * sz
        H_A += kron_n(ops)
    # two-site XX and ZZ terms fully inside A
    for i in range(len(A_sites)-1):
        ai = A_sites[i]
        aj = A_sites[i+1]
        # map to local indices
        ii = i
        jj = i+1
        H_A += -Jx * (two_local_bond(sx, ii, sx, jj, len(A_sites)) + two_local_bond(sy, ii, sy, jj, len(A_sites)))
        H_A += -jx_commuting * two_local_bond(sz, ii, sz, jj, len(A_sites))
    return H_A

def two_local_bond(opA, i, opB, j, LA):
    ops = [np.eye(2)] * LA
    ops[i] = opA
    ops[j] = opB
    return kron_n(ops)

# -----------------------
# Build system & run
# -----------------------
CUT_START = cut_start
CUT_END = cut_start + ELL - 1
A_sites = list(range(CUT_START, CUT_START + ELL))
B_sites = [i for i in range(L) if i not in A_sites]

def run_once(eps_nc=0.0, verbose=True, make_plot=True):
    # Hamiltonian build according to modes
    if use_commuting_baseline:
        H = build_total_H(L, 0.0, 0.0, commuting_jzz, eps_nc=0.0)
    else:
        H = build_total_H(L, Jx, hz, jx_commuting, eps_nc=eps_nc)

    # Boundary decoupling diagnostic (Hilbert–Schmidt coefficient of boundary term)
    i_b = CUT_END
    if 0 <= i_b < L-1:
        B = (-Jx) * (two_site_bond(sx, i_b, sx, i_b+1, L) + two_site_bond(sy, i_b, sy, i_b+1, L)) + (-jx_commuting) * two_site_bond(sz, i_b, sz, i_b+1, L)
        num = np.trace(H.conj().T @ B)
        den = np.trace(B.conj().T @ B)
        hs_coeff = 0.0 if abs(den) < 1e-15 else float(np.real(num/den))
        # Frobenius norm of boundary component present in H
        frob = float(np.sqrt(np.real(np.trace((hs_coeff*B).conj().T @ (hs_coeff*B))))) if den != 0 else 0.0
    else:
        hs_coeff = 0.0; frob = 0.0

    # Initial state vector for baseline; build density matrix
    local_use_open = use_open_system
    if use_commuting_baseline or decoupled_cut:
        # Enforce product Z-basis and closed evolution for controls
        psi0 = np.array([1.0, 0.0], dtype=complex)
        for _ in range(L-1):
            psi0 = np.kron(psi0, np.array([1.0, 0.0], dtype=complex))
        local_use_open = False
    else:
        rng = np.random.default_rng(seed)
        def rand_qubit():
            u, v = rng.random(2)
            theta = np.arccos(1 - 2*u)
            phi = 2*np.pi*v
            return np.array([np.cos(theta/2.0), np.exp(1j*phi)*np.sin(theta/2.0)], dtype=complex)
        psi0 = rand_qubit()
        for _ in range(L-1):
            psi0 = np.kron(psi0, rand_qubit())
    rho = np.outer(psi0, psi0.conj())

    # Precompute evolution operator per step (exact)
    Udt = expm(-1j * H * dt)

    # Local generator on A (approximate from model terms)
    H_A = local_H_on_A(L, H, A_sites)
    opnorm = float(np.linalg.norm(H_A, 2))
    denom = 2.0 * max(opnorm, 1e-12)

    # Sites to dephase: boundary across cut
    dephase_sites = [CUT_END, CUT_END + 1]

    times = []
    sqrtF = []
    SA = []
    IA = []
    dH_series = []  # ΔH_A on reduced ρ_A for energy triad
    psi_hat_series = []  # intensive triad Ψ̂ per Eq. (12)

    for n in range(steps):
        t = n * dt
        rho = Udt @ rho @ Udt.conj().T
        if local_use_open:
            rho = apply_local_dephasing(rho, dephase_sites, L, gamma, dt)

        rhoA = reduced_dm(rho, A_sites, L)
        rhoB = reduced_dm(rho, B_sites, L)
        sA = vN_entropy_dm(rhoA)
        sB = vN_entropy_dm(rhoB)
        sAB = vN_entropy_dm(rho)
        IAB = max(0.0, sA + sB - sAB)

        Fq = sld_qfi(rhoA, H_A)
        sqrtF_raw = float(np.sqrt(max(Fq, 0.0)))
        sqrtF_val = sqrtF_raw / denom  # normalized by 2||H_A||
        sqrtF_val = float(min(1.0, max(0.0, sqrtF_val)))

        # ΔH_A on reduced ρ_A (energy leg approximation)
        EA = H_A
        mu = float(np.real(np.trace(rhoA @ EA)))
        mu2 = float(np.real(np.trace(rhoA @ (EA @ EA))))
        dH = float(np.sqrt(max(0.0, mu2 - mu*mu)))

        # Intensive Ψ̂ per Eq. (12): (ħ√F_Q)/(2||H_A||) * (S_A/log d_A) * (I/(2 log d_min))
        hbar = 1.0
        dA = 2**len(A_sites)
        dbar = 2**(L - len(A_sites))
        dmin = min(dA, dbar)
        psi_hat = (hbar * (sqrtF_raw / (2.0 * max(opnorm, 1e-12)))) * (sA / max(np.log(dA), 1e-12)) * (IAB / (2.0 * max(np.log(dmin), 1e-12)))

        times.append(t); sqrtF.append(sqrtF_val); SA.append(sA); IA.append(IAB)
        dH_series.append(dH); psi_hat_series.append(float(psi_hat))

    times = np.array(times); sqrtF = np.array(sqrtF); SA = np.array(SA); IA = np.array(IA)
    dH_series = np.array(dH_series); psi_hat_series = np.array(psi_hat_series)

    # Percentile thresholds (per leg) + positive floor
    thr_F = eps_F; thr_S = eps_S; thr_I = eps_I
    if use_percentiles:
        thr_F = float(np.quantile(sqrtF, quantiles.get('F', quantile_q)))
        thr_S = float(np.quantile(SA,     quantiles.get('S', quantile_q)))
        thr_I = float(np.quantile(IA,     quantiles.get('I', quantile_q)))
    thr_F = max(thr_F, 1e-6); thr_S = max(thr_S, 1e-6); thr_I = max(thr_I, 1e-6)

    K_F = sqrtF >= thr_F; K_S = SA >= thr_S; K_I = IA >= thr_I; K_all = K_F & K_S & K_I

    def wilson_ci(k, n, z=1.96):
        if n <= 0: return (0.0, 0.0)
        p = k / n; denomW = 1 + z*z/n
        center = (p + z*z/(2*n)) / denomW
        half = z * np.sqrt((p*(1-p) + z*z/(4*n))/n) / denomW
        return float(max(0.0, center - half)), float(min(1.0, center + half))

    n = steps
    dF = K_F.mean(); dS = K_S.mean(); dI = K_I.mean(); dAll = K_all.mean()
    ciF = wilson_ci(int(K_F.sum()), n); ciS = wilson_ci(int(K_S.sum()), n); ciI = wilson_ci(int(K_I.sum()), n); ciAll = wilson_ci(int(K_all.sum()), n)

    # Multi-window lower density proxy (min of block-wise fractions)
    def window_min_density(indicator):
        m = max(1, num_windows)
        size = n // m
        if size == 0:
            return float(indicator.mean())
        vals = []
        for i in range(m):
            a = i*size; b = (i+1)*size if i < m-1 else n
            if b > a:
                vals.append(float(np.mean(indicator[a:b])))
        return float(min(vals)) if vals else float(indicator.mean())

    dAll_window_min = window_min_density(K_all)

    # Lagged overlap (+/- 1 step tolerance) and pairwise correlations
    ker = np.array([1, 1, 1], dtype=int)
    def widen(x):
        y = np.convolve(x.astype(int), ker, mode='same')
        return y > 0
    K_all_lag = widen(K_F) & widen(K_S) & widen(K_I)
    dAll_lag = float(np.mean(K_all_lag))

    def corr01(a, b):
        a = a.astype(float); b = b.astype(float)
        if a.std() == 0.0 or b.std() == 0.0: return 0.0
        return float(np.corrcoef(a, b)[0, 1])
    c_FS = corr01(K_F, K_S); c_SI = corr01(K_S, K_I); c_IF = corr01(K_I, K_F)

    # Canonical structural thresholds (Elephant Eqs. 308 & 343)
    # Energy triad (LHS): time-average of ΔH_A * S_A * I(A:Ā)
    energy_triad_lhs = float(np.mean(dH_series * SA * IA)) if dH_series.size else 0.0
    # Intensive triad average <Ψ̂>
    psi_hat_mean = float(np.mean(psi_hat_series)) if psi_hat_series.size else 0.0
    # For single cut: H_min = ||H_A||_op, d_min = min(d_A, d_Ā)
    dA = 2**len(A_sites); dbar = 2**(L - len(A_sites)); dmin = min(dA, dbar)
    rhs_capped = 2.0 * (np.log(dmin)**2) * opnorm * psi_hat_mean
    # Speed-based: convert normalized eps_F to raw via sqrtF_raw scaling factor (2||H_A||/ħ)
    hbar = 1.0
    epsF_raw = (thr_F * 2.0 * max(opnorm, 1e-12)) / hbar
    rhs_speed = 0.5 * hbar * dAll * epsF_raw * thr_S * thr_I

    if verbose:
        mode = []
        mode.append('open' if local_use_open else 'closed')
        if use_commuting_baseline:
            mode.append('commuting')
        elif decoupled_cut:
            mode.append('decoupled')
        else:
            mode.append('noncommuting')
        mode_str = ','.join(mode)
        base = f"[mode={mode_str}, seed={seed}, eps_nc={eps_nc}]"
        print(f"{base} thresholds: F>={thr_F:.4f}, S>={thr_S:.4f}, I>={thr_I:.4f}")
        print(f"{base} densities:  D(F)={dF:.3f} [{ciF[0]:.3f},{ciF[1]:.3f}]  D(S)={dS:.3f} [{ciS[0]:.3f},{ciS[1]:.3f}]  D(I)={dI:.3f} [{ciI[0]:.3f},{ciI[1]:.3f}]  D(all)={dAll:.3f} [{ciAll[0]:.3f},{ciAll[1]:.3f}]  min-window D(all)={dAll_window_min:.3f}")
        print(f"  D(all) with +/-1-step tolerance = {dAll_lag:.3f}")
        print(f"  corr(KF,KS)={c_FS:.3f}  corr(KS,KI)={c_SI:.3f}  corr(KI,KF)={c_IF:.3f}")
        print(f"  boundary HS coefficient (should be ~0 if decoupled): {hs_coeff:.3e}; boundary Fro norm: {frob:.3e}")
        avg_product = float(np.mean(sqrtF * SA * IA))
        print(f"  avg < sqrtF_norm * S * I >  = {avg_product:.6f}")
        lb = (thr_F * thr_S * thr_I) * dAll
        print(f"  lower-bound estimate      = {lb:.6f}  (thr_F*thr_S*thr_I * D(all))")
        # Canonical thresholds summary
        print("\n=== Canonical Structural Thresholds (Elephant) ===")
        print(f"Energy triad LHS ⟨ΔH_A·S_A·I⟩      : {energy_triad_lhs:.6f}")
        print(f"RHS_capped 2(log d_min)^2||H||⟨Ψ̂⟩ : {rhs_capped:.6f}  (d_min={dmin}, ||H_A||={opnorm:.6f}, ⟨Ψ̂⟩={psi_hat_mean:.6f})")
        print(f"RHS_speed (ħ/2)·δ·ε_F·ε_S·ε_I     : {rhs_speed:.6f}  (δ={dAll:.3f}, ε_F_raw~{epsF_raw:.6f}, ε_S={thr_S:.6f}, ε_I={thr_I:.6f})")

    # Save plots and NPZ
    tag_parts = [f"L{L}", f"A{ELL}", f"seed{seed}", ('open' if local_use_open else 'closed')]
    if use_commuting_baseline: tag_parts.append('comm')
    elif decoupled_cut: tag_parts.append('decpl')
    else: tag_parts.append(f"eps{eps_nc}")
    tag = '_'.join(tag_parts)

    # Assert decoupled/closed control yields S=I=0
    if assert_controls and decoupled_cut and (not local_use_open):
        if (np.max(SA) > 1e-10) or (np.max(IA) > 1e-10):
            raise AssertionError(f"Decoupled-closed control violated: max S_A={np.max(SA):.3e}, max I={np.max(IA):.3e}")

    png_path = None
    if make_plot and plt is not None:
        fig, ax = plt.subplots(5, 1, figsize=(8, 10), sharex=True)
        ax[0].plot(times, sqrtF, 'b-', label='sqrt(F_Q)/(2||H_A||)'); ax[0].axhline(thr_F, color='b', ls='--', alpha=0.5); ax[0].set_ylabel('sqrtF_norm')
        ax[1].plot(times, SA, 'g-', label='S_A (nats)'); ax[1].axhline(thr_S, color='g', ls='--', alpha=0.5); ax[1].set_ylabel('S_A')
        ax[2].plot(times, IA, 'm-', label='I(A:Ā) (nats)'); ax[2].axhline(thr_I, color='m', ls='--', alpha=0.5); ax[2].set_ylabel('I')
        ax[3].plot(times, K_F.astype(float), 'c-', label='1_{K_F}')
        ax[3].plot(times, K_S.astype(float), 'y-', label='1_{K_S}')
        ax[3].plot(times, K_I.astype(float), 'r-', label='1_{K_I}')
        ax[3].plot(times, K_all.astype(float), 'k-', label='1_{K_all}')
        ax[3].set_ylabel('Indicators')
        ax[4].plot(times, K_all.astype(float), 'k-', label='K_all'); ax[4].set_xlabel('t'); ax[4].set_ylabel('K_all')
        for a in ax:
            a.grid(True, alpha=0.3); a.legend(loc='best')
        fig.suptitle(f"Triad intersection: {mode_str}")
        plt.tight_layout()
        png_path = os.path.join(out_dir, f"triad_density_{tag}.png")
        plt.savefig(png_path, dpi=150); plt.close(fig)

    if save_npz:
        npz_path = os.path.join(out_dir, f"triad_density_{tag}.npz")
        np.savez(npz_path,
                 times=times, sqrtF=sqrtF, S=SA, I=IA,
                 K_F=K_F.astype(np.uint8), K_S=K_S.astype(np.uint8), K_I=K_I.astype(np.uint8), K_all=K_all.astype(np.uint8),
                 K_all_lag=K_all_lag.astype(np.uint8),
                 thr_F=thr_F, thr_S=thr_S, thr_I=thr_I,
                 energy_triad_lhs=energy_triad_lhs, rhs_capped=rhs_capped, rhs_speed=rhs_speed,
                 mode=mode_str, seed=seed,
                 hs_coeff=hs_coeff, frob=frob)

    return {
        "times": times, "sqrtF": sqrtF, "S": SA, "I": IA,
        "dens_F": dF, "dens_S": dS, "dens_I": dI, "dens_all": dAll,
        "ci_F": ciF, "ci_S": ciS, "ci_I": ciI, "ci_all": ciAll,
        "thresholds": (thr_F, thr_S, thr_I),
        "energy_triad_lhs": energy_triad_lhs,
        "rhs_capped": rhs_capped,
        "rhs_speed": rhs_speed,
        "plot": png_path
    }

def main():
    print("=== Triad Intersection Density Test (closed pure state) ===")
    print(f"L={L}, A=[{CUT_START}..{CUT_END}], |A|={ELL}, Jx={Jx}, hz={hz}, jx_commuting={jx_commuting}")
    print(f"dt={dt}, steps={steps}, thresholds: eps_F={eps_F}, eps_S={eps_S}, eps_I={eps_I}")
    print("Note: For pure global states, I(A:Ā)=2*S_A (Nielsen & Chuang, 2010).")

    if do_eps_sweep:
        sweep = []
        for e in eps_sweep_vals:
            res = run_once(eps_nc=e, verbose=True, make_plot=True)
            sweep.append((e, res["dens_all"]))
        # Sweep plot
        if plt is not None:
            es = np.array([x for x, _ in sweep])
            ds = np.array([y for _, y in sweep])
            plt.figure(figsize=(5,4))
            plt.plot(es, ds, 'o-', label='D(all) vs ε')
            plt.xlabel('ε (cross-cut perturbation)'); plt.ylabel('Empirical lower density D(all)')
            plt.grid(True, alpha=0.3); plt.legend(); plt.tight_layout()
            sw_png = os.path.join(out_dir, f"triad_density_sweep_L{L}_A{ELL}.png")
            plt.savefig(sw_png, dpi=150); plt.close()
            print(f"Saved sweep plot: {sw_png}")
    else:
        res = run_once(eps_nc=0.0, verbose=True, make_plot=True)
        if res.get("plot"):
            print(f"Saved time-series plot: {res['plot']}")

if __name__ == "__main__":
    main()
