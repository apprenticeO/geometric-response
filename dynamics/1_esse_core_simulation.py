# esse_refactor.py
# Full refactor implementing PDF-consistent Π(t), admissibility, and boundary-aware H_A
from __future__ import annotations
import numpy as np
import time
from typing import List, Tuple, Dict
import os

from qutip import (
    Qobj, basis, tensor, qeye,
    sigmax, sigmay, sigmaz,
    entropy_vn, sesolve, expect, ket2dm
)
import math
import warnings
import matplotlib.pyplot as plt

class ESSESimulation:
    """
    Phase: validation/simulation
    Inputs: N, hbar, hx, hz, J, J2, partition, commuting flag, t_max, n_points, n_trials
    Outputs: trial dictionaries containing Π(t), ⟨Π⟩_T, sliding liminf proxy, S_A(t), ΔE_A(t), commutator norms

    Paper alignment (Convention A):
    - Π_triad_rate(t) = Σ_A sqrt(F_Q(ρ_A; H_A_small)) · S_A · I(A:Ā) (rate units).
      F_Q is SLD-QFI (dimensionless), so sqrt(F_Q) has units 1/time for time parameterizations.
    - Energy product for capped bound: per sample we can compute Σ_A ΔE_A · S_A · I(A:Ā) (energy units).
    - Generators: H_A_small (A-restricted on-site/internal terms, conservative for QFI) vs H_A_full (boundary-aware via 1/2-splitting) for energy variance and commutator norms.
    - Intensive triad normalization uses global 2·log d_min in the I-leg denominator and H_min = min_A ||H_A||_op to match the theorem's (log d_min)^2 H_min ⟨Ψ̂⟩ form.
    """
    def __init__(self, N=6, hbar=1.0, hx=1.0, hz=0.5, J=1.0, J2=0.3, seed: int | None = 42):
        self.N = N
        self.hbar = hbar
        self.hx = hx
        self.hz = hz
        self.J = J
        self.J2 = J2
        self.rng = np.random.default_rng(seed)
        self.seed = seed

        # Single-qubit ops
        self.sx, self.sy, self.sz = sigmax(), sigmay(), sigmaz()
        self.id2 = qeye(2)

    # ---------- Utility: tensor embedding ----------
    def local_op(self, op: Qobj, i: int) -> Qobj:
        ops = [self.id2]*self.N
        ops[i] = op
        return tensor(ops)

    def two_site_op(self, op_i: Qobj, i: int, op_j: Qobj, j: int) -> Qobj:
        if i == j:
            raise ValueError("i and j must differ for two-site operator")
        ops = [self.id2]*self.N
        ops[i] = op_i
        ops[j] = op_j
        return tensor(ops)

    # ---------- Hamiltonians ----------
    def build_hamiltonian(self, commuting: bool = False) -> Tuple[Qobj, str]:
        H = 0
        # Nearest-neighbor Ising ZZ
        for i in range(self.N-1):
            H += -self.J * self.two_site_op(self.sz, i, self.sz, i+1)
        # Next-nearest neighbor ZZ (breaks integrability)
        if not commuting and self.N >= 3:
            for i in range(self.N-2):
                H += -self.J2 * self.two_site_op(self.sz, i, self.sz, i+2)
        # Fields
        if commuting:
            for i in range(self.N):
                H += self.hz * self.local_op(self.sz, i)
            label = "Commuting: -J Σ σzσz + hz Σ σz"
        else:
            for i in range(self.N):
                H += self.hx * self.local_op(self.sx, i)
            for i in range(self.N):
                H += self.hz * (i % 2 - 0.5) * self.local_op(self.sz, i)  # staggered
            label = "Non-commuting: -J Σ σzσz - J2 Σ σzσz(next) + hx Σ σx + staggered hz"
        return H, label

    def build_local_energy_ops(self, commuting: bool = False) -> List[Qobj]:
        """
        Construct full-space operators H_A^{(i)} that represent the local generator for site i
        including boundary couplings via symmetric splitting:
          H_i = on-site(i) + 1/2 Σ_{j neigh} J_{ij} σz_i σz_j + 1/2 Σ_{k next} J2_{ik} σz_i σz_k
        This allows computing ΔH_i from the full state without tracing tricks.
        """
        H_loc_full = [0 for _ in range(self.N)]
        # On-site parts
        for i in range(self.N):
            if commuting:
                H_loc_full[i] += self.hz * self.local_op(self.sz, i)
            else:
                H_loc_full[i] += self.hx * self.local_op(self.sx, i)
                H_loc_full[i] += self.hz * (i % 2 - 0.5) * self.local_op(self.sz, i)
        # Nearest neighbors (split half to each endpoint)
        for i in range(self.N-1):
            term = -self.J * self.two_site_op(self.sz, i, self.sz, i+1)
            H_loc_full[i]     += 0.5 * term
            H_loc_full[i+1]   += 0.5 * term
        # Next-nearest (if present)
        if not commuting and self.N >= 3:
            for i in range(self.N-2):
                term = -self.J2 * self.two_site_op(self.sz, i, self.sz, i+2)
                H_loc_full[i]     += 0.5 * term
                H_loc_full[i+2]   += 0.5 * term
        return H_loc_full

    # ---------- Per-A local operator (small space) for gap floor ----------
    def build_HA_small(self, A: Tuple[int,...], commuting: bool = False) -> Qobj:
        """Build H_A restricted to subsystem A (2^{|A|}×2^{|A|}), including on-site and internal couplings.
        Boundary couplings are not representable entirely inside A; we omit them (conservative floor)."""
        m = len(A)
        # map original site index -> position in A
        pos = {site: idx for idx, site in enumerate(sorted(A))}
        I2 = qeye(2)
        def loc(op: Qobj, i_pos: int) -> Qobj:
            ops = [I2]*m
            ops[i_pos] = op
            return tensor(ops)
        def two(op_i: Qobj, i_pos: int, op_j: Qobj, j_pos: int) -> Qobj:
            ops = [I2]*m
            ops[i_pos] = op_i
            ops[j_pos] = op_j
            return tensor(ops)
        H = 0
        # On-site
        for s in A:
            i = pos[s]
            if commuting:
                H += self.hz * loc(self.sz, i)
            else:
                H += self.hx * loc(self.sx, i)
                H += self.hz * (s % 2 - 0.5) * loc(self.sz, i)
        # Internal nearest-neighbor couplings within A
        for s in A:
            t = s+1
            if t in pos:
                H += -self.J * two(self.sz, pos[s], self.sz, pos[t])
        # Internal next-nearest couplings (noncomm)
        if not commuting:
            for s in A:
                t = s+2
                if t in pos:
                    H += -self.J2 * two(self.sz, pos[s], self.sz, pos[t])
        return H

    # ---------- QFI (SLD) for unitary encoding on subsystem A ----------
    def qfi_spectral_unitary(self, rho_A: Qobj, H_loc_A: Qobj) -> float:
        evals, evecs = rho_A.eigenstates()
        lam = np.real(np.array(evals, dtype=float))
        F = 0.0
        d = len(lam)
        for i in range(d):
            for j in range(d):
                denom = lam[i] + lam[j]
                if denom <= 1e-15:
                    continue
                Hij_q = (evecs[i].dag() * H_loc_A * evecs[j])
                try:
                    Hij = complex(Hij_q.full()[0, 0])
                except Exception:
                    Hij = complex(Hij_q)
                diff = lam[i] - lam[j]
                F += 2.0 * (diff * diff / denom) * (abs(Hij) ** 2)
        return float(np.real(F))

    # ---------- Trace-distance to product state on full space ----------
    def product_trace_distance(self, rho_full: Qobj, A: Tuple[int,...]) -> float:
        N = self.N
        A_sorted = tuple(sorted(A))
        bar = tuple([i for i in range(N) if i not in A_sorted])
        try:
            rho_A = rho_full.ptrace(list(A_sorted))
            rho_bar = rho_full.ptrace(list(bar))
            sigma = tensor(rho_A, rho_bar)
            # Reorder to original subsystem order [0..N-1]
            concat = list(A_sorted) + list(bar)
            perm = [concat.index(i) for i in range(N)]
            sigma_full = sigma.permute(perm)
            assert sigma_full.dims == rho_full.dims, "dims mismatch after permute"
            diff = (rho_full - sigma_full)
            d = float(diff.norm('tr'))  # trace norm
            return d
        except Exception:
            warnings.warn("Trace norm failed; falling back to Frobenius norm lower bound.")
            try:
                rho_A = rho_full.ptrace(list(A_sorted))
                rho_bar = rho_full.ptrace(list(bar))
                sigma = tensor(rho_A, rho_bar)
                concat = list(A_sorted) + list(bar)
                perm = [concat.index(i) for i in range(N)]
                sigma_full = sigma.permute(perm)
                diff = (rho_full - sigma_full)
                d_fro = float(diff.norm('fro'))
                print("[warn] Frobenius fallback used in product_trace_distance")
                return d_fro
            except Exception:
                return 0.0

    # ---------- Partitions ----------
    def default_partition(self) -> List[Tuple[int]]:
        return [(i,) for i in range(self.N)]

    # ---------- Initial states ----------
    def haar_random_pure(self) -> Qobj:
        dim = 2**self.N
        v = self.rng.normal(size=dim) + 1j*self.rng.normal(size=dim)
        v /= np.linalg.norm(v)
        return Qobj(v, dims=[[2]*self.N, [1]*self.N]).unit()

    def plus_product(self) -> Qobj:
        plus = (basis(2,0) + basis(2,1)).unit()
        return tensor([plus]*self.N)

    # ---------- Admissibility ----------
    def level_spacing_ratio(self, H: Qobj) -> float:
        e = np.sort(H.eigenenergies())
        s = np.diff(e)
        # Filter out tiny/zero spacings to avoid divide-by-zero
        eps = 1e-12
        s = s[np.abs(s) > eps]
        if s.size < 2:
            return float('nan')
        a, b = s[:-1], s[1:]
        denom = np.maximum(np.abs(a), np.abs(b))
        mask = denom > eps
        ratios = np.full_like(denom, np.nan, dtype=float)
        ratios[mask] = np.minimum(np.abs(a[mask]), np.abs(b[mask])) / denom[mask]
        return float(np.nanmean(ratios))

    # ---------- Core simulation ----------
    def run(self, commuting: bool = False, t_max: float = 20.0, n_points: int = 400, n_trials: int = 3,
            partition: List[Tuple[int]] | None = None, window: int = 50,
            epsilons: Tuple[float,...] = (1e-8, 1e-6, 1e-4), visualize: bool = False,
            save_fig: bool = True, outdir: str = "plots/core", persist_w: int = 3) -> List[Dict]:
        if partition is None:
            partition = self.default_partition()
        print(f"Running {'commuting' if commuting else 'non-commuting'} simulation, N={self.N}")
        H, label = self.build_hamiltonian(commuting)
        H_loc_full = self.build_local_energy_ops(commuting)

        # Admissibility diagnostics
        r = self.level_spacing_ratio(H)
        r_robust = None
        r_note = None
        if np.isnan(r):
            r_note = "degenerate/too-few-gaps; diagnostic only; ESSE metrics unaffected"
            print("Level-spacing r: n/a (degenerate/too-few-gaps). ESSE metrics are unaffected.")
            # Robust estimate with tiny diagonal jitter (for reporting only)
            try:
                eps = 1e-9
                jitter = 0
                for i in range(self.N):
                    jitter += float(self.rng.normal()) * self.local_op(self.sz, i)
                H_eps = H + eps * jitter
                r_eps = self.level_spacing_ratio(H_eps)
                if not np.isnan(r_eps):
                    r_robust = float(r_eps)
                    print(f"Level-spacing r_robust ≈ {r_robust:.3f} (ε-jitter diagnostic)")
            except Exception:
                pass
        else:
            print(f"Level-spacing r ≈ {r:.3f} (GOE ~0.53)")

        tlist = np.linspace(0.0, t_max, n_points)
        trials: List[Dict] = []

        for trial in range(n_trials):
            psi0 = self.haar_random_pure() if trial > 0 else self.plus_product()
            print(f"Trial {trial+1}/{n_trials}: state={'Haar' if trial>0 else '|+>^⊗N'}")

            t0 = time.time()
            result = sesolve(H, psi0, tlist)
            print(f"Time evolution: {time.time()-t0:.2f}s, {len(tlist)} steps")

            Pi_t = np.zeros_like(tlist, dtype=float)
            Pi_tilde = np.zeros_like(tlist, dtype=float)
            time_avg = np.zeros_like(tlist, dtype=float)
            liminf_proxy = np.zeros_like(tlist, dtype=float)

            S_blocks = np.zeros((len(partition), len(tlist)), dtype=float)
            dH_blocks = np.zeros((len(partition), len(tlist)), dtype=float)
            comm_norms = np.zeros((len(partition), len(tlist)), dtype=float)
            d_blocks = np.zeros((len(partition), len(tlist)), dtype=float)  # trace distance to product

            # Precompute per-partition operators
            HA_full_by_A = []
            HA_small_by_A = []
            for A in partition:
                Hsum = 0
                for i in A:
                    Hsum += H_loc_full[i]
                HA_full_by_A.append(Hsum)
                HA_small_by_A.append(self.build_HA_small(A, commuting=commuting))

            # Precompute norms and logs for intensive triad
            opnorms = []
            for op in HA_full_by_A:
                try:
                    vals = np.real(np.asarray(op.eigenenergies(), dtype=float))
                    opnorms.append(float(np.max(np.abs(vals))))
                except Exception:
                    # Fallback to induced 1-norm as a safe upper bound
                    opnorms.append(float(op.norm('one')))
            log_dA_list = [len(A)*math.log(2.0) for A in partition]
            log_dmin_list = [math.log(min(2**len(A), 2**(self.N-len(A)))) for A in partition]
            Hmin = float(np.min(opnorms)) if len(opnorms) else 1.0
            log_dmin_global = float(np.min(log_dmin_list)) if len(log_dmin_list) else math.log(2.0)

            # Series containers
            energy_prod_t = np.zeros_like(tlist, dtype=float)
            triad_hat_t = np.zeros_like(tlist, dtype=float)
            sqrtF_blocks = np.zeros((len(partition), len(tlist)), dtype=float)
            I_blocks = np.zeros((len(partition), len(tlist)), dtype=float)

            for tidx, ket in enumerate(result.states):
                rho = ket2dm(ket)
                triad_sum = 0.0
                energy_prod_sum = 0.0
                triad_hat_sum = 0.0
                # For each block A
                for a_idx, A in enumerate(partition):
                    # Reduced state and entropy
                    rho_A = rho.ptrace(list(A))
                    S_A = entropy_vn(rho_A, base=np.e)
                    S_blocks[a_idx, tidx] = S_A

                    # Fluctuations ΔH_A computed on FULL state directly
                    exp_H = expect(HA_full_by_A[a_idx], rho)
                    exp_H2 = expect(HA_full_by_A[a_idx]*HA_full_by_A[a_idx], rho)
                    var = max(0.0, float(np.real(exp_H2 - exp_H**2)))
                    dH = np.sqrt(var)
                    dH_blocks[a_idx, tidx] = dH

                    # Commutator norm ‖[H_A, ρ]‖_F (full-space proxy for admissibility condition)
                    comm = HA_full_by_A[a_idx]*rho - rho*HA_full_by_A[a_idx]
                    comm_norms[a_idx, tidx] = float(comm.norm('fro'))

                    # Trace distance to product state on full space
                    d_blocks[a_idx, tidx] = self.product_trace_distance(rho, A)

                    # Elephant triad contribution
                    H_A_small = HA_small_by_A[a_idx]
                    F = self.qfi_spectral_unitary(rho_A, H_A_small)
                    sqrtF = float(np.sqrt(max(F, 0.0)))
                    sqrtF_blocks[a_idx, tidx] = sqrtF
                    # Mutual information I(A:Ā)
                    rho_bar = rho.ptrace([i for i in range(self.N) if i not in A])
                    S_bar = entropy_vn(rho_bar, base=np.e)
                    S_AB = entropy_vn(rho, base=np.e)
                    I_AB = float(S_A + S_bar - S_AB)
                    I_blocks[a_idx, tidx] = I_AB
                    Psi_A = sqrtF * S_A * I_AB  # Remove ℏ factor - matches LaTeX Eq. 78
                    triad_sum += Psi_A
                    energy_prod_sum += dH * S_A * I_AB

                    # Intensive triad contribution for this A
                    denom_op = max(opnorms[a_idx], 1e-12)
                    denom_S = max(log_dA_list[a_idx], 1e-12)
                    denom_I = max(2.0*log_dmin_global, 1e-12)
                    triad_hat_sum += (self.hbar*sqrtF)/(2.0*denom_op) * (S_A/denom_S) * (I_AB/denom_I)

                # System-wide series
                triad_sum = max(triad_sum, 0.0)
                Pi_t[tidx] = triad_sum
                Pi_tilde[tidx] = triad_sum / max(1, len(partition))
                energy_prod_t[tidx] = energy_prod_sum
                triad_hat_t[tidx] = triad_hat_sum / max(1, len(partition))
                time_avg[tidx] = Pi_t[:tidx+1].mean()
                w0 = max(0, tidx - window + 1)
                liminf_proxy[tidx] = np.min(time_avg[w0:tidx+1])

            # Per-A metrics: δ_A and τ_A^4 (sliding liminf of time-averaged d^4)
            c_lin = 0.5
            deltaA = {e: [] for e in epsilons}
            tau4_list = []
            v0_list = []
            for a_idx, A in enumerate(partition):
                # v0(A) from boundary-split full-space H_A_full gap (used in good-time filter)
                H_A_full = HA_full_by_A[a_idx]
                evals = np.sort(np.real(np.asarray(H_A_full.eigenenergies(), dtype=float)))
                gaps = np.diff(evals)
                gaps = gaps[gaps > 1e-12]
                deltaE = float(np.min(gaps)) if gaps.size>0 else 0.0
                v0 = 0.25*deltaE
                v0_list.append(v0)

                # δ_A at each epsilon: frames where comm_norm>e AND ΔH_A >= v0(A)
                good_var = dH_blocks[a_idx, :] >= v0
                for e in epsilons:
                    good_comm = comm_norms[a_idx, :] > e
                    frac = float(np.mean(np.logical_and(good_comm, good_var)))
                    deltaA[e].append(frac)

                # τ_A^4 via sliding liminf of average of d^4
                d4 = d_blocks[a_idx, :]**4
                cums = np.cumsum(d4)
                avg = cums/(np.arange(len(d4))+1)
                w = max(1, len(d4)//8)
                lim = np.array([np.min(avg[max(0,i-w+1):i+1]) for i in range(len(d4))])
                tau4 = float(lim[-1])
                tau4_list.append(tau4)

            # Intensive triad average and RHS estimates
            triad_hat_mean = float(np.mean(triad_hat_t))
            # Floors (proxy via 10th percentile); delta via persistence window across all A
            try:
                epsF = float(np.quantile(sqrtF_blocks.flatten(), 0.10))
                epsS = float(np.quantile(S_blocks.flatten(), 0.10))
                epsI = float(np.quantile(I_blocks.flatten(), 0.10))
                goodF = (sqrtF_blocks >= epsF)
                goodS = (S_blocks >= epsS)
                goodI = (I_blocks >= epsI)
                ok_t = np.all(goodF & goodS & goodI, axis=0)  # all A good at time t
                if persist_w > 1 and ok_t.size >= persist_w:
                    ok_erosion = np.array([np.all(ok_t[t:t+persist_w]) for t in range(ok_t.size - persist_w + 1)])
                    delta_frac = float(np.mean(ok_erosion))
                else:
                    delta_frac = float(np.mean(ok_t))
            except Exception:
                epsF = epsS = epsI = 0.0
                delta_frac = 0.0
            capped_rhs_energy = 2.0*(log_dmin_global**2)*Hmin*triad_hat_mean
            speed_rhs_energy = 0.5*self.hbar*delta_frac*epsF*epsS*epsI

            # Leg means and commutator mean
            avg_sqrtF = float(np.mean(sqrtF_blocks))
            avg_S = float(np.mean(S_blocks))
            avg_I = float(np.mean(I_blocks))
            comm_mean = float(np.mean(comm_norms))

            # Trial aggregation
            trials.append(dict(
                tlist=tlist, Pi_t=Pi_t, Pi_tilde=Pi_tilde, time_avg=time_avg, liminf_proxy=liminf_proxy,
                S_blocks=S_blocks, dH_blocks=dH_blocks, comm_norms=comm_norms,
                d_blocks=d_blocks, sqrtF_blocks=sqrtF_blocks, I_blocks=I_blocks,
                energy_prod_t=energy_prod_t, energy_prod_mean=float(energy_prod_t.mean()),
                triad_hat_t=triad_hat_t, triad_hat_mean=triad_hat_mean,
                capped_rhs_energy=capped_rhs_energy, speed_rhs_energy=speed_rhs_energy,
                floors=dict(epsF=epsF, epsS=epsS, epsI=epsI, delta=delta_frac, proxy='q10', persist_w=persist_w),
                leg_means=dict(sqrtF=avg_sqrtF, S=avg_S, I=avg_I, Psi_hat=triad_hat_mean),
                deltaA=deltaA, tau4_list=tau4_list, v0_list=v0_list,
                label=label, r_stat=r, r_robust=r_robust, r_note=r_note, psi0_type=('Haar' if trial>0 else 'plus'),
                seed=self.seed,
                comm_mean=comm_mean
            ))

            if visualize:
                try:
                    os.makedirs(outdir, exist_ok=True)
                    # Improve style
                    try:
                        plt.style.use('seaborn-v0_8-colorblind')
                    except Exception:
                        pass
                    fig, axs = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
                    # Π(t) and running average (rate units) + energy product (secondary axis)
                    axs[0,0].plot(tlist, Pi_t, label='Π_rate(t)', color='#1f77b4', lw=2)
                    axs[0,0].plot(tlist, time_avg, label='⟨Π_rate⟩(t)', color='#ff7f0e', lw=2)
                    axs[0,0].plot(tlist, liminf_proxy, label='liminf proxy', color='#2ca02c', lw=2)
                    axs[0,0].set_title('ESSE activity (rate)')
                    axs[0,0].set_xlabel('time')
                    axs[0,0].set_ylabel('rate')
                    ax2 = axs[0,0].twinx()
                    ax2.plot(tlist, energy_prod_t, label='Σ ΔE·S·I (energy)', color='#d62728', lw=1.5, ls='--', alpha=0.8)
                    ax2.set_ylabel('energy')
                    # Merge legends
                    lines1, labs1 = axs[0,0].get_legend_handles_labels()
                    lines2, labs2 = ax2.get_legend_handles_labels()
                    axs[0,0].legend(lines1+lines2, labs1+labs2, loc='upper left', frameon=True)
                    axs[0,0].grid(True, alpha=0.3)
                    # τ_A^4
                    axs[0,1].bar(range(len(partition)), tau4_list, color='#ff7f0e')
                    axs[0,1].set_title('τ_A^4 (liminf of avg d^4)')
                    axs[0,1].set_xlabel('A index'); axs[0,1].set_ylabel('τ_A^4')
                    axs[0,1].grid(True, axis='y', alpha=0.2)
                    # δ_A at mid epsilon
                    e_plot = epsilons[1] if len(epsilons) > 1 else epsilons[0]
                    axs[1,0].bar(range(len(partition)), deltaA[e_plot], color='#2ca02c')
                    axs[1,0].set_title(f'δ_A at ε={e_plot}')
                    axs[1,0].set_xlabel('A index'); axs[1,0].set_ylabel('fraction')
                    axs[1,0].set_ylim(0, 1.05)
                    axs[1,0].grid(True, axis='y', alpha=0.2)
                    # Add fourth subplot for triad components
                    axs[1,1].plot(tlist, sqrtF_blocks.mean(axis=0), label='⟨√F_Q⟩', color='#1f77b4', lw=2)
                    axs[1,1].plot(tlist, S_blocks.mean(axis=0), label='⟨S_A⟩', color='#ff7f0e', lw=2)
                    axs[1,1].plot(tlist, I_blocks.mean(axis=0), label='⟨I(A:Ā)⟩', color='#2ca02c', lw=2)
                    axs[1,1].set_title('Triad Components')
                    axs[1,1].set_xlabel('time'); axs[1,1].set_ylabel('value')
                    axs[1,1].legend(frameon=True)
                    axs[1,1].grid(True, alpha=0.3)
                    tag = 'comm' if commuting else 'noncomm'
                    fig.suptitle(f"{label} | seed={self.seed}")
                    if save_fig:
                        fname = os.path.join(outdir, f"A1_{tag}_N{self.N}.png")
                        fig.savefig(fname, dpi=200, bbox_inches='tight')
                        print(f"Saved figure: {fname}")
                    plt.show()
                except Exception as _:
                    pass
        return trials


if __name__ == "__main__":
    sim = ESSESimulation(N=6)
    trials_nc = sim.run(commuting=False, t_max=12, n_points=200, n_trials=1, visualize=True)
    trials_c  = sim.run(commuting=True,  t_max=12, n_points=200, n_trials=1, visualize=True)

    # Minimal textual summary
    for name, trials in [('Non-comm', trials_nc), ('Comm', trials_c)]:
        finals = [float(tr['time_avg'][-1]) for tr in trials]
        eprod = [float(tr.get('energy_prod_mean', 0.0)) for tr in trials]
        caps = [float(tr.get('capped_rhs_energy', 0.0)) for tr in trials]
        spds = [float(tr.get('speed_rhs_energy', 0.0)) for tr in trials]
        print(f"{name}: <Π_rate> = {np.mean(finals):.6f} ± {np.std(finals):.6f} [units: rate]")
        print(f"{name}: <Σ_A ΔE·S·I> = {np.mean(eprod):.6f} ± {np.std(eprod):.6f} [units: energy]")
        print(f"{name}: RHS_capped≈{np.mean(caps):.6f}, RHS_speed≈{np.mean(spds):.6f}")
        print(f"{name}: seed={trials[0].get('seed')}")
        lm = trials[0].get('leg_means', {})
        fl = trials[0].get('floors', {})
        cm = trials[0].get('comm_mean', float('nan'))
        if lm:
            print(f"{name}: legs⟨·⟩ sqrtF={lm.get('sqrtF', float('nan')):.3f}, S={lm.get('S', float('nan')):.3f}, I={lm.get('I', float('nan')):.3f}, Ψ̂={lm.get('Psi_hat', float('nan')):.3f}")
        if fl:
            print(f"{name}: floors(proxy={fl.get('proxy')}, q10) epsF={fl.get('epsF'):.3f}, epsS={fl.get('epsS'):.3f}, epsI={fl.get('epsI'):.3f}, δ_persist(w={fl.get('persist_w')} )={fl.get('delta'):.3f}")
        print(f"{name}: mean ||[H_A,ρ]||_F over (t,A) ≈ {cm:.3e}")
