#!/usr/bin/env python3
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional, List


@dataclass
class GKResults:
    tau_g: float
    tau_bound: Optional[float]
    l1_abs: float
    var_pi: float


def detrend_mean(x: np.ndarray) -> np.ndarray:
    return x - np.mean(x)


def normalized_acf(x: np.ndarray, dt: float, max_lag: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Return lags (s) and normalized ACF C(τ) for τ>=0 using FFT."""
    n = len(x)
    if max_lag is None:
        max_lag = n - 1
    x0 = detrend_mean(x)
    # next power of two for FFT
    nfft = 1 << (2 * n - 1).bit_length()
    fx = np.fft.rfft(x0, n=nfft)
    sxx = fx * np.conj(fx)
    acf_full = np.fft.irfft(sxx, n=nfft)[:n]
    acf_full /= acf_full[0]
    lags = dt * np.arange(0, max_lag + 1)
    return lags, acf_full[: max_lag + 1]


def gk_tau_from_acf(lags: np.ndarray, acf: np.ndarray) -> float:
    """Trapezoidal integral of normalized ACF for τ>=0."""
    return np.trapz(acf, lags)


def short_time_bound(pi: np.ndarray, dt: float) -> Optional[float]:
    """τ_bound = Var[Π] / ⟨(dΠ/dt)^2⟩ if derivative makes sense."""
    if len(pi) < 3:
        return None
    dpi_dt = np.gradient(pi, dt)
    var_pi = np.var(pi)
    mean_d2 = np.mean(dpi_dt ** 2)
    if mean_d2 <= 0:
        return None
    return var_pi / mean_d2


def one_sided_psd_zero(acf: np.ndarray, dt: float) -> float:
    """S(0) = 2 ∫_0^∞ C(τ) dτ ≈ 2 * dt * sum C for τ>=0."""
    return 2.0 * np.trapz(acf, dx=dt)


def welch_cross_spectrum(
    x: np.ndarray,
    y: np.ndarray,
    fs: float,
    nperseg: int = 256,
    noverlap: Optional[int] = None,
    nfft: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Cross-spectrum via Welch (Hann window). Returns f (Hz) and S_xy(f).
    Supports optional zero-padding via nfft to densify low-ω bins.
    """
    from scipy.signal import csd, get_window

    if noverlap is None:
        noverlap = nperseg // 2
    f, pxy = csd(
        x,
        y,
        fs=fs,
        window=get_window('hann', nperseg),
        nperseg=nperseg,
        noverlap=noverlap,
        nfft=nfft,
        return_onesided=True,
    )
    return f, pxy


def small_omega_mask(omega: np.ndarray, tau: float, max_frac: float = 0.2) -> np.ndarray:
    return np.abs(omega) * tau <= max_frac


def tau_kk_from_chi(
    omega: np.ndarray,
    chi_im: np.ndarray,
    tau_guess: Optional[float] = None,
    min_points: int = 10,
    max_frac_start: Optional[float] = None,
) -> float:
    """Adaptive robust slope proxy τ_KK.
    Default band: |ω| τ ≤ 0.2; if too few points, widen progressively.
    Returns NaN if no valid band is found.
    """
    # sanitize
    omega = np.asarray(omega)
    chi_im = np.asarray(chi_im)
    finite = np.isfinite(omega) & np.isfinite(chi_im)
    omega = omega[finite]
    chi_im = chi_im[finite]
    if omega.size == 0:
        return np.nan

    if tau_guess is None or not np.isfinite(tau_guess) or tau_guess <= 0:
        # choose a scale from omega grid (5th percentile of 1/|ω|)
        nonzero = omega[np.nonzero(omega)]
        if nonzero.size == 0:
            return np.nan
        tau_guess = float(np.percentile(1.0 / np.abs(nonzero), 5))

    # candidate bands: |ω| τ ≤ thresholds OR ω < k/τ
    thresholds = [0.2, 0.3, 0.5, 1.0]
    if isinstance(max_frac_start, (int, float)) and max_frac_start > 0:
        # Prepend the requested band and ensure monotone increase
        thresholds = sorted({float(max_frac_start)} | set(thresholds))
    factors = [1.0, 2.0, 5.0]
    for th in thresholds:
        m = (np.abs(omega) * tau_guess <= th) & (omega != 0)
        if np.count_nonzero(m) >= min_points:
            vals = -chi_im[m] / omega[m]
            return float(np.median(vals))
    for k in factors:
        m = (omega > 0) & (omega < (k / tau_guess))
        if np.count_nonzero(m) >= min_points:
            vals = -chi_im[m] / omega[m]
            return float(np.median(vals))
    # fallback: first N positive-frequency points
    pos = np.where(omega > 0)[0]
    if pos.size >= min_points:
        idx = pos[:min_points]
        vals = -chi_im[idx] / omega[idx]
        return float(np.median(vals))
    return np.nan


def window_indices(n: int, win: int, step: int) -> List[Tuple[int, int]]:
    idxs: List[Tuple[int, int]] = []
    start = 0
    while start + win <= n:
        idxs.append((start, start + win))
        start += step
    return idxs


def bootstrap_tau_g(
    pi: np.ndarray,
    dt: float,
    window_s: float = 20.0,
    step_s: Optional[float] = None,
    acf_cut_zero: bool = True,
) -> Tuple[float, float, float, int]:
    """
    Sliding-window GK estimates; returns (mean, low95, high95, n_windows).
    Uses normalized ACF and integrates to first zero crossing if acf_cut_zero.
    """
    n = len(pi)
    win = max(2, int(round(window_s / dt)))
    step = max(1, int(round((step_s if step_s is not None else (window_s / 2.0)) / dt)))
    idxs = window_indices(n, win, step)
    vals: List[float] = []
    for a, b in idxs:
        seg = pi[a:b]
        if seg.size < 3 or np.allclose(np.var(seg), 0.0):
            continue
        lags, acf = normalized_acf(seg, dt)
        if acf_cut_zero:
            nz = np.where(acf <= 0)[0]
            if nz.size > 0:
                cut_idx = int(nz[0])
                lags = lags[: cut_idx + 1]
                acf = acf[: cut_idx + 1]
        vals.append(gk_tau_from_acf(lags, acf))
    if not vals:
        return (float("nan"), float("nan"), float("nan"), 0)
    arr = np.asarray(vals, dtype=float)
    mean = float(np.mean(arr))
    low = float(np.percentile(arr, 2.5))
    high = float(np.percentile(arr, 97.5))
    return (mean, low, high, len(vals))


