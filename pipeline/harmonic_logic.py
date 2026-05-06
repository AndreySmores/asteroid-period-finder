import numpy as np


def find_peak_clusters(frequencies, powers, n_peaks=20, cluster_width=0.05):
    """Group top N peaks into clusters, return best representative per cluster"""
    top_idx = np.argsort(powers)[-n_peaks:][::-1]
    top_freqs = frequencies[top_idx]
    top_powers = powers[top_idx]
    
    clusters = []
    used = np.zeros(len(top_freqs), dtype=bool)
    
    for i in range(len(top_freqs)):
        if used[i]:
            continue
        mask = np.abs(top_freqs - top_freqs[i]) < cluster_width
        cluster_freqs = top_freqs[mask]
        cluster_powers = top_powers[mask]
        used[mask] = True
        
        best_idx = np.argmax(cluster_powers)
        clusters.append((cluster_freqs[best_idx], cluster_powers[best_idx]))
    
    return clusters


def fourier_model(t, frot, coeffs):
    """Evaluate second order Fourier series at times t"""
    C, A1, B1, A2, B2 = coeffs
    return (C
            + A1 * np.sin(2 * np.pi * frot * t)
            + B1 * np.cos(2 * np.pi * frot * t)
            + A2 * np.sin(4 * np.pi * frot * t)
            + B2 * np.cos(4 * np.pi * frot * t))


def fit_fourier(t, y, frot):
    """Fit second order Fourier series via linear least squares at a given frot"""
    A = np.column_stack([
        np.ones(len(t)),
        np.sin(2 * np.pi * frot * t),
        np.cos(2 * np.pi * frot * t),
        np.sin(4 * np.pi * frot * t),
        np.cos(4 * np.pi * frot * t)
    ])
    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    y_pred = A @ coeffs
    chi2_reduced = np.sum((y - y_pred) ** 2) / (len(t) - len(coeffs))
    return coeffs, chi2_reduced

def fit_single(t, y, frot):
    """Fit single harmonic (first order Fourier series)"""
    A = np.column_stack([
        np.ones(len(t)),
        np.sin(2 * np.pi * frot * t),
        np.cos(2 * np.pi * frot * t)
    ])
    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    y_pred = A @ coeffs
    chi2_reduced = np.sum((y - y_pred) ** 2) / (len(t) - 3)
    return coeffs, chi2_reduced

def refine_period(t, y, frot_candidate, window=0.1, n_steps=10000):
    """
    Refine frot around a candidate by grid searching +-window
    and picking the frequency with lowest chi2.
    frot_candidate is treated directly as frot — no halving.
    """
    f_lo = frot_candidate * (1 - window)
    f_hi = frot_candidate * (1 + window)
    fine_freqs = np.linspace(f_lo, f_hi, n_steps)

    chi2s = [fit_fourier(t, y, f)[1] for f in fine_freqs]

    best_frot = fine_freqs[np.argmin(chi2s)]
    best_coeffs, best_chi2 = fit_fourier(t, y, best_frot)

    return best_frot, best_coeffs, best_chi2


def compare_models(t, y, frot):
    """
    Compare single vs double harmonic model at frot.
    Returns (chi2_single, chi2_double).
    """
    A_single = np.column_stack([
        np.ones(len(t)),
        np.sin(2 * np.pi * frot * t),
        np.cos(2 * np.pi * frot * t)
    ])
    coeffs_single, _, _, _ = np.linalg.lstsq(A_single, y, rcond=None)
    y_pred_single = A_single @ coeffs_single
    chi2_single = np.sum((y - y_pred_single) ** 2) / (len(t) - 3)

    _, chi2_double = fit_fourier(t, y, frot)

    return chi2_single, chi2_double


def harmonic_logic(t, y, frequencies, powers, n_peaks=20, cluster_width=0.05, window=0.1):
    """
    For each LS cluster frequency f_ls:
      - Interpret as single peaked: frot = f_ls, refine around f_ls
      - Interpret as double peaked: frot = f_ls/2, refine around f_ls/2
    Pick the (interpretation, cluster) pair with lowest chi2.
    is_double_peaked determined by second harmonic amplitude ratio.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    clusters = find_peak_clusters(frequencies, powers, n_peaks, cluster_width)
    print(f"Found {len(clusters)} clusters, refining each...")

    best_chi2 = np.inf
    best_frot = None
    best_coeffs = None

    for f_ls, power in clusters:
        # f_ls as frot, single harmonic
        coeffs_s1, chi2_s1 = fit_single(t, y, f_ls)
        
        # f_ls as frot, double harmonic
        coeffs_s2, chi2_s2 = fit_fourier(t, y, f_ls)
        
        # f_ls/2 as frot, single harmonic
        coeffs_d1, chi2_d1 = fit_single(t, y, f_ls / 2)
        
        # f_ls/2 as frot, double harmonic
        coeffs_d2, chi2_d2 = fit_fourier(t, y, f_ls / 2)

        candidates = [
            (f_ls,     coeffs_s1, chi2_s1, False),  # single peaked, 1st order
            (f_ls,     coeffs_s2, chi2_s2, True),   # single peaked, 2nd order
            (f_ls/2,   coeffs_d1, chi2_d1, False),  # double peaked, 1st order
            (f_ls/2,   coeffs_d2, chi2_d2, True),   # double peaked, 2nd order
        ]

        for frot, coeffs, chi2, is_double in candidates:
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_frot = frot
                best_coeffs = coeffs
                best_is_double = is_double

    A1, B1 = best_coeffs[1], best_coeffs[2]
    A2, B2 = best_coeffs[3], best_coeffs[4]
    fundamental_amplitude = np.sqrt(A1**2 + B1**2)
    harmonic_amplitude = np.sqrt(A2**2 + B2**2)
    is_double_peaked = (harmonic_amplitude / fundamental_amplitude) > 0.1

    period_hours = 24 / best_frot
    chi2_single, chi2_double = compare_models(t, y, best_frot)

    print(f"\nBest rotation frequency: {best_frot:.6f} c/d")
    print(f"Best period: {period_hours:.4f} hours")
    print(f"Chi2 single: {chi2_single:.6f}")
    print(f"Chi2 double: {chi2_double:.6f}")
    print(f"Harmonic amplitude ratio: {harmonic_amplitude/fundamental_amplitude:.4f}")
    print(f"Double peaked: {is_double_peaked}")

    return best_frot, period_hours, best_coeffs, bool(is_double_peaked)