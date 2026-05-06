import numpy as np


def find_peak_clusters(frequencies, powers, n_peaks=20, cluster_width=0.05):
    """Group top N peaks into clusters, other methods methods 
    will refine the best frequnecy with any one bin"""
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


def amplitude_ratio(coeffs):
    """
    Compute the ratio of 2nd harmonic amplitude to fundamental amplitude
    from a set of 2nd-order Fourier coefficients [C, A1, B1, A2, B2].
    Both amplitudes are computed in quadrature: sqrt(A_k^2 + B_k^2).
    Returns inf if the fundamental amplitude is zero. 

    There is no science decsisons based on this ratio ~usually~ since the condition 
    defined in harmonic_logic is nearly always true. 
    """
    A1, B1 = coeffs[1], coeffs[2]
    A2, B2 = coeffs[3], coeffs[4]
    fund = np.sqrt(A1**2 + B1**2)
    harm = np.sqrt(A2**2 + B2**2)
    return harm / fund if fund > 0 else np.inf


def refine_period(t, y, frot_candidate, window=0.1, n_steps=10000):
    """
    Refine frot around a candidate by grid searching +-window
    and picking the frequency with lowest chi2.
    frot_candidate is treated directly as frot — no halving.

    Number of steps is set at 10000 for high precion.
    TODO: Add a user flag that is passed through to this method.
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
    A_single = np.column_stack([np.ones(len(t)),np.sin(2 * np.pi * frot * t),np.cos(2 * np.pi * frot * t)])
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
    
    Note: t is NOT shifted internally — caller is responsible for consistency, works as is, but
    if things are shifted in one method, they MUST be shifted across all scripts.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    clusters = find_peak_clusters(frequencies, powers, n_peaks, cluster_width)
    print(f"Found {len(clusters)} clusters, refining each...")

    best_chi2 = np.inf
    best_frot = None
    best_coeffs = None

    for f_ls, power in clusters:
        # Interpretation 1: f_ls is the rotation frequency (single peaked)
        frot_s, coeffs_s, chi2_s = refine_period(t, y, f_ls, window=window)
        ratio_s = amplitude_ratio(coeffs_s)

        # Interpretation 2: f_ls is the double frequency (double peaked)
        frot_d, coeffs_d, chi2_d = refine_period(t, y, f_ls / 2, window=window)
        ratio_d = amplitude_ratio(coeffs_d)

        print(f"  f_ls={f_ls:.4f} c/d | "
              f"single: frot={frot_s:.4f} (P={24/frot_s:.3f} hr) chi2={chi2_s:.6f} R2={ratio_s:.4f} | "
              f"double: frot={frot_d:.4f} (P={24/frot_d:.3f} hr) chi2={chi2_d:.6f} R2={ratio_d:.4f}")

        for frot, coeffs, chi2 in [(frot_s, coeffs_s, chi2_s),
                                    (frot_d, coeffs_d, chi2_d)]:
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_frot = frot
                best_coeffs = coeffs

    best_ratio = amplitude_ratio(best_coeffs)
    is_double_peaked = best_ratio > 0.1

    period_hours = 24 / best_frot
    chi2_single, chi2_double = compare_models(t, y, best_frot)

    print(f"\nBest rotation frequency: {best_frot:.6f} c/d")
    print(f"Best period: {period_hours:.4f} hours")
    print(f"Chi2 single: {chi2_single:.6f}")
    print(f"Chi2 double: {chi2_double:.6f}")
    print(f"Harmonic amplitude ratio: {best_ratio:.4f}")
    print(f"Double peaked: {is_double_peaked}")

    return best_frot, period_hours, best_coeffs, bool(is_double_peaked)