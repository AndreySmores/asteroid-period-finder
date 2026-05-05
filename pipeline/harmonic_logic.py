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
    """Fourier series at times t given coefficients and frequency"""
    C, A1, B1, A2, B2 = coeffs
    return (C + A1 * np.sin(2 * np.pi * frot * t) + B1 * np.cos(2 * np.pi * frot * t)+ A2 * np.sin(2 * np.pi * 2 * frot * t) + B2 * np.cos(2 * np.pi * 2 * frot * t))


def fit_fourier(t, y, frot):
    """Fit second order Fourier series via linear least squares"""
    A = np.column_stack([
        np.ones(len(t)),
        np.sin(2 * np.pi * frot * t),
        np.cos(2 * np.pi * frot * t),
        np.sin(2 * np.pi * 2 * frot * t),
        np.cos(2 * np.pi * 2 * frot * t)
    ])
    
    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    y_pred = A @ coeffs
    chi2_reduced = np.sum((y - y_pred) ** 2) / (len(t) - len(coeffs))
    
    return coeffs, chi2_reduced


def refine_period(t, y, f_candidate, window=0.1, n_steps=1000):
    """
    Refine period around candidate frequency by minimizing chi2.
    Always treats f_candidate as the double frequency (frot = f/2).
    The single vs double decision is made later in harmonic_logic.
    """
    f_lo = f_candidate * (1 - window)
    f_hi = f_candidate * (1 + window)
    fine_freqs = np.linspace(f_lo, f_hi, n_steps)

    chi2s = []
    for f in fine_freqs:
        _, chi2_r = fit_fourier(t, y, f / 2)
        chi2s.append(chi2_r)

    best_f = fine_freqs[np.argmin(chi2s)]
    best_frot = best_f / 2
    best_coeffs, best_chi2 = fit_fourier(t, y, best_frot)

    return best_frot, best_coeffs, best_chi2


def compare_models(t, y, frot):
    """
    Compare single vs double harmonic model using reduced chi2.
    Single: y = C + A1*sin(2pi*f*t) + B1*cos(2pi*f*t)
    Double: y = C + A1*sin + B1*cos + A2*sin(4pi*f*t) + B2*cos(4pi*f*t)
    """
    A_single = np.column_stack([np.ones(len(t)), np.sin(2 * np.pi * frot * t), np.cos(2 * np.pi * frot * t)])
    coeffs_single, _, _, _ = np.linalg.lstsq(A_single, y, rcond=None)
    y_pred_single = A_single @ coeffs_single
    chi2_single = np.sum((y - y_pred_single) ** 2) / (len(t) - 3)

    _, chi2_double = fit_fourier(t, y, frot)

    return chi2_single, chi2_double

def harmonic_logic(t, y, frequencies, powers, n_peaks=20, cluster_width=0.05, window=0.1):
    t = np.asarray(t)
    y = np.asarray(y)

    clusters = find_peak_clusters(frequencies, powers, n_peaks, cluster_width)
    print(f"Found {len(clusters)} clusters, refining each...")

    best_chi2 = np.inf
    best_frot = None
    best_coeffs = None
    best_is_double = False

    for f_cluster, power in clusters:
        # Refine around candidate frequency
        frot_refined, coeffs_refined, chi2_refined = refine_period(t, y, f_cluster, window=window)

        # Try refined frequency as fundamental
        coeffs_fund, chi2_fund = fit_fourier(t, y, frot_refined * 2)  # undo halving

        # Try refined frequency as double (keep as is)
        coeffs_double, chi2_double_fit = fit_fourier(t, y, frot_refined)

        # Pick whichever interpretation gives better chi2
        if chi2_fund < chi2_double_fit:
            best_frot_cluster = frot_refined * 2  # undo halving
            best_coeffs_cluster = coeffs_fund
            best_chi2_cluster = chi2_fund
            is_double = False
        else:
            best_frot_cluster = frot_refined
            best_coeffs_cluster = coeffs_double
            best_chi2_cluster = chi2_double_fit
            is_double = True

        print(f"  Cluster f={f_cluster:.4f} c/d → frot={best_frot_cluster:.4f} c/d "
              f"(P={24/best_frot_cluster:.3f} hr), chi2={best_chi2_cluster:.6f}, double={is_double}")

        if best_chi2_cluster < best_chi2:
            best_chi2 = best_chi2_cluster
            best_frot = best_frot_cluster
            best_coeffs = best_coeffs_cluster
            best_is_double = is_double

    # compare_models is now just for reporting
    chi2_single, chi2_double = compare_models(t, y, best_frot)
    period_hours = 24 / best_frot

    print(f"\nBest rotation frequency: {best_frot:.6f} c/d")
    print(f"Best period: {period_hours:.4f} hours")
    print(f"Chi2 single harmonic: {chi2_single:.6f}")
    print(f"Chi2 double harmonic: {chi2_double:.6f}")
    print(f"Double peaked: {best_is_double}")

    return best_frot, period_hours, best_coeffs, best_is_double