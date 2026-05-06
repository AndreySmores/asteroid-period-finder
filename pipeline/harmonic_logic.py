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
    """Evaluate Fourier series of arbitrary order at times t"""
    result = np.full(len(t), coeffs[0])  # C term
    n_harmonics = (len(coeffs) - 1) // 2
    for n in range(1, n_harmonics + 1):
        a = coeffs[2*n - 1]  # sin coefficient
        b = coeffs[2*n]      # cos coefficient
        result += a * np.sin(2 * np.pi * n * frot * t)
        result += b * np.cos(2 * np.pi * n * frot * t)
    return result


def fit_nth_order(t, y, frot, order):
    """
    Fit nth order Fourier series via linear least squares.
    Returns coeffs padded to length 2*max_order+1 for uniform amplitude_ratio,
    and reduced chi2.
    order: number of harmonics to include (1=single, 2=double, 3=triple, 4=quadruple)
    """
    cols = [np.ones(len(t))]
    for n in range(1, order + 1):
        cols.append(np.sin(2 * np.pi * n * frot * t))
        cols.append(np.cos(2 * np.pi * n * frot * t))
    
    A = np.column_stack(cols)
    coeffs, _, _, _ = np.linalg.lstsq(A, y, rcond=None)
    y_pred = A @ coeffs
    n_params = 1 + 2 * order
    chi2_reduced = np.sum((y - y_pred) ** 2) / (len(t) - n_params)

    # Pad coeffs to length 9 (up to 4th order: C + 4*2 = 9)
    # so amplitude_ratio always works on coeffs[1], coeffs[2], coeffs[3], coeffs[4]
    coeffs_padded = np.zeros(13)
    coeffs_padded[:len(coeffs)] = coeffs

    return coeffs_padded, chi2_reduced


def fit_fourier(t, y, frot):
    """Fit second order Fourier series — kept for compatibility"""
    coeffs, chi2 = fit_nth_order(t, y, frot, order=2)
    return coeffs, chi2


def fit_single(t, y, frot):
    """Fit first order Fourier series"""
    coeffs, chi2 = fit_nth_order(t, y, frot, order=1)
    return coeffs, chi2


def amplitude_ratio(coeffs):
    """Second to first harmonic amplitude ratio"""
    fundamental = np.sqrt(coeffs[1]**2 + coeffs[2]**2)
    harmonic = np.sqrt(coeffs[3]**2 + coeffs[4]**2)
    if fundamental < 1e-10:
        return 0.0
    return harmonic / fundamental


def refine_period(t, y, frot_candidate, window=0.1, n_steps=1000):
    """Refine frot around candidate by minimizing chi2 over fine grid"""
    f_lo = frot_candidate * (1 - window)
    f_hi = frot_candidate * (1 + window)
    fine_freqs = np.linspace(f_lo, f_hi, n_steps)
    chi2s = [fit_nth_order(t, y, f, order=6)[1] for f in fine_freqs]
    best_frot = fine_freqs[np.argmin(chi2s)]
    best_coeffs, best_chi2 = fit_nth_order(t, y, best_frot, order=6)
    return best_frot, best_coeffs, best_chi2


def compare_models(t, y, frot):
    """Compare single vs double harmonic at frot"""
    _, chi2_single = fit_nth_order(t, y, frot, order=1)
    _, chi2_double = fit_nth_order(t, y, frot, order=2)
    return chi2_single, chi2_double


def harmonic_logic(t, y, frequencies, powers, n_peaks=20, cluster_width=0.05, window=0.1):
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)

    clusters = find_peak_clusters(frequencies, powers, n_peaks, cluster_width)

    print("=" * 80)
    print(f"HARMONIC ANALYSIS — {len(clusters)} clusters found")
    print("=" * 80)

    best_chi2 = np.inf
    best_frot = None
    best_coeffs = None

    for cluster_idx, (f_ls, power) in enumerate(clusters):
        print(f"\nCluster {cluster_idx+1}: f_ls = {f_ls:.6f} c/d  "
              f"(P_ls = {24/f_ls:.4f} hr)  LS power = {power:.2f}")
        print(f"  {'Interpretation':<30} {'frot (c/d)':<14} {'Period (hr)':<14} "
              f"{'chi2':<12} {'ratio':<10} {'order'}")
        print(f"  {'-'*90}")

        candidates = []

        # --- f_ls as frot, all orders ---
        for order in [1, 2, 3, 4, 5, 6]:
            c, chi2 = fit_nth_order(t, y, f_ls, order=order)
            r = amplitude_ratio(c)
            candidates.append((f_ls, c, chi2, r, 'f_ls', f'{order}th'))

        # --- f_ls/2 refined as frot, all orders ---
        frot_half_refined, _, _ = refine_period(t, y, f_ls / 2, window=window)
        for order in [1, 2, 3, 4, 5, 6]:
            c, chi2 = fit_nth_order(t, y, frot_half_refined, order=order)
            r = amplitude_ratio(c)
            candidates.append((frot_half_refined, c, chi2, r, 'f_ls/2 refined', f'{order}th'))

        for frot, coeffs, chi2, ratio, interp, order in candidates:
            marker = ' ◄ best' if chi2 < best_chi2 else ''
            print(f"  {interp+' '+order:<30} {frot:<14.6f} {24/frot:<14.4f} "
                  f"{chi2:<12.6f} {ratio:<10.4f}{marker}")
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_frot = frot
                best_coeffs = coeffs


    # keep the coeffs from whichever order won
    # don't refit to 2nd order — use the winning model
    ratio = amplitude_ratio(best_coeffs)
    is_double_peaked = ratio > 0.1
    period_hours = 24 / best_frot
    chi2_single, chi2_double = compare_models(t, y, best_frot)

    print(f"\n{'=' * 80}")
    print(f"FINAL RESULT")
    print(f"  Rotation frequency : {best_frot:.6f} c/d")
    print(f"  Rotation period    : {period_hours:.4f} hours")
    print(f"  Chi2 single        : {chi2_single:.6f}")
    print(f"  Chi2 double        : {chi2_double:.6f}")
    print(f"  Amplitude ratio    : {ratio:.4f}")
    print(f"  Double peaked      : {is_double_peaked}")
    print("=" * 80)

    return best_frot, period_hours, best_coeffs, bool(is_double_peaked)