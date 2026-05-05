import numpy as np

def compute_tau(f, t):
    """Compute time offset tau from equation 5"""
    numerator = np.sum(np.sin(4 * np.pi * f * t))
    denominator = np.sum(np.cos(4 * np.pi * f * t))
    return np.arctan2(numerator, denominator) / (4 * np.pi * f)

def compute_power(f, t, y):
    """Compute Lomb-Scargle power at frequency f from equation. Reelveant equation reference in Wiki"""
    y_mean = np.mean(y)
    y_var = np.var(y)
    tau = compute_tau(f, t)

    cos_term = np.cos(2 * np.pi * f * (t - tau))
    sin_term = np.sin(2 * np.pi * f * (t - tau))

    y_centered = y - y_mean

    power = (1 / (2 * y_var)) * (
        (np.sum(y_centered * cos_term) ** 2) / np.sum(cos_term ** 2) +
        (np.sum(y_centered * sin_term) ** 2) / np.sum(sin_term ** 2)
    )

    return power

def lomb_scargle(t, y, fmin=0.024, fmax=12.0):
    """
    Compute Lomb-Scargle periodogram.
    
    t: time array (days)
    y: magnitude array
    fmin: minimum frequency (cycles/day), default 1/1000hr
    fmax: maximum frequency (cycles/day), default 1/2hr
    
    Returns: frequencies, powers, fpeak
    """
    T = np.minimum(t.max() - t.min(), 365) # Longest known asteroid rotation period is ~200 days, for multi year observations, we know want an overly inflated T
    df = 1 / (10 * T)  # frequency grid spacing from spec

    frequencies = np.arange(fmin, fmax, df)
    powers = np.array([compute_power(f, t, y) for f in frequencies])

    fpeak = frequencies[np.argmax(powers)]

    return frequencies, powers, fpeak