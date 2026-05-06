import numpy as np
import matplotlib.pyplot as plt
import json
from pathlib import Path


def phase_fold(t, period_hours):
    """Fold times onto phase using period in hours"""
    period_days = period_hours / 24
    phase = (t % period_days) / period_days
    return phase


def plot_results(t, y, frequencies, powers, period_hours, best_frot,
                 best_coeffs, is_double_peaked, asteroid_name=None,
                 save_path=None):
    """
    Dual panel figure:
    Top: Lomb-Scargle periodogram with LS peak marked
    Bottom: Phase folded lightcurve with Fourier model
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    title = f"Asteroid {asteroid_name}" if asteroid_name else "Asteroid"
    fig.suptitle(title, fontsize=14)

    # --- Top panel: Periodogram ---
    periods_hr = 24 / frequencies

    # LS peak is at 2*frot for double peaked, frot for single
    ls_peak_hours = 24 / (best_frot * 2) if is_double_peaked else 24 / best_frot

    ax1.plot(periods_hr, powers, 'b-', linewidth=0.5, alpha=0.7)
    ax1.axvline(ls_peak_hours, color='red', linestyle='--', linewidth=0.5,
                label=f'LS peak = {ls_peak_hours:.3f} hr → P_rot = {period_hours:.3f} hr')
    ax1.set_xlabel('Period (hours)')
    ax1.set_ylabel('Lomb-Scargle Power')
    ax1.set_title('Periodogram')
    ax1.set_xlim(50, 0)  # inverted: long periods on left, short on right
    ax1.legend()

    # --- Bottom panel: Phase folded lightcurve ---
    phase = phase_fold(t, period_hours)

    # Sort by phase for model line
    sort_idx = np.argsort(phase)
    phase_sorted = phase[sort_idx]

    # Evaluate Fourier model over sorted phase
    t_model = phase_sorted * (period_hours / 24)
    from pipeline.harmonic_logic import fourier_model
    y_model = fourier_model(t_model, best_frot, best_coeffs)

    ax2.scatter(phase, y, s=2, alpha=0.5, color='blue', label='Data')
    ax2.plot(phase_sorted, y_model, 'r-', linewidth=1.5, label='Fourier model')
    ax2.set_xlabel('Phase')
    ax2.set_ylabel('Magnitude')
    ax2.set_title(f'Phase Folded Lightcurve (P = {period_hours:.3f} hr, '
                  f'{"double" if is_double_peaked else "single"} peaked)')
    ax2.invert_yaxis()
    ax2.legend()

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    plt.show()
    return fig


def export_json(asteroid_name, period_hours, best_frot, best_coeffs,
                is_double_peaked, chi2_single, chi2_double,
                phase_slope=None, save_path=None):
    """Export JSON summary of results."""
    summary = {
        "asteroid": asteroid_name,
        "rotation_period_hours": round(period_hours, 4),
        "rotation_frequency_cpd": round(best_frot, 6),
        "is_double_peaked": is_double_peaked,
        "fourier_coefficients": {
            "C":  round(float(best_coeffs[0]), 6),
            "A1": round(float(best_coeffs[1]), 6),
            "B1": round(float(best_coeffs[2]), 6),
            "A2": round(float(best_coeffs[3]), 6),
            "B2": round(float(best_coeffs[4]), 6),
        },
        "chi2_single_harmonic": round(chi2_single, 6),
        "chi2_double_harmonic": round(chi2_double, 6),
        "chi2_ratio": round(chi2_single / chi2_double, 4),
        "phase_slope_beta": round(phase_slope, 4) if phase_slope is not None else None,
    }

    if save_path:
        Path(save_path).write_text(json.dumps(summary, indent=2))
        print(f"JSON summary saved to {save_path}")

    return summary

def generate_output(t, y, frequencies, powers, period_hours, best_frot,
                    best_coeffs, is_double_peaked, asteroid_name=None,
                    output_dir=None):
    """
    Main entry point — generates both figure and JSON.
    """
    from pipeline.harmonic_logic import compare_models
    chi2_single, chi2_double = compare_models(
        np.asarray(t), np.asarray(y), best_frot
    )

    fig_path = None
    json_path = None

    if output_dir:
        name = asteroid_name or "asteroid"
        fig_path = f"{output_dir}/{name}_results.png"
        json_path = f"{output_dir}/{name}_results.json"

    fig = plot_results(
        np.asarray(t), np.asarray(y),
        frequencies, powers,
        period_hours, best_frot, best_coeffs,
        is_double_peaked, asteroid_name,
        save_path=fig_path
    )

    summary = {
        "asteroid": str(asteroid_name),
        "rotation_period_hours": float(round(period_hours, 4)),
        "rotation_frequency_cpd": float(round(best_frot, 6)),
        "is_double_peaked": bool(is_double_peaked),
        "fourier_coefficients": {
            "C":  float(round(float(best_coeffs[0]), 6)),
            "A1": float(round(float(best_coeffs[1]), 6)),
            "B1": float(round(float(best_coeffs[2]), 6)),
            "A2": float(round(float(best_coeffs[3]), 6)),
            "B2": float(round(float(best_coeffs[4]), 6)),
        },
        "chi2_single_harmonic": float(round(chi2_single, 6)),
        "chi2_double_harmonic": float(round(chi2_double, 6)),
        "chi2_ratio": float(round(chi2_single / chi2_double, 4))
        }

    return fig, summary