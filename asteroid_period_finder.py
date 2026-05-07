"""
period_finder.py — Asteroid Rotation Period Pipeline
Entry point for the asteroid-period-finder pipeline.

Run from the repo root:
    python period_finder.py --file data/ALCDEF/ALCDEF_4_Vesta_....txt
    python period_finder.py --folder data/ALCDEF --ext txt --file_id Vesta
"""

import argparse
import sys
import numpy as np
from pathlib import Path

from pipeline import data_handler, fourier_engine, harmonic_logic, output_gen



def parse_args():
    parser = argparse.ArgumentParser(
        description="Asteroid rotation period finder using Lomb-Scargle + Fourier fitting.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file, asteroid ID from metadata
  python period_finder.py --file data/ALCDEF/ALCDEF_4_Vesta_20260505.txt

  # Single file, override asteroid ID
  python period_finder.py --file data/ALCDEF/ALCDEF_4_Vesta_20260505.txt --asteroid-id 4

  # Folder of txt files filtered by name
  python period_finder.py --folder data/ALCDEF --ext txt --file_id Vesta --asteroid-id 4

  # Custom output directory, no sigma filter, custom G slope
  python period_finder.py --file data/ALCDEF/ALCDEF_4_Vesta.txt \\
      --output-dir results/ --no-sigma-filter --G 0.25

  # Skip JSON output
  python period_finder.py --file data/ALCDEF/ALCDEF_4_Vesta.txt --no-json
        """
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--file', metavar='PATH', help='Path to a single lightcurve file (.txt or .tab).')
    input_group.add_argument('--folder', metavar='PATH', help='Path to a folder containing lightcurve files.' )

    parser.add_argument('--ext', metavar='EXTENSION', default='txt', help='File extension to look for when using --folder (default: txt). Do not include the dot!!!')
    parser.add_argument('--file_id', metavar='ID_STRING', dest='id_string', default=None, help='Optional filter string: only files whose names contain this string will be loaded (used with --folder).')

    parser.add_argument('--asteroid-id', metavar='ID', default=None,
        help=('Asteroid ID or name for JPL Horizons queries (e.g. "4" or "Vesta"). '
            'If not provided, the pipeline will attempt to read it from file metadata. '
            'Required for .tab files, which carry no metadata.'))

    parser.add_argument('--G', type=float, default=0.15, metavar='SLOPE',
        help='HG phase function slope parameter G for phase angle correction (default: 0.15).')

    parser.add_argument('--no-sigma-filter', action='store_true', help='Disable the running sigma clip outlier filter (enabled by default).')
    parser.add_argument('--output-dir', metavar='PATH', default='outputs', help='Directory where output files are saved (default: outputs/).')
    parser.add_argument('--output-name', metavar='NAME', default=None, 
        help=('Base name for output files (e.g. "Vesta" → Vesta_results.pdf and Vesta_results.json). '
            'Defaults to the asteroid name/ID from metadata or --asteroid-id.'))
    parser.add_argument('--no-json', action='store_true', help='Skip writing the JSON summary file (on by default).')

    parser.add_argument('--fmin', type=float, default=0.024, metavar='FREQ',help='Minimum search frequency in cycles/day (default: 0.024, i.e. ~1000 hr (~40 day) period).')
    parser.add_argument('--fmax',type=float, default=24.0, metavar='FREQ', help='Maximum search frequency in cycles/day (default: 24.0, i.e. 1 hr period). You must change this for fast asteroids, ie. ~100m fast rotators.')
    return parser.parse_args()


def load_data(args):
    """Load lightcurve data from a single file or a folder."""
    if args.file:
        print(f"Loading single file: {args.file}")
        data, metadata, status = data_handler.load_all_single(args.file)
    else:
        print(f"Loading folder: {args.folder}  extension=.{args.ext}  filter={args.id_string!r}")
        data, metadata, status = data_handler.load_all_multifile(args.folder, args.ext, args.id_string)

    if status < 0 or data is None:
        print("ERROR: Failed to load any lightcurve data. Check your input path and file format.")
        sys.exit(1)

    print(f"Loaded {len(data)} lightcurve segment(s).")
    return data, metadata


def resolve_asteroid_name(args, metadata):
    """
    Determine the asteroid name/ID to use throughout the run.
    Priority: --asteroid-id flag then metadata, else: error.
    """
    if args.asteroid_id is not None:
        return str(args.asteroid_id)

    # Try metadata from first segment
    name = data_handler.fetch_asteroid_id(metadata[0])
    if name is not None:
        print(f"Asteroid ID read from metadata: {name}")
        return str(name)

    print("ERROR: Could not determine asteroid ID from metadata and --asteroid-id was not provided.\n"
        "       Please rerun with --asteroid-id <id> (e.g. --asteroid-id 4 or --asteroid-id Vesta).")
    sys.exit(1)


def resolve_output_name(args, asteroid_id):
    """Determine the base name used for output files."""
    if args.output_name is not None:
        return args.output_name
    return asteroid_id


def main():
    args = parse_args()
    data, metadata = load_data(args) #Load the data, metadata

    asteroid_id = resolve_asteroid_name(args, metadata)
    output_name = resolve_output_name(args, asteroid_id)

    print(f"\nProcessing lightcurves for asteroid: {asteroid_id}")
    sigma_filter = not args.no_sigma_filter
    combined = data_handler.process_lightcurve(data, metadata, asteroid_id=asteroid_id, sigma_filter=sigma_filter)

    if 'jd_corrected' in combined.columns:
        t = combined['jd_corrected'].values
    else:
        t = combined[0].values

    y = combined[1].values

    print(f"Combined dataset: {len(t)} data points after processing.")
    print(f"\nRunning Lomb-Scargle periodgram  (fmin={args.fmin}, fmax={args.fmax}) ...")
    frequencies, powers, fpeak = fourier_engine.lomb_scargle(t, y, fmin=args.fmin, fmax=args.fmax)
    print(f"LS peak frequency: {fpeak:.6f} c/d  ({24/fpeak:.3f} hr)")

    print("\nRunning harmonic analysis ...")
    best_frot, period_hours, best_coeffs, is_double = harmonic_logic.harmonic_logic(t, y, frequencies, powers)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fig_path  = str(output_dir / f"{output_name}_results.pdf")
    json_path = str(output_dir / f"{output_name}_results.json") if not args.no_json else None

    print(f"\nGenerating outputs → {output_dir}/")

    chi2_single, chi2_double = harmonic_logic.compare_models(t, y, best_frot)

    output_gen.plot_results(t, y, frequencies, powers, period_hours, best_frot, best_coeffs, is_double, asteroid_name=output_name, save_path=fig_path)

    if json_path:
        output_gen.export_json(asteroid_name=output_name, period_hours=period_hours, best_frot=best_frot, best_coeffs=best_coeffs, is_double_peaked=is_double, chi2_single=chi2_single, chi2_double=chi2_double, save_path=json_path)

    # summary
    print("\n" + "="*50)
    print(f"  Asteroid         : {output_name}")
    print(f"  Rotation period  : {period_hours:.4f} hours")
    print(f"  Frequency        : {best_frot:.6f} cycles/day")
    print(f"  Double-peaked    : {is_double}")
    print(f"  chi2 (single)    : {chi2_single:.6f}")
    print(f"  chi2 (double)    : {chi2_double:.6f}")
    print(f"  PDF saved to     : {fig_path}")
    if json_path:
        print(f"  JSON saved to    : {json_path}")

if __name__ == '__main__':
    main()