# Asteroid Period Finder

A pipeline for determining asteroid rotation periods from photometric lightcurve data. It fetches ephemeris data from JPL Horizons, applies phase angle and light travel time corrections, runs a Lomb-Scargle periodogram, and fits a Fourier model to determine the best rotation period.

---

## Directory Structure

```
asteroid-period-finder/
├── period_finder.py       # Main entry point — run this
├── pipeline/
│   ├── data_handler.py    # Loading, corrections, sigma filtering
│   ├── fourier_engine.py  # Lomb-Scargle periodogram
│   ├── harmonic_logic.py  # Period refinement and harmonic selection
│   ├── io.py              # File readers (.txt, .tab, .csv, .json)
│   └── output_gen.py      # Plot and JSON export
├── data/
│   └── ALCDEF/            # Put your lightcurve files here
├── outputs/               # Default output location (PDF + JSON)
├── tests/
│   └── module_tests.ipynb
└── readme.md
```

---

## Installation

Python 3.10+ is required. Install dependencies with:

```bash
pip install numpy pandas matplotlib astropy astroquery sbpy
```

All dependencies:

| Package      | Purpose                                      |
|--------------|----------------------------------------------|
| `numpy`      | Numerical computation                        |
| `pandas`     | Dataframe handling for lightcurve data       |
| `matplotlib` | Plotting periodogram and phase-folded curve  |
| `astropy`    | Units and astropy table handling             |
| `astroquery` | JPL Horizons ephemeris queries               |
| `sbpy`       | HG phase function (photometric correction)   |

---

## Basic Usage

Always run `period_finder.py` from the **repo root** so that `pipeline/` imports resolve correctly.

### Single file

```bash
python period_finder.py --file data/ALCDEF/ALCDEF_4_Vesta_20260505_185525.txt
```

### Single file with explicit asteroid ID

Needed when the file metadata does not contain `OBJECTNUMBER` or `OBJECTNAME` (e.g. `.tab` files).

```bash
python period_finder.py --file data/LCDB/asteroid.tab --asteroid-id 4
```

### Folder of files

Loads all files in a folder matching the given extension. Use `--file_id` to filter by filename substring.

```bash
python period_finder.py --folder data/ALCDEF --ext txt --file_id Vesta --asteroid-id 4
```

### Custom output location

```bash
python period_finder.py --file data/ALCDEF/ALCDEF_4_Vesta.txt --output-dir results/
```

### Skip JSON, disable sigma filter, set custom G slope

```bash
python period_finder.py --file data/ALCDEF/ALCDEF_4_Vesta.txt \
    --no-json --no-sigma-filter --G 0.25
```

---

## All Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--file PATH` | — | Path to a single lightcurve file. Mutually exclusive with `--folder`. |
| `--folder PATH` | — | Path to a folder of lightcurve files. Mutually exclusive with `--file`. |
| `--ext EXT` | `txt` | File extension to scan for when using `--folder`. |
| `--file_id STRING` | None | Filename substring filter when using `--folder` (e.g. `Vesta`). |
| `--asteroid-id ID` | from metadata | Asteroid number or name for JPL Horizons (e.g. `4` or `Vesta`). Required if metadata does not contain it. |
| `--G SLOPE` | `0.15` | HG phase function slope parameter for phase angle correction. |
| `--no-sigma-filter` | off | Disables the running sigma-clip outlier filter. |
| `--output-dir PATH` | `outputs/` | Directory where output PDF and JSON are written. Created if it does not exist. |
| `--output-name NAME` | asteroid ID | Base name for output files (e.g. `Vesta` → `Vesta_results.pdf`). |
| `--no-json` | off | Skips writing the JSON summary file. |
| `--fmin FREQ` | `0.024` | Minimum search frequency in cycles/day (~1000 hr period). |
| `--fmax FREQ` | `12.0` | Maximum search frequency in cycles/day (~2 hr period). |

---

## Outputs

By default, two files are written to `outputs/` (or your `--output-dir`):

**`<name>_results.pdf`** — A two-panel figure:
- Top: Lomb-Scargle periodogram with the LS peak marked
- Bottom: Phase-folded lightcurve with the best-fit Fourier model overlaid

**`<name>_results.json`** — A summary including rotation period, frequency, Fourier coefficients, chi² values for single and double harmonic models, and whether the lightcurve is double-peaked. Pass `--no-json` to skip this.

---

## Supported File Formats

### `.txt` — ALCDEF format (preferred)

Files downloaded from [alcdef.org](https://alcdef.org). These are structured text files with embedded metadata blocks so the pipeline knows which corrections have already been applied and does not double-correct.

A single `.txt` file may contain multiple observation sessions, each with its own metadata block.

An important note about the ALCDEF files: They compile observations across multiple observation periods from different observers who often choose to handle their data in different ways. The most notable implication of this is the difference in zero points across observations. As such, we have implemented a rough normalizer, but ideally, all data should be inputted with the same zero points. In future, this will be supported through .csv user reduced light curves. 

### `.tab` — LCDB format

Tab-separated files from the Asteroid Lightcurve Data Base (LCDB). These files contain no embedded metadata. The pipeline assumes they are pre-calibrated differential magnitudes in MJD with no light travel time correction applied. Because there is no `OBJECTNUMBER` or `OBJECTNAME` in the file, **`--asteroid-id` is required** when using `.tab` files.

### `.csv` — *(not yet implemented)*

The file reader has a placeholder for CSV support. It is not functional in the current version.

### `.json` — *(not yet implemented)*

The file reader has a placeholder for JSON support. It is not functional in the current version.

---

## Notes

- The pipeline queries JPL Horizons for ephemeris data (RA, Dec, heliocentric distance, geocentric distance, phase angle). *An internet connection is required.*
- Ephemeris queries are chunked to avoid server timeouts. A short delay is introduced between chunks to respect rate limits.
- Observer location defaults to geocenter (`500`). Higher-accuracy topocentric queries are not yet wired up.
- The Lomb-Scargle implementation follows the Scargle (1982) formulation and is intentionally written from scratch rather than wrapping `astropy.timeseries.LombScargle`, for pipeline-specific control over frequency gridding and period refinement.
- Period refinement tests each LS cluster peak as both the fundamental and the double frequency, selecting whichever interpretation minimises reduced chi².
