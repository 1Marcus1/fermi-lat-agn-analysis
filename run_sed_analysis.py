"""
Fermi-LAT source-finding, spectral fit, and SED analysis for a single source,
built on top of `fermipy.gtanalysis.GTAnalysis`.

Workflow: run an initial optimization, free nearby sources and the diffuse
components, fit the region of interest (ROI), search for and remove
under-significant sources found in the residual/TS maps, fit the target
source's spectral energy distribution (SED), and localize it.

Usage
-----
    python run_sed_analysis.py

Expects a fermipy `config.yaml` in the working directory.
"""

from pathlib import Path

import numpy as np
from fermipy.gtanalysis import GTAnalysis

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CONFIG_PATH = "config.yaml"
SOURCE_NAME = "MKN501"

# Label used to tag this run's output directory and result files (e.g. the
# month or observation period being analyzed).
RUN_LABEL = "run1"
OUTPUT_DIR = f"{SOURCE_NAME.lower()}_fit_{RUN_LABEL}"

# Source-finding / cleanup settings
FREE_SOURCES_RADIUS_DEG = 10.0
NEW_SOURCE_TS_THRESHOLD = 5.0
NEW_SOURCE_MIN_SEPARATION_DEG = 0.5
DELETE_SOURCES_MAX_TS = 9.0  # sources with TS in [-inf, this] are removed


def build_gta(config_path: str, output_dir: str) -> GTAnalysis:
    """Construct a `GTAnalysis` instance for this run.

    Parameters
    ----------
    config_path : str
        Path to the fermipy `config.yaml` file.
    output_dir : str
        Directory fermipy will write its output (plots, ROI files, logs) to.

    Returns
    -------
    GTAnalysis
        A configured (but not yet set up) analysis object.
    """
    return GTAnalysis(config_path, logging={"verbosity": 3}, fileio={"outdir": output_dir})


def fit_region_of_interest(gta: GTAnalysis) -> None:
    """Run the initial ROI optimization and fit, freeing normalizations.

    Frees the normalization of all point sources within
    `FREE_SOURCES_RADIUS_DEG` of the ROI center, plus the Galactic and
    isotropic diffuse components, then performs a full likelihood fit.

    Parameters
    ----------
    gta : GTAnalysis
        Analysis object, already set up (`gta.setup()` already called).
    """
    gta.optimize()
    gta.free_sources(distance=FREE_SOURCES_RADIUS_DEG, pars="norm")
    gta.free_source("galdiff", pars="norm")
    gta.free_source("isodiff", pars="norm")
    gta.fit()


def find_and_remove_spurious_sources(gta: GTAnalysis) -> None:
    """Search for additional sources in the residual/TS maps and prune weak ones.

    Produces residual and TS maps, searches for new sources above
    `NEW_SOURCE_TS_THRESHOLD`, refreshes the TS map, then deletes any source
    with `TS <= DELETE_SOURCES_MAX_TS` (including any spurious sources just
    added by `find_sources`).

    Parameters
    ----------
    gta : GTAnalysis
        Analysis object, already fit at least once.
    """
    gta.residmap(make_plots=True)
    gta.tsmap(make_plots=True)
    gta.find_sources(sqrt_ts_threshold=NEW_SOURCE_TS_THRESHOLD, min_separation=NEW_SOURCE_MIN_SEPARATION_DEG)
    gta.tsmap(make_plots=True)
    gta.delete_sources(minmax_ts=[-np.inf, DELETE_SOURCES_MAX_TS])


def run_sed_and_localization(gta: GTAnalysis, source_name: str, run_label: str) -> dict:
    """Fit the target source's SED and localize it, saving results to disk.

    Parameters
    ----------
    gta : GTAnalysis
        Analysis object, with the ROI already cleaned up
        (`find_and_remove_spurious_sources` already called).
    source_name : str
        Name of the target source, as known to fermipy/the source catalog.
    run_label : str
        Label used to tag the output `sed_<label>.txt` /
        `spectral_index_<label>.txt` files.

    Returns
    -------
    dict
        The raw SED result dictionary returned by `gta.sed`.
    """
    gta.free_source(source_name)
    gta.write_roi("fit0", make_plots=True)

    gta.sed(source_name, make_plots=True)
    sed = gta.sed(source_name, outfile="sed.fits")

    print("SED")
    print("dnde:        ", sed["dnde"])
    print("e2dnde:      ", sed["e2dnde"])
    print("e2dnde_err:  ", sed["e2dnde_err"])
    print("spectral index")
    print("alpha:       ", sed["param_values"][1])
    print("alpha_err:   ", sed["param_errors"][1])
    print("beta:        ", sed["param_values"][2])
    print("beta_err:    ", sed["param_errors"][2])

    np.savetxt(
        f"sed_{run_label}.txt",
        np.c_[sed["dnde"], sed["e2dnde"], sed["e2dnde_err"]],
        delimiter=";",
        header="dnde;e2dnde;e2dnde_err",
    )
    np.savetxt(
        f"spectral_index_{run_label}.txt",
        np.c_[sed["param_values"][1], sed["param_errors"][1], sed["param_values"][2], sed["param_errors"][2]],
        delimiter=";",
        header="alpha;alpha_err;beta;beta_err",
    )

    gta.localize(source_name)
    return sed


def main() -> None:
    gta = build_gta(CONFIG_PATH, OUTPUT_DIR)
    gta.setup()

    fit_region_of_interest(gta)
    find_and_remove_spurious_sources(gta)
    run_sed_and_localization(gta, SOURCE_NAME, RUN_LABEL)


if __name__ == "__main__":
    main()
