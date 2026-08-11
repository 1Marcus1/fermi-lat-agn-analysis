"""
Fermi-LAT light curve extraction for a single source, built on top of
`fermipy.gtanalysis.GTAnalysis`.

Splits the time range defined in `config.yaml` into equal-width bins (see
`fermi_lat_utils.compute_time_bins`), runs an independent ROI fit + SED in
each bin, collects the energy flux and its uncertainty from each bin's ROI
output, and saves/plots the resulting light curve.

Usage
-----
    python run_lightcurve_analysis.py

Expects a fermipy `config.yaml` (with a `selection.tmin`/`selection.tmax`
time range) in the working directory.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from fermipy.gtanalysis import GTAnalysis

from fermi_lat_utils import bin_centroids, build_output_names, compute_time_bins, load_time_selection

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CONFIG_PATH = "config.yaml"
SOURCE_NAME_4FGL = "4FGL J1653.8+3945"  # name as known to fermipy/the source catalog
SOURCE_LABEL = "J1653.8+3945"  # short label used for output directory/file names

BIN_LENGTH_SECONDS = 86400.0  # 1 day
RUN_LABEL = "run1"  # tags the output lightcurve_<label>.txt / .png files

FREE_SOURCES_RADIUS_DEG = 10.0


def run_one_bin(config_path: str, t_start: float, t_end: float, output_dir: str, roi_name: str) -> None:
    """Run a full ROI fit and SED for a single time bin.

    Parameters
    ----------
    config_path : str
        Path to the fermipy `config.yaml` file.
    t_start, t_end : float
        Start and end of this bin's time selection (MET seconds).
    output_dir : str
        Directory fermipy will write this bin's output to.
    roi_name : str
        Base filename (without extension) fermipy will use for the saved ROI.
    """
    gta = GTAnalysis(
        config_path,
        selection={"tmin": t_start, "tmax": t_end},
        logging={"verbosity": 3},
        fileio={"outdir": output_dir},
    )
    gta.setup()
    gta.optimize()
    gta.free_sources(distance=FREE_SOURCES_RADIUS_DEG, pars="norm")
    gta.free_source("galdiff", pars="norm")
    gta.free_source("isodiff", pars="norm")
    gta.fit()
    gta.sed(SOURCE_NAME_4FGL)
    gta.write_roi(roi_name)


def collect_energy_flux(directories: list[str], roi_names: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Read the fitted energy flux and its uncertainty back from each bin's saved ROI.

    Parameters
    ----------
    directories : list[str]
        Per-bin output directories, as returned by `fermi_lat_utils.build_output_names`.
    roi_names : list[str]
        Per-bin ROI base filenames, as returned by `fermi_lat_utils.build_output_names`.

    Returns
    -------
    eflux, eflux_err : np.ndarray
        Energy flux and uncertainty for the target source, one value per bin.
    """
    home = os.getcwd()
    eflux, eflux_err = [], []
    try:
        for directory, roi_name in zip(directories, roi_names):
            os.chdir(directory)
            try:
                results = Table.read(roi_name + ".fits")
                eflux.append(results["eflux"][0])
                eflux_err.append(results["eflux_err"][0])
            finally:
                os.chdir(home)
    finally:
        os.chdir(home)
    return np.array(eflux), np.array(eflux_err)


def plot_lightcurve(
    t_centroid: np.ndarray, eflux: np.ndarray, eflux_err: np.ndarray, bin_length_seconds: float, output_path: str
) -> None:
    """Plot the energy flux light curve with symmetric time/flux error bars."""
    fig, ax = plt.subplots()
    ax.errorbar(
        t_centroid, eflux, xerr=bin_length_seconds / 2, yerr=eflux_err, fmt="o", ecolor="k",
    )
    ax.set_xlim(t_centroid[0] - bin_length_seconds, t_centroid[-1] + bin_length_seconds)
    ax.set_yscale("log")
    ax.grid(linestyle="--", linewidth=0.5)
    ax.set_xlabel("Time (MET)")
    ax.set_ylabel(r"Energy flux")
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def main() -> None:
    t_min, t_max = load_time_selection(CONFIG_PATH)
    bin_edges = compute_time_bins(t_min, t_max, BIN_LENGTH_SECONDS)
    bin_starts, bin_ends = bin_edges[:-1], bin_edges[1:]
    t_centroid = bin_centroids(bin_edges, BIN_LENGTH_SECONDS)
    directories, roi_names = build_output_names(SOURCE_LABEL, bin_starts)

    for t_start, t_end, directory, roi_name in zip(bin_starts, bin_ends, directories, roi_names):
        run_one_bin(CONFIG_PATH, t_start, t_end, directory, roi_name)

    eflux, eflux_err = collect_energy_flux(directories, roi_names)

    print("Light curve")
    print("x:     ", t_centroid)
    print("x_err: ", BIN_LENGTH_SECONDS / 2)
    print("y:     ", eflux)
    print("y_err: ", eflux_err)

    np.savetxt(
        f"lightcurve_{RUN_LABEL}.txt",
        np.c_[t_centroid, eflux, eflux_err],
        delimiter=";",
        header="time;flux;flux_error",
    )
    plot_lightcurve(t_centroid, eflux, eflux_err, BIN_LENGTH_SECONDS, f"lightcurve_{RUN_LABEL}.png")


if __name__ == "__main__":
    main()
