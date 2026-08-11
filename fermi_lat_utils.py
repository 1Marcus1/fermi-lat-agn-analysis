"""
Shared utilities for Fermi-LAT AGN analyses built on top of fermipy.

These helpers are used by both `run_sed_analysis.py` and
`run_lightcurve_analysis.py`: reading the time selection from a fermipy
`config.yaml`, and building equal-width time bins for a light curve.
"""

from pathlib import Path

import math

import numpy as np
import yaml


def load_time_selection(config_path: str) -> tuple[float, float]:
    """Read the `tmin`/`tmax` time selection from a fermipy config YAML file.

    Reads the actual `selection: {tmin, tmax}` keys via a YAML parser, rather
    than assuming they sit on fixed line numbers of the file (the original
    implementation read `lista[12]`/`lista[13]` directly, which silently
    breaks if the config file is reordered or reformatted).

    Parameters
    ----------
    config_path : str
        Path to the fermipy `config.yaml` file.

    Returns
    -------
    t_min, t_max : float
        Mission Elapsed Time (MET) bounds of the time selection, in seconds.

    Raises
    ------
    KeyError
        If the config file does not have a `selection.tmin`/`selection.tmax` entry.
    """
    with open(config_path) as fh:
        config = yaml.safe_load(fh)

    try:
        selection = config["selection"]
        return float(selection["tmin"]), float(selection["tmax"])
    except KeyError as exc:
        raise KeyError(
            f"Could not find 'selection.tmin'/'selection.tmax' in {config_path!r}. "
            "Check that this is a valid fermipy config file."
        ) from exc


def compute_time_bins(t_min: float, t_max: float, bin_length_seconds: float = 86400.0) -> np.ndarray:
    """Compute equal-width time bin edges between `t_min` and `t_max`.

    The number of bins is `ceil((t_max - t_min) / bin_length_seconds)`; if the
    total time range is not an exact multiple of `bin_length_seconds`, the
    *last* bin absorbs the remainder and is therefore wider than the others
    (matching the original implementation's behavior, which is equivalent to
    this `ceil` even though it wasn't expressed that way explicitly).

    Parameters
    ----------
    t_min, t_max : float
        Start and end of the time range (MET seconds).
    bin_length_seconds : float, default 86400.0 (1 day)
        Nominal bin width, in seconds.

    Returns
    -------
    np.ndarray
        Bin edges, length `n_bins + 1`, where
        `n_bins = ceil((t_max - t_min) / bin_length_seconds)`. The final
        edge is always exactly `t_max`.

    Raises
    ------
    ValueError
        If the time range is shorter than one bin.
    """
    n_bins_float = (t_max - t_min) / bin_length_seconds
    n_bins = math.ceil(n_bins_float)
    if n_bins < 1:
        raise ValueError(
            f"Time range ({t_max - t_min:.0f} s) is shorter than one bin "
            f"({bin_length_seconds:.0f} s); no bins to compute."
        )
    edges = t_min + bin_length_seconds * np.arange(n_bins)
    return np.append(edges, t_max)


def bin_centroids(bin_edges: np.ndarray, bin_length_seconds: float = 86400.0) -> np.ndarray:
    """Compute the centroid of each bin, given its edges.

    Every bin's centroid is offset from its start edge by half of the
    *nominal* bin width, `bin_length_seconds`, matching the original
    implementation. Note this makes the centroid of the last (possibly
    wider, see `compute_time_bins`) bin only approximate.

    Parameters
    ----------
    bin_edges : np.ndarray
        Bin edges, as returned by `compute_time_bins` (length `n_bins + 1`).
    bin_length_seconds : float, default 86400.0
        Nominal bin width used to offset each start edge to its centroid.

    Returns
    -------
    np.ndarray
        Centroid of each of the `n_bins` bins.
    """
    bin_starts = bin_edges[:-1]
    return bin_starts + 0.5 * bin_length_seconds


def build_output_names(source_label: str, bin_starts: np.ndarray) -> tuple[list[str], list[str]]:
    """Build per-bin output directory and ROI file base names.

    Parameters
    ----------
    source_label : str
        Short label identifying the source, used as a filename/directory prefix
        (e.g. "J1653.8+3945").
    bin_starts : np.ndarray
        Start time of each bin (MET seconds), typically `compute_time_bins(...)[:-1]`.

    Returns
    -------
    directories, roi_names : list[str]
        One directory name and one ROI base name per bin, in the same order
        as `bin_starts`.
    """
    directories = [f"{source_label}_{t}" for t in bin_starts]
    roi_names = [f"{source_label}_fit_{t}" for t in bin_starts]
    return directories, roi_names
