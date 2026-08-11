# Fermi-LAT AGN Analysis

Fermi-LAT source-finding, spectral (SED), and light curve analysis for AGN,
built on top of [fermipy](https://fermipy.readthedocs.io/). Includes a full
region-of-interest (ROI) fit and SED extraction for a single time window
(`run_sed_analysis.py`) and a binned light curve extraction that re-fits the
ROI independently in each time bin (`run_lightcurve_analysis.py`).

## Overview

Both scripts follow the standard fermipy workflow: load a region of interest
around the target source, run an initial optimization, free the normalization
of nearby sources and the diffuse backgrounds, perform a full likelihood fit,
and extract science products (a spectral energy distribution, a light curve
point, a source position). `run_sed_analysis.py` additionally searches the
residual/TS maps for missed sources and prunes spurious low-significance ones
before the final fit.

## Repository structure

```
fermi-lat-agn-analysis/
├── fermi_lat_utils.py           # Shared utilities: config parsing, time binning
├── run_sed_analysis.py          # Single-window ROI fit + SED + localization
├── run_lightcurve_analysis.py   # Time-binned light curve extraction
├── requirements.txt
├── config.yaml                  # fermipy configuration (see "Data" below)
├── ft1.fits, ft2.fits            # Fermi-LAT photon & spacecraft files (not tracked)
├── gll_iem_v07.fits              # Galactic diffuse emission model (not tracked)
├── iso_P8R3_SOURCE_V3_v1.txt     # Isotropic diffuse template (not tracked)
├── gll_psc_v*.fit                # 4FGL source catalog (not tracked)
└── README.md
```

**`fermi_lat_utils.py`** reads the `tmin`/`tmax` time selection from
`config.yaml` (via a proper YAML parser, not by assuming fixed line numbers)
and computes the equal-width time bins used by the light curve script.

**`run_sed_analysis.py`** fits the ROI once over the full time range defined
in `config.yaml`, cleans up the source model (residual/TS maps, new-source
search, pruning of low-TS sources), and extracts the target source's SED and
best-fit position.

**`run_lightcurve_analysis.py`** splits `config.yaml`'s time range into
equal-width bins (1 day by default), runs an independent ROI fit + SED in
each bin, and collects the energy flux into a light curve.

## Installation

`fermipy` additionally requires the **Fermi Science Tools** (or the
`fermitools` conda package) to be installed and initialized in your
environment — `pip install` alone does not provide the underlying
Fermi-LAT instrument response functions and low-level analysis tools that
fermipy wraps. See the
[fermipy installation guide](https://fermipy.readthedocs.io/en/latest/install.html)
if you don't already have a working Fermi Science Tools setup.

## Data

Both scripts expect a fermipy `config.yaml` in the working directory, along
with the data files it references. None of this data is included in this
repository (it's source- and time-range-specific, and some files are large);
here's how to get it for a given source and time window.

**1. Query the Fermi-LAT photon and spacecraft data.**
Go to the FSSC
[LAT Data Query page](https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi)
and submit a query centered on your target's RA/Dec, with the energy range,
time range (in MET seconds — this is what `config.yaml`'s `selection.tmin`/
`tmax` should match), and search radius you want to analyze. The query
returns a photon file (`..._PH*.fits`, referred to as `evfile` in
`config.yaml`) and a spacecraft file (`..._SC*.fits`, `scfile`). Fermipy can
also submit this query for you automatically if `evfile`/`scfile` aren't
already present — see the
[fermipy data selection docs](https://fermipy.readthedocs.io/en/latest/config.html#data).

**2. Download the diffuse background models.**
From the FSSC
[background models page](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html),
download the current Galactic diffuse emission model
(`gll_iem_v07.fits` or the current recommended version) and the isotropic
spectral template matching your event class/type
(`iso_P8R3_SOURCE_V3_v1.txt` for `SOURCE` class, Pass 8).

**3. Download the source catalog.**
From the FSSC
[LAT catalogs page](https://fermi.gsc.nasa.gov/ssc/data/access/lat/), download
the current 4FGL catalog FITS file (`gll_psc_v*.fit`), used by fermipy to
seed the ROI with known nearby sources.

**4. Place everything in the same folder.**
Put `config.yaml`, the photon/spacecraft files, the diffuse models, and the
catalog file together in the directory you'll run the scripts from (or point
to them with absolute paths in `config.yaml`), then reference them in
`config.yaml`, e.g.:

```yaml
data:
  evfile : ft1.fits
  scfile : ft2.fits

selection:
  ra      : <target RA, deg>
  dec     : <target Dec, deg>
  radius  : 15
  emin    : 100
  emax    : 300000
  tmin    : <start, MET seconds>
  tmax    : <end, MET seconds>
  zmax    : 90
  evclass : 128
  evtype  : 3

model:
  galdiff  : gll_iem_v07.fits
  isodiff  : iso_P8R3_SOURCE_V3_v1.txt
  catalogs : ['gll_psc_v*.fit']
```

See the [fermipy configuration reference](https://fermipy.readthedocs.io/en/latest/config.html)
for the full set of options (binning, ROI width, etc.).

## Usage

```bash
python run_sed_analysis.py            # single-window ROI fit + SED
python run_lightcurve_analysis.py     # time-binned light curve
```

Edit the configuration block at the top of each script (source name,
`config.yaml` path, run label, bin length for the light curve) before
running. Both scripts write their output (ROI files, plots, SED/light curve
text files) into the working directory and per-bin subdirectories.

## References

- [Fermipy documentation](https://fermipy.readthedocs.io/)
- [Fermi Science Support Center (FSSC)](https://fermi.gsc.nasa.gov/ssc/)
- Abdollahi, S., et al. (2022). *Incremental Fermi Large Area Telescope
  Fourth Source Catalog.* ApJS, 260, 53 (4FGL-DR3).
  [arXiv:2201.11184](https://arxiv.org/abs/2201.11184)
