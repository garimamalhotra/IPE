# Standalone IPE with a neutral file reader

This branch builds IPE as a standalone model (no WAM/coupling) and adds the
option to drive it with a neutral atmosphere read from NetCDF files instead of
the built-in MSIS/HWM empirical models. Three neutral sources are supported from
the same executable:

1. **MSIS-HWM** (empirical) — the default.
2. **GSM-format files** on the 90×91 geographic grid.
3. **High-resolution files** (e.g. FV3WAM, 384×190) interpolated directly to the
   apex grid.

The reader selects the legacy (90×91) or direct-to-apex path automatically from
the file dimensions, and interpolates linearly in time between hourly frames.

## Requirements

An MPI Fortran toolchain plus these libraries (typically from a module or spack
environment):

- `mpif90` / `mpicc`
- COMIO
- ESMF (`esmf.mk`)
- netCDF-Fortran (`nf-config` in `PATH`) and netCDF-C
- HDF5, (parallel-)netCDF

## Build

```bash
export ESMFMKFILE=/path/to/esmf.mk
export COMIO_PREFIX=/path/to/comio/install     # or have comio-config in PATH
./build_standalone.sh
```

The script configures with coupling disabled and adds the netCDF-Fortran link
flags from `nf-config` (the reader calls the `nf90_*` API, which the base build
does not link). The executable is installed to `install/bin/ipe.x`.

## Selecting the neutral source

The neutral source is chosen in the `IPE.inp` namelist:

**MSIS-HWM** — omit `&NeutralFileIO`, or:
```fortran
&NeutralFileIO
  READ_GSM_NEUTRALS = F
/
```

**File-based (GSM or high-resolution)** — point `GSM_NEUTRALS_DIR` at a directory
of `gsm.YYYYMMDD_HHMMSS.nc` files:
```fortran
&NeutralFileIO
  READ_GSM_NEUTRALS     = T,
  GSM_NEUTRALS_DIR      = '/path/to/neutral_files/',
  NEUTRAL_INTERP_METHOD = 'linear',
  VERBOSE_DIAG_LOCAL    = F,
  NEUTRAL_T_SCALE  = 1.0,
  NEUTRAL_O_SCALE  = 1.0,
  NEUTRAL_O2_SCALE = 1.0,
  NEUTRAL_N2_SCALE = 1.0
/
```

The same block works for both the 90×91 GSM grid and higher-resolution sources;
the reader detects which path to use from the file's `nlon`/`nlat`. The optional
`NEUTRAL_*_SCALE` factors multiply the corresponding neutral field (default 1.0).
`VERBOSE_DIAG_LOCAL = T` enables the diagnostic output.

## Preparing FV3WAM files

FV3WAM history files (`atmf*.nc`) are converted to the `gsm.*.nc` format with
`scripts/fv3wam_to_gsm.py`:

```bash
python scripts/fv3wam_to_gsm.py -i <atmf_dir> -o <gsm_dir> --nproc 8
```

Point `GSM_NEUTRALS_DIR` at `<gsm_dir>`. The vertical levels are fixed from the
first file (or `--ref`) so every output has the same number of levels, which the
reader requires. Needs `h5py`, `hdf5plugin`, `netCDF4`, `numpy`.

## Input file format

Each hourly `gsm.YYYYMMDD_HHMMSS.nc` file holds, on a `(lon, lat, level)` grid:

| variable | meaning |
|---|---|
| `temp_neutral` | neutral temperature (K) |
| `O_Density`, `O2_Density`, `N2_Density` | number densities (m⁻³) |
| `eastward_wind_neutral`, `northward_wind_neutral`, `upward_wind_neutral` | winds (m/s) |
| `wam_height_levels` | geometric height of each level (m) |
| `grid_xt`, `grid_yt`, `pfull` | longitude, latitude, level coordinates |

Grid, initial condition, and solar-forcing files are provided separately (they
are large and not part of the repository).

## Output

IPE writes the full apex-grid state to `IPE_State.apex.<timestamp>.h5` at
`FILE_OUTPUT_FREQUENCY`. Each file is a complete state and can be used directly
as the initial condition for a restart.
