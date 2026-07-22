#!/usr/bin/env python3
"""
Convert FV3WAM atmf*.nc history files to the gsm-format NetCDF files read by
the IPE neutral file reader. Native 384x190; height coordinate from delz+hgtsfc,
number densities from mixing ratios, latitude flipped S->N, output named
gsm.YYYYMMDD_HH0000.nc.

  python fv3wam_to_gsm.py -i <atmf_dir> -o <gsm_dir> [--nproc N]
"""

# Levels whose global-mean height is below this (meters) are dropped.
DEFAULT_MIN_HEIGHT = 80000.0

import os
import sys
import glob
import argparse
import numpy as np

# Must import hdf5plugin before h5py to register zstd codec
import hdf5plugin
import h5py
from netCDF4 import Dataset
from datetime import datetime, timedelta

# Physical constants
G       = 9.80665       # m/s^2
R_DRY   = 287.058       # J/(kg*K)  specific gas constant for dry air
NA      = 6.02214076e23 # Avogadro's number (mol^-1)

# Molecular masses (kg/mol)
M_O    = 0.015999       # atomic oxygen
M_O2   = 0.031998       # molecular oxygen
M_N2   = 0.028014       # molecular nitrogen
M_H2O  = 0.018015       # water
M_O3   = 0.047997       # ozone


def read_fv3_file(filepath):
    """Read one FV3WAM atmf file (h5py handles the zstd compression)."""
    f = h5py.File(filepath, 'r')

    # Parse ISO time string
    time_bytes = f['time_iso'][0]  # array of single bytes
    time_str = b''.join(time_bytes).decode().strip('\x00').strip()
    file_dt = datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ')

    data = {
        'datetime':  file_dt,
        'grid_xt':   f['grid_xt'][:],            # (384,)  longitude 0..359.06
        'grid_yt':   f['grid_yt'][:],            # (190,)  latitude N->S
        'pfull':     f['pfull'][:] * 100.0,       # (196,)  convert mb -> Pa
        'hgtsfc':    f['hgtsfc'][0, :, :],       # (190, 384)  surface height m
        'delz':      f['delz'][0, :, :, :],      # (196, 190, 384) layer thickness (negative)
        'tmp':       f['tmp'][0, :, :, :],        # (196, 190, 384)
        'ugrd':      f['ugrd'][0, :, :, :],
        'vgrd':      f['vgrd'][0, :, :, :],
        'dzdt':      f['dzdt'][0, :, :, :],       # Pa/s
        'spfo':      f['spfo'][0, :, :, :],       # O  mass mixing ratio
        'spfo2':     f['spfo2'][0, :, :, :],      # O2 mass mixing ratio
        'spfh':      f['spfh'][0, :, :, :],       # specific humidity
        'o3mr':      f['o3mr'][0, :, :, :],       # ozone mixing ratio
    }
    f.close()
    return data


def compute_heights(delz, hgtsfc):
    """Mid-layer geometric heights (m); level 0 = model top."""
    nlev, nlat, nlon = delz.shape

    # Interface heights: nlev+1 interfaces, index nlev = surface
    z_ifc = np.empty((nlev + 1, nlat, nlon), dtype=np.float64)
    z_ifc[nlev, :, :] = hgtsfc  # bottom interface = surface

    # Accumulate upward from surface (k=nlev-1 is bottom layer, k=0 is top layer)
    for k in range(nlev - 1, -1, -1):
        z_ifc[k, :, :] = z_ifc[k + 1, :, :] + np.abs(delz[k, :, :])

    # Mid-layer heights
    z_mid = 0.5 * (z_ifc[:-1, :, :] + z_ifc[1:, :, :])

    return z_mid.astype(np.float32)


def compute_densities(pfull, tmp, spfo, spfo2, spfh, o3mr):
    """Mass mixing ratios -> number densities (m^-3); returns n_O, n_O2, n_N2."""
    # Broadcast pressure to 3D
    p3d = pfull[:, np.newaxis, np.newaxis]  # (nlev, 1, 1)

    # N2 as residual: 1 - O - O2 - H2O - O3
    mmr_n2 = np.clip(1.0 - spfo - spfo2 - spfh - o3mr, 0.0, None)

    # Compute mean molecular weight from composition (varies with altitude)
    # M_mean = 1 / sum(w_i / M_i) where w_i are mass mixing ratios
    sum_wi_Mi = spfo / M_O + spfo2 / M_O2 + mmr_n2 / M_N2 + spfh / M_H2O + o3mr / M_O3
    M_mean = 1.0 / sum_wi_Mi  # kg/mol

    # Air density from ideal gas law with composition-dependent gas constant
    # R_specific = R_universal / M_mean;  rho = p / (R_specific * T)
    R_UNIV = 8.314462  # J/(mol*K)
    rho_air = p3d * M_mean / (R_UNIV * tmp)  # kg/m^3

    # Number density = (MMR * rho_air / M) * NA
    n_O  = (spfo  * rho_air / M_O)  * NA
    n_O2 = (spfo2 * rho_air / M_O2) * NA
    n_N2 = (mmr_n2 * rho_air / M_N2) * NA

    return n_O.astype(np.float32), n_O2.astype(np.float32), n_N2.astype(np.float32)


def compute_vertical_wind(dzdt, pfull, tmp, spfo, spfo2, spfh, o3mr):
    """omega (dzdt, Pa/s) -> vertical wind (m/s), clipped to +/-200 m/s."""
    p3d = pfull[:, np.newaxis, np.newaxis]

    # Compute mean molecular weight from composition
    mmr_n2 = np.clip(1.0 - spfo - spfo2 - spfh - o3mr, 0.0, None)
    sum_wi_Mi = spfo / M_O + spfo2 / M_O2 + mmr_n2 / M_N2 + spfh / M_H2O + o3mr / M_O3
    M_mean = 1.0 / sum_wi_Mi  # kg/mol

    R_UNIV = 8.314462  # J/(mol*K)
    rho_air = np.maximum(p3d * M_mean / (R_UNIV * tmp), 1.0e-12)
    w = -dzdt / (rho_air * G)
    w = np.clip(w, -200.0, 200.0)  # physical limit for thermospheric winds
    return w.astype(np.float32)


def filter_levels(z_mid, min_height):
    """Indices of levels with mean height >= min_height, sorted ascending."""
    # Mean height per level
    zmean = np.mean(z_mid, axis=(1, 2))  # (nlev,)

    # Levels above threshold
    mask = zmean >= min_height

    # Level indices, sorted by ascending mean height
    valid_idx = np.where(mask)[0]
    # Sort by mean height ascending (FV3 is top-down, so reversed)
    valid_idx = valid_idx[np.argsort(zmean[valid_idx])]

    return valid_idx


def flip_lat(arr, lat_axis):
    """Flip array along latitude axis so latitude goes S->N."""
    return np.flip(arr, axis=lat_axis)


def write_gsm_netcdf(filepath, nlon, nlat, nlev, heights, temp,
                     n_O, n_O2, n_N2, u, v, w, grid_xt, grid_yt, pfull_subset):
    """Write one gsm-format NetCDF file; arrays are (nlon, nlat, nlev), lat S->N."""
    nc = Dataset(filepath, 'w', format='NETCDF4')

    # Dimensions (GSM convention: x01=lon, x02=lat, x03=levels)
    nc.createDimension('x01', nlon)
    nc.createDimension('x02', nlat)
    nc.createDimension('x03', nlev)

    def write_var(name, data, units=''):
        # GSM convention: variable dims are (x03, x02, x01) so that
        # Fortran (column-major) reads them as (nlon, nlat, nlev)
        v = nc.createVariable(name, 'f4', ('x03', 'x02', 'x01'), zlib=True)
        v[:] = np.transpose(data, (2, 1, 0))  # (nlon,nlat,nlev) -> (nlev,nlat,nlon)
        if units:
            v.units = units

    write_var('wam_height_levels', heights, 'm')
    write_var('temp_neutral', temp, 'K')
    write_var('O_Density', n_O, 'm-3')
    write_var('O2_Density', n_O2, 'm-3')
    write_var('N2_Density', n_N2, 'm-3')
    write_var('eastward_wind_neutral', u, 'm/s')
    write_var('northward_wind_neutral', v, 'm/s')
    write_var('upward_wind_neutral', w, 'm/s')

    # Store coordinate arrays as 1D variables for FileReader to use
    vlon = nc.createVariable('grid_xt', 'f8', ('x01',))
    vlon[:] = grid_xt
    vlon.units = 'degrees_east'

    vlat = nc.createVariable('grid_yt', 'f8', ('x02',))
    vlat[:] = grid_yt
    vlat.units = 'degrees_north'

    vpres = nc.createVariable('pfull', 'f4', ('x03',))
    vpres[:] = pfull_subset
    vpres.units = 'Pa'

    nc.source = 'fv3wam_to_gsm.py'
    nc.close()


def process_one_file(input_path, output_dir, min_height, valid_idx_ref=None):
    """Process a single FV3WAM file. Returns (output_path, valid_idx) or (None, None) on error."""
    print(f'  Reading {os.path.basename(input_path)}...', flush=True)
    data = read_fv3_file(input_path)

    # Compute heights
    z_mid = compute_heights(data['delz'], data['hgtsfc'])

    # Determine which levels to keep
    if valid_idx_ref is not None:
        valid_idx = valid_idx_ref
    else:
        valid_idx = filter_levels(z_mid, min_height)
        print(f'    Keeping {len(valid_idx)} of {z_mid.shape[0]} levels '
              f'(>= {min_height/1000:.0f} km)')

    nlev_out = len(valid_idx)
    if nlev_out == 0:
        print(f'    WARNING: No levels above {min_height/1000:.0f} km, skipping')
        return None, valid_idx

    # Subset and sort by ascending height
    z_sub   = z_mid[valid_idx]      # (nlev_out, nlat, nlon)
    tmp_sub = data['tmp'][valid_idx]
    u_sub   = data['ugrd'][valid_idx]
    v_sub   = data['vgrd'][valid_idx]
    pfull_sub = data['pfull'][valid_idx]

    # Compute derived fields on subsetted levels
    n_O, n_O2, n_N2 = compute_densities(
        data['pfull'][valid_idx], tmp_sub,
        data['spfo'][valid_idx], data['spfo2'][valid_idx],
        data['spfh'][valid_idx], data['o3mr'][valid_idx])

    w_sub = compute_vertical_wind(
        data['dzdt'][valid_idx], data['pfull'][valid_idx], tmp_sub,
        data['spfo'][valid_idx], data['spfo2'][valid_idx],
        data['spfh'][valid_idx], data['o3mr'][valid_idx])

    # Flip latitude from N->S to S->N (axis 1 is lat in (lev, lat, lon))
    z_sub   = flip_lat(z_sub, 1)
    tmp_sub = flip_lat(tmp_sub, 1)
    n_O     = flip_lat(n_O, 1)
    n_O2    = flip_lat(n_O2, 1)
    n_N2    = flip_lat(n_N2, 1)
    u_sub   = flip_lat(u_sub, 1)
    v_sub   = flip_lat(v_sub, 1)
    w_sub   = flip_lat(w_sub, 1)

    # Transpose from (lev, lat, lon) to (lon, lat, lev) for GSM format
    z_out   = np.transpose(z_sub, (2, 1, 0))
    tmp_out = np.transpose(tmp_sub, (2, 1, 0))
    nO_out  = np.transpose(n_O, (2, 1, 0))
    nO2_out = np.transpose(n_O2, (2, 1, 0))
    nN2_out = np.transpose(n_N2, (2, 1, 0))
    u_out   = np.transpose(u_sub, (2, 1, 0))
    v_out   = np.transpose(v_sub, (2, 1, 0))
    w_out   = np.transpose(w_sub, (2, 1, 0))

    nlon = z_out.shape[0]
    nlat = z_out.shape[1]

    # Latitude array flipped to S->N
    grid_yt_sn = np.flip(data['grid_yt'])

    # Output filename: gsm.YYYYMMDD_HH0000.nc
    dt = data['datetime']
    outname = f'gsm.{dt.strftime("%Y%m%d_%H")}0000.nc'
    outpath = os.path.join(output_dir, outname)

    write_gsm_netcdf(outpath, nlon, nlat, nlev_out,
                     z_out, tmp_out, nO_out, nO2_out, nN2_out,
                     u_out, v_out, w_out,
                     data['grid_xt'], grid_yt_sn, pfull_sub)

    print(f'    -> {outname}  ({dt.isoformat()})  '
          f'nlon={nlon} nlat={nlat} nlev={nlev_out}')
    return outpath, valid_idx


def output_name_for(filepath):
    """Read only the timestamp and return the gsm.* output filename."""
    f = h5py.File(filepath, 'r')
    time_str = b''.join(f['time_iso'][0]).decode().strip('\x00').strip()
    f.close()
    dt = datetime.strptime(time_str, '%Y-%m-%dT%H:%M:%SZ')
    return f'gsm.{dt.strftime("%Y%m%d_%H")}0000.nc'


# worker state (set by the pool initializer so valid_idx is not re-pickled per task)
_VALID_IDX = None
_OUTDIR = None
_MINH = None

def _worker_init(valid_idx, outdir, min_height):
    global _VALID_IDX, _OUTDIR, _MINH
    _VALID_IDX, _OUTDIR, _MINH = valid_idx, outdir, min_height

def _convert(fpath):
    try:
        outpath, _ = process_one_file(fpath, _OUTDIR, _MINH, _VALID_IDX)
        return ('ok' if outpath else 'nolev', fpath)
    except Exception as e:
        return (f'ERR {type(e).__name__}: {e}', fpath)


def main():
    ap = argparse.ArgumentParser(
        description='Convert FV3WAM atmf*.nc history files to the gsm-format '
                    'NetCDF files read by the IPE neutral file reader.')
    ap.add_argument('-i', '--input', required=True,
                    help='directory containing FV3WAM atmf*.nc files')
    ap.add_argument('-o', '--output', required=True,
                    help='output directory for gsm.*.nc files')
    ap.add_argument('--pattern', default='atmf*.nc',
                    help='glob for input files (default: atmf*.nc)')
    ap.add_argument('--min-height', type=float, default=DEFAULT_MIN_HEIGHT,
                    help='drop levels below this height in meters (default: 80000)')
    ap.add_argument('--ref', default=None,
                    help='file that fixes the level selection for all outputs '
                         '(default: first file), keeping nlev uniform')
    ap.add_argument('--nproc', type=int, default=1, help='parallel workers (default: 1)')
    ap.add_argument('--force', action='store_true',
                    help='reconvert even if the output file already exists')
    args = ap.parse_args()
    os.makedirs(args.output, exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.input, args.pattern)))
    if not files:
        raise SystemExit(f'no files matching {args.pattern} in {args.input}')

    # Fix the level selection once so every output has the same number of levels.
    ref = args.ref or files[0]
    d = read_fv3_file(ref)
    valid_idx = filter_levels(compute_heights(d['delz'], d['hgtsfc']), args.min_height)
    print(f'Reference {os.path.basename(ref)}: keeping {len(valid_idx)} levels '
          f'(>= {args.min_height/1000:.0f} km)')

    # Skip files already converted (unless --force).
    todo = []
    for fp in files:
        if not args.force:
            try:
                if os.path.exists(os.path.join(args.output, output_name_for(fp))):
                    continue
            except Exception:
                pass
        todo.append(fp)
    print(f'{len(todo)} of {len(files)} files to convert; nproc={args.nproc}')

    n_ok = n_err = 0
    if args.nproc > 1:
        from multiprocessing import Pool
        with Pool(args.nproc, initializer=_worker_init,
                  initargs=(valid_idx, args.output, args.min_height)) as pool:
            for i, (status, fp) in enumerate(pool.imap_unordered(_convert, todo, chunksize=4), 1):
                if status == 'ok':
                    n_ok += 1
                elif status.startswith('ERR'):
                    n_err += 1
                    print(f'  {os.path.basename(fp)}: {status}')
                if i % 200 == 0:
                    print(f'  ...{i}/{len(todo)}  ok={n_ok} err={n_err}')
    else:
        _worker_init(valid_idx, args.output, args.min_height)
        for fp in todo:
            status, _ = _convert(fp)
            if status == 'ok':
                n_ok += 1
            elif status.startswith('ERR'):
                n_err += 1
                print(f'  {os.path.basename(fp)}: {status}')

    print(f'Done: {n_ok} converted, {n_err} errors -> {args.output}')


if __name__ == '__main__':
    main()
