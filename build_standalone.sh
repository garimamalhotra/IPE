#!/bin/bash
#===========================================================================
# build_standalone.sh -- build standalone IPE (no coupling), with the
# neutral file reader. Netcdf link flags are derived portably from nf-config,
# so there are no hardcoded library paths.
#
# Requirements in your environment (module load / spack / etc.):
#   - an MPI Fortran compiler (mpif90) and mpicc
#   - COMIO            (set COMIO_PREFIX, or have comio-config in PATH)
#   - ESMF             (set ESMFMKFILE to esmf.mk)
#   - netCDF-Fortran   (nf-config in PATH)  and netCDF-C, HDF5, (p)netCDF
#
# Usage:
#   export ESMFMKFILE=/path/to/esmf.mk
#   export COMIO_PREFIX=/path/to/comio/install      # or ensure comio-config is in PATH
#   ./build_standalone.sh
#
# Binary is installed to ./install/bin/ipe.x
#===========================================================================
set -e
cd "$(dirname "$0")"

# --- sanity checks ---
command -v mpif90    >/dev/null || { echo "ERROR: mpif90 not in PATH"; exit 1; }
command -v nf-config >/dev/null || { echo "ERROR: nf-config not in PATH (netCDF-Fortran)"; exit 1; }
: "${ESMFMKFILE:?ERROR: set ESMFMKFILE to your esmf.mk}"
COMIO_ARG=""
if [ -n "${COMIO_PREFIX:-}" ]; then
  COMIO_ARG="--with-comio=${COMIO_PREFIX}"
elif command -v comio-config >/dev/null; then
  COMIO_ARG="--with-comio=$(dirname "$(dirname "$(command -v comio-config)")")"
else
  echo "ERROR: set COMIO_PREFIX or put comio-config in PATH"; exit 1
fi

PREFIX=$PWD/install
FCFLAGS_OPT="-O2 -fp-model precise -ftz -fast-transcendentals -no-prec-div -no-prec-sqrt -align array64byte -align sequence -march=core-avx2"

# --- generate configure if needed ---
[ -x ./configure ] || autoreconf -i

# --- configure (standalone: coupling disabled) ---
./configure --disable-coupling \
  --prefix="$PREFIX" --datarootdir="$PREFIX" --libdir="$PREFIX" \
  $COMIO_ARG \
  FC=mpif90 F77=mpif90 CC=mpicc \
  "FCFLAGS=$FCFLAGS_OPT"

# --- add netCDF to the link line (the reader uses nf90_*; the base build only
#     links comio/hdf5/pnetcdf). netCDF-Fortran and netCDF-C often live in
#     separate prefixes, so add both library dirs (fortran from nf-config,
#     C from pkg-config / nc-config) — no hardcoded paths. ---
NF_LIBDIR=$(nf-config --prefix)/lib
NC_LIBDIR=$(pkg-config --variable=libdir netcdf 2>/dev/null)
[ -z "$NC_LIBDIR" ] && NC_LIBDIR=$(nc-config --libdir 2>/dev/null)
: "${NC_LIBDIR:?ERROR: could not locate the netCDF-C library dir (need pkg-config or nc-config)}"
sed -i "s|-lmpi_intel -lpnetcdf|-lmpi_intel -lmpifort_intel -L${NF_LIBDIR} -L${NC_LIBDIR} -lnetcdff -lnetcdf -lpnetcdf|" src/Makefile

# --- build ipelib first (dependency ordering for the dynamo/ subdir), then all ---
( cd src/ipelib && make -j1 libipe_a-IPE_Precision.o libipe_a-IPE_Constants_Dictionary.o && make -j4 )
make -j4
make install

echo
echo "Done: $PREFIX/bin/ipe.x"
