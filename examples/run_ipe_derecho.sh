#!/bin/bash
#PBS -N ipe_standalone
#PBS -A YOUR_ACCOUNT
#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=128:mpiprocs=40:mem=230GB
#PBS -q main
#PBS -j oe

# Example run script for Derecho (PBS + spack-stack): standalone IPE with the
# neutral file reader. Adjust the account, paths, and dates for your setup.

# --- User-configurable settings ---
CDATE=2021022200
FHMAX=8                           # Forecast hours
DELTIM_IPE=360                    # IPE timestep (seconds) - 2x default for speed
IPEFREQ=3600                      # Output frequency (seconds)
MSIS_TIME_STEP=900                # MSIS/HWM call frequency
NPROCIPE=40                       # All MPI ranks for IPE

# Solar/geomagnetic forcing (real-time indices to match coupled run)
INPUT_PARAMETERS=realtime
DCOM_PATH=/glade/work/akubaryk/dcom
FIX_F107=120.0
FIX_KP=3.0

# GSM neutral file reader settings
READ_GSM_NEUTRALS=T
GSM_NEUTRALS_DIR='/glade/derecho/scratch/garimam/ptmp/WAM62_coupled_2021_sparta_realsolar_ipeop_fv3test/'
NEUTRAL_INTERP_METHOD='linear'    # 'linear' or 'msis_ratio'

# --- Paths (adjust for your setup) ---
IPEDIR=${IPEDIR:-$HOME/IPE}                    # this repo
IPEEXEC=${IPEEXEC:-$IPEDIR/install/bin/ipe.x}
IPE_IC_DIR=${IPE_IC_DIR:?grid + fix files: IPE_Grid.h5, wei96*, *.dat, ...}
BASEDIR=${BASEDIR:?GSMWAM-IPE tree (modulefiles + input-parameter scripts)}
PARAMETER_PATH=${PARAMETER_PATH:-/glade/work/akubaryk/noscrub/WAM-IPE_INPUT_PARAMETERS}
CONDA_PYTHON=${CONDA_PYTHON:-python}
POSTSCRIPT=$IPEDIR/scripts/post.py

RUNDIR=${RUNDIR:-/glade/derecho/scratch/$USER/ipe_run}
ROTDIR=$RUNDIR

echo "============================================================"
echo " Standalone IPE with GSM Reader"
echo " CDATE: ${CDATE}  FHMAX: ${FHMAX}h"
echo " GSM files: ${GSM_NEUTRALS_DIR}"
echo " Output: ${RUNDIR}"
echo "============================================================"

# --- Load modules (must use ncarenv/23.09 to match compiled libraries) ---
module --force purge
module load ncarenv/23.09
module use $BASEDIR/modulefiles/derecho.intel
module load wam-ipe
module list

# --- Memory/stack settings ---
ulimit -s unlimited
export KMP_STACKSIZE=512m

# --- Derived variables ---
INI_YEAR=${CDATE:0:4}
INI_MONTH=${CDATE:4:2}
INI_DAY=${CDATE:6:2}
INI_HOUR=${CDATE:8:2}
CIPEDATE="${CDATE}00"
START_UT_SEC=$((10#$INI_HOUR * 3600))
END_TIME=$((START_UT_SEC + FHMAX * 3600))

# F10.7/Kp parameters (matching standalone IPE)
F107_KP_INTERVAL=60
F107_KP_SIZE=4800
F107_KP_SKIP_SIZE=$((36*3600/$F107_KP_INTERVAL))
F107_KP_DATA_SIZE=56
F107_KP_READ_IN_START=0
F107_KP_REALTIME_INTERVAL=60

# Operational parameters
COLFAC=1.3
OFFSET1_DEG=5.0
OFFSET2_DEG=20.0
POTENTIAL_MODEL=2
HPEQ=0.0
TRANSPORT_HIGHLAT_LP=30
PERP_TRANSPORT_MAX_LP=151
VERTICAL_WIND_LIMIT=100.0

# --- Setup run directory ---
# Clean only model output, preserve plots and other files
mkdir -p $RUNDIR $RUNDIR/netcdf $RUNDIR/output
rm -f $RUNDIR/output/IPE_State.apex.*.h5
rm -f $RUNDIR/netcdf/*.nc4
rm -f $RUNDIR/run_output.log
rm -f $RUNDIR/IPE.inp $RUNDIR/ipe.x $RUNDIR/wam_input_f107_kp.txt
cd $RUNDIR

# Log output
exec > >(tee -a $RUNDIR/run_output.log) 2>&1

# Copy executable
if [ ! -f "$IPEEXEC" ]; then
    echo "ERROR: ipe.x not found at $IPEEXEC"
    echo "Run build_standalone.sh first."
    exit 1
fi
cp $IPEEXEC $RUNDIR/

# Link IPE fix files
ln -sf ${IPE_IC_DIR}/IPE_Grid.h5 .
ln -sf ${IPE_IC_DIR}/wei96* .
ln -sf ${IPE_IC_DIR}/*.dat .
ln -sf ${IPE_IC_DIR}/*.bin .
ln -sf ${IPE_IC_DIR}/ionprof .
ln -sf ${IPE_IC_DIR}/tiros_spectra .

# Link netcdf input files
for f in ${IPE_IC_DIR}/*.nc; do
    [ -f "$f" ] && ln -sf $f .
done

# real file (not a symlink) so IPE can overwrite it
# Prefer coupled run's IPE state (Feb 22 00Z), then IPE_FIX
COUPLED_IC=${GSM_NEUTRALS_DIR}/IPE_State.apex.${CIPEDATE}.h5
FIX_IC=${IPE_IC_DIR}/IPE_State.apex.${CIPEDATE}.h5
JAN_IC=${IPE_IC_DIR}/IPE_State.apex.202101100000.h5

IC_FOUND=""
for candidate in "$COUPLED_IC" "$FIX_IC" "$JAN_IC"; do
    if [ -f "$candidate" ] && [ $(stat -c%s "$candidate") -gt 1000000 ]; then
        cp "$candidate" IPE_State.apex.${CIPEDATE}.h5
        echo "IC: $candidate (copied as ${CIPEDATE})"
        IC_FOUND=1
        break
    fi
done

if [ -z "$IC_FOUND" ]; then
    echo "ERROR: No valid IC found"
    exit 1
fi

# Setup F10.7/Kp input
SCRIPTSDIR=$BASEDIR/scripts
START_36H=$($(which ndate) -36 ${CDATE})
if [ "$INPUT_PARAMETERS" = "realtime" ] && [ -d "$DCOM_PATH" ]; then
    echo "Using real-time F10.7/Kp from $DCOM_PATH"
    python $SCRIPTSDIR/interpolate_input_parameters/parse_realtime.py \
        -s ${START_36H}00 \
        -d $(((36 + FHMAX) * 60)) \
        -p $DCOM_PATH 2>&1 || echo "WARNING: parse_realtime.py failed"
    # Update sizes based on generated file
    # Compute skip from actual file start time to CDATE
    LEN_F107=$(wc -l < wam_input_f107_kp.txt)
    F107_KP_SIZE=$((LEN_F107 - 5))
    F107_KP_DATA_SIZE=$F107_KP_SIZE
    F107_KP_INTERVAL=60
    # Extract first data timestamp and compute minutes to CDATE
    FILE_START=$(sed -n '6p' wam_input_f107_kp.txt | awk '{print $1}')
    SKIP_MINS=$($CONDA_PYTHON -c "
from datetime import datetime
t0 = datetime.strptime('${FILE_START}', '%Y-%m-%dT%H:%M:%SZ')
t1 = datetime.strptime('${CDATE}', '%Y%m%d%H')
print(int((t1-t0).total_seconds()/60))
")
    F107_KP_SKIP_SIZE=${SKIP_MINS}
    F107_KP_READ_IN_START=0
    # Validate: need SKIP + FHMAX*60 minutes of data after file start
    NEEDED=$((SKIP_MINS + FHMAX * 60))
    if [ $F107_KP_SIZE -lt $NEEDED ]; then
        echo "ERROR: F10.7/Kp file too short: have $F107_KP_SIZE minutes, need $NEEDED (SKIP=$SKIP_MINS + forecast=$((FHMAX*60)))"
        echo "  File covers: $FILE_START to $(tail -1 wam_input_f107_kp.txt | awk '{print $1}')"
        exit 1
    fi
    echo "F10.7/Kp: $LEN_F107 lines, SIZE=$F107_KP_SIZE, SKIP=$F107_KP_SKIP_SIZE, post-CDATE=$((F107_KP_SIZE - SKIP_MINS)) min (need $((FHMAX*60)))"
else
    echo "Using fixed F10.7=${FIX_F107}, Kp=${FIX_KP}"
    cat > temp_fix << EOF
${FIX_F107}
${FIX_KP}
EOF
    if [ -d "$PARAMETER_PATH" ]; then
        python $SCRIPTSDIR/interpolate_input_parameters/interpolate_input_parameters.py \
            -d $((36 + FHMAX)) \
            -s ${START_36H} \
            -p $PARAMETER_PATH \
            -m $INPUT_PARAMETERS \
            -f temp_fix 2>&1 || echo "WARNING: F10.7/Kp interpolation failed"
    fi
fi

# --- Create IPE.inp namelist ---
cat > IPE.inp << EOF
&SPACEMANAGEMENT
 GRID_FILE = 'IPE_Grid.h5',
/
&TIMESTEPPING
 TIME_STEP   = ${DELTIM_IPE}D0,
 START_TIME  = ${START_UT_SEC}.0D0,
 END_TIME    = ${END_TIME}.0D0,
 MSIS_TIME_STEP   = ${MSIS_TIME_STEP}.D0,
 INITIAL_TIMESTAMP = ${CIPEDATE}
/
&FORCING
 SOLAR_FORCING_TIME_STEP = ${F107_KP_INTERVAL}.0D0,
 F107_KP_SIZE            = ${F107_KP_SIZE},
 F107_KP_INTERVAL        = ${F107_KP_INTERVAL},
 F107_KP_SKIP_SIZE       = ${F107_KP_SKIP_SIZE},
 F107_KP_DATA_SIZE       = ${F107_KP_DATA_SIZE},
 F107_KP_READ_IN_START   = ${F107_KP_READ_IN_START},
 F107_KP_FILE            = './wam_input_f107_kp.txt',
 F107_KP_REALTIME_INTERVAL = ${F107_KP_REALTIME_INTERVAL}
/
&FILEIO
 READ_APEX_NEUTRALS        = F,
 WRITE_APEX_NEUTRALS       = T,
 WRITE_GEOGRAPHIC_NEUTRALS = T,
 FILE_OUTPUT_FREQUENCY     = ${IPEFREQ}.0D0,
 FILE_PREFIX               = 'output/IPE_State.apex.'
/
&ipecap
  mesh_height_min = 0.
  mesh_height_max = 2000.
  mesh_write      = 0
  mesh_write_file = 'ipemesh'
  mesh_fill       = 1
/
&ELDYN
  DYNAMO_EFIELD          = T
/
&OPERATIONAL
  COLFAC                 = ${COLFAC}
  OFFSET1_DEG            = ${OFFSET1_DEG}
  OFFSET2_DEG            = ${OFFSET2_DEG}
  POTENTIAL_MODEL        = ${POTENTIAL_MODEL}
  HPEQ                   = ${HPEQ}
  TRANSPORT_HIGHLAT_LP   = ${TRANSPORT_HIGHLAT_LP}
  PERP_TRANSPORT_MAX_LP  = ${PERP_TRANSPORT_MAX_LP}
  VERTICAL_WIND_LIMIT    = ${VERTICAL_WIND_LIMIT}
/
&NeutralFileIO
  READ_GSM_NEUTRALS     = ${READ_GSM_NEUTRALS},
  GSM_NEUTRALS_DIR      = '${GSM_NEUTRALS_DIR}',
  NEUTRAL_INTERP_METHOD = '${NEUTRAL_INTERP_METHOD}',
  VERBOSE_DIAG_LOCAL    = F
/
EOF

echo "=== IPE.inp ==="
cat IPE.inp
echo "=== Run directory contents ==="
ls -la
echo "=== Starting standalone IPE with GSM Reader ==="

# --- Run ---
mpiexec -n $NPROCIPE ./ipe.x 2>&1
IPE_EXIT=$?

echo "=== IPE finished with exit code $IPE_EXIT ==="

if [ $IPE_EXIT -ne 0 ] && [ $IPE_EXIT -ne 134 ]; then
    echo "ERROR: IPE failed with exit code $IPE_EXIT"
    exit $IPE_EXIT
fi

# --- Post-process: HDF5 -> NetCDF ---
echo "=== Post-processing (HDF5 -> NetCDF) ==="
$CONDA_PYTHON $POSTSCRIPT \
    -g $RUNDIR/IPE_Grid.h5 \
    -i $RUNDIR/output \
    -o $RUNDIR/netcdf 2>&1

echo "=== Done ==="
echo "HDF5 output: $RUNDIR/output/IPE_State.apex.*.h5"
echo "NetCDF output: $RUNDIR/netcdf/"
ls $RUNDIR/netcdf/*.nc4 2>/dev/null | wc -l
echo "netcdf files produced"
