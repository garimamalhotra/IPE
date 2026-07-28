MODULE IPE_Neutrals_FileReader
!
! Direct GSM/NetCDF neutral file reader for IPE.
! Reads hourly GSM-format neutral atmosphere files and interpolates
! vertically and temporally onto IPE's geographic grid or directly
! onto the apex grid (for high-resolution source data).
!
! GSM grid (90x91) matches IPE geographic grid exactly.
! For larger source grids (e.g. FV3WAM 384x190), data is interpolated
! directly to apex grid points, bypassing the 90x91 geographic grid
! bottleneck and preserving the full source resolution.
!
  USE IPE_Precision
  USE IPE_Constants_Dictionary
  USE IPE_Time_Class
  USE IPE_MPI_Layer_Class
  USE ipe_error_module
  USE netcdf
  USE physics_msis        ! gtd7
  USE utils_constants, ONLY : msis_dp => dp
  USE IPE_Forcing_Class
  USE IPE_Grid_Class

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: NFIELDS = 7
  INTEGER, PARAMETER :: FLD_TEMP   = 1
  INTEGER, PARAMETER :: FLD_O      = 2
  INTEGER, PARAMETER :: FLD_O2     = 3
  INTEGER, PARAMETER :: FLD_N2     = 4
  INTEGER, PARAMETER :: FLD_UEAST  = 5
  INTEGER, PARAMETER :: FLD_UNORTH = 6
  INTEGER, PARAMETER :: FLD_UUP    = 7

  ! MSIS field indices (for msis buffers)
  INTEGER, PARAMETER :: NFIELDS_MSIS = 4
  INTEGER, PARAMETER :: MFLD_TEMP = 1, MFLD_O = 2, MFLD_O2 = 3, MFLD_N2 = 4

  TYPE, PUBLIC :: IPE_FileReader

    LOGICAL :: initialized = .FALSE.
    CHARACTER(512) :: file_dir

    ! Source grid dimensions
    INTEGER :: nlon_src, nlat_src, nlev_src

    ! Two time-level buffers on IPE geographic grid (legacy path)
    ! Shape: (NFIELDS, nlon_geo, nlat_geo, nheights_geo)
    REAL(prec), ALLOCATABLE :: geo_t0(:,:,:,:)
    REAL(prec), ALLOCATABLE :: geo_t1(:,:,:,:)

    ! Time tags (hours from start of simulation day)
    INTEGER :: hour_t0 = -1
    INTEGER :: hour_t1 = -1

    ! Reference date for constructing file names
    INTEGER :: ref_year, ref_month, ref_day

    ! Interpolation method: 'linear' or 'msis_ratio'
    CHARACTER(32) :: interp_method = 'linear'

    ! Horizontal interpolation (for source grids != IPE geo grid)
    LOGICAL :: needs_horiz_interp = .FALSE.
    REAL(prec), ALLOCATABLE :: src_lons(:)    ! (nlon_src) source longitudes, deg
    REAL(prec), ALLOCATABLE :: src_lats(:)    ! (nlat_src) source latitudes, deg (S->N)
    INTEGER, ALLOCATABLE    :: ilon_lo(:)     ! (nlon_geo) lower lon index in source
    REAL(prec), ALLOCATABLE :: wlon(:)        ! (nlon_geo) lon interp weight
    INTEGER, ALLOCATABLE    :: jlat_lo(:)     ! (nlat_geo) lower lat index in source
    REAL(prec), ALLOCATABLE :: wlat(:)        ! (nlat_geo) lat interp weight

    ! Direct-to-apex interpolation (bypasses 90x91 geo grid)
    LOGICAL :: direct_apex_interp = .FALSE.

    ! Per-apex-point precomputed horizontal interpolation weights in source grid
    INTEGER :: n_local_apex = 0
    INTEGER,    ALLOCATABLE :: apex_kp_arr(:)     ! (n_local_apex) kp index
    INTEGER,    ALLOCATABLE :: apex_lp_arr(:)     ! (n_local_apex) lp index
    INTEGER,    ALLOCATABLE :: apex_mp_arr(:)     ! (n_local_apex) mp index
    INTEGER,    ALLOCATABLE :: apex_ilon_lo(:)    ! (n_local_apex) lower lon in source
    INTEGER,    ALLOCATABLE :: apex_jlat_lo(:)    ! (n_local_apex) lower lat in source
    REAL(prec), ALLOCATABLE :: apex_wlon_arr(:)   ! (n_local_apex) lon weight
    REAL(prec), ALLOCATABLE :: apex_wlat_arr(:)   ! (n_local_apex) lat weight
    REAL(prec), ALLOCATABLE :: apex_alt_arr(:)    ! (n_local_apex) altitude in km

    ! Two time-level buffers on apex grid (direct path)
    ! Shape: (NFIELDS, nFluxTube, NLP, mp_low:mp_high)
    REAL(prec), ALLOCATABLE :: apex_t0(:,:,:,:)
    REAL(prec), ALLOCATABLE :: apex_t1(:,:,:,:)
    INTEGER :: apex_hour_t0 = -1
    INTEGER :: apex_hour_t1 = -1
    ! Grid dimensions needed for apex path
    INTEGER :: apex_nFluxTube = 0, apex_NLP = 0
    INTEGER :: apex_mp_low = 0, apex_mp_high = 0

    ! MSIS buffers for ratio-based interpolation (on geographic grid)
    ! Shape: (NFIELDS_MSIS, nlon_geo, nlat_geo, nheights_geo)
    REAL(prec), ALLOCATABLE :: msis_t0(:,:,:,:)
    REAL(prec), ALLOCATABLE :: msis_t1(:,:,:,:)
    INTEGER :: msis_hour_t0 = -1
    INTEGER :: msis_hour_t1 = -1

  CONTAINS
    PROCEDURE :: Init     => Init_FileReader
    PROCEDURE :: Finalize => Finalize_FileReader
    PROCEDURE :: Update   => Update_FileReader
    PROCEDURE :: InitApexWeights
    PROCEDURE :: UpdateApex => Update_FileReader_Apex
    PROCEDURE, PRIVATE :: ReadAndInterp
    PROCEDURE, PRIVATE :: ReadAndInterpToApex
    PROCEDURE, PRIVATE :: VertInterp_ToGeoGrid
    PROCEDURE, PRIVATE :: BuildGSMFileName
    PROCEDURE, PRIVATE :: Compute_MSIS_GeoGrid
    PROCEDURE, PRIVATE :: ComputeHorizWeights
    PROCEDURE, PRIVATE :: HorizInterp
  END TYPE IPE_FileReader

  PUBLIC :: NFIELDS, FLD_TEMP, FLD_O, FLD_O2, FLD_N2
  PUBLIC :: FLD_UEAST, FLD_UNORTH, FLD_UUP

CONTAINS

  !---------------------------------------------------------------------------
  ! Init_FileReader: Initialize the file reader, probe first file for dims
  !---------------------------------------------------------------------------
  SUBROUTINE Init_FileReader( reader, file_dir, ref_year, ref_month, ref_day, &
                              mpi_layer, interp_method, rc )

    CLASS(IPE_FileReader), INTENT(inout) :: reader
    CHARACTER(*),          INTENT(in)    :: file_dir
    INTEGER,               INTENT(in)    :: ref_year, ref_month, ref_day
    TYPE(IPE_MPI_Layer),   INTENT(in)    :: mpi_layer
    CHARACTER(*), OPTIONAL, INTENT(in)   :: interp_method
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    ! Local
    INTEGER :: localrc, ncid, dimid
    INTEGER :: nlon, nlat, nlev
    CHARACTER(512) :: filename

    IF (PRESENT(rc)) rc = IPE_SUCCESS

    reader % file_dir  = TRIM(file_dir)
    reader % ref_year  = ref_year
    reader % ref_month = ref_month
    reader % ref_day   = ref_day

    ! Probe the first GSM file for dimensions (rank 0 only)
    IF ( mpi_layer % rank_id == 0 ) THEN
      CALL reader % BuildGSMFileName( ref_year, ref_month, ref_day, 0, filename )
      WRITE(6,*) 'IPE_FileReader: probing ', TRIM(filename)

      localrc = nf90_open( TRIM(filename), NF90_NOWRITE, ncid )
      IF ( localrc /= NF90_NOERR ) THEN
        WRITE(6,*) 'IPE_FileReader ERROR: cannot open ', TRIM(filename)
        WRITE(6,*) '  ', TRIM(nf90_strerror(localrc))
        IF (PRESENT(rc)) rc = IPE_FAILURE
        RETURN
      ENDIF

      localrc = nf90_inq_dimid( ncid, 'x01', dimid )
      localrc = nf90_inquire_dimension( ncid, dimid, len=nlon )
      localrc = nf90_inq_dimid( ncid, 'x02', dimid )
      localrc = nf90_inquire_dimension( ncid, dimid, len=nlat )
      localrc = nf90_inq_dimid( ncid, 'x03', dimid )
      localrc = nf90_inquire_dimension( ncid, dimid, len=nlev )

      localrc = nf90_close( ncid )

      WRITE(6,*) 'IPE_FileReader: GSM grid = ', nlon, ' x ', nlat, ' x ', nlev
    ENDIF

    ! Broadcast dimensions
#ifdef HAVE_MPI
    CALL MPI_BCAST( nlon, 1, MPI_INTEGER, 0, mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( nlat, 1, MPI_INTEGER, 0, mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( nlev, 1, MPI_INTEGER, 0, mpi_layer % mpi_communicator, localrc )
#endif

    reader % nlon_src = nlon
    reader % nlat_src = nlat
    reader % nlev_src = nlev

    ! Check if source grid differs from IPE geographic grid
    IF ( nlon /= nlon_geo .OR. nlat /= nlat_geo ) THEN
      reader % needs_horiz_interp = .TRUE.
      IF ( mpi_layer % rank_id == 0 ) THEN
        WRITE(6,*) 'IPE_FileReader: source grid (', nlon, 'x', nlat, &
                   ') differs from IPE geo grid (', nlon_geo, 'x', nlat_geo, ')'
        WRITE(6,*) '  Will use direct-to-apex interpolation (after InitApexWeights)'
      ENDIF
      ! Read source grid coordinates (needed for apex weight computation)
      CALL reader % ComputeHorizWeights( mpi_layer )
    ENDIF

    ! Allocate geographic grid buffers (always needed for legacy/fallback path)
    ALLOCATE( reader % geo_t0(NFIELDS, nlon_geo, nlat_geo, nheights_geo), &
              reader % geo_t1(NFIELDS, nlon_geo, nlat_geo, nheights_geo) )
    reader % geo_t0 = 0.0_prec
    reader % geo_t1 = 0.0_prec
    reader % hour_t0 = -1
    reader % hour_t1 = -1

    reader % initialized = .TRUE.

    ! Store interpolation method
    IF ( PRESENT(interp_method) ) THEN
      reader % interp_method = TRIM(interp_method)
    ELSE
      reader % interp_method = 'linear'
    ENDIF

    ! Allocate MSIS buffers if using msis_ratio interpolation
    IF ( TRIM(reader % interp_method) == 'msis_ratio' ) THEN
      ALLOCATE( reader % msis_t0(NFIELDS_MSIS, nlon_geo, nlat_geo, nheights_geo), &
                reader % msis_t1(NFIELDS_MSIS, nlon_geo, nlat_geo, nheights_geo) )
      reader % msis_t0 = 0.0_prec
      reader % msis_t1 = 0.0_prec
      reader % msis_hour_t0 = -1
      reader % msis_hour_t1 = -1
      IF ( mpi_layer % rank_id == 0 ) THEN
        WRITE(6,*) 'IPE_FileReader: MSIS-ratio interpolation enabled'
      ENDIF
    ENDIF

    IF ( mpi_layer % rank_id == 0 ) THEN
      WRITE(6,*) 'IPE_FileReader: initialized successfully, method=', TRIM(reader % interp_method)
    ENDIF

  END SUBROUTINE Init_FileReader


  !---------------------------------------------------------------------------
  ! InitApexWeights: Precompute interpolation weights from source grid
  ! to each local apex grid point. Called from IPE_Neutrals_Class after
  ! the grid is available.
  !---------------------------------------------------------------------------
  SUBROUTINE InitApexWeights( reader, grid, mpi_layer, rc )

    CLASS(IPE_FileReader), INTENT(inout) :: reader
    TYPE(IPE_Grid),        INTENT(in)    :: grid
    TYPE(IPE_MPI_Layer),   INTENT(in)    :: mpi_layer
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    ! Local
    INTEGER :: kp, lp, mp, n, count
    REAL(prec) :: geo_lon, geo_lat, geo_alt
    REAL(prec) :: dlon, frac_lon
    INTEGER :: jlo
    REAL(prec) :: wlat_val

    IF (PRESENT(rc)) rc = IPE_SUCCESS

    ! Store apex grid dimensions
    reader % apex_nFluxTube = grid % nFluxTube
    reader % apex_NLP       = grid % NLP
    reader % apex_mp_low    = grid % mp_low
    reader % apex_mp_high   = grid % mp_high

    ! Count local apex points
    count = 0
    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        DO kp = 1, grid % flux_tube_max(lp)
          count = count + 1
        ENDDO
      ENDDO
    ENDDO
    reader % n_local_apex = count

    IF ( mpi_layer % rank_id == 0 ) THEN
      WRITE(6,*) 'IPE_FileReader: InitApexWeights: n_local_apex=', count
      WRITE(6,*) '  apex grid: nFluxTube=', grid % nFluxTube, &
                 ' NLP=', grid % NLP, ' NMP=', grid % NMP
      WRITE(6,*) '  local mp range:', grid % mp_low, '-', grid % mp_high
    ENDIF

    ! Allocate arrays
    ALLOCATE( reader % apex_kp_arr(count),   &
              reader % apex_lp_arr(count),   &
              reader % apex_mp_arr(count),   &
              reader % apex_ilon_lo(count),  &
              reader % apex_jlat_lo(count),  &
              reader % apex_wlon_arr(count), &
              reader % apex_wlat_arr(count), &
              reader % apex_alt_arr(count) )

    ! Source longitude spacing (assumed uniform)
    dlon = reader % src_lons(2) - reader % src_lons(1)

    ! Compute weights for each apex point
    n = 0
    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        DO kp = 1, grid % flux_tube_max(lp)
          n = n + 1

          reader % apex_kp_arr(n) = kp
          reader % apex_lp_arr(n) = lp
          reader % apex_mp_arr(n) = mp

          ! Geographic coordinates of this apex point
          geo_lon = rtd * grid % longitude(kp, lp, mp)       ! degrees (0-360)
          geo_lat = 90.0_prec - rtd * grid % colatitude(kp, lp, mp)  ! degrees (-90 to 90)
          geo_alt = m_to_km * grid % altitude(kp, lp)        ! km

          reader % apex_alt_arr(n) = geo_alt

          ! --- Longitude weight ---
          ! Ensure lon is in [0, 360)
          geo_lon = MODULO(geo_lon, 360.0_prec)
          frac_lon = geo_lon / dlon
          reader % apex_ilon_lo(n) = INT(frac_lon) + 1    ! 1-based
          IF ( reader % apex_ilon_lo(n) < 1 ) reader % apex_ilon_lo(n) = reader % nlon_src
          IF ( reader % apex_ilon_lo(n) > reader % nlon_src ) reader % apex_ilon_lo(n) = 1
          reader % apex_wlon_arr(n) = frac_lon - REAL(INT(frac_lon), prec)

          ! --- Latitude weight ---
          IF ( geo_lat <= reader % src_lats(1) ) THEN
            reader % apex_jlat_lo(n) = 1
            reader % apex_wlat_arr(n) = 0.0_prec
          ELSEIF ( geo_lat >= reader % src_lats(reader % nlat_src) ) THEN
            reader % apex_jlat_lo(n) = reader % nlat_src - 1
            reader % apex_wlat_arr(n) = 1.0_prec
          ELSE
            CALL BinarySearchLat( reader % src_lats, reader % nlat_src, &
                                  geo_lat, jlo, wlat_val )
            reader % apex_jlat_lo(n) = jlo
            reader % apex_wlat_arr(n) = wlat_val
          ENDIF

        ENDDO
      ENDDO
    ENDDO

    ! Allocate apex time-level buffers
    ALLOCATE( reader % apex_t0(NFIELDS, grid % nFluxTube, grid % NLP, &
                                grid % mp_low : grid % mp_high), &
              reader % apex_t1(NFIELDS, grid % nFluxTube, grid % NLP, &
                                grid % mp_low : grid % mp_high) )
    reader % apex_t0 = 0.0_prec
    reader % apex_t1 = 0.0_prec
    reader % apex_hour_t0 = -1
    reader % apex_hour_t1 = -1

    reader % direct_apex_interp = .TRUE.

    IF ( mpi_layer % rank_id == 0 ) THEN
      WRITE(6,*) 'IPE_FileReader: direct-to-apex interpolation enabled'
    ENDIF

  END SUBROUTINE InitApexWeights


  !---------------------------------------------------------------------------
  ! Finalize_FileReader: Clean up allocated memory
  !---------------------------------------------------------------------------
  SUBROUTINE Finalize_FileReader( reader )

    CLASS(IPE_FileReader), INTENT(inout) :: reader

    IF ( ALLOCATED( reader % geo_t0 ) ) DEALLOCATE( reader % geo_t0 )
    IF ( ALLOCATED( reader % geo_t1 ) ) DEALLOCATE( reader % geo_t1 )
    IF ( ALLOCATED( reader % msis_t0 ) ) DEALLOCATE( reader % msis_t0 )
    IF ( ALLOCATED( reader % msis_t1 ) ) DEALLOCATE( reader % msis_t1 )
    IF ( ALLOCATED( reader % src_lons ) ) DEALLOCATE( reader % src_lons )
    IF ( ALLOCATED( reader % src_lats ) ) DEALLOCATE( reader % src_lats )
    IF ( ALLOCATED( reader % ilon_lo ) ) DEALLOCATE( reader % ilon_lo )
    IF ( ALLOCATED( reader % wlon ) )    DEALLOCATE( reader % wlon )
    IF ( ALLOCATED( reader % jlat_lo ) ) DEALLOCATE( reader % jlat_lo )
    IF ( ALLOCATED( reader % wlat ) )    DEALLOCATE( reader % wlat )
    IF ( ALLOCATED( reader % apex_kp_arr ) )   DEALLOCATE( reader % apex_kp_arr )
    IF ( ALLOCATED( reader % apex_lp_arr ) )   DEALLOCATE( reader % apex_lp_arr )
    IF ( ALLOCATED( reader % apex_mp_arr ) )   DEALLOCATE( reader % apex_mp_arr )
    IF ( ALLOCATED( reader % apex_ilon_lo ) )  DEALLOCATE( reader % apex_ilon_lo )
    IF ( ALLOCATED( reader % apex_jlat_lo ) )  DEALLOCATE( reader % apex_jlat_lo )
    IF ( ALLOCATED( reader % apex_wlon_arr ) ) DEALLOCATE( reader % apex_wlon_arr )
    IF ( ALLOCATED( reader % apex_wlat_arr ) ) DEALLOCATE( reader % apex_wlat_arr )
    IF ( ALLOCATED( reader % apex_alt_arr ) )  DEALLOCATE( reader % apex_alt_arr )
    IF ( ALLOCATED( reader % apex_t0 ) ) DEALLOCATE( reader % apex_t0 )
    IF ( ALLOCATED( reader % apex_t1 ) ) DEALLOCATE( reader % apex_t1 )
    reader % initialized = .FALSE.

  END SUBROUTINE Finalize_FileReader


  !---------------------------------------------------------------------------
  ! Update_FileReader: Legacy entry point (geo grid output).
  ! Used when source grid matches IPE geo grid (90x91 GSM files).
  !---------------------------------------------------------------------------
  SUBROUTINE Update_FileReader( reader, time, mpi_layer, &
       geo_temperature, geo_oxygen, geo_mol_oxygen, geo_mol_nitrogen, &
       geo_velocity, altitude_geo, forcing, rc )

    CLASS(IPE_FileReader),  INTENT(inout) :: reader
    TYPE(IPE_Time),         INTENT(in)    :: time
    TYPE(IPE_MPI_Layer),    INTENT(in)    :: mpi_layer
    REAL(prec), INTENT(out) :: geo_temperature(:,:,:)     ! (nlon_geo, nlat_geo, nheights_geo)
    REAL(prec), INTENT(out) :: geo_oxygen(:,:,:)
    REAL(prec), INTENT(out) :: geo_mol_oxygen(:,:,:)
    REAL(prec), INTENT(out) :: geo_mol_nitrogen(:,:,:)
    REAL(prec), INTENT(out) :: geo_velocity(:,:,:,:)       ! (3, nlon_geo, nlat_geo, nheights_geo)
    REAL(prec), INTENT(in)  :: altitude_geo(:)             ! (nheights_geo) in km
    TYPE(IPE_Forcing), OPTIONAL, INTENT(in) :: forcing
    INTEGER, OPTIONAL,      INTENT(out) :: rc

    ! Local
    INTEGER :: hour_before, hour_after
    INTEGER :: localrc
    REAL(prec) :: weight
    REAL(prec) :: current_hour_frac
    INTEGER :: i, j, k
    REAL(prec) :: val_t0, val_t1
    LOGICAL :: use_msis_ratio

    ! MSIS-ratio locals
    REAL(prec), ALLOCATABLE :: msis_now(:,:,:,:)
    REAL(prec) :: ratio_t0, ratio_t1, ratio_interp
    REAL(prec) :: diff_t0, diff_t1, diff_interp
    REAL(prec) :: msis_val_t0, msis_val_t1, msis_val_now
    REAL(prec), PARAMETER :: RATIO_MIN = 0.01_prec
    REAL(prec), PARAMETER :: RATIO_MAX = 100.0_prec

    IF (PRESENT(rc)) rc = IPE_SUCCESS

    use_msis_ratio = ( TRIM(reader % interp_method) == 'msis_ratio' )

    ! Determine current fractional hour from elapsed time
    ! time % utime is seconds from midnight UTC
    current_hour_frac = REAL(time % utime, prec) / 3600.0_prec

    ! Bracketing hours
    hour_before = INT(current_hour_frac)
    hour_after  = hour_before + 1
    weight = current_hour_frac - REAL(hour_before, prec)

    ! At exact hour boundaries, no interpolation needed
    IF ( ABS(weight) < 1.0e-10_prec ) THEN
      hour_after = hour_before
      weight = 0.0_prec
    ENDIF

    ! Load t0 buffer if needed
    IF ( reader % hour_t0 /= hour_before ) THEN
      ! Check if t1 buffer already has what we need
      IF ( reader % hour_t1 == hour_before ) THEN
        ! Shift: t0 <- t1
        reader % geo_t0 = reader % geo_t1
        reader % hour_t0 = hour_before
        ! Also shift MSIS buffers
        IF ( use_msis_ratio ) THEN
          reader % msis_t0 = reader % msis_t1
          reader % msis_hour_t0 = reader % msis_hour_t1
        ENDIF
      ELSE
        ! Read new file for hour_before
        CALL reader % ReadAndInterp( hour_before, mpi_layer, altitude_geo, &
             reader % geo_t0, localrc )
        IF ( ipe_error_check( localrc, msg="FileReader: failed to read t0 file", rc=rc ) ) RETURN
        reader % hour_t0 = hour_before
        ! Compute MSIS at this hour boundary
        IF ( use_msis_ratio .AND. reader % msis_hour_t0 /= hour_before ) THEN
          CALL reader % Compute_MSIS_GeoGrid( &
               hour_before * 3600.0_prec, time % day_of_year, &
               forcing, altitude_geo, mpi_layer, reader % msis_t0, localrc )
          IF ( ipe_error_check( localrc, msg="FileReader: MSIS t0 failed", rc=rc ) ) RETURN
          reader % msis_hour_t0 = hour_before
        ENDIF
      ENDIF
    ENDIF

    ! Load t1 buffer if needed (and different from t0)
    IF ( hour_after /= hour_before .AND. reader % hour_t1 /= hour_after ) THEN
      CALL reader % ReadAndInterp( hour_after, mpi_layer, altitude_geo, &
           reader % geo_t1, localrc )
      IF ( ipe_error_check( localrc, msg="FileReader: failed to read t1 file", rc=rc ) ) RETURN
      reader % hour_t1 = hour_after
      ! Compute MSIS at this hour boundary
      IF ( use_msis_ratio .AND. reader % msis_hour_t1 /= hour_after ) THEN
        CALL reader % Compute_MSIS_GeoGrid( &
             hour_after * 3600.0_prec, time % day_of_year, &
             forcing, altitude_geo, mpi_layer, reader % msis_t1, localrc )
        IF ( ipe_error_check( localrc, msg="FileReader: MSIS t1 failed", rc=rc ) ) RETURN
        reader % msis_hour_t1 = hour_after
      ENDIF
    ENDIF

    ! --- Temporal interpolation: fill output arrays ---

    IF ( hour_after == hour_before .OR. ABS(weight) < 1.0e-10_prec ) THEN
      ! Exact hour boundary, no interpolation
      DO k = 1, nheights_geo
        DO j = 1, nlat_geo
          DO i = 1, nlon_geo
            geo_temperature(i,j,k)     = reader % geo_t0(FLD_TEMP,i,j,k)
            geo_oxygen(i,j,k)          = reader % geo_t0(FLD_O,i,j,k)
            geo_mol_oxygen(i,j,k)      = reader % geo_t0(FLD_O2,i,j,k)
            geo_mol_nitrogen(i,j,k)    = reader % geo_t0(FLD_N2,i,j,k)
            geo_velocity(1,i,j,k)      = reader % geo_t0(FLD_UEAST,i,j,k)
            geo_velocity(2,i,j,k)      = reader % geo_t0(FLD_UNORTH,i,j,k)
            geo_velocity(3,i,j,k)      = reader % geo_t0(FLD_UUP,i,j,k)
          ENDDO
        ENDDO
      ENDDO

    ELSEIF ( use_msis_ratio .AND. PRESENT(forcing) ) THEN
      ! MSIS-ratio interpolation
      ALLOCATE( msis_now(NFIELDS_MSIS, nlon_geo, nlat_geo, nheights_geo) )
      CALL reader % Compute_MSIS_GeoGrid( &
           REAL(time % utime, prec), time % day_of_year, &
           forcing, altitude_geo, mpi_layer, msis_now, localrc )
      IF ( ipe_error_check( localrc, msg="FileReader: MSIS now failed", rc=rc ) ) THEN
        DEALLOCATE( msis_now )
        RETURN
      ENDIF

      DO k = 1, nheights_geo
        DO j = 1, nlat_geo
          DO i = 1, nlon_geo
            ! Oxygen: ratio interpolation
            msis_val_t0  = MAX(reader % msis_t0(MFLD_O,i,j,k), 1.0e-30_prec)
            msis_val_t1  = MAX(reader % msis_t1(MFLD_O,i,j,k), 1.0e-30_prec)
            msis_val_now = MAX(msis_now(MFLD_O,i,j,k), 1.0e-30_prec)
            ratio_t0 = MIN(MAX(reader % geo_t0(FLD_O,i,j,k) / msis_val_t0, RATIO_MIN), RATIO_MAX)
            ratio_t1 = MIN(MAX(reader % geo_t1(FLD_O,i,j,k) / msis_val_t1, RATIO_MIN), RATIO_MAX)
            ratio_interp = EXP( (1.0_prec - weight)*LOG(ratio_t0) + weight*LOG(ratio_t1) )
            geo_oxygen(i,j,k) = msis_val_now * ratio_interp

            ! O2
            msis_val_t0  = MAX(reader % msis_t0(MFLD_O2,i,j,k), 1.0e-30_prec)
            msis_val_t1  = MAX(reader % msis_t1(MFLD_O2,i,j,k), 1.0e-30_prec)
            msis_val_now = MAX(msis_now(MFLD_O2,i,j,k), 1.0e-30_prec)
            ratio_t0 = MIN(MAX(reader % geo_t0(FLD_O2,i,j,k) / msis_val_t0, RATIO_MIN), RATIO_MAX)
            ratio_t1 = MIN(MAX(reader % geo_t1(FLD_O2,i,j,k) / msis_val_t1, RATIO_MIN), RATIO_MAX)
            ratio_interp = EXP( (1.0_prec - weight)*LOG(ratio_t0) + weight*LOG(ratio_t1) )
            geo_mol_oxygen(i,j,k) = msis_val_now * ratio_interp

            ! N2
            msis_val_t0  = MAX(reader % msis_t0(MFLD_N2,i,j,k), 1.0e-30_prec)
            msis_val_t1  = MAX(reader % msis_t1(MFLD_N2,i,j,k), 1.0e-30_prec)
            msis_val_now = MAX(msis_now(MFLD_N2,i,j,k), 1.0e-30_prec)
            ratio_t0 = MIN(MAX(reader % geo_t0(FLD_N2,i,j,k) / msis_val_t0, RATIO_MIN), RATIO_MAX)
            ratio_t1 = MIN(MAX(reader % geo_t1(FLD_N2,i,j,k) / msis_val_t1, RATIO_MIN), RATIO_MAX)
            ratio_interp = EXP( (1.0_prec - weight)*LOG(ratio_t0) + weight*LOG(ratio_t1) )
            geo_mol_nitrogen(i,j,k) = msis_val_now * ratio_interp

            ! Temperature: difference-based interpolation
            diff_t0 = reader % geo_t0(FLD_TEMP,i,j,k) - reader % msis_t0(MFLD_TEMP,i,j,k)
            diff_t1 = reader % geo_t1(FLD_TEMP,i,j,k) - reader % msis_t1(MFLD_TEMP,i,j,k)
            diff_interp = (1.0_prec - weight) * diff_t0 + weight * diff_t1
            geo_temperature(i,j,k) = msis_now(MFLD_TEMP,i,j,k) + diff_interp

            ! Winds: simple linear
            geo_velocity(1,i,j,k) = (1.0_prec - weight) * reader % geo_t0(FLD_UEAST,i,j,k) &
                                   + weight * reader % geo_t1(FLD_UEAST,i,j,k)
            geo_velocity(2,i,j,k) = (1.0_prec - weight) * reader % geo_t0(FLD_UNORTH,i,j,k) &
                                   + weight * reader % geo_t1(FLD_UNORTH,i,j,k)
            geo_velocity(3,i,j,k) = (1.0_prec - weight) * reader % geo_t0(FLD_UUP,i,j,k) &
                                   + weight * reader % geo_t1(FLD_UUP,i,j,k)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE( msis_now )

    ELSE
      ! Simple temporal interpolation
      DO k = 1, nheights_geo
        DO j = 1, nlat_geo
          DO i = 1, nlon_geo
            ! Temperature: linear
            geo_temperature(i,j,k) = (1.0_prec - weight) * reader % geo_t0(FLD_TEMP,i,j,k) &
                                    + weight * reader % geo_t1(FLD_TEMP,i,j,k)

            ! Densities: log-linear
            val_t0 = MAX(reader % geo_t0(FLD_O,i,j,k), 1.0e-30_prec)
            val_t1 = MAX(reader % geo_t1(FLD_O,i,j,k), 1.0e-30_prec)
            geo_oxygen(i,j,k) = EXP( (1.0_prec - weight)*LOG(val_t0) + weight*LOG(val_t1) )

            val_t0 = MAX(reader % geo_t0(FLD_O2,i,j,k), 1.0e-30_prec)
            val_t1 = MAX(reader % geo_t1(FLD_O2,i,j,k), 1.0e-30_prec)
            geo_mol_oxygen(i,j,k) = EXP( (1.0_prec - weight)*LOG(val_t0) + weight*LOG(val_t1) )

            val_t0 = MAX(reader % geo_t0(FLD_N2,i,j,k), 1.0e-30_prec)
            val_t1 = MAX(reader % geo_t1(FLD_N2,i,j,k), 1.0e-30_prec)
            geo_mol_nitrogen(i,j,k) = EXP( (1.0_prec - weight)*LOG(val_t0) + weight*LOG(val_t1) )

            ! Winds: linear
            geo_velocity(1,i,j,k) = (1.0_prec - weight) * reader % geo_t0(FLD_UEAST,i,j,k) &
                                   + weight * reader % geo_t1(FLD_UEAST,i,j,k)
            geo_velocity(2,i,j,k) = (1.0_prec - weight) * reader % geo_t0(FLD_UNORTH,i,j,k) &
                                   + weight * reader % geo_t1(FLD_UNORTH,i,j,k)
            geo_velocity(3,i,j,k) = (1.0_prec - weight) * reader % geo_t0(FLD_UUP,i,j,k) &
                                   + weight * reader % geo_t1(FLD_UUP,i,j,k)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE Update_FileReader


  !---------------------------------------------------------------------------
  ! Update_FileReader_Apex: Direct-to-apex entry point.
  ! Reads source files and interpolates directly to apex grid arrays,
  ! bypassing the 90x91 geographic grid.
  !---------------------------------------------------------------------------
  SUBROUTINE Update_FileReader_Apex( reader, time, mpi_layer, &
       apex_temperature, apex_oxygen, apex_mol_oxygen, apex_mol_nitrogen, &
       apex_velocity_geo, rc )

    CLASS(IPE_FileReader),  INTENT(inout) :: reader
    TYPE(IPE_Time),         INTENT(in)    :: time
    TYPE(IPE_MPI_Layer),    INTENT(in)    :: mpi_layer
    REAL(prec), INTENT(out) :: apex_temperature(:,:,:)   ! (nFluxTube, NLP, mp_low:mp_high)
    REAL(prec), INTENT(out) :: apex_oxygen(:,:,:)
    REAL(prec), INTENT(out) :: apex_mol_oxygen(:,:,:)
    REAL(prec), INTENT(out) :: apex_mol_nitrogen(:,:,:)
    REAL(prec), INTENT(out) :: apex_velocity_geo(:,:,:,:) ! (3, nFluxTube, NLP, mp_low:mp_high)
    INTEGER, OPTIONAL,      INTENT(out)   :: rc

    ! Local
    INTEGER :: hour_before, hour_after
    INTEGER :: localrc
    REAL(prec) :: weight, current_hour_frac
    INTEGER :: n, kp, lp, mp, mp_loc
    REAL(prec) :: val_t0, val_t1

    IF (PRESENT(rc)) rc = IPE_SUCCESS

    ! Determine current fractional hour
    current_hour_frac = REAL(time % utime, prec) / 3600.0_prec
    hour_before = INT(current_hour_frac)
    hour_after  = hour_before + 1
    weight = current_hour_frac - REAL(hour_before, prec)

    IF ( ABS(weight) < 1.0e-10_prec ) THEN
      hour_after = hour_before
      weight = 0.0_prec
    ENDIF

    ! Load apex_t0 buffer if needed
    IF ( reader % apex_hour_t0 /= hour_before ) THEN
      IF ( reader % apex_hour_t1 == hour_before ) THEN
        reader % apex_t0 = reader % apex_t1
        reader % apex_hour_t0 = hour_before
      ELSE
        CALL reader % ReadAndInterpToApex( hour_before, mpi_layer, &
             reader % apex_t0, localrc )
        IF ( ipe_error_check( localrc, msg="FileReader: failed to read apex t0", rc=rc ) ) RETURN
        reader % apex_hour_t0 = hour_before
      ENDIF
    ENDIF

    ! Load apex_t1 buffer if needed
    IF ( hour_after /= hour_before .AND. reader % apex_hour_t1 /= hour_after ) THEN
      CALL reader % ReadAndInterpToApex( hour_after, mpi_layer, &
           reader % apex_t1, localrc )
      IF ( ipe_error_check( localrc, msg="FileReader: failed to read apex t1", rc=rc ) ) RETURN
      reader % apex_hour_t1 = hour_after
    ENDIF

    ! --- Temporal interpolation to apex output arrays ---
    IF ( hour_after == hour_before .OR. ABS(weight) < 1.0e-10_prec ) THEN
      ! Exact hour boundary
      DO n = 1, reader % n_local_apex
        kp = reader % apex_kp_arr(n)
        lp = reader % apex_lp_arr(n)
        mp = reader % apex_mp_arr(n)
        ! Output arrays are assumed-shape (lower bound=1), apex_t0 is type member (retains mp_low:mp_high)
        mp_loc = mp - reader % apex_mp_low + 1
        apex_temperature(kp,lp,mp_loc)    = reader % apex_t0(FLD_TEMP,kp,lp,mp)
        apex_oxygen(kp,lp,mp_loc)         = reader % apex_t0(FLD_O,kp,lp,mp)
        apex_mol_oxygen(kp,lp,mp_loc)     = reader % apex_t0(FLD_O2,kp,lp,mp)
        apex_mol_nitrogen(kp,lp,mp_loc)   = reader % apex_t0(FLD_N2,kp,lp,mp)
        apex_velocity_geo(1,kp,lp,mp_loc) = reader % apex_t0(FLD_UEAST,kp,lp,mp)
        apex_velocity_geo(2,kp,lp,mp_loc) = reader % apex_t0(FLD_UNORTH,kp,lp,mp)
        apex_velocity_geo(3,kp,lp,mp_loc) = 0.0_prec  ! vertical wind zeroed (see ReadAndInterpToApex)
      ENDDO
    ELSE
      ! Linear temporal interpolation (log-linear for densities)
      DO n = 1, reader % n_local_apex
        kp = reader % apex_kp_arr(n)
        lp = reader % apex_lp_arr(n)
        mp = reader % apex_mp_arr(n)
        mp_loc = mp - reader % apex_mp_low + 1

        ! Temperature: linear
        apex_temperature(kp,lp,mp_loc) = (1.0_prec - weight) * reader % apex_t0(FLD_TEMP,kp,lp,mp) &
                                    + weight * reader % apex_t1(FLD_TEMP,kp,lp,mp)

        ! Densities: log-linear
        val_t0 = MAX(reader % apex_t0(FLD_O,kp,lp,mp), 1.0e-30_prec)
        val_t1 = MAX(reader % apex_t1(FLD_O,kp,lp,mp), 1.0e-30_prec)
        apex_oxygen(kp,lp,mp_loc) = EXP( (1.0_prec - weight)*LOG(val_t0) + weight*LOG(val_t1) )

        val_t0 = MAX(reader % apex_t0(FLD_O2,kp,lp,mp), 1.0e-30_prec)
        val_t1 = MAX(reader % apex_t1(FLD_O2,kp,lp,mp), 1.0e-30_prec)
        apex_mol_oxygen(kp,lp,mp_loc) = EXP( (1.0_prec - weight)*LOG(val_t0) + weight*LOG(val_t1) )

        val_t0 = MAX(reader % apex_t0(FLD_N2,kp,lp,mp), 1.0e-30_prec)
        val_t1 = MAX(reader % apex_t1(FLD_N2,kp,lp,mp), 1.0e-30_prec)
        apex_mol_nitrogen(kp,lp,mp_loc) = EXP( (1.0_prec - weight)*LOG(val_t0) + weight*LOG(val_t1) )

        ! Winds: linear (horizontal only; vertical wind zeroed)
        apex_velocity_geo(1,kp,lp,mp_loc) = (1.0_prec - weight) * reader % apex_t0(FLD_UEAST,kp,lp,mp) &
                                        + weight * reader % apex_t1(FLD_UEAST,kp,lp,mp)
        apex_velocity_geo(2,kp,lp,mp_loc) = (1.0_prec - weight) * reader % apex_t0(FLD_UNORTH,kp,lp,mp) &
                                        + weight * reader % apex_t1(FLD_UNORTH,kp,lp,mp)
        apex_velocity_geo(3,kp,lp,mp_loc) = 0.0_prec
      ENDDO
    ENDIF

  END SUBROUTINE Update_FileReader_Apex


  !---------------------------------------------------------------------------
  ! ReadAndInterp: Read a GSM file and interpolate to IPE geographic grid.
  ! Legacy path for 90x91 source grids.
  !---------------------------------------------------------------------------
  SUBROUTINE ReadAndInterp( reader, hour, mpi_layer, altitude_geo, geo_buf, rc )

    CLASS(IPE_FileReader), INTENT(in)    :: reader
    INTEGER,               INTENT(in)    :: hour
    TYPE(IPE_MPI_Layer),   INTENT(in)    :: mpi_layer
    REAL(prec),            INTENT(in)    :: altitude_geo(:)  ! (nheights_geo) in km
    REAL(prec),            INTENT(out)   :: geo_buf(:,:,:,:) ! (NFIELDS, nlon_geo, nlat_geo, nheights_geo)
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    ! Local
    INTEGER :: localrc
    INTEGER :: nlev, nlat, nlon
    CHARACTER(512) :: filename
    REAL(sp), ALLOCATABLE :: heights_src(:,:,:)  ! (nlon, nlat, nlev) float - NetCDF Fortran order
    REAL(sp), ALLOCATABLE :: field_src(:,:,:)    ! (nlon, nlat, nlev) float
    REAL(sp), ALLOCATABLE :: temp_src(:,:,:)     ! (nlon, nlat, nlev) temperature for density extrap
    INTEGER :: ncid
    INTEGER :: file_year, file_month, file_day, file_hour
    INTEGER :: bcast_size

    ! Buffers for horizontal interpolation path
    REAL(prec), ALLOCATABLE :: vert_buf(:,:,:)    ! (nlon_src, nlat_src, nheights_geo) after vert interp

    IF (PRESENT(rc)) rc = IPE_SUCCESS

    nlev = reader % nlev_src
    nlat = reader % nlat_src
    nlon = reader % nlon_src

    ! Determine file date/hour (handle day rollover)
    file_year  = reader % ref_year
    file_month = reader % ref_month
    file_day   = reader % ref_day
    file_hour  = hour
    ! normalize hour -> day, then day -> month/year (handles month/year rollover)
    DO WHILE ( file_hour >= 24 )
      file_hour = file_hour - 24
      file_day  = file_day + 1
    END DO
    DO WHILE ( file_day > DaysInMonth( file_year, file_month ) )
      file_day   = file_day - DaysInMonth( file_year, file_month )
      file_month = file_month + 1
      IF ( file_month > 12 ) THEN
        file_month = 1
        file_year  = file_year + 1
      END IF
    END DO

    CALL reader % BuildGSMFileName( file_year, file_month, file_day, file_hour, filename )

    IF ( verbose_diag .AND. mpi_layer % rank_id == 0 ) THEN
      WRITE(6,*) 'IPE_FileReader: reading ', TRIM(filename)
    ENDIF

    ALLOCATE( heights_src(nlon, nlat, nlev), &
              field_src(nlon, nlat, nlev),   &
              temp_src(nlon, nlat, nlev) )

    ! Rank 0 reads the file
    IF ( mpi_layer % rank_id == 0 ) THEN

      localrc = nf90_open( TRIM(filename), NF90_NOWRITE, ncid )
      IF ( localrc /= NF90_NOERR ) THEN
        WRITE(6,*) 'IPE_FileReader ERROR: cannot open ', TRIM(filename)
        IF (PRESENT(rc)) rc = IPE_FAILURE
        DEALLOCATE( heights_src, field_src, temp_src )
        RETURN
      ENDIF

      ! Read heights
      CALL ReadVar3D( ncid, 'wam_height_levels', heights_src, localrc )
      IF (localrc /= IPE_SUCCESS) THEN
        IF (PRESENT(rc)) rc = IPE_FAILURE
        DEALLOCATE( heights_src, field_src, temp_src )
        RETURN
      ENDIF
      CALL ReplaceFillValues( heights_src, nlon, nlat, nlev )

      IF ( reader % needs_horiz_interp ) THEN
        ! === Horizontal interpolation path ===
        ! Vert interp at source resolution, then horiz interp to geo grid
        ALLOCATE( vert_buf(nlon, nlat, nheights_geo) )

        ! Temperature
        CALL ReadVar3D( ncid, 'temp_neutral', temp_src, localrc )
        CALL ReplaceFillValues( temp_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, temp_src, altitude_geo, &
             vert_buf, .FALSE. )
        CALL reader % HorizInterp( vert_buf, geo_buf(FLD_TEMP,:,:,:) )

        ! O Density
        CALL ReadVar3D( ncid, 'O_Density', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             vert_buf, .TRUE., species_mass=O_mass, src_temp=temp_src )
        CALL reader % HorizInterp( vert_buf, geo_buf(FLD_O,:,:,:) )

        ! O2 Density
        CALL ReadVar3D( ncid, 'O2_Density', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             vert_buf, .TRUE., species_mass=O2_mass, src_temp=temp_src )
        CALL reader % HorizInterp( vert_buf, geo_buf(FLD_O2,:,:,:) )

        ! N2 Density
        CALL ReadVar3D( ncid, 'N2_Density', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             vert_buf, .TRUE., species_mass=N2_mass, src_temp=temp_src )
        CALL reader % HorizInterp( vert_buf, geo_buf(FLD_N2,:,:,:) )

        ! Eastward wind
        CALL ReadVar3D( ncid, 'eastward_wind_neutral', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             vert_buf, .FALSE. )
        CALL reader % HorizInterp( vert_buf, geo_buf(FLD_UEAST,:,:,:) )

        ! Northward wind
        CALL ReadVar3D( ncid, 'northward_wind_neutral', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             vert_buf, .FALSE. )
        CALL reader % HorizInterp( vert_buf, geo_buf(FLD_UNORTH,:,:,:) )

        ! Upward wind
        CALL ReadVar3D( ncid, 'upward_wind_neutral', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             vert_buf, .FALSE. )
        CALL reader % HorizInterp( vert_buf, geo_buf(FLD_UUP,:,:,:) )

        DEALLOCATE( vert_buf )

      ELSE
        ! === Direct path (source matches geo grid) ===

        ! Read temperature first (needed for density extrapolation above model top)
        CALL ReadVar3D( ncid, 'temp_neutral', temp_src, localrc )
        CALL ReplaceFillValues( temp_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, temp_src, altitude_geo, &
             geo_buf(FLD_TEMP,:,:,:), .FALSE. )

        ! O Density (pass temperature for scale-height extrapolation)
        CALL ReadVar3D( ncid, 'O_Density', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             geo_buf(FLD_O,:,:,:), .TRUE., species_mass=O_mass, src_temp=temp_src )

        ! O2 Density
        CALL ReadVar3D( ncid, 'O2_Density', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             geo_buf(FLD_O2,:,:,:), .TRUE., species_mass=O2_mass, src_temp=temp_src )

        ! N2 Density
        CALL ReadVar3D( ncid, 'N2_Density', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             geo_buf(FLD_N2,:,:,:), .TRUE., species_mass=N2_mass, src_temp=temp_src )

        ! Eastward wind
        CALL ReadVar3D( ncid, 'eastward_wind_neutral', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             geo_buf(FLD_UEAST,:,:,:), .FALSE. )

        ! Northward wind
        CALL ReadVar3D( ncid, 'northward_wind_neutral', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             geo_buf(FLD_UNORTH,:,:,:), .FALSE. )

        ! Upward wind
        CALL ReadVar3D( ncid, 'upward_wind_neutral', field_src, localrc )
        CALL ReplaceFillValues( field_src, nlon, nlat, nlev )
        CALL reader % VertInterp_ToGeoGrid( heights_src, field_src, altitude_geo, &
             geo_buf(FLD_UUP,:,:,:), .FALSE. )

      ENDIF  ! needs_horiz_interp

      localrc = nf90_close( ncid )

    ENDIF

    DEALLOCATE( heights_src, field_src, temp_src )

    ! Broadcast the interpolated geographic grid to all ranks
#ifdef HAVE_MPI
    bcast_size = NFIELDS * nheights_geo * nlat_geo * nlon_geo
    CALL MPI_BCAST( geo_buf, bcast_size, mpi_layer % mpi_prec, 0, &
                    mpi_layer % mpi_communicator, localrc )
#endif

  END SUBROUTINE ReadAndInterp


  !---------------------------------------------------------------------------
  ! ReadAndInterpToApex: Read a source file and interpolate directly to
  ! each local apex grid point. All ranks participate: rank 0 reads the
  ! file and broadcasts raw source data, then each rank interpolates
  ! its own apex points.
  !---------------------------------------------------------------------------
  SUBROUTINE ReadAndInterpToApex( reader, hour, mpi_layer, apex_buf, rc )

    CLASS(IPE_FileReader), INTENT(in)    :: reader
    INTEGER,               INTENT(in)    :: hour
    TYPE(IPE_MPI_Layer),   INTENT(in)    :: mpi_layer
    REAL(prec),            INTENT(out)   :: apex_buf(:,:,:,:)  ! (NFIELDS, nFluxTube, NLP, mp_low:mp_high)
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    ! Local
    INTEGER :: localrc, ncid
    INTEGER :: nlev, nlat, nlon
    CHARACTER(512) :: filename
    INTEGER :: file_year, file_month, file_day, file_hour
    INTEGER :: bcast_size

    ! Source data arrays (all ranks hold a copy after broadcast)
    REAL(sp), ALLOCATABLE :: heights_src(:,:,:)  ! (nlon, nlat, nlev)
    REAL(sp), ALLOCATABLE :: temp_src(:,:,:)     ! (nlon, nlat, nlev)
    REAL(sp), ALLOCATABLE :: field_o(:,:,:)
    REAL(sp), ALLOCATABLE :: field_o2(:,:,:)
    REAL(sp), ALLOCATABLE :: field_n2(:,:,:)
    REAL(sp), ALLOCATABLE :: field_ue(:,:,:)
    REAL(sp), ALLOCATABLE :: field_un(:,:,:)
    REAL(sp), ALLOCATABLE :: field_uw(:,:,:)

    ! Per-apex-point interpolation
    INTEGER :: n, kp, lp, mp, mp_loc
    INTEGER :: i0, i1, j0, j1
    REAL(prec) :: wx, wy
    REAL(prec) :: alt_km
    REAL(prec) :: v00, v10, v01, v11

    IF (PRESENT(rc)) rc = IPE_SUCCESS

    nlev = reader % nlev_src
    nlat = reader % nlat_src
    nlon = reader % nlon_src

    ! Determine file date/hour
    file_year  = reader % ref_year
    file_month = reader % ref_month
    file_day   = reader % ref_day
    file_hour  = hour
    ! normalize hour -> day, then day -> month/year (handles month/year rollover)
    DO WHILE ( file_hour >= 24 )
      file_hour = file_hour - 24
      file_day  = file_day + 1
    END DO
    DO WHILE ( file_day > DaysInMonth( file_year, file_month ) )
      file_day   = file_day - DaysInMonth( file_year, file_month )
      file_month = file_month + 1
      IF ( file_month > 12 ) THEN
        file_month = 1
        file_year  = file_year + 1
      END IF
    END DO

    CALL reader % BuildGSMFileName( file_year, file_month, file_day, file_hour, filename )

    ! Allocate source data on all ranks
    ALLOCATE( heights_src(nlon, nlat, nlev), &
              temp_src(nlon, nlat, nlev),    &
              field_o(nlon, nlat, nlev),     &
              field_o2(nlon, nlat, nlev),    &
              field_n2(nlon, nlat, nlev),    &
              field_ue(nlon, nlat, nlev),    &
              field_un(nlon, nlat, nlev),    &
              field_uw(nlon, nlat, nlev) )

    ! Rank 0 reads all fields from file
    IF ( mpi_layer % rank_id == 0 ) THEN
      IF ( verbose_diag ) WRITE(6,*) 'IPE_FileReader (apex): reading ', TRIM(filename)

      localrc = nf90_open( TRIM(filename), NF90_NOWRITE, ncid )
      IF ( localrc /= NF90_NOERR ) THEN
        WRITE(6,*) 'IPE_FileReader ERROR: cannot open ', TRIM(filename)
        WRITE(6,*) '  ', TRIM(nf90_strerror(localrc))
        IF (PRESENT(rc)) rc = IPE_FAILURE
        ! Set heights to negative to signal error to all ranks after broadcast
        heights_src = -1.0
      ELSE
        CALL ReadVar3D( ncid, 'wam_height_levels', heights_src, localrc )
        CALL ReplaceFillValues( heights_src, nlon, nlat, nlev )
        CALL ReadVar3D( ncid, 'temp_neutral', temp_src, localrc )
        CALL ReplaceFillValues( temp_src, nlon, nlat, nlev )
        CALL ReadVar3D( ncid, 'O_Density', field_o, localrc )
        CALL ReplaceFillValues( field_o, nlon, nlat, nlev )
        CALL ReadVar3D( ncid, 'O2_Density', field_o2, localrc )
        CALL ReplaceFillValues( field_o2, nlon, nlat, nlev )
        CALL ReadVar3D( ncid, 'N2_Density', field_n2, localrc )
        CALL ReplaceFillValues( field_n2, nlon, nlat, nlev )
        CALL ReadVar3D( ncid, 'eastward_wind_neutral', field_ue, localrc )
        CALL ReplaceFillValues( field_ue, nlon, nlat, nlev )
        CALL ReadVar3D( ncid, 'northward_wind_neutral', field_un, localrc )
        CALL ReplaceFillValues( field_un, nlon, nlat, nlev )
        CALL ReadVar3D( ncid, 'upward_wind_neutral', field_uw, localrc )
        CALL ReplaceFillValues( field_uw, nlon, nlat, nlev )

        localrc = nf90_close( ncid )

        ! Smooth source fields to reduce tube-to-tube noise.
        ! NOTE: Do NOT smooth heights_src — heights define the vertical
        ! coordinate for each column. Smoothing them corrupts the
        ! height-field relationship and produces incorrect vertical interp.
        ! No smoothing — match coupled run behavior (ESMF direct regrid)

        IF ( verbose_diag ) WRITE(6,*) 'IPE_FileReader (apex): file read (no smoothing)'
      ENDIF
    ENDIF

    ! Broadcast all source data to all ranks
#ifdef HAVE_MPI
    bcast_size = nlon * nlat * nlev
    CALL MPI_BCAST( heights_src, bcast_size, MPI_REAL, 0, &
                    mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( temp_src, bcast_size, MPI_REAL, 0, &
                    mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( field_o, bcast_size, MPI_REAL, 0, &
                    mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( field_o2, bcast_size, MPI_REAL, 0, &
                    mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( field_n2, bcast_size, MPI_REAL, 0, &
                    mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( field_ue, bcast_size, MPI_REAL, 0, &
                    mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( field_un, bcast_size, MPI_REAL, 0, &
                    mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( field_uw, bcast_size, MPI_REAL, 0, &
                    mpi_layer % mpi_communicator, localrc )
#endif

    ! Each rank interpolates its own local apex points
    DO n = 1, reader % n_local_apex
      kp = reader % apex_kp_arr(n)
      lp = reader % apex_lp_arr(n)
      mp = reader % apex_mp_arr(n)
      ! apex_buf is assumed-shape so lower bound is 1, not mp_low
      mp_loc = mp - reader % apex_mp_low + 1

      i0 = reader % apex_ilon_lo(n)
      i1 = MOD(i0, nlon) + 1  ! wrap-around
      j0 = reader % apex_jlat_lo(n)
      j1 = MIN(j0 + 1, nlat)
      wx = reader % apex_wlon_arr(n)
      wy = reader % apex_wlat_arr(n)
      alt_km = reader % apex_alt_arr(n)

      ! Temperature: vert interp at 4 surrounding columns, then bilinear
      v00 = VertInterpColumn( heights_src(i0,j0,:), temp_src(i0,j0,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i0,j0,:) )
      v10 = VertInterpColumn( heights_src(i1,j0,:), temp_src(i1,j0,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i1,j0,:) )
      v01 = VertInterpColumn( heights_src(i0,j1,:), temp_src(i0,j1,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i0,j1,:) )
      v11 = VertInterpColumn( heights_src(i1,j1,:), temp_src(i1,j1,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i1,j1,:) )
      apex_buf(FLD_TEMP,kp,lp,mp_loc) = BilinearVal(v00, v10, v01, v11, wx, wy)

      ! O density (log-space bilinear to match legacy trilinear behavior)
      v00 = VertInterpColumn( heights_src(i0,j0,:), field_o(i0,j0,:), nlev, alt_km, &
            .TRUE., O_mass, temp_src(i0,j0,:) )
      v10 = VertInterpColumn( heights_src(i1,j0,:), field_o(i1,j0,:), nlev, alt_km, &
            .TRUE., O_mass, temp_src(i1,j0,:) )
      v01 = VertInterpColumn( heights_src(i0,j1,:), field_o(i0,j1,:), nlev, alt_km, &
            .TRUE., O_mass, temp_src(i0,j1,:) )
      v11 = VertInterpColumn( heights_src(i1,j1,:), field_o(i1,j1,:), nlev, alt_km, &
            .TRUE., O_mass, temp_src(i1,j1,:) )
      apex_buf(FLD_O,kp,lp,mp_loc) = BilinearValLog(v00, v10, v01, v11, wx, wy)

      ! O2 density (log-space bilinear)
      v00 = VertInterpColumn( heights_src(i0,j0,:), field_o2(i0,j0,:), nlev, alt_km, &
            .TRUE., O2_mass, temp_src(i0,j0,:) )
      v10 = VertInterpColumn( heights_src(i1,j0,:), field_o2(i1,j0,:), nlev, alt_km, &
            .TRUE., O2_mass, temp_src(i1,j0,:) )
      v01 = VertInterpColumn( heights_src(i0,j1,:), field_o2(i0,j1,:), nlev, alt_km, &
            .TRUE., O2_mass, temp_src(i0,j1,:) )
      v11 = VertInterpColumn( heights_src(i1,j1,:), field_o2(i1,j1,:), nlev, alt_km, &
            .TRUE., O2_mass, temp_src(i1,j1,:) )
      apex_buf(FLD_O2,kp,lp,mp_loc) = BilinearValLog(v00, v10, v01, v11, wx, wy)

      ! N2 density (log-space bilinear)
      v00 = VertInterpColumn( heights_src(i0,j0,:), field_n2(i0,j0,:), nlev, alt_km, &
            .TRUE., N2_mass, temp_src(i0,j0,:) )
      v10 = VertInterpColumn( heights_src(i1,j0,:), field_n2(i1,j0,:), nlev, alt_km, &
            .TRUE., N2_mass, temp_src(i1,j0,:) )
      v01 = VertInterpColumn( heights_src(i0,j1,:), field_n2(i0,j1,:), nlev, alt_km, &
            .TRUE., N2_mass, temp_src(i0,j1,:) )
      v11 = VertInterpColumn( heights_src(i1,j1,:), field_n2(i1,j1,:), nlev, alt_km, &
            .TRUE., N2_mass, temp_src(i1,j1,:) )
      apex_buf(FLD_N2,kp,lp,mp_loc) = BilinearValLog(v00, v10, v01, v11, wx, wy)

      ! Eastward wind
      v00 = VertInterpColumn( heights_src(i0,j0,:), field_ue(i0,j0,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i0,j0,:) )
      v10 = VertInterpColumn( heights_src(i1,j0,:), field_ue(i1,j0,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i1,j0,:) )
      v01 = VertInterpColumn( heights_src(i0,j1,:), field_ue(i0,j1,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i0,j1,:) )
      v11 = VertInterpColumn( heights_src(i1,j1,:), field_ue(i1,j1,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i1,j1,:) )
      apex_buf(FLD_UEAST,kp,lp,mp_loc) = BilinearVal(v00, v10, v01, v11, wx, wy)

      ! Northward wind
      v00 = VertInterpColumn( heights_src(i0,j0,:), field_un(i0,j0,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i0,j0,:) )
      v10 = VertInterpColumn( heights_src(i1,j0,:), field_un(i1,j0,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i1,j0,:) )
      v01 = VertInterpColumn( heights_src(i0,j1,:), field_un(i0,j1,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i0,j1,:) )
      v11 = VertInterpColumn( heights_src(i1,j1,:), field_un(i1,j1,:), nlev, alt_km, &
            .FALSE., O_mass, temp_src(i1,j1,:) )
      apex_buf(FLD_UNORTH,kp,lp,mp_loc) = BilinearVal(v00, v10, v01, v11, wx, wy)

      ! Upward wind — set to zero (FV3WAM dzdt->w conversion is unreliable
      ! in the thermosphere due to very low rho, producing ±200 m/s noise.
      ! This matches HWM14 behavior which also sets vertical wind to zero.)
      apex_buf(FLD_UUP,kp,lp,mp_loc) = 0.0_prec

    ENDDO  ! n (local apex points)

    DEALLOCATE( heights_src, temp_src, field_o, field_o2, field_n2, &
                field_ue, field_un, field_uw )

  END SUBROUTINE ReadAndInterpToApex


  !---------------------------------------------------------------------------
  ! VertInterpColumn: Vertically interpolate one column to a single target
  ! altitude. Same logic as VertInterp_ToGeoGrid but for a single point.
  !---------------------------------------------------------------------------
  FUNCTION VertInterpColumn( src_heights, src_field, nlev, target_alt_km, &
       is_density, species_mass, src_temp ) RESULT(val)

    REAL(sp),   INTENT(in) :: src_heights(:)  ! (nlev) in meters
    REAL(sp),   INTENT(in) :: src_field(:)    ! (nlev)
    INTEGER,    INTENT(in) :: nlev
    REAL(prec), INTENT(in) :: target_alt_km   ! target altitude in km
    LOGICAL,    INTENT(in) :: is_density
    REAL(prec), INTENT(in) :: species_mass    ! molecular mass in AMU
    REAL(sp),   INTENT(in) :: src_temp(:)     ! (nlev) temperature for density extrap

    REAL(prec) :: val

    ! Local
    INTEGER :: kk
    REAL(prec) :: z_target, z_lo, z_hi, f_lo, f_hi, frac
    REAL(prec) :: z_top, f_top, T_top, dz

    z_target = target_alt_km  ! km

    ! Find top of source data for this column
    z_top = REAL(src_heights(nlev), prec) * m_to_km

    IF ( z_target <= REAL(src_heights(1), prec) * m_to_km ) THEN
      ! Below source data: use lowest level
      val = REAL(src_field(1), prec)

    ELSEIF ( z_target >= z_top ) THEN
      ! Above source data: extrapolate
      f_top = REAL(src_field(nlev), prec)

      IF ( is_density ) THEN
        T_top = MAX( REAL(src_temp(nlev), prec), 200.0_prec )
        dz = (z_target - z_top) * km_to_m  ! convert to meters
        val = f_top * EXP( -dz * species_mass * AMU * G0 / (kBoltz * T_top) )
        val = MAX(val, 1.0e-30_prec)
      ELSE
        val = f_top  ! isothermal / constant above top
      ENDIF

    ELSE
      ! Within source range: find bracketing levels and interpolate
      val = REAL(src_field(1), prec)  ! fallback
      DO kk = 1, nlev - 1
        z_lo = REAL(src_heights(kk  ), prec) * m_to_km
        z_hi = REAL(src_heights(kk+1), prec) * m_to_km
        IF ( z_target >= z_lo .AND. z_target < z_hi ) THEN
          frac = (z_target - z_lo) / (z_hi - z_lo)
          f_lo = REAL(src_field(kk  ), prec)
          f_hi = REAL(src_field(kk+1), prec)

          IF ( is_density ) THEN
            f_lo = MAX(f_lo, 1.0e-30_prec)
            f_hi = MAX(f_hi, 1.0e-30_prec)
            val = EXP( (1.0_prec - frac)*LOG(f_lo) + frac*LOG(f_hi) )
          ELSE
            val = (1.0_prec - frac) * f_lo + frac * f_hi
          ENDIF
          EXIT
        ENDIF
      ENDDO
    ENDIF

  END FUNCTION VertInterpColumn


  !---------------------------------------------------------------------------
  ! BilinearVal: Simple bilinear interpolation from 4 corner values.
  !---------------------------------------------------------------------------
  FUNCTION BilinearVal( v00, v10, v01, v11, wx, wy ) RESULT(val)

    REAL(prec), INTENT(in) :: v00, v10, v01, v11, wx, wy
    REAL(prec) :: val

    val = (1.0_prec - wx) * (1.0_prec - wy) * v00 &
        + wx              * (1.0_prec - wy) * v10 &
        + (1.0_prec - wx) * wy              * v01 &
        + wx              * wy              * v11

  END FUNCTION BilinearVal


  !---------------------------------------------------------------------------
  ! BilinearValLog: Bilinear interpolation in log-space for density fields.
  ! Takes 4 positive corner values, interpolates LOG(v), returns EXP(result).
  !---------------------------------------------------------------------------
  FUNCTION BilinearValLog( v00, v10, v01, v11, wx, wy ) RESULT(val)

    REAL(prec), INTENT(in) :: v00, v10, v01, v11, wx, wy
    REAL(prec) :: val
    REAL(prec) :: l00, l10, l01, l11

    l00 = LOG(MAX(v00, 1.0e-30_prec))
    l10 = LOG(MAX(v10, 1.0e-30_prec))
    l01 = LOG(MAX(v01, 1.0e-30_prec))
    l11 = LOG(MAX(v11, 1.0e-30_prec))

    val = EXP( (1.0_prec - wx) * (1.0_prec - wy) * l00 &
             + wx              * (1.0_prec - wy) * l10 &
             + (1.0_prec - wx) * wy              * l01 &
             + wx              * wy              * l11 )

  END FUNCTION BilinearValLog


  !---------------------------------------------------------------------------
  ! VertInterp_ToGeoGrid: Vertically interpolate one field from source levels
  ! to IPE geographic height grid for all columns.
  !---------------------------------------------------------------------------
  SUBROUTINE VertInterp_ToGeoGrid( reader, src_heights, src_field, altitude_geo, &
       dst_field, is_density, species_mass, src_temp )

    CLASS(IPE_FileReader), INTENT(in)  :: reader
    REAL(sp),              INTENT(in)  :: src_heights(:,:,:)  ! (nlon, nlat, nlev) in meters
    REAL(sp),              INTENT(in)  :: src_field(:,:,:)    ! (nlon, nlat, nlev)
    REAL(prec),            INTENT(in)  :: altitude_geo(:)     ! (nheights_geo) in km
    REAL(prec),            INTENT(out) :: dst_field(:,:,:)    ! (nlon, nlat, nheights_geo)
    LOGICAL,               INTENT(in)  :: is_density
    REAL(prec), OPTIONAL,  INTENT(in)  :: species_mass        ! molecular mass in AMU
    REAL(sp),   OPTIONAL,  INTENT(in)  :: src_temp(:,:,:)     ! (nlon, nlat, nlev) temperature

    ! Local
    INTEGER :: i, j, k, kk
    INTEGER :: nlev, nlat, nlon
    REAL(prec) :: z_target, z_lo, z_hi, f_lo, f_hi, frac
    REAL(prec) :: z_top, f_top, T_top, dz
    REAL(prec) :: mol_mass

    nlev = reader % nlev_src
    nlat = reader % nlat_src
    nlon = reader % nlon_src

    ! Default molecular mass for density extrapolation
    mol_mass = O_mass
    IF ( PRESENT(species_mass) ) mol_mass = species_mass

    DO i = 1, nlon
      DO j = 1, nlat

        ! Find top of source data for this column
        z_top = REAL(src_heights(i, j, nlev), prec) * m_to_km  ! convert m to km

        DO k = 1, nheights_geo

          z_target = altitude_geo(k)  ! km

          ! Find bracketing source levels
          ! Source heights increase with level index (bottom to top)
          IF ( z_target <= REAL(src_heights(i, j, 1), prec) * m_to_km ) THEN
            ! Below source data: use lowest level
            dst_field(i, j, k) = REAL(src_field(i, j, 1), prec)

          ELSEIF ( z_target >= z_top ) THEN
            ! Above source data: extrapolate
            f_top = REAL(src_field(i, j, nlev), prec)

            IF ( is_density ) THEN
              ! Scale-height extrapolation using temperature at model top
              IF ( PRESENT(src_temp) ) THEN
                T_top = MAX( REAL(src_temp(i, j, nlev), prec), 200.0_prec )
              ELSE
                T_top = 1000.0_prec  ! fallback exospheric temperature
              ENDIF
              ! n(z) = n(z_top) * exp(-m*g*dz / (kB*T))
              dz = (z_target - z_top) * km_to_m  ! convert to meters
              dst_field(i, j, k) = f_top * EXP( -dz * mol_mass * AMU * G0 / &
                                   (kBoltz * T_top) )
              dst_field(i, j, k) = MAX(dst_field(i, j, k), 1.0e-30_prec)
            ELSE
              ! Temperature/winds: hold constant above model top (isothermal)
              dst_field(i, j, k) = f_top
            ENDIF

          ELSE
            ! Within source range: interpolate
            DO kk = 1, nlev - 1
              z_lo = REAL(src_heights(i, j, kk  ), prec) * m_to_km
              z_hi = REAL(src_heights(i, j, kk+1), prec) * m_to_km
              IF ( z_target >= z_lo .AND. z_target < z_hi ) THEN
                frac = (z_target - z_lo) / (z_hi - z_lo)
                f_lo = REAL(src_field(i, j, kk  ), prec)
                f_hi = REAL(src_field(i, j, kk+1), prec)

                IF ( is_density ) THEN
                  ! Log-linear interpolation for densities
                  f_lo = MAX(f_lo, 1.0e-30_prec)
                  f_hi = MAX(f_hi, 1.0e-30_prec)
                  dst_field(i, j, k) = EXP( (1.0_prec - frac)*LOG(f_lo) + frac*LOG(f_hi) )
                ELSE
                  ! Linear interpolation for temperature and winds
                  dst_field(i, j, k) = (1.0_prec - frac) * f_lo + frac * f_hi
                ENDIF
                EXIT
              ENDIF
            ENDDO
          ENDIF

        ENDDO  ! k (heights)
      ENDDO    ! j (lat)
    ENDDO      ! i (lon)

  END SUBROUTINE VertInterp_ToGeoGrid


  !---------------------------------------------------------------------------
  ! BuildGSMFileName: Construct GSM filename from date/hour
  ! Format: <dir>/gsm.YYYYMMDD_HH0000.nc
  !---------------------------------------------------------------------------
  SUBROUTINE BuildGSMFileName( reader, year, month, day, hour, filename )

    CLASS(IPE_FileReader), INTENT(in)  :: reader
    INTEGER,               INTENT(in)  :: year, month, day, hour
    CHARACTER(*),          INTENT(out) :: filename

    CHARACTER(8)  :: date_str
    CHARACTER(6)  :: time_str

    WRITE(date_str, '(I4.4,I2.2,I2.2)') year, month, day
    WRITE(time_str, '(I2.2,A)') hour, '0000'

    filename = TRIM(reader % file_dir) // '/gsm.' // date_str // '_' // time_str // '.nc'

  END SUBROUTINE BuildGSMFileName

  !---------------------------------------------------------------------------
  ! DaysInMonth: days in a given month, Gregorian leap years included.
  !---------------------------------------------------------------------------
  PURE INTEGER FUNCTION DaysInMonth( year, month )
    INTEGER, INTENT(in) :: year, month
    INTEGER, PARAMETER :: dim(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    DaysInMonth = dim(month)
    IF ( month == 2 .AND. MOD(year,4) == 0 .AND. &
         ( MOD(year,100) /= 0 .OR. MOD(year,400) == 0 ) ) DaysInMonth = 29
  END FUNCTION DaysInMonth


  !---------------------------------------------------------------------------
  ! ReadVar3D: Helper to read a 3D float variable from NetCDF
  !---------------------------------------------------------------------------
  SUBROUTINE ReadVar3D( ncid, varname, data, rc )

    INTEGER,      INTENT(in)  :: ncid
    CHARACTER(*), INTENT(in)  :: varname
    REAL(sp),     INTENT(out) :: data(:,:,:)
    INTEGER,      INTENT(out) :: rc

    INTEGER :: varid, localrc

    rc = IPE_SUCCESS

    localrc = nf90_inq_varid( ncid, varname, varid )
    IF ( localrc /= NF90_NOERR ) THEN
      WRITE(6,*) 'IPE_FileReader ERROR: variable not found: ', TRIM(varname)
      WRITE(6,*) '  ', TRIM(nf90_strerror(localrc))
      rc = IPE_FAILURE
      RETURN
    ENDIF

    localrc = nf90_get_var( ncid, varid, data )
    IF ( localrc /= NF90_NOERR ) THEN
      WRITE(6,*) 'IPE_FileReader ERROR: cannot read variable: ', TRIM(varname)
      WRITE(6,*) '  ', TRIM(nf90_strerror(localrc))
      rc = IPE_FAILURE
      RETURN
    ENDIF

  END SUBROUTINE ReadVar3D


  !---------------------------------------------------------------------------
  ! ReplaceFillValues: Replace fill values (-99999) with nearest valid
  ! neighbor along latitude dimension. Forward pass (south to north),
  ! then backward pass (north to south) to cover both poles.
  ! Same logic as dataWamCap.F90.
  !---------------------------------------------------------------------------
  SUBROUTINE ReplaceFillValues( data, nlon, nlat, nlev )

    REAL(sp), INTENT(inout) :: data(:,:,:)  ! (nlon, nlat, nlev)
    INTEGER,  INTENT(in)    :: nlon, nlat, nlev

    INTEGER :: i, j, k
    REAL(sp), PARAMETER :: fill_thresh = -99990.0

    ! Forward pass: south to north
    DO k = 1, nlev
      DO j = 2, nlat
        DO i = 1, nlon
          IF ( data(i,j,k) < fill_thresh .AND. &
               data(i,j-1,k) > fill_thresh ) THEN
            data(i,j,k) = data(i,j-1,k)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! Backward pass: north to south
    DO k = 1, nlev
      DO j = nlat-1, 1, -1
        DO i = 1, nlon
          IF ( data(i,j,k) < fill_thresh .AND. &
               data(i,j+1,k) > fill_thresh ) THEN
            data(i,j,k) = data(i,j+1,k)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE ReplaceFillValues


  !---------------------------------------------------------------------------
  ! Compute_MSIS_GeoGrid: Compute MSIS climatological values on the IPE
  ! geographic grid at a given UT time. Only rank 0 computes; result is
  ! broadcast to all ranks.
  !
  ! Output msis_buf has shape (NFIELDS_MSIS, nlon_geo, nlat_geo, nheights_geo)
  ! Fields: MFLD_TEMP (K), MFLD_O, MFLD_O2, MFLD_N2 (cm^-3)
  !---------------------------------------------------------------------------
  SUBROUTINE Compute_MSIS_GeoGrid( reader, utime_sec, day_of_year, &
       forcing, altitude_geo, mpi_layer, msis_buf, rc )

    CLASS(IPE_FileReader), INTENT(in)    :: reader
    REAL(prec),            INTENT(in)    :: utime_sec     ! UT seconds
    INTEGER,               INTENT(in)    :: day_of_year
    TYPE(IPE_Forcing),     INTENT(in)    :: forcing
    REAL(prec),            INTENT(in)    :: altitude_geo(:)  ! (nheights_geo) in km
    TYPE(IPE_MPI_Layer),   INTENT(in)    :: mpi_layer
    REAL(prec),            INTENT(out)   :: msis_buf(:,:,:,:) ! (NFIELDS_MSIS, nlon_geo, nlat_geo, nheights_geo)
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    ! Local
    INTEGER :: i, j, k, localrc
    INTEGER(4) :: iyd
    REAL(msis_dp) :: msis_sec, msis_alt, msis_lat, msis_lon, msis_stl
    REAL(msis_dp) :: msis_f107a, msis_f107d
    REAL(msis_dp), DIMENSION(7) :: msis_ap
    REAL(msis_dp), DIMENSION(9) :: densities
    REAL(msis_dp), DIMENSION(2) :: temperatures
    INTEGER, PARAMETER :: msis_mass = 48
    REAL(prec), DIMENSION(7) :: AP_arr
    REAL(prec) :: geo_lat, geo_lon
    INTEGER :: bcast_size

    IF (PRESENT(rc)) rc = IPE_SUCCESS

    IF ( mpi_layer % rank_id == 0 ) THEN

      ! Get forcing parameters
      AP_arr    = forcing % GetAP()
      msis_ap   = REAL(AP_arr, KIND=msis_dp)
      msis_f107a = REAL(forcing % f107_81day_avg(forcing % current_index), KIND=msis_dp)
      msis_f107d = REAL(forcing % f107(forcing % current_index), KIND=msis_dp)
      msis_sec   = REAL(utime_sec, KIND=msis_dp)

      iyd = 99000 + day_of_year

      DO k = 1, nheights_geo
        msis_alt = REAL(altitude_geo(k), KIND=msis_dp)
        DO j = 1, nlat_geo
          ! lat: -90 to 90 at 2 deg spacing
          geo_lat = -90.0_prec + REAL(j - 1, prec) * 2.0_prec
          msis_lat = REAL(geo_lat, KIND=msis_dp)
          DO i = 1, nlon_geo
            ! lon: 0 to 356 at 4 deg spacing
            geo_lon = REAL(i - 1, prec) * 4.0_prec
            msis_lon = REAL(geo_lon, KIND=msis_dp)
            msis_stl = REAL(utime_sec / 3600.0_prec + geo_lon / 15.0_prec, KIND=msis_dp)

            densities    = 0.0_msis_dp
            temperatures = 0.0_msis_dp

            CALL gtd7( iyd, msis_sec, msis_alt, msis_lat, msis_lon, &
                       msis_stl, msis_f107a, msis_f107d, msis_ap, &
                       msis_mass, densities, temperatures )

            msis_buf(MFLD_TEMP, i, j, k) = REAL(temperatures(2), prec)
            msis_buf(MFLD_O,    i, j, k) = REAL(densities(2), prec)  ! O in cm^-3
            msis_buf(MFLD_O2,   i, j, k) = REAL(densities(4), prec)  ! O2 in cm^-3
            msis_buf(MFLD_N2,   i, j, k) = REAL(densities(3), prec)  ! N2 in cm^-3
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    ! Broadcast to all ranks
#ifdef HAVE_MPI
    bcast_size = NFIELDS_MSIS * nlon_geo * nlat_geo * nheights_geo
    CALL MPI_BCAST( msis_buf, bcast_size, mpi_layer % mpi_prec, 0, &
                    mpi_layer % mpi_communicator, localrc )
#endif

  END SUBROUTINE Compute_MSIS_GeoGrid


  !---------------------------------------------------------------------------
  ! ComputeHorizWeights: Read source grid coordinates from the first GSM file
  ! and precompute bilinear interpolation weights to the IPE geographic grid.
  ! Called once during Init when source grid != IPE geo grid.
  !---------------------------------------------------------------------------
  SUBROUTINE ComputeHorizWeights( reader, mpi_layer )

    CLASS(IPE_FileReader), INTENT(inout) :: reader
    TYPE(IPE_MPI_Layer),   INTENT(in)    :: mpi_layer

    ! Local
    INTEGER :: i, j, nlon_s, nlat_s
    INTEGER :: localrc, ncid, varid, dimid
    REAL(prec) :: target_lon, target_lat, dlon
    REAL(prec) :: frac_lon
    CHARACTER(512) :: filename

    nlon_s = reader % nlon_src
    nlat_s = reader % nlat_src

    ALLOCATE( reader % src_lons(nlon_s), reader % src_lats(nlat_s) )
    ALLOCATE( reader % ilon_lo(nlon_geo), reader % wlon(nlon_geo) )
    ALLOCATE( reader % jlat_lo(nlat_geo), reader % wlat(nlat_geo) )

    ! Rank 0 reads the coordinate arrays from the first file
    IF ( mpi_layer % rank_id == 0 ) THEN

      CALL reader % BuildGSMFileName( reader % ref_year, reader % ref_month, &
                                       reader % ref_day, 0, filename )
      localrc = nf90_open( TRIM(filename), NF90_NOWRITE, ncid )

      ! Read grid_xt (longitudes)
      localrc = nf90_inq_varid( ncid, 'grid_xt', varid )
      IF ( localrc == NF90_NOERR ) THEN
        localrc = nf90_get_var( ncid, varid, reader % src_lons )
      ELSE
        ! Fallback: assume uniform longitude grid
        dlon = 360.0_prec / REAL(nlon_s, prec)
        DO i = 1, nlon_s
          reader % src_lons(i) = REAL(i - 1, prec) * dlon
        ENDDO
      ENDIF

      ! Read grid_yt (latitudes, should be S->N after pre-processor)
      localrc = nf90_inq_varid( ncid, 'grid_yt', varid )
      IF ( localrc == NF90_NOERR ) THEN
        localrc = nf90_get_var( ncid, varid, reader % src_lats )
      ELSE
        ! Fallback: assume uniform latitude grid
        DO j = 1, nlat_s
          reader % src_lats(j) = -90.0_prec + REAL(j - 1, prec) * (180.0_prec / REAL(nlat_s - 1, prec))
        ENDDO
      ENDIF

      localrc = nf90_close( ncid )

      WRITE(6,*) 'IPE_FileReader: src lon range = ', &
                  reader % src_lons(1), ' to ', reader % src_lons(nlon_s)
      WRITE(6,*) 'IPE_FileReader: src lat range = ', &
                  reader % src_lats(1), ' to ', reader % src_lats(nlat_s)

    ENDIF

    ! Broadcast coordinate arrays to all ranks
#ifdef HAVE_MPI
    CALL MPI_BCAST( reader % src_lons, nlon_s, mpi_layer % mpi_prec, 0, &
                    mpi_layer % mpi_communicator, localrc )
    CALL MPI_BCAST( reader % src_lats, nlat_s, mpi_layer % mpi_prec, 0, &
                    mpi_layer % mpi_communicator, localrc )
#endif

    ! --- Compute longitude weights ---
    ! Source longitude is assumed uniform. IPE geo lon: 0, 4, 8, ..., 356 degrees
    dlon = reader % src_lons(2) - reader % src_lons(1)
    DO i = 1, nlon_geo
      target_lon = REAL(i - 1, prec) * 4.0_prec  ! IPE geo lon
      frac_lon = target_lon / dlon
      reader % ilon_lo(i) = INT(frac_lon) + 1    ! 1-based
      ! Handle wrap-around
      IF ( reader % ilon_lo(i) < 1 ) reader % ilon_lo(i) = nlon_s
      IF ( reader % ilon_lo(i) > nlon_s ) reader % ilon_lo(i) = 1
      reader % wlon(i) = frac_lon - REAL(INT(frac_lon), prec)
    ENDDO

    ! --- Compute latitude weights ---
    ! Source lats may be non-uniform (Gaussian). IPE geo lat: -90, -88, ..., 90 degrees.
    ! Source lats are S->N (ascending) after pre-processor.
    DO j = 1, nlat_geo
      target_lat = -90.0_prec + REAL(j - 1, prec) * 2.0_prec  ! IPE geo lat

      ! Clamp to source range
      IF ( target_lat <= reader % src_lats(1) ) THEN
        reader % jlat_lo(j) = 1
        reader % wlat(j)    = 0.0_prec
      ELSEIF ( target_lat >= reader % src_lats(nlat_s) ) THEN
        reader % jlat_lo(j) = nlat_s - 1
        reader % wlat(j)    = 1.0_prec
      ELSE
        ! Binary search for bracketing indices
        CALL BinarySearchLat( reader % src_lats, nlat_s, target_lat, &
                              reader % jlat_lo(j), reader % wlat(j) )
      ENDIF
    ENDDO

    IF ( mpi_layer % rank_id == 0 ) THEN
      WRITE(6,*) 'IPE_FileReader: horizontal interp weights computed'
    ENDIF

  END SUBROUTINE ComputeHorizWeights


  !---------------------------------------------------------------------------
  ! BinarySearchLat: Find bracketing index and weight for target in ascending
  ! latitude array. Returns jlo such that lats(jlo) <= target < lats(jlo+1).
  !---------------------------------------------------------------------------
  SUBROUTINE BinarySearchLat( lats, n, target, jlo, weight )

    REAL(prec), INTENT(in)  :: lats(:)
    INTEGER,    INTENT(in)  :: n
    REAL(prec), INTENT(in)  :: target
    INTEGER,    INTENT(out) :: jlo
    REAL(prec), INTENT(out) :: weight

    INTEGER :: lo, hi, mid

    lo = 1
    hi = n

    DO WHILE ( hi - lo > 1 )
      mid = (lo + hi) / 2
      IF ( lats(mid) <= target ) THEN
        lo = mid
      ELSE
        hi = mid
      ENDIF
    ENDDO

    jlo = lo
    weight = (target - lats(lo)) / (lats(lo + 1) - lats(lo))
    weight = MAX(0.0_prec, MIN(1.0_prec, weight))

  END SUBROUTINE BinarySearchLat


  !---------------------------------------------------------------------------
  ! HorizInterp: Bilinear horizontal interpolation from source grid to IPE
  ! geographic grid, using precomputed weights. Operates on a single 3D field.
  !---------------------------------------------------------------------------
  SUBROUTINE HorizInterp( reader, src, dst )

    CLASS(IPE_FileReader), INTENT(in)  :: reader
    REAL(prec),            INTENT(in)  :: src(:,:,:)  ! (nlon_src, nlat_src, nheights_geo)
    REAL(prec),            INTENT(out) :: dst(:,:,:)  ! (nlon_geo, nlat_geo, nheights_geo)

    INTEGER :: i, j, k
    INTEGER :: i0, i1, j0, j1
    REAL(prec) :: wx, wy
    INTEGER :: nlon_s

    nlon_s = reader % nlon_src

    DO k = 1, nheights_geo
      DO j = 1, nlat_geo
        j0 = reader % jlat_lo(j)
        j1 = j0 + 1
        IF ( j1 > reader % nlat_src ) j1 = reader % nlat_src
        wy = reader % wlat(j)
        DO i = 1, nlon_geo
          i0 = reader % ilon_lo(i)
          i1 = MOD(i0, nlon_s) + 1    ! wrap-around for longitude
          wx = reader % wlon(i)
          dst(i,j,k) = (1.0_prec - wx) * (1.0_prec - wy) * src(i0, j0, k) &
                      + wx              * (1.0_prec - wy) * src(i1, j0, k) &
                      + (1.0_prec - wx) * wy              * src(i0, j1, k) &
                      + wx              * wy              * src(i1, j1, k)
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE HorizInterp


  !---------------------------------------------------------------------------
  ! SmoothSourceField3D: Apply boxcar smoothing to a 3D source field.
  ! Smooths in lon and lat dimensions with periodic lon boundaries.
  ! hw_lon, hw_lat are half-widths: a hw of 2 means a 5-point average.
  ! On a ~1 deg grid, hw_lon=2, hw_lat=1 gives ~5 deg x ~3 deg smoothing,
  ! comparable to the 90x91 (4 deg x 2 deg) geo grid resolution.
  !---------------------------------------------------------------------------
  SUBROUTINE SmoothSourceField3D( field, nlon, nlat, nlev, hw_lon, hw_lat )

    REAL(sp),  INTENT(inout) :: field(:,:,:)  ! (nlon, nlat, nlev)
    INTEGER,   INTENT(in)    :: nlon, nlat, nlev
    INTEGER,   INTENT(in)    :: hw_lon   ! half-width in longitude (grid cells)
    INTEGER,   INTENT(in)    :: hw_lat   ! half-width in latitude (grid cells)

    ! Local
    REAL(sp), ALLOCATABLE :: tmp(:,:)  ! scratch for one level
    INTEGER :: i, j, k, di, dj, ii, jj, cnt
    REAL(sp) :: total

    IF ( hw_lon <= 0 .AND. hw_lat <= 0 ) RETURN

    ALLOCATE( tmp(nlon, nlat) )

    DO k = 1, nlev
      ! Smooth this level
      DO j = 1, nlat
        DO i = 1, nlon
          total = 0.0
          cnt   = 0
          DO dj = -hw_lat, hw_lat
            jj = j + dj
            IF ( jj < 1 .OR. jj > nlat ) CYCLE
            DO di = -hw_lon, hw_lon
              ii = MOD(i - 1 + di + nlon, nlon) + 1   ! periodic lon
              total = total + field(ii, jj, k)
              cnt   = cnt + 1
            ENDDO
          ENDDO
          tmp(i,j) = total / REAL(cnt)
        ENDDO
      ENDDO
      field(:,:,k) = tmp(:,:)
    ENDDO

    DEALLOCATE( tmp )

  END SUBROUTINE SmoothSourceField3D


END MODULE IPE_Neutrals_FileReader
