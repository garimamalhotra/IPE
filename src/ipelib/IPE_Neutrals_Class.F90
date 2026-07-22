MODULE IPE_Neutrals_Class

  USE IPE_Precision
  USE IPE_Constants_Dictionary
  USE IPE_Grid_Class
  USE IPE_Time_Class
  USE IPE_Forcing_Class
  USE IPE_MPI_Layer_Class
  USE IPE_Model_Parameters_Class
  USE IPE_Neutrals_FileReader
  USE ipe_error_module

  ! MSIS
  USE physics_msis ! gtd7
  USE utils_constants, ONLY : msis_dp => dp

  IMPLICIT NONE

  TYPE IPE_Neutrals

    INTEGER :: nFluxTube, NLP, NMP
    INTEGER :: mp_low, mp_high

    REAL(prec), POINTER :: helium(:,:,:)             ! he
    REAL(prec), POINTER :: oxygen(:,:,:)             ! on
    REAL(prec), POINTER :: molecular_oxygen(:,:,:)   ! o2n
    REAL(prec), POINTER :: molecular_nitrogen(:,:,:) ! n2n
    REAL(prec), POINTER :: nitrogen(:,:,:)           ! n4s
    REAL(prec), POINTER :: hydrogen(:,:,:)           ! hn
    REAL(prec), POINTER :: temperature(:,:,:)
    REAL(prec), POINTER :: temperature_inf(:,:,:)
    REAL(prec), POINTER :: velocity_geographic(:,:,:,:)
    REAL(prec), POINTER :: velocity_apex(:,:,:,:)

    ! GSM/NetCDF file reader for direct neutral input
    TYPE(IPE_FileReader) :: file_reader

    ! Interpolated fields
    REAL(prec), ALLOCATABLE :: geo_helium(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_oxygen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_molecular_oxygen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_molecular_nitrogen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_nitrogen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_hydrogen(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_temperature(:,:,:)
    REAL(prec), ALLOCATABLE :: geo_velocity(:,:,:,:)

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Neutrals
      PROCEDURE :: Trash => Trash_IPE_Neutrals

      PROCEDURE :: Update => Update_IPE_Neutrals

      ! PRIVATE Routines
      PROCEDURE, PRIVATE :: IPE_Neutrals_Empirical
      PROCEDURE, PRIVATE :: IPE_Neutrals_Extrapolate
      PROCEDURE, PRIVATE :: IPE_Neutrals_SetTempInf
      PROCEDURE, PRIVATE :: Geographic_to_Apex_Velocity
      PROCEDURE, PRIVATE :: Geographic_to_Apex_Neutrals

  END TYPE IPE_Neutrals


  CHARACTER(250), PARAMETER      :: hwm_path     = './'

  INTEGER,    PARAMETER, PRIVATE :: N_heights    = 72
  INTEGER,    PARAMETER, PRIVATE :: N_Latitudes  = 19
  INTEGER,    PARAMETER, PRIVATE :: N_Longitudes = 36
  REAL(prec), PARAMETER, PRIVATE :: small_power  = -3.0_prec
  REAL(prec), PARAMETER, PRIVATE :: small_number = 1.0e-03_prec
  REAL(prec), PARAMETER, PRIVATE :: min_density  = 1.0e-12_prec


CONTAINS

  SUBROUTINE Build_IPE_Neutrals( neutrals, nFluxTube, NLP, NMP, mp_low, mp_high, rc )
    IMPLICIT NONE
    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    INTEGER,               INTENT(in)    :: nFluxTube
    INTEGER,               INTENT(in)    :: NLP
    INTEGER,               INTENT(in)    :: NMP
    INTEGER,               INTENT(in)    :: mp_low, mp_high
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    neutrals % nFluxTube = nFluxTube
    neutrals % NLP       = NLP
    neutrals % NMP       = NMP
    neutrals % mp_low    = mp_low
    neutrals % mp_high   = mp_high

    ALLOCATE( neutrals % helium(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % oxygen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % molecular_oxygen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % molecular_nitrogen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % nitrogen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % hydrogen(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % temperature(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % temperature_inf(1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % velocity_geographic(1:3,1:nFluxTube,1:NLP,mp_low:mp_high), &
              neutrals % velocity_apex(1:3,nFluxTube,1:NLP,mp_low:mp_high), &
              stat = stat )
    IF ( ipe_alloc_check( stat, msg="Failed to allocate neutrals internal arrays", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    neutrals % helium              = 0.0_prec
    neutrals % oxygen              = 0.0_prec
    neutrals % molecular_oxygen    = 0.0_prec
    neutrals % molecular_nitrogen  = 0.0_prec
    neutrals % nitrogen            = 0.0_prec
    neutrals % hydrogen            = 0.0_prec
    neutrals % temperature         = 0.0_prec
    neutrals % temperature_inf     = 0.0_prec
    neutrals % velocity_geographic = 0.0_prec
    neutrals % velocity_apex       = 0.0_prec

    ALLOCATE( neutrals % geo_helium(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_oxygen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_molecular_oxygen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_molecular_nitrogen(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_nitrogen(1:nlon_geo,1:nlat_geo,1:nheights_geo),&
              neutrals % geo_hydrogen(1:nlon_geo,1:nlat_geo,1:nheights_geo),&
              neutrals % geo_temperature(1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              neutrals % geo_velocity(1:3,1:nlon_geo,1:nlat_geo,1:nheights_geo), &
              stat = stat )
    IF ( ipe_alloc_check( stat, msg="Failed to allocate neutrals internal geo arrays", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    neutrals % geo_helium             = 0.0_prec
    neutrals % geo_oxygen             = 0.0_prec
    neutrals % geo_molecular_oxygen   = 0.0_prec
    neutrals % geo_molecular_nitrogen = 0.0_prec
    neutrals % geo_nitrogen           = 0.0_prec
    neutrals % geo_hydrogen           = 0.0_prec
    neutrals % geo_temperature        = 0.0_prec
    neutrals % geo_velocity           = 0.0_prec

  END SUBROUTINE Build_IPE_Neutrals


  SUBROUTINE Trash_IPE_Neutrals( neutrals, rc )

    IMPLICIT NONE

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    INTEGER :: stat

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    IF ( ASSOCIATED( neutrals % helium              ) ) &
         DEALLOCATE( neutrals % helium              , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % hydrogen            ) ) &
         DEALLOCATE( neutrals % hydrogen            , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % molecular_nitrogen  ) ) &
         DEALLOCATE( neutrals % molecular_nitrogen  , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % molecular_oxygen    ) ) &
         DEALLOCATE( neutrals % molecular_oxygen    , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % nitrogen            ) ) &
         DEALLOCATE( neutrals % nitrogen            , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % oxygen              ) ) &
         DEALLOCATE( neutrals % oxygen              , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % temperature         ) ) &
         DEALLOCATE( neutrals % temperature         , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % temperature_inf     ) ) &
         DEALLOCATE( neutrals % temperature_inf     , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % velocity_apex       ) ) &
         DEALLOCATE( neutrals % velocity_apex       , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    IF ( ASSOCIATED( neutrals % velocity_geographic ) ) &
         DEALLOCATE( neutrals % velocity_geographic , stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN


    DEALLOCATE( neutrals % geo_helium, &
                 neutrals % geo_oxygen, &
                 neutrals % geo_molecular_oxygen, &
                 neutrals % geo_molecular_nitrogen, &
                 neutrals % geo_nitrogen, &
                 neutrals % geo_hydrogen, &
                 neutrals % geo_temperature, &
                 neutrals % geo_velocity, &
                 stat=stat )
    IF ( ipe_dealloc_check( stat, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    ! Clean up file reader if it was initialized
    IF ( neutrals % file_reader % initialized ) THEN
      CALL neutrals % file_reader % Finalize( )
    ENDIF

  END SUBROUTINE Trash_IPE_Neutrals


  SUBROUTINE Update_IPE_Neutrals( neutrals, params, grid, time, forcing, mpi_layer, vertical_wind_limit, rc )

    CLASS( IPE_Neutrals         ), INTENT(inout) :: neutrals
    TYPE ( IPE_Model_Parameters ), INTENT(in   ) :: params
    TYPE ( IPE_Grid             ), INTENT(in   ) :: grid
    TYPE ( IPE_Time             ), INTENT(in   ) :: time
    TYPE ( IPE_Forcing          ), INTENT(in   ) :: forcing
    TYPE ( IPE_MPI_Layer        ), INTENT(in   ) :: mpi_layer
    REAL ( prec ),                 INTENT(in   ) :: vertical_wind_limit

    INTEGER, OPTIONAL,             INTENT(out  ) :: rc

    ! Local
    LOGICAL :: msis_switch
    INTEGER :: localrc

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    ! --- MSIS for empirical neutral atmosphere ---
    ! Run MSIS FIRST so it fills He, H, N (which FileReader doesn't provide).
    ! When read_gsm_neutrals is enabled, FileReader runs AFTER and overwrites
    ! T, O, O2, N2, winds — preventing MSIS from discarding file-based neutrals.
    msis_switch = mod(time % elapsed_sec,params % msis_time_step) == 0.0

    IF ( msis_switch .and. (time % elapsed_sec > 0._prec .or. .NOT. params % read_apex_neutrals) ) THEN
      IF( mpi_layer % rank_id == 0 ) write(6,*) 'Calling MSIS ', int(time % elapsed_sec / 60), ' Mins UT'
      CALL neutrals % IPE_Neutrals_Empirical( grid, time, forcing, rc=localrc )
      IF ( ipe_error_check( localrc, msg="call to IPE_Neutrals_Empirical failed", rc=rc ) ) RETURN
    ENDIF

    ! --- Read neutrals from GSM files if enabled ---
    ! This runs AFTER MSIS so that file-based T, O, O2, N2, winds override
    ! the MSIS/HWM values. MSIS-provided He, H, N are preserved.
    IF ( params % read_gsm_neutrals ) THEN

      ! Initialize file reader on first call
      IF ( .NOT. neutrals % file_reader % initialized ) THEN
        CALL neutrals % file_reader % Init( &
             TRIM(params % gsm_neutrals_dir), &
             time % year, time % month, time % day, &
             mpi_layer, &
             interp_method=TRIM(params % neutral_interp_method), &
             rc=localrc )
        IF ( ipe_error_check( localrc, msg="FileReader Init failed", rc=rc ) ) RETURN

        ! For high-res source grids, enable direct-to-apex interpolation
        IF ( neutrals % file_reader % needs_horiz_interp ) THEN
          CALL neutrals % file_reader % InitApexWeights( grid, mpi_layer, localrc )
          IF ( ipe_error_check( localrc, msg="FileReader InitApexWeights failed", rc=rc ) ) RETURN
        ENDIF
      ENDIF

      IF ( neutrals % file_reader % direct_apex_interp ) THEN
        ! === Direct-to-apex path (high-res source, e.g. FV3WAM 384x190) ===
        ! Interpolate source data directly to apex grid points,
        ! bypassing the 90x91 geographic grid bottleneck
        CALL neutrals % file_reader % UpdateApex( time, mpi_layer, &
             neutrals % temperature, &
             neutrals % oxygen, &
             neutrals % molecular_oxygen, &
             neutrals % molecular_nitrogen, &
             neutrals % velocity_geographic, &
             rc=localrc )
        IF ( ipe_error_check( localrc, msg="FileReader UpdateApex failed", rc=rc ) ) RETURN

        ! Diagnostic: check apex array ranges (rank 0 only)
        IF ( mpi_layer % rank_id == 0 .AND. verbose_diag ) THEN
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' Apex(direct) T    min/max: ', &
            MINVAL(neutrals % temperature(:,:,grid%mp_low:grid%mp_high)), ' / ', &
            MAXVAL(neutrals % temperature(:,:,grid%mp_low:grid%mp_high))
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' Apex(direct) O    min/max: ', &
            MINVAL(neutrals % oxygen(:,:,grid%mp_low:grid%mp_high)), ' / ', &
            MAXVAL(neutrals % oxygen(:,:,grid%mp_low:grid%mp_high))
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' Apex(direct) O2   min/max: ', &
            MINVAL(neutrals % molecular_oxygen(:,:,grid%mp_low:grid%mp_high)), ' / ', &
            MAXVAL(neutrals % molecular_oxygen(:,:,grid%mp_low:grid%mp_high))
        ENDIF

        ! Apply min_density floor matching IPE_Neutrals_Extrapolate (1e-12).
        ! VertInterpColumn uses 1e-30 floor, but for plasmaspheric kp (above
        ! ~2000 km) constant-G0 scale-height gives densities << 1e-12.
        ! FLIP requires densities >= min_density to avoid NaN in rate equations.
        WHERE (neutrals % oxygen         (1:grid%nFluxTube, 1:grid%NLP, grid%mp_low:grid%mp_high) < min_density)
          neutrals % oxygen         (1:grid%nFluxTube, 1:grid%NLP, grid%mp_low:grid%mp_high) = min_density
        END WHERE
        WHERE (neutrals % molecular_oxygen (1:grid%nFluxTube, 1:grid%NLP, grid%mp_low:grid%mp_high) < min_density)
          neutrals % molecular_oxygen (1:grid%nFluxTube, 1:grid%NLP, grid%mp_low:grid%mp_high) = min_density
        END WHERE
        WHERE (neutrals % molecular_nitrogen(1:grid%nFluxTube, 1:grid%NLP, grid%mp_low:grid%mp_high) < min_density)
          neutrals % molecular_nitrogen(1:grid%nFluxTube, 1:grid%NLP, grid%mp_low:grid%mp_high) = min_density
        END WHERE

        ! Skip Geographic_to_Apex_Neutrals — data is already on apex grid

      ELSE
        ! === Legacy path (90x91 GSM files) ===
        ! Read and interpolate GSM data to geographic grid
        CALL neutrals % file_reader % Update( time, mpi_layer, &
             neutrals % geo_temperature, &
             neutrals % geo_oxygen, &
             neutrals % geo_molecular_oxygen, &
             neutrals % geo_molecular_nitrogen, &
             neutrals % geo_velocity, &
             grid % altitude_geo, &
             forcing=forcing, &
             rc=localrc )
        IF ( ipe_error_check( localrc, msg="FileReader Update failed", rc=rc ) ) RETURN

        ! Diagnostic: check geo array ranges (rank 0 only)
        IF ( mpi_layer % rank_id == 0 .AND. verbose_diag ) THEN
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' FileReader geo_T    min/max: ', &
            MINVAL(neutrals % geo_temperature), ' / ', MAXVAL(neutrals % geo_temperature)
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' FileReader geo_O    min/max: ', &
            MINVAL(neutrals % geo_oxygen), ' / ', MAXVAL(neutrals % geo_oxygen)
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' FileReader geo_O2   min/max: ', &
            MINVAL(neutrals % geo_molecular_oxygen), ' / ', MAXVAL(neutrals % geo_molecular_oxygen)
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' FileReader geo_N2   min/max: ', &
            MINVAL(neutrals % geo_molecular_nitrogen), ' / ', MAXVAL(neutrals % geo_molecular_nitrogen)
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' FileReader geo_Ue   min/max: ', &
            MINVAL(neutrals % geo_velocity(1,:,:,:)), ' / ', MAXVAL(neutrals % geo_velocity(1,:,:,:))
        ENDIF

        ! Interpolate geographic grid -> apex grid
        CALL neutrals % Geographic_to_Apex_Neutrals( grid )

        ! Diagnostic: check apex array ranges after interpolation (rank 0 only)
        IF ( mpi_layer % rank_id == 0 .AND. verbose_diag ) THEN
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' Apex T    min/max: ', &
            MINVAL(neutrals % temperature(:,:,grid%mp_low:grid%mp_high)), ' / ', &
            MAXVAL(neutrals % temperature(:,:,grid%mp_low:grid%mp_high))
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' Apex O    min/max: ', &
            MINVAL(neutrals % oxygen(:,:,grid%mp_low:grid%mp_high)), ' / ', &
            MAXVAL(neutrals % oxygen(:,:,grid%mp_low:grid%mp_high))
          WRITE(6,'(A,ES12.4,A,ES12.4)') ' Apex O2   min/max: ', &
            MINVAL(neutrals % molecular_oxygen(:,:,grid%mp_low:grid%mp_high)), ' / ', &
            MAXVAL(neutrals % molecular_oxygen(:,:,grid%mp_low:grid%mp_high))
        ENDIF

      ENDIF  ! direct_apex_interp

      ! Apply neutral scaling factors (from NeutralFileIO namelist)
      IF ( params % neutral_T_scale /= 1.0_prec ) &
        neutrals % temperature(:,:,grid%mp_low:grid%mp_high) = &
          neutrals % temperature(:,:,grid%mp_low:grid%mp_high) * params % neutral_T_scale
      IF ( params % neutral_O_scale /= 1.0_prec ) &
        neutrals % oxygen(:,:,grid%mp_low:grid%mp_high) = &
          neutrals % oxygen(:,:,grid%mp_low:grid%mp_high) * params % neutral_O_scale
      IF ( params % neutral_O2_scale /= 1.0_prec ) &
        neutrals % molecular_oxygen(:,:,grid%mp_low:grid%mp_high) = &
          neutrals % molecular_oxygen(:,:,grid%mp_low:grid%mp_high) * params % neutral_O2_scale
      IF ( params % neutral_N2_scale /= 1.0_prec ) &
        neutrals % molecular_nitrogen(:,:,grid%mp_low:grid%mp_high) = &
          neutrals % molecular_nitrogen(:,:,grid%mp_low:grid%mp_high) * params % neutral_N2_scale

      IF ( verbose_diag .AND. mpi_layer % rank_id == 0 .AND. &
           mod(time % elapsed_sec, 3600.0_prec) == 0.0_prec ) &
        WRITE(6,'(A,4F8.3)') ' Neutral scaling T/O/O2/N2: ', &
          params % neutral_T_scale, params % neutral_O_scale, &
          params % neutral_O2_scale, params % neutral_N2_scale

    ENDIF

    IF ( params % read_gsm_neutrals .AND. &
         neutrals % file_reader % direct_apex_interp ) THEN
      ! Direct-apex path: FileReader already handles above-top extrapolation
      ! via VertInterpColumn (isothermal for T/winds, scale-height for densities).
      ! Only compute temperature_inf; skip overwriting T/O/O2/N2/winds.
      CALL neutrals % IPE_Neutrals_SetTempInf( grid )
    ELSE
      CALL neutrals % IPE_Neutrals_Extrapolate( grid, forcing )
    ENDIF

    CALL neutrals % Geographic_to_Apex_Velocity( grid, vertical_wind_limit )

  END SUBROUTINE Update_IPE_Neutrals


  SUBROUTINE IPE_Neutrals_Extrapolate( neutrals, grid, forcing )

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE ( IPE_Grid     ), INTENT(in   ) :: grid
    TYPE ( IPE_Forcing  ), INTENT(in   ) :: forcing

    ! Local
    INTEGER    :: kp, kpp, kpStart, kpStep, kpStop, lp, mp
    REAL(prec) :: r, w
    REAL(prec) :: t( grid % mp_low:grid % mp_high )

    REAL(prec), PARAMETER :: fac = 2.0_prec * G0 / GSCON

    ! -- NOTE: requires correct northern_top_index and southern_top_index arrays

    DO lp = 1, grid % NLP

      ! loop over hemispheres, starting with northern one
      kpStart = grid % northern_top_index(lp)
      kpStop  = grid % flux_tube_midpoint(lp)

      DO kpStep = 1, -1, -2

        DO kp = kpStart + kpStep, kpStop, kpStep

          ! extend temperature
          neutrals % temperature(kp,lp,grid % mp_low:grid % mp_high) = &
            neutrals % temperature(kpStart,lp,grid % mp_low:grid % mp_high)

          ! extend winds
          neutrals % velocity_geographic(:,kp,lp,grid % mp_low:grid % mp_high) = &
            neutrals % velocity_geographic(:,kpStart,lp,grid % mp_low:grid % mp_high)

          r = 0.0_prec
          t = 0.0_prec
          DO kpp = kp - kpStep, kp, kpStep
            w = 1.0_prec + grid % altitude(kpp, lp) / earth_radius
            t = t + w * w * neutrals % temperature(kpp,lp,:)
            r = grid % altitude(kpp, lp) - r
          ENDDO
          t = -fac * r / t

          kpp = kp - kpStep
          ! extrapolating atomic oxygen
          neutrals % oxygen(kp,lp,grid % mp_low:grid % mp_high) = &
            max( min_density, neutrals % oxygen(kpp,lp,grid % mp_low:grid % mp_high) * exp( O_mass * t ) )
          ! extrapolating molecular oxygen
          neutrals % molecular_oxygen(kp,lp,grid % mp_low:grid % mp_high) = &
            max( min_density, neutrals % molecular_oxygen(kpp,lp,grid % mp_low:grid % mp_high) * exp( O2_mass * t ) )
          ! extrapolating molecular nitrogen
          neutrals % molecular_nitrogen(kp,lp,grid % mp_low:grid % mp_high) = &
            max( min_density, neutrals % molecular_nitrogen(kpp,lp,grid % mp_low:grid % mp_high) * exp( N2_mass * t ) )

        ENDDO
        ! preparing for the southern hemisphere
        kpStart = grid % southern_top_index(lp)
        kpStop  = kpStop + 1
      ENDDO
    ENDDO

    ! set exospheric temperature
    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        neutrals % temperature_inf(1:grid % flux_tube_midpoint(lp),lp,mp) = &
          neutrals % temperature(grid % northern_top_index(lp),lp,mp)
        neutrals % temperature_inf(grid % flux_tube_midpoint(lp)+1:grid % flux_tube_max(lp),lp,mp) = &
          neutrals % temperature(grid % southern_top_index(lp),lp,mp)
      ENDDO
    ENDDO

  END SUBROUTINE IPE_Neutrals_Extrapolate


  !---------------------------------------------------------------------------
  ! IPE_Neutrals_SetTempInf: Set exospheric temperature only.
  ! Used when FileReader direct-apex path provides all other neutrals.
  !---------------------------------------------------------------------------
  SUBROUTINE IPE_Neutrals_SetTempInf( neutrals, grid )

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ),      INTENT(in)    :: grid

    INTEGER :: lp, mp

    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        neutrals % temperature_inf(1:grid % flux_tube_midpoint(lp),lp,mp) = &
          neutrals % temperature(grid % northern_top_index(lp),lp,mp)
        neutrals % temperature_inf(grid % flux_tube_midpoint(lp)+1:grid % flux_tube_max(lp),lp,mp) = &
          neutrals % temperature(grid % southern_top_index(lp),lp,mp)
      ENDDO
    ENDDO

  END SUBROUTINE IPE_Neutrals_SetTempInf


  SUBROUTINE IPE_Neutrals_Empirical( neutrals, grid, time, forcing, rc )
  !
  ! Usage :
  !
  !   CALL neutrals % Update( grid, utime, year, day, f107d, f107a, ap )
  ! ================================================================================================== !

    IMPLICIT NONE

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ),      INTENT(in)    :: grid
    TYPE( IPE_Time ),      INTENT(in)    :: time
    TYPE( IPE_Forcing ),   INTENT(in)    :: forcing
    INTEGER, OPTIONAL,     INTENT(out)   :: rc

    ! Local
    INTEGER    :: localrc
    INTEGER    :: kp, lp, mp, day
    REAL(prec) :: geo_alt, geo_lat, geo_lon, utime
    REAL(prec), DIMENSION(7) :: AP

    INTEGER(4)            :: iyd
    REAL(4)               :: hwm_sec, hwm_f107d, hwm_f107a, hwm_alt, hwm_lat, hwm_lon
    REAL(4), DIMENSION(2) :: hwm_ap, w

    REAL(msis_dp)               :: msis_alt, msis_f107d, msis_f107a, msis_lat, msis_lon, msis_sec, msis_stl
    REAL(msis_dp), DIMENSION(7) :: msis_ap
    REAL(msis_dp), DIMENSION(2) :: temperatures
    REAL(msis_dp), DIMENSION(9) :: densities

    INTEGER, PARAMETER    :: msis_mass = 48


    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    AP = forcing % GetAP( )

    w         = 0.0
    hwm_ap    = REAL(AP(1:2), KIND=4)
    hwm_f107a = REAL(forcing % f107_81day_avg( forcing % current_index ), KIND=4)
    hwm_f107d = REAL(forcing % f107( forcing % current_index ),           KIND=4)

    msis_ap      = REAL(AP,    KIND=msis_dp)
    msis_f107a   = REAL(forcing % f107_81day_avg( forcing % current_index ), KIND=msis_dp)
    msis_f107d   = REAL(forcing % f107( forcing % current_index ),           KIND=msis_dp)
    densities    = 0.0_msis_dp
    temperatures = 0.0_msis_dp

    iyd = 99000 + time % day_of_year    ! Input, year and day as yyddd

    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        DO kp = 1, grid % flux_tube_max(lp)
          geo_lon = rtd * grid % longitude(kp,lp,mp)
          geo_lat = 90.0 - rtd * grid % colatitude(kp,lp,mp)
          geo_alt = m_to_km * grid % altitude(kp,lp)

          ! -- horizontal & vertical wind
          iF ( .NOT. forcing % coupled ) THEN

            hwm_sec = REAL(time % utime, KIND=4)
            hwm_alt = geo_alt
            hwm_lat = geo_lat
            hwm_lon = geo_lon
            w       = 0.0

            call hwm14( iyd,       &    ! Input, year and day as yyddd
                        hwm_sec,   &    ! Input, universal time ( sec )
                        hwm_alt,   &    ! Input, altitude ( km )
                        hwm_lat,   &    ! Input, geodetic latitude ( degrees )
                        hwm_lon,   &    ! Input, geodetic longitude ( degrees )
                        0.0,       &    ! Input, local apparent solar time ( hrs )[ not used ]
                        hwm_f107a, &    ! Input, 3 month average of f10.7 flux [ not used ]
                        hwm_f107d, &    ! Input, daily average of f10.7 flux for the previous day [ not used ]
                        hwm_ap,    &    ! Input, magnetic index ( daily ), current 3hr ap index
                        hwm_path,  &    ! Input, default datafile path
                        w,         &    ! Ouput, neutral wind velocity meridional-northwards(1) and zonal-eastwards(2) components
                        localrc )       ! Ouput, return code
            IF ( ipe_error_check( localrc, msg="call to hwm14 failed", &
              line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

            neutrals % velocity_geographic(1,kp,lp,mp) = w(2)
            neutrals % velocity_geographic(2,kp,lp,mp) = w(1)
            neutrals % velocity_geographic(3,kp,lp,mp) = 0.0_prec

          ENDIF

          ! -- composition & temperature
          msis_sec     = REAL(time % utime, KIND=msis_dp)
          msis_alt     = geo_alt
          msis_lat     = geo_lat
          msis_lon     = geo_lon
          msis_stl     = REAL(time % utime / 3600.0_prec + geo_lon / 15.0_prec, KIND=msis_dp)
          densities    = 0.0_msis_dp
          temperatures = 0.0_msis_dp

          call gtd7( iyd,         &    ! Input, year and day as yyddd
                     msis_sec,    &    ! Input, universal time ( sec )
                     msis_alt,    &    ! Input, altitude ( km )
                     msis_lat,    &    ! Input, geodetic latitude ( degrees )
                     msis_lon,    &    ! Input, geodetic longitude ( degrees )
                     msis_stl,    &    ! Input, local apparent solar time ( hrs )
                     msis_f107a,  &    ! Input, 3 month average of f10.7 flux
                     msis_f107d,  &    ! Input, daily average of f10.7 flux for the previous day
                     msis_ap,     &    ! Input, magnetic index ( daily ), current, 3,6,9hrs prior 3hr ap index, 12-33 hr prior ap average, 36-57 hr prior ap average
                     msis_mass,   &    ! Mass number ( see src/msis/physics_msis.f90 for more details )
                     densities,   &    ! Ouput, neutral densities in cm-3
                     temperatures )    ! Output, exospheric temperature and temperature at altitude

          IF ( .NOT. forcing % coupled ) THEN

            neutrals % temperature_inf(kp,lp,mp) = temperatures(1)
            neutrals % temperature(kp,lp,mp)     = temperatures(2)

            neutrals % oxygen(kp,lp,mp)             = cm_3_to_m_3 * densities(2)
            neutrals % molecular_nitrogen(kp,lp,mp) = cm_3_to_m_3 * densities(3)
            neutrals % molecular_oxygen(kp,lp,mp)   = cm_3_to_m_3 * densities(4)

          ENDIF

          neutrals % helium(kp,lp,mp)   = cm_3_to_m_3 * densities(1)
          neutrals % hydrogen(kp,lp,mp) = cm_3_to_m_3 * densities(7)
          neutrals % nitrogen(kp,lp,mp) = cm_3_to_m_3 * densities(8)

        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE IPE_Neutrals_Empirical


  SUBROUTINE Geographic_to_Apex_Velocity( neutrals, grid, vertical_wind_limit )

    IMPLICIT NONE

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ),      INTENT(in)    :: grid
    REAL(prec), INTENT(in) :: vertical_wind_limit

    ! Local
    INTEGER    :: i_D_vec, kp, lp, mp
    REAL(prec) :: dotprod
    REAL(prec) :: neutral_vertical_velocity

    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        DO kp = 1, grid % flux_tube_max(lp)

          DO i_D_vec = 1, 3 ! D1, D2 and D3

             dotprod = grid % apex_d_vectors(1,i_D_vec,kp,lp,mp)*grid % apex_d_vectors(1,i_D_vec,kp,lp,mp) + &
                       grid % apex_d_vectors(2,i_D_vec,kp,lp,mp)*grid % apex_d_vectors(2,i_D_vec,kp,lp,mp) + &
                       grid % apex_d_vectors(3,i_D_vec,kp,lp,mp)*grid % apex_d_vectors(3,i_D_vec,kp,lp,mp)

             IF ( dotprod > 0.0_prec ) THEN

               if (neutrals % velocity_geographic(3,kp,lp,mp).gt.vertical_wind_limit) then
                 neutral_vertical_velocity = vertical_wind_limit
               else if (neutrals % velocity_geographic(3,kp,lp,mp).lt.0.0 - vertical_wind_limit) then
                 neutral_vertical_velocity = 0.0 - vertical_wind_limit
               else
                 neutral_vertical_velocity = neutrals % velocity_geographic(3,kp,lp,mp)
               endif

               neutrals % velocity_apex(i_D_vec,kp,lp,mp) = &
                 ( grid % apex_d_vectors(1,i_D_vec,kp,lp,mp)*neutrals % velocity_geographic(1,kp,lp,mp) + &
                   grid % apex_d_vectors(2,i_D_vec,kp,lp,mp)*neutrals % velocity_geographic(2,kp,lp,mp) + &
                   grid % apex_d_vectors(3,i_D_vec,kp,lp,mp)*neutral_vertical_velocity ) &
                 / SQRT( dotprod )

             ELSE

                neutrals % velocity_apex(i_D_vec,kp,lp,mp) = 0.0_prec

             END IF

           ENDDO

           ! dbg20110131 : the midpoint values become NaN otherwise because
           ! of inappropriate D1/3 values...
           IF ( lp >= 1 .AND. lp <= 6 )THEN
             IF( kp == grid % flux_tube_midpoint(lp) )THEN
               neutrals % velocity_apex(1:3,kp,lp,mp) = neutrals % velocity_apex(1:3,kp-1,lp,mp)
             ENDIF
           ENDIF

        END DO
      END DO
    END DO

  END SUBROUTINE Geographic_to_Apex_Velocity


  !---------------------------------------------------------------------------
  ! Geographic_to_Apex_Neutrals: Interpolate neutral scalar fields and winds
  ! from the regular geographic grid (geo_*) to apex grid points using
  ! trilinear interpolation.
  !---------------------------------------------------------------------------
  SUBROUTINE Geographic_to_Apex_Neutrals( neutrals, grid )

    IMPLICIT NONE

    CLASS( IPE_Neutrals ), INTENT(inout) :: neutrals
    TYPE( IPE_Grid ),      INTENT(in)    :: grid

    ! Local
    INTEGER    :: kp, lp, mp
    REAL(prec) :: geo_lon, geo_lat, geo_alt
    REAL(prec) :: fi, fj, fk
    INTEGER    :: i0, j0, k0, i1, j1, k1
    REAL(prec) :: wi, wj, wk
    REAL(prec) :: dlon, dlat, dalt
    REAL(prec) :: c000, c100, c010, c110, c001, c101, c011, c111
    REAL(prec) :: log000, log100, log010, log110, log001, log101, log011, log111

    ! Geographic grid spacing
    dlon = 360.0_prec / REAL(nlon_geo, prec)   ! 4 degrees
    dlat = 180.0_prec / REAL(nlat_geo - 1, prec)  ! 2 degrees
    dalt = 5.0_prec  ! km

    DO mp = grid % mp_low, grid % mp_high
      DO lp = 1, grid % NLP
        DO kp = 1, grid % flux_tube_max(lp)

          ! Get geographic coordinates of this apex grid point
          geo_lon = rtd * grid % longitude(kp,lp,mp)          ! degrees (0-360)
          geo_lat = 90.0_prec - rtd * grid % colatitude(kp,lp,mp)  ! degrees (-90 to 90)
          geo_alt = m_to_km * grid % altitude(kp,lp)          ! km

          ! Fractional indices into regular geographic grid
          ! lon: 0 to 356 at 4 deg spacing
          fi = geo_lon / dlon
          ! lat: -90 to 90 at 2 deg spacing
          fj = (geo_lat + 90.0_prec) / dlat
          ! alt: 90 to 1000 at 5 km spacing — clamp to grid range
          ! (out-of-range points get nearest boundary value;
          !  IPE_Neutrals_Extrapolate will properly extend along field lines)
          fk = (MAX(90.0_prec, MIN(geo_alt, 1000.0_prec)) - 90.0_prec) / dalt

          ! Integer indices (1-based)
          i0 = INT(fi) + 1
          j0 = INT(fj) + 1
          k0 = INT(fk) + 1

          ! Weights (fractional part)
          wi = fi - REAL(i0 - 1, prec)
          wj = fj - REAL(j0 - 1, prec)
          wk = fk - REAL(k0 - 1, prec)

          ! Neighbor indices with bounds checking
          i1 = MOD(i0, nlon_geo) + 1  ! wrap longitude
          j0 = MAX(1, MIN(j0, nlat_geo - 1))
          j1 = j0 + 1
          k0 = MAX(1, MIN(k0, nheights_geo - 1))
          k1 = k0 + 1

          ! Clamp i0 to valid range
          i0 = MAX(1, MIN(i0, nlon_geo))

          ! --- Temperature: trilinear interpolation ---
          c000 = neutrals % geo_temperature(i0,j0,k0)
          c100 = neutrals % geo_temperature(i1,j0,k0)
          c010 = neutrals % geo_temperature(i0,j1,k0)
          c110 = neutrals % geo_temperature(i1,j1,k0)
          c001 = neutrals % geo_temperature(i0,j0,k1)
          c101 = neutrals % geo_temperature(i1,j0,k1)
          c011 = neutrals % geo_temperature(i0,j1,k1)
          c111 = neutrals % geo_temperature(i1,j1,k1)
          neutrals % temperature(kp,lp,mp) = &
            (1.0_prec-wk)*((1.0_prec-wj)*((1.0_prec-wi)*c000 + wi*c100) &
                          +           wj *((1.0_prec-wi)*c010 + wi*c110)) &
           +          wk *((1.0_prec-wj)*((1.0_prec-wi)*c001 + wi*c101) &
                          +           wj *((1.0_prec-wi)*c011 + wi*c111))

          ! --- Oxygen: trilinear in log-space ---
          log000 = LOG(MAX(neutrals % geo_oxygen(i0,j0,k0), 1.0e-30_prec))
          log100 = LOG(MAX(neutrals % geo_oxygen(i1,j0,k0), 1.0e-30_prec))
          log010 = LOG(MAX(neutrals % geo_oxygen(i0,j1,k0), 1.0e-30_prec))
          log110 = LOG(MAX(neutrals % geo_oxygen(i1,j1,k0), 1.0e-30_prec))
          log001 = LOG(MAX(neutrals % geo_oxygen(i0,j0,k1), 1.0e-30_prec))
          log101 = LOG(MAX(neutrals % geo_oxygen(i1,j0,k1), 1.0e-30_prec))
          log011 = LOG(MAX(neutrals % geo_oxygen(i0,j1,k1), 1.0e-30_prec))
          log111 = LOG(MAX(neutrals % geo_oxygen(i1,j1,k1), 1.0e-30_prec))
          neutrals % oxygen(kp,lp,mp) = EXP( &
            (1.0_prec-wk)*((1.0_prec-wj)*((1.0_prec-wi)*log000 + wi*log100) &
                          +           wj *((1.0_prec-wi)*log010 + wi*log110)) &
           +          wk *((1.0_prec-wj)*((1.0_prec-wi)*log001 + wi*log101) &
                          +           wj *((1.0_prec-wi)*log011 + wi*log111)) )

          ! --- Molecular oxygen: trilinear in log-space ---
          log000 = LOG(MAX(neutrals % geo_molecular_oxygen(i0,j0,k0), 1.0e-30_prec))
          log100 = LOG(MAX(neutrals % geo_molecular_oxygen(i1,j0,k0), 1.0e-30_prec))
          log010 = LOG(MAX(neutrals % geo_molecular_oxygen(i0,j1,k0), 1.0e-30_prec))
          log110 = LOG(MAX(neutrals % geo_molecular_oxygen(i1,j1,k0), 1.0e-30_prec))
          log001 = LOG(MAX(neutrals % geo_molecular_oxygen(i0,j0,k1), 1.0e-30_prec))
          log101 = LOG(MAX(neutrals % geo_molecular_oxygen(i1,j0,k1), 1.0e-30_prec))
          log011 = LOG(MAX(neutrals % geo_molecular_oxygen(i0,j1,k1), 1.0e-30_prec))
          log111 = LOG(MAX(neutrals % geo_molecular_oxygen(i1,j1,k1), 1.0e-30_prec))
          neutrals % molecular_oxygen(kp,lp,mp) = EXP( &
            (1.0_prec-wk)*((1.0_prec-wj)*((1.0_prec-wi)*log000 + wi*log100) &
                          +           wj *((1.0_prec-wi)*log010 + wi*log110)) &
           +          wk *((1.0_prec-wj)*((1.0_prec-wi)*log001 + wi*log101) &
                          +           wj *((1.0_prec-wi)*log011 + wi*log111)) )

          ! --- Molecular nitrogen: trilinear in log-space ---
          log000 = LOG(MAX(neutrals % geo_molecular_nitrogen(i0,j0,k0), 1.0e-30_prec))
          log100 = LOG(MAX(neutrals % geo_molecular_nitrogen(i1,j0,k0), 1.0e-30_prec))
          log010 = LOG(MAX(neutrals % geo_molecular_nitrogen(i0,j1,k0), 1.0e-30_prec))
          log110 = LOG(MAX(neutrals % geo_molecular_nitrogen(i1,j1,k0), 1.0e-30_prec))
          log001 = LOG(MAX(neutrals % geo_molecular_nitrogen(i0,j0,k1), 1.0e-30_prec))
          log101 = LOG(MAX(neutrals % geo_molecular_nitrogen(i1,j0,k1), 1.0e-30_prec))
          log011 = LOG(MAX(neutrals % geo_molecular_nitrogen(i0,j1,k1), 1.0e-30_prec))
          log111 = LOG(MAX(neutrals % geo_molecular_nitrogen(i1,j1,k1), 1.0e-30_prec))
          neutrals % molecular_nitrogen(kp,lp,mp) = EXP( &
            (1.0_prec-wk)*((1.0_prec-wj)*((1.0_prec-wi)*log000 + wi*log100) &
                          +           wj *((1.0_prec-wi)*log010 + wi*log110)) &
           +          wk *((1.0_prec-wj)*((1.0_prec-wi)*log001 + wi*log101) &
                          +           wj *((1.0_prec-wi)*log011 + wi*log111)) )

          ! --- Winds: trilinear interpolation (3 components) ---
          ! Eastward
          c000 = neutrals % geo_velocity(1,i0,j0,k0)
          c100 = neutrals % geo_velocity(1,i1,j0,k0)
          c010 = neutrals % geo_velocity(1,i0,j1,k0)
          c110 = neutrals % geo_velocity(1,i1,j1,k0)
          c001 = neutrals % geo_velocity(1,i0,j0,k1)
          c101 = neutrals % geo_velocity(1,i1,j0,k1)
          c011 = neutrals % geo_velocity(1,i0,j1,k1)
          c111 = neutrals % geo_velocity(1,i1,j1,k1)
          neutrals % velocity_geographic(1,kp,lp,mp) = &
            (1.0_prec-wk)*((1.0_prec-wj)*((1.0_prec-wi)*c000 + wi*c100) &
                          +           wj *((1.0_prec-wi)*c010 + wi*c110)) &
           +          wk *((1.0_prec-wj)*((1.0_prec-wi)*c001 + wi*c101) &
                          +           wj *((1.0_prec-wi)*c011 + wi*c111))

          ! Northward
          c000 = neutrals % geo_velocity(2,i0,j0,k0)
          c100 = neutrals % geo_velocity(2,i1,j0,k0)
          c010 = neutrals % geo_velocity(2,i0,j1,k0)
          c110 = neutrals % geo_velocity(2,i1,j1,k0)
          c001 = neutrals % geo_velocity(2,i0,j0,k1)
          c101 = neutrals % geo_velocity(2,i1,j0,k1)
          c011 = neutrals % geo_velocity(2,i0,j1,k1)
          c111 = neutrals % geo_velocity(2,i1,j1,k1)
          neutrals % velocity_geographic(2,kp,lp,mp) = &
            (1.0_prec-wk)*((1.0_prec-wj)*((1.0_prec-wi)*c000 + wi*c100) &
                          +           wj *((1.0_prec-wi)*c010 + wi*c110)) &
           +          wk *((1.0_prec-wj)*((1.0_prec-wi)*c001 + wi*c101) &
                          +           wj *((1.0_prec-wi)*c011 + wi*c111))

          ! Upward
          c000 = neutrals % geo_velocity(3,i0,j0,k0)
          c100 = neutrals % geo_velocity(3,i1,j0,k0)
          c010 = neutrals % geo_velocity(3,i0,j1,k0)
          c110 = neutrals % geo_velocity(3,i1,j1,k0)
          c001 = neutrals % geo_velocity(3,i0,j0,k1)
          c101 = neutrals % geo_velocity(3,i1,j0,k1)
          c011 = neutrals % geo_velocity(3,i0,j1,k1)
          c111 = neutrals % geo_velocity(3,i1,j1,k1)
          neutrals % velocity_geographic(3,kp,lp,mp) = &
            (1.0_prec-wk)*((1.0_prec-wj)*((1.0_prec-wi)*c000 + wi*c100) &
                          +           wj *((1.0_prec-wi)*c010 + wi*c110)) &
           +          wk *((1.0_prec-wj)*((1.0_prec-wi)*c001 + wi*c101) &
                          +           wj *((1.0_prec-wi)*c011 + wi*c111))

        ENDDO  ! kp
      ENDDO    ! lp
    ENDDO      ! mp

  END SUBROUTINE Geographic_to_Apex_Neutrals


END MODULE IPE_Neutrals_Class
