MODULE IPE_Model_Class

  USE IPE_Precision
  USE IPE_Model_Parameters_Class
  USE IPE_Time_Class
  USE IPE_Grid_Class
  USE IPE_Neutrals_Class
  USE IPE_Forcing_Class
  USE IPE_Plasma_Class
  USE IPE_Electrodynamics_Class
  USE IPE_MPI_Layer_Class
  USE ipe_error_module

  USE COMIO

  USE IPE_Constants_Dictionary, ONLY: verbose_diag   ! gates diagnostic output

  IMPLICIT NONE


  ! The IPE_Model serves as a wrapper for all of the underlying attributes.
  ! This class should be used to orchestrate model setup, updating, and
  ! breakdown. At this level, we define the API for interacting with the
  ! deeper attributes within IPE.

  TYPE IPE_Model

    TYPE( IPE_Time )              :: time_tracker
    TYPE( IPE_Model_Parameters )  :: parameters
    TYPE( IPE_Grid )              :: grid
    TYPE( IPE_Forcing )           :: forcing
    TYPE( IPE_Neutrals )          :: neutrals
    TYPE( IPE_Plasma )            :: plasma
    TYPE( IPE_Electrodynamics )   :: eldyn
    TYPE( IPE_MPI_Layer )         :: mpi_layer
    CLASS( COMIO_T ), ALLOCATABLE :: io

    CONTAINS

      PROCEDURE :: Build => Build_IPE_Model
      PROCEDURE :: Trash => Trash_IPE_Model

      PROCEDURE :: Update     => Update_IPE_Model
      PROCEDURE :: Initialize => Initialize_IPE_Model
      PROCEDURE :: Write      => Write_IPE_State
      PROCEDURE :: Read       => Read_IPE_State

  END TYPE IPE_Model

CONTAINS

  SUBROUTINE Build_IPE_Model( ipe, mpi_comm, rc )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    INTEGER, OPTIONAL,  INTENT(in)    :: mpi_comm
    INTEGER, OPTIONAL,  INTENT(out)   :: rc

    ! Local
    INTEGER :: localrc

    IF (PRESENT(rc)) rc = IPE_SUCCESS

    CALL ipe % mpi_layer % Initialize( comm = mpi_comm )

    CALL ipe % parameters % Build( ipe % mpi_layer, rc=localrc )
    IF ( ipe_error_check( localrc, line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    ! Initialize I/O
    IF ( ipe % mpi_layer % enabled ) THEN
      call COMIO_Create(ipe % io, COMIO_FMT_HDF5, &
                        comm=ipe % mpi_layer % mpi_communicator, &
                        info=ipe % mpi_layer % mpi_info, &
                        rc=localrc)
      IF ( ipe_error_check( localrc, &
        msg="Failed to initialize I/O layer", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    ELSE
      call COMIO_Create(ipe % io, COMIO_FMT_HDF5, &
                        rc=localrc)
      IF ( ipe_error_check( localrc, &
        msg="Failed to initialize I/O layer", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    END IF

    ! Initialize the clock
    CALL ipe % time_tracker % Build( ipe % parameters % initial_timestamp, ipe % mpi_layer % rank_id )

    ! ////// Forcing ////// !

    CALL ipe % forcing % Build( ipe % parameters, rc=localrc )
    IF ( ipe_error_check( localrc, msg="Failed to initialize model forcing", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    ! ////// grid ////// !

    CALL ipe % grid % Build( ipe % io, ipe % mpi_layer, ipe % parameters, "IPE_Grid.h5", rc=localrc )
    IF ( ipe_error_check( localrc, msg="Failed to initialize model grid", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    ! ////// neutrals ////// !

    CALL ipe % neutrals % Build( nFluxtube = ipe % grid % nFluxTube, &
                                 NLP       = ipe % grid % NLP, &
                                 NMP       = ipe % grid % NMP, &
                                 mp_low    = ipe % mpi_layer % mp_low, &
                                 mp_high   = ipe % mpi_layer % mp_high, &
                                 rc        = localrc )
    IF ( ipe_error_check( localrc, msg="Failed to initialize neutrals", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    ! ////// electric field /////// !

    CALL ipe % eldyn % Build( nFluxTube = ipe % grid % nFluxTube,&
                              NLP       = ipe % grid % NLP, &
                              NMP       = ipe % grid % NMP, &
                              dynamo    = ipe % parameters % dynamo_efield, &
                              mp_low    = ipe % mpi_layer % mp_low, &
                              mp_high   = ipe % mpi_layer % mp_high, &
                              halo      = ipe % mpi_layer % mp_halo_size, &
                              rc        = localrc )
    IF ( ipe_error_check( localrc, msg="Failed to initialize electrodynamics", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN


    ! ////// plasma ////// !

    CALL ipe % plasma % Build( nFluxTube = ipe % grid % nFluxTube, &
                               NLP       = ipe % grid % NLP, &
                               NMP       = ipe % grid % NMP, &
                               mp_low    = ipe % mpi_layer % mp_low, &
                               mp_high   = ipe % mpi_layer % mp_high, &
                               halo      = ipe % mpi_layer % mp_halo_size, &
                               rc        = localrc )
    IF ( ipe_error_check( localrc, msg="Failed to initialize plasma", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

  END SUBROUTINE Build_IPE_Model


  SUBROUTINE Trash_IPE_Model( ipe, rc )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    INTEGER, OPTIONAL,  INTENT(out)   :: rc

    INTEGER :: localrc

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    CALL ipe % forcing  % Trash( rc=localrc )
    IF ( ipe_error_check( localrc, msg="Failed to shutdown forcing", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    CALL ipe % grid     % Trash( rc=localrc )
    IF ( ipe_error_check( localrc, msg="Failed to shutdown grid", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    CALL ipe % neutrals % Trash( rc=localrc )
    IF ( ipe_error_check( localrc, msg="Failed to shutdown neutrals", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    CALL ipe % plasma   % Trash( rc=localrc )
    IF ( ipe_error_check( localrc, msg="Failed to shutdown plasma", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN
    CALL ipe % eldyn    % Trash( rc=localrc )
    IF ( ipe_error_check( localrc, msg="Failed to shutdown electrodynamics", &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    CALL ipe % io       % Shutdown( )

    CALL ipe % mpi_layer % Finalize( )

  END SUBROUTINE Trash_IPE_Model


  SUBROUTINE Update_IPE_Model( ipe, t0, t1, rc )

    IMPLICIT NONE

    CLASS( IPE_Model ),   INTENT(inout) :: ipe
    REAL(prec), OPTIONAL, INTENT(in)    :: t0, t1
    INTEGER, OPTIONAL,    INTENT(out)   :: rc

    ! Local
    INTEGER    :: i, nSteps, localrc

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

    IF( PRESENT( t0 ) .AND. PRESENT(t1) )THEN

      nSteps = INT((t1-t0)/ipe % parameters % time_step)

    ELSE

      nSteps = 1

    ENDIF

    DO i = 1, nSteps

      call ipe % forcing % Update_Current_Index( ipe % parameters, &
                                                 ipe % time_tracker % elapsed_sec, &
                                                 rc = localrc )
      IF ( ipe_error_check( localrc, msg="Failed to update driver index", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      CALL ipe % neutrals % Update( ipe % parameters, &
                                    ipe % grid, &
                                    ipe % time_tracker, &
                                    ipe % forcing, &
                                    ipe % mpi_layer, &
                                    ipe % parameters % vertical_wind_limit, &
                                    rc = localrc )
      IF ( ipe_error_check( localrc, msg="Failed to update neutrals", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      ! NaN check on neutrals
      IF ( verbose_diag ) THEN
      block
        integer :: n_nan_T, n_nan_O, n_nan_O2, n_nan_N2, n_nan_vel
        n_nan_T = count(ipe%neutrals%temperature(:,:,ipe%grid%mp_low:ipe%grid%mp_high) /= &
                        ipe%neutrals%temperature(:,:,ipe%grid%mp_low:ipe%grid%mp_high))
        n_nan_O = count(ipe%neutrals%oxygen(:,:,ipe%grid%mp_low:ipe%grid%mp_high) /= &
                        ipe%neutrals%oxygen(:,:,ipe%grid%mp_low:ipe%grid%mp_high))
        n_nan_O2 = count(ipe%neutrals%molecular_oxygen(:,:,ipe%grid%mp_low:ipe%grid%mp_high) /= &
                         ipe%neutrals%molecular_oxygen(:,:,ipe%grid%mp_low:ipe%grid%mp_high))
        n_nan_N2 = count(ipe%neutrals%molecular_nitrogen(:,:,ipe%grid%mp_low:ipe%grid%mp_high) /= &
                         ipe%neutrals%molecular_nitrogen(:,:,ipe%grid%mp_low:ipe%grid%mp_high))
        n_nan_vel = count(ipe%neutrals%velocity_geographic(:,:,:,ipe%grid%mp_low:ipe%grid%mp_high) /= &
                          ipe%neutrals%velocity_geographic(:,:,:,ipe%grid%mp_low:ipe%grid%mp_high))
        if (n_nan_T + n_nan_O + n_nan_O2 + n_nan_N2 + n_nan_vel > 0) then
          write(6,'(A,I4,A,5I8)') 'NaN AFTER NEUTRALS rank=', ipe%mpi_layer%rank_id, &
            ' T/O/O2/N2/vel: ', n_nan_T, n_nan_O, n_nan_O2, n_nan_N2, n_nan_vel
        endif
        if (ipe%mpi_layer%rank_id == 0) then
          write(6,'(A,2ES12.4)') 'DIAG neutrals T range: ', &
            minval(ipe%neutrals%temperature(:,:,ipe%grid%mp_low:ipe%grid%mp_high)), &
            maxval(ipe%neutrals%temperature(:,:,ipe%grid%mp_low:ipe%grid%mp_high))
          write(6,'(A,2ES12.4)') 'DIAG neutrals O range: ', &
            minval(ipe%neutrals%oxygen(:,:,ipe%grid%mp_low:ipe%grid%mp_high)), &
            maxval(ipe%neutrals%oxygen(:,:,ipe%grid%mp_low:ipe%grid%mp_high))
          ! geographic wind magnitude (1)=east, (2)=north, (3)=up
          write(6,'(A,3ES12.4)') 'DIAG WIND max|vg(E/N/Up)|: ', &
            maxval(abs(ipe%neutrals%velocity_geographic(1,:,:,ipe%grid%mp_low:ipe%grid%mp_high))), &
            maxval(abs(ipe%neutrals%velocity_geographic(2,:,:,ipe%grid%mp_low:ipe%grid%mp_high))), &
            maxval(abs(ipe%neutrals%velocity_geographic(3,:,:,ipe%grid%mp_low:ipe%grid%mp_high)))
        endif
      end block
      ENDIF

      ! EIA migration: O+ and neutrals at lp=90 over three altitudes
      IF ( verbose_diag ) THEN
      block
        integer :: diag_lp, dmp, dkp
        integer :: kp_lo, kp_mid, kp_hi
        real(prec) :: op_lo, op_mid, op_hi
        real(prec) :: neut_O_val, neut_T_val
        real(prec) :: alt_lo, alt_mid, alt_hi
        real(prec) :: geo_lon_deg
        diag_lp = 90
        ! low (~200km), mid (~350km), high (~500km) footpoint
        kp_lo  = max(1, ipe%grid%flux_tube_max(diag_lp) / 6)
        kp_mid = max(1, ipe%grid%flux_tube_max(diag_lp) / 4)
        kp_hi  = max(1, ipe%grid%flux_tube_max(diag_lp) / 3)
        alt_lo  = ipe%grid%altitude(kp_lo, diag_lp) * 1.0e-3_prec
        alt_mid = ipe%grid%altitude(kp_mid, diag_lp) * 1.0e-3_prec
        alt_hi  = ipe%grid%altitude(kp_hi, diag_lp) * 1.0e-3_prec
        do dmp = ipe%grid%mp_low, ipe%grid%mp_high
          op_lo  = ipe%plasma%ion_densities(1, kp_lo,  diag_lp, dmp)
          op_mid = ipe%plasma%ion_densities(1, kp_mid, diag_lp, dmp)
          op_hi  = ipe%plasma%ion_densities(1, kp_hi,  diag_lp, dmp)
          neut_O_val = ipe%neutrals%oxygen(kp_mid, diag_lp, dmp)
          neut_T_val = ipe%neutrals%temperature(kp_mid, diag_lp, dmp)
          geo_lon_deg = 57.2957795_prec * ipe%grid%longitude(kp_mid, diag_lp, dmp)
          write(6,'(A,I4,A,I3,A,3ES11.3,A,ES11.3,A,F7.1,A,F6.1)') &
            'EIA_DIAG PRE rank=', ipe%mpi_layer%rank_id, &
            ' mp=', dmp, &
            ' O+(lo/mid/hi)=', op_lo, op_mid, op_hi, &
            ' nO=', neut_O_val, ' nT=', neut_T_val, &
            ' glon=', geo_lon_deg
        enddo
        if (ipe%mpi_layer%rank_id == 0) then
          write(6,'(A,3F7.1)') 'EIA_DIAG altitudes(km) lo/mid/hi=', alt_lo, alt_mid, alt_hi
        endif
      end block
      ENDIF

      CALL ipe % eldyn % Update( ipe % grid, &
                                 ipe % forcing, &
                                 ipe % time_tracker, &
                                 ipe % plasma, &
                                 ipe % parameters % offset1_deg, &
                                 ipe % parameters % offset2_deg, &
                                 ipe % parameters % potential_model, &
                                 ipe % mpi_layer, &
                                 rc = localrc )
      IF ( ipe_error_check( localrc, msg="Failed to update electrodynamics", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      ! NaN check on electrodynamics
      IF ( verbose_diag ) THEN
      block
        integer :: n_nan_exb
        n_nan_exb = count(ipe%eldyn%v_ExB_apex(:,:,ipe%grid%mp_low:ipe%grid%mp_high) /= &
                          ipe%eldyn%v_ExB_apex(:,:,ipe%grid%mp_low:ipe%grid%mp_high))
        if (n_nan_exb > 0) then
          write(6,'(A,I4,A,I8)') 'NaN AFTER ELDYN rank=', ipe%mpi_layer%rank_id, &
            ' ExB: ', n_nan_exb
        endif
        if (ipe%mpi_layer%rank_id == 0) then
          write(6,'(A,2ES12.4)') 'DIAG ExB(1) range: ', &
            minval(ipe%eldyn%v_ExB_apex(1,:,ipe%grid%mp_low:ipe%grid%mp_high)), &
            maxval(ipe%eldyn%v_ExB_apex(1,:,ipe%grid%mp_low:ipe%grid%mp_high))
        endif
      end block
      ENDIF

      ! EIA ExB at the crest and near the equator.
      ! v_ExB_apex: (1)=zonal, (2)=meridional; geographic ExB gives vertical drift.
      IF ( verbose_diag ) THEN
      block
        integer :: dmp, diag_lp_eia, diag_lp_eq, kp_eq
        real(prec) :: exb_z90, exb_m90, exb_z_eq, exb_m_eq
        real(prec) :: exb_geo_e, exb_geo_n, exb_geo_up
        diag_lp_eia = 90   ! EIA crest (~19 deg mlat)
        diag_lp_eq  = 150  ! near equator (~6 deg mlat)
        ! kp at ~300 km footpoint for equatorial lp
        kp_eq = max(1, ipe%grid%flux_tube_max(diag_lp_eq) / 4)
        do dmp = ipe%grid%mp_low, ipe%grid%mp_high
          exb_z90 = ipe%eldyn%v_ExB_apex(1, diag_lp_eia, dmp)
          exb_m90 = ipe%eldyn%v_ExB_apex(2, diag_lp_eia, dmp)
          exb_z_eq = ipe%eldyn%v_ExB_apex(1, diag_lp_eq, dmp)
          exb_m_eq = ipe%eldyn%v_ExB_apex(2, diag_lp_eq, dmp)
          exb_geo_e  = ipe%eldyn%v_ExB_geographic(1, kp_eq, diag_lp_eq, dmp)
          exb_geo_n  = ipe%eldyn%v_ExB_geographic(2, kp_eq, diag_lp_eq, dmp)
          exb_geo_up = ipe%eldyn%v_ExB_geographic(3, kp_eq, diag_lp_eq, dmp)
          write(6,'(A,I4,A,I3,A,2F8.1,A,2F8.1,A,3F8.1)') &
            'EIA_DIAG ExB rank=', ipe%mpi_layer%rank_id, &
            ' mp=', dmp, &
            ' lp90_z/m=', exb_z90, exb_m90, &
            ' eq_z/m=', exb_z_eq, exb_m_eq, &
            ' eq_geo_E/N/Up=', exb_geo_e, exb_geo_n, exb_geo_up
        enddo
      end block
      ENDIF

      CALL ipe % plasma % Update( ipe % grid, &
                                  ipe % neutrals, &
                                  ipe % forcing, &
                                  ipe % time_tracker, &
                                  ipe % mpi_layer, &
                                  ipe % eldyn % v_ExB_apex, &
                                  ipe % parameters % time_step, &
                                  ipe % parameters % colfac, &
                                  ipe % parameters % hpeq, &
                                  ipe % parameters % transport_highlat_lp, &
                                  ipe % parameters % perp_transport_max_lp, &
                                  rc = localrc )
      IF ( ipe_error_check( localrc, msg="Failed to update plasma", &
        line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

      ! NaN check on plasma
      IF ( verbose_diag ) THEN
      block
        integer :: n_nan_ion, n_nan_te, n_nan_ti
        n_nan_ion = count(ipe%plasma%ion_densities(:,:,:,ipe%plasma%mp_low:ipe%plasma%mp_high) /= &
                          ipe%plasma%ion_densities(:,:,:,ipe%plasma%mp_low:ipe%plasma%mp_high))
        n_nan_te = count(ipe%plasma%electron_temperature(:,:,ipe%plasma%mp_low:ipe%plasma%mp_high) /= &
                         ipe%plasma%electron_temperature(:,:,ipe%plasma%mp_low:ipe%plasma%mp_high))
        n_nan_ti = count(ipe%plasma%ion_temperature(:,:,ipe%plasma%mp_low:ipe%plasma%mp_high) /= &
                         ipe%plasma%ion_temperature(:,:,ipe%plasma%mp_low:ipe%plasma%mp_high))
        if (n_nan_ion + n_nan_te + n_nan_ti > 0) then
          write(6,'(A,I4,A,3I8)') 'NaN AFTER PLASMA rank=', ipe%mpi_layer%rank_id, &
            ' ion/Te/Ti: ', n_nan_ion, n_nan_te, n_nan_ti
        endif
        ! report the (kp,lp,mp) where O+ exceeds a physical ceiling (~1e13)
        block
          integer :: opl(3)
          real(prec) :: opmax
          opmax = maxval(ipe%plasma%ion_densities(1,:,:,ipe%plasma%mp_low:ipe%plasma%mp_high))
          if (opmax > 1.0e13_prec) then
            opl = maxloc(ipe%plasma%ion_densities(1,:,:,ipe%plasma%mp_low:ipe%plasma%mp_high))
            write(6,'(A,I4,A,ES11.3,A,I4,A,I4,A,I3)') &
              'DIAG OPLOC rank=', ipe%mpi_layer%rank_id, ' O+max=', opmax, &
              ' kp=', opl(1), ' lp=', opl(2), ' mp=', opl(3)+ipe%plasma%mp_low-1
          endif
        end block
        if (ipe%mpi_layer%rank_id == 0) then
          write(6,'(A,2ES12.4)') 'DIAG O+ range: ', &
            minval(ipe%plasma%ion_densities(1,:,:,ipe%plasma%mp_low:ipe%plasma%mp_high)), &
            maxval(ipe%plasma%ion_densities(1,:,:,ipe%plasma%mp_low:ipe%plasma%mp_high))
          write(6,'(A,2ES12.4)') 'DIAG Te range: ', &
            minval(ipe%plasma%electron_temperature(:,:,ipe%plasma%mp_low:ipe%plasma%mp_high)), &
            maxval(ipe%plasma%electron_temperature(:,:,ipe%plasma%mp_low:ipe%plasma%mp_high))
        endif
      end block
      ENDIF

      ! EIA migration: O+ at lp=90 after plasma transport
      IF ( verbose_diag ) THEN
      block
        integer :: diag_lp, dmp
        integer :: kp_lo, kp_mid, kp_hi
        real(prec) :: op_lo, op_mid, op_hi
        diag_lp = 90
        kp_lo  = max(1, ipe%grid%flux_tube_max(diag_lp) / 6)
        kp_mid = max(1, ipe%grid%flux_tube_max(diag_lp) / 4)
        kp_hi  = max(1, ipe%grid%flux_tube_max(diag_lp) / 3)
        do dmp = ipe%grid%mp_low, ipe%grid%mp_high
          op_lo  = ipe%plasma%ion_densities(1, kp_lo,  diag_lp, dmp)
          op_mid = ipe%plasma%ion_densities(1, kp_mid, diag_lp, dmp)
          op_hi  = ipe%plasma%ion_densities(1, kp_hi,  diag_lp, dmp)
          write(6,'(A,I4,A,I3,A,3ES11.3)') &
            'EIA_DIAG POST rank=', ipe%mpi_layer%rank_id, &
            ' mp=', dmp, &
            ' O+(lo/mid/hi)=', op_lo, op_mid, op_hi
        enddo
      end block
      ENDIF

      ! Update the timer
      CALL ipe % time_tracker % Increment( ipe % parameters % time_step )

    ENDDO

  END SUBROUTINE Update_IPE_Model

  SUBROUTINE Initialize_IPE_Model( ipe, filename, rc )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER, OPTIONAL,  INTENT(out)   :: rc

    ! Local
    INTEGER :: localrc

    ! Begin
    IF (PRESENT(rc)) rc = IPE_SUCCESS

    CALL ipe % Read( filename, rc=localrc )
    IF ( ipe_error_check( localrc, msg="Failed to read file "//filename, &
      line=__LINE__, file=__FILE__, rc=rc ) ) RETURN

    ! Field-line integrals (and the conductivities they populate) are computed
    ! inside the time loop by Calculate_Field_Line_Integrals once MSIS has filled
    ! the neutrals; do NOT call it here with zero-valued neutrals.

  END SUBROUTINE Initialize_IPE_Model


  SUBROUTINE Read_IPE_State( ipe, filename, rc )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    CHARACTER(*),       INTENT(in)    :: filename
    INTEGER, OPTIONAL,  INTENT(out)   :: rc

    ! Local
    INTEGER, PARAMETER :: num_ion_densities = 9
    CHARACTER(LEN=*), DIMENSION(num_ion_densities), PARAMETER :: ion_densities = &
      (/ &
        "/apex/o_plus_density   ", &
        "/apex/h_plus_density   ", &
        "/apex/he_plus_density  ", &
        "/apex/n_plus_density   ", &
        "/apex/no_plus_density  ", &
        "/apex/o2_plus_density  ", &
        "/apex/n2_plus_density  ", &
        "/apex/o_plus_2D_density", &
        "/apex/o_plus_2P_density"  &
      /)

    INTEGER, PARAMETER :: num_ion_velocities = 3
    CHARACTER(LEN=*), DIMENSION(num_ion_velocities), PARAMETER :: ion_velocities = &
      (/ &
        "/apex/o_plus_velocity ", &
        "/apex/h_plus_velocity ", &
        "/apex/he_plus_velocity"  &
      /)

    INTEGER, PARAMETER :: num_plasma_datasets = 2
    CHARACTER(LEN=*), DIMENSION(num_plasma_datasets), PARAMETER :: plasma_datasets = &
      (/ &
        "/apex/ion_temperature     ", &
        "/apex/electron_temperature"  &
      /)

    INTEGER, PARAMETER :: num_apex_velocities = 3
    CHARACTER(LEN=*), DIMENSION(num_apex_velocities), PARAMETER :: apex_velocities = &
      (/ &
        "/apex/neutral_apex1_velocity", &
        "/apex/neutral_apex2_velocity", &
        "/apex/neutral_apex3_velocity"  &
      /)


    INTEGER, PARAMETER :: num_geo_datasets = 3
    CHARACTER(LEN=*), DIMENSION(num_geo_datasets), PARAMETER :: geo_datasets = &
      (/ &
        "/apex/neutral_geographic_velocity1", &
        "/apex/neutral_geographic_velocity2", &
        "/apex/neutral_geographic_velocity3"  &
      /)

    INTEGER :: item

    ! Begin
    IF ( PRESENT( rc ) ) rc = IPE_FAILURE

    IF( ipe % mpi_layer % rank_id == 0 )THEN
      PRINT *, '  Reading output file : '//TRIM(filename)
    ENDIF

    ! Open HDF5 input file
    CALL ipe % io % open(filename, "r")
    IF (ipe % io % err % check(msg="Failed to open file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    ! Setup common data decomposition for datasets
    CALL ipe % io % domain( (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % grid % NMP /), &
      (/ 1, 1, ipe % mpi_layer % mp_low /), &
      (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % mpi_layer % mp_high - ipe % mpi_layer % mp_low + 1 /) )
    IF (ipe % io % err % check(msg="Failed to setup I/O data decomposition", &
      file=__FILE__, line=__LINE__)) RETURN

    ! Read individual datasets
    ! -- Ion densities
    DO item = 1, num_ion_densities
      CALL ipe % io % read(ion_densities(item), &
        ipe % plasma % ion_densities(item,:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset "//ion_densities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Ion velocities
    DO item = 1, num_ion_velocities
      CALL ipe % io % read(ion_velocities(item), &
        ipe % plasma % ion_velocities(item,:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset "//ion_velocities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Plasma temperatures
    DO item = 1, num_plasma_datasets
      SELECT CASE (item)
      CASE (1)
        CALL ipe % io % read(plasma_datasets(item), &
          ipe % plasma % ion_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Failed to read dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      CASE (2)
        CALL ipe % io % read(plasma_datasets(item), &
          ipe % plasma % electron_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Failed to read dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END SELECT
    END DO

    ! -- Compute plasma electron density from ion densities
    ipe % plasma % electron_density2(:,:,ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high) = &
      SUM(ipe % plasma % ion_densities(1:9, :, :, ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high), dim=1)

    ! -- Read neutral datasets if requested
    IF( ipe % parameters % read_apex_neutrals )THEN
      CALL ipe % io % read("/apex/o_density", &
        ipe % neutrals % oxygen(:,:,                &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/o_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/h_density", &
        ipe % neutrals % hydrogen(:,:,              &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/h_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/he_density", &
        ipe % neutrals % helium(:,:,                 &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/he_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/n_density", &
        ipe % neutrals % nitrogen(:,:,              &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/n_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/o2_density", &
        ipe % neutrals % molecular_oxygen(:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/o2_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/n2_density", &
        ipe % neutrals % molecular_nitrogen(:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/n2_density", &
        file=__FILE__, line=__LINE__)) RETURN

      CALL ipe % io % read("/apex/neutral_temperature", &
        ipe % neutrals % temperature(:,:,                     &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Failed to read dataset /apex/neutral_temperature", &
        file=__FILE__, line=__LINE__)) RETURN

      ! Read neutral velocities on apex grid
      DO item = 1, num_apex_velocities
        CALL ipe % io % read(apex_velocities(item), &
          ipe % neutrals % velocity_apex(item,:,:,         &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Failed to read dataset "//apex_velocities(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Read neutral velocities on geographic grid, if requested
    IF( ipe % parameters % read_geographic_neutrals )THEN
      DO item = 1, num_geo_datasets
        CALL ipe % io % read(geo_datasets(item), &
          ipe % neutrals % velocity_geographic(item,:,:, &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Failed to read dataset "//geo_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Close HDF5 file
    CALL ipe % io % close()
    IF (ipe % io % err % check(msg="Failed to close file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

  END SUBROUTINE Read_IPE_State

  SUBROUTINE Write_IPE_State( ipe, rc )

    IMPLICIT NONE

    CLASS( IPE_Model ), INTENT(inout) :: ipe
    INTEGER, OPTIONAL,  INTENT(out)   :: rc

    ! Local
    CHARACTER(215)  :: filename
    INTEGER, PARAMETER :: num_groups = 1
    CHARACTER(LEN=*), DIMENSION(num_groups),   PARAMETER :: groups = (/ "apex" /)

    INTEGER, PARAMETER :: num_ion_densities = 9
    CHARACTER(LEN=*), DIMENSION(num_ion_densities), PARAMETER :: ion_densities = &
      (/ &
        "/apex/o_plus_density   ", &
        "/apex/h_plus_density   ", &
        "/apex/he_plus_density  ", &
        "/apex/n_plus_density   ", &
        "/apex/no_plus_density  ", &
        "/apex/o2_plus_density  ", &
        "/apex/n2_plus_density  ", &
        "/apex/o_plus_2D_density", &
        "/apex/o_plus_2P_density"  &
      /)

    INTEGER, PARAMETER :: num_ion_velocities = 3
    CHARACTER(LEN=*), DIMENSION(num_ion_velocities), PARAMETER :: ion_velocities = &
      (/ &
        "/apex/o_plus_velocity ", &
        "/apex/h_plus_velocity ", &
        "/apex/he_plus_velocity"  &
      /)

    INTEGER, PARAMETER :: num_plasma_datasets = 2
    CHARACTER(LEN=*), DIMENSION(num_plasma_datasets), PARAMETER :: plasma_datasets = &
      (/ &
        "/apex/ion_temperature     ", &
        "/apex/electron_temperature"  &
      /)

    INTEGER, PARAMETER :: num_apex_velocities = 3
    CHARACTER(LEN=*), DIMENSION(num_apex_velocities), PARAMETER :: apex_velocities = &
      (/ &
        "/apex/neutral_apex1_velocity", &
        "/apex/neutral_apex2_velocity", &
        "/apex/neutral_apex3_velocity"  &
      /)

    INTEGER, PARAMETER :: num_geo_datasets = 3
    CHARACTER(LEN=*), DIMENSION(num_geo_datasets), PARAMETER :: geo_datasets = &
      (/ &
        "/apex/neutral_geographic_velocity1", &
        "/apex/neutral_geographic_velocity2", &
        "/apex/neutral_geographic_velocity3"  &
      /)

    INTEGER, PARAMETER :: num_exb_datasets = 3
    CHARACTER(LEN=*), DIMENSION(num_exb_datasets), PARAMETER :: exb_datasets = &
      (/ &
        "/apex/eastward_exb_velocity ", &
        "/apex/northward_exb_velocity", &
        "/apex/upward_exb_velocity   "  &
      /)

    INTEGER, PARAMETER :: num_conductivities = 8
    CHARACTER(LEN=*), DIMENSION(num_conductivities), PARAMETER :: conductivities = &
      (/ &
        "/apex/integral513", &
        "/apex/integral514", &
        "/apex/integral517", &
        "/apex/integral518", &
        "/apex/integral519", &
        "/apex/integral520", &
        "/apex/integral1  ", &
        "/apex/integral2  "  &
      /)

    INTEGER, PARAMETER :: num_elef = 2
    CHARACTER(LEN=*), DIMENSION(num_elef), PARAMETER :: electric_field = &
      (/ &
        "/apex/Ed1", &
        "/apex/Ed2"  &
      /)

    INTEGER :: item
    CHARACTER(LEN=28) :: dset_name

    ! Begin

    IF ( PRESENT( rc ) ) rc = IPE_FAILURE

    filename = TRIM(ipe % parameters % file_prefix) // &
               ipe % time_tracker % DateStamp ( ) // &
               ipe % parameters % file_extension

    IF( ipe % mpi_layer % rank_id == 0 )THEN
      PRINT *, '  Writing output file : '//TRIM(filename)
    ENDIF

    ! Create HDF5 input file
    CALL ipe % io % open(filename, "c")
    IF (ipe % io % err % check(msg="Unable to open file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    ! Setup common data decomposition for datasets
    CALL ipe % io % domain( (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % grid % NMP /), &
      (/ 1, 1, ipe % mpi_layer % mp_low /), &
      (/ ipe % grid % nFluxtube, ipe % grid % NLP, ipe % mpi_layer % mp_high - ipe % mpi_layer % mp_low + 1 /) )
    IF (ipe % io % err % check(msg="Failed to setup I/O data decomposition", &
      file=__FILE__, line=__LINE__)) RETURN

    ! Write individual datasets
    ! -- Ion densities
    DO item = 1, num_ion_densities
      CALL ipe % io % write(ion_densities(item), &
        ipe % plasma % ion_densities(item,:,:,       &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//ion_densities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Ion velocities
    DO item = 1, num_ion_velocities
      CALL ipe % io % write(ion_velocities(item), &
        ipe % plasma % ion_velocities(item,:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//ion_velocities(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO
    ! -- Plasma temperatures
    DO item = 1, num_plasma_datasets
      SELECT CASE (item)
      CASE (1)
        CALL ipe % io % write(plasma_datasets(item), &
          ipe % plasma % ion_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      CASE (2)
        CALL ipe % io % write(plasma_datasets(item), &
          ipe % plasma % electron_temperature(:,:,            &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//plasma_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END SELECT
    END DO

    ! -- Read neutral datasets if requested
    IF( ipe % parameters % write_apex_neutrals )THEN

      dset_name = "/apex/o_density"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % oxygen(:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/h_density"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % hydrogen(:,:,      &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/he_density"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % helium(:,:,        &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/n_density"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % nitrogen(:,:,      &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/o2_density"
      CALL ipe % io % write(dset_name,    &
        ipe % neutrals % molecular_oxygen(:,:, &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/n2_density"
      CALL ipe % io % write(dset_name,      &
        ipe % neutrals % molecular_nitrogen(:,:, &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      dset_name = "/apex/neutral_temperature"
      CALL ipe % io % write(dset_name, &
        ipe % neutrals % temperature(:,:,   &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//dset_name, &
        file=__FILE__, line=__LINE__)) RETURN

      ! Read neutral velocities on apex grid
      DO item = 1, num_apex_velocities
        CALL ipe % io % write(apex_velocities(item), &
          ipe % neutrals % velocity_apex(item,:,:,         &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//apex_velocities(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Read neutral velocities on geographic grid, if requested
    IF( ipe % parameters % write_geographic_neutrals )THEN
      DO item = 1, num_geo_datasets
        CALL ipe % io % write(geo_datasets(item), &
          ipe % neutrals % velocity_geographic(item,:,:, &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//geo_datasets(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO
    END IF

    ! Write ExB drift velocities on geographic grid
    DO item = 1, num_exb_datasets
      CALL ipe % io % write(exb_datasets(item), &
        ipe % eldyn % v_ExB_geographic(item,:,:, &
        ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
      IF (ipe % io % err % check(msg="Unable to write dataset "//exb_datasets(item), &
        file=__FILE__, line=__LINE__)) RETURN
    END DO

    ! -- Field-line-integrated conductances and electric field, if requested
    IF( ipe % parameters % write_conductivities )THEN

      ! Conductivities decomposition (hemisphere, NLP, NMP)
      CALL ipe % io % domain( (/ 2, ipe % grid % NLP, ipe % grid % NMP /), &
        (/ 1, 1, ipe % mpi_layer % mp_low /), &
        (/ 2, ipe % grid % NLP, ipe % mpi_layer % mp_high - ipe % mpi_layer % mp_low + 1 /) )
      IF (ipe % io % err % check(msg="Failed to setup I/O data decomposition", &
        file=__FILE__, line=__LINE__)) RETURN

      DO item = 1, num_conductivities
        CALL ipe % io % write(conductivities(item), &
          ipe % plasma % conductivities(item,:,:, &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//conductivities(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO

      ! Electric field decomposition (NLP, NMP)
      CALL ipe % io % domain( (/ ipe % grid % NLP, ipe % grid % NMP /), &
        (/ 1, ipe % mpi_layer % mp_low /), &
        (/ ipe % grid % NLP, ipe % mpi_layer % mp_high - ipe % mpi_layer % mp_low + 1 /) )
      IF (ipe % io % err % check(msg="Failed to setup I/O data decomposition", &
        file=__FILE__, line=__LINE__)) RETURN

      DO item = 1, num_elef
        CALL ipe % io % write(electric_field(item), &
          ipe % eldyn % electric_field(item,:, &
          ipe % mpi_layer % mp_low:ipe % mpi_layer % mp_high))
        IF (ipe % io % err % check(msg="Unable to write dataset "//electric_field(item), &
          file=__FILE__, line=__LINE__)) RETURN
      END DO

    END IF

    ! Close HDF5 file
    CALL ipe % io % close()
    IF (ipe % io % err % check(msg="Unable to close file "//filename, &
      file=__FILE__, line=__LINE__)) RETURN

    IF ( PRESENT( rc ) ) rc = IPE_SUCCESS

  END SUBROUTINE Write_IPE_State


END MODULE IPE_Model_Class
