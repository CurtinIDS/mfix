!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_IC_DEM                                            !
!  Author:   R.Garg                                   Date: 11-Mar-14  !
!                                                                      !
!  Purpose: check the initial conditions input section for DEM model   !
!     - calculate the number of particles needed to initialize the      !
!        DEM model                                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_IC_DEM

! Global variables
!---------------------------------------------------------------------//
! Runtime Flag: Generate initial particle configuration.
      USE discretelement, only : gener_part_config
! Simulation dimension (2D/3D)
      USE discretelement, only: DIMN
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! runtime flag to specify if the simulation uses both discrete/continuum
! solids
      USE discretelement, only: DES_CONTINUUM_HYBRID

! direction wise spans of the domain and grid spacing in each direction
      Use geometry, only: zlength
! Use the error manager for posting error messages.
      use error_manager

      use physprop, only: mmax, d_p0
      use toleranc

      implicit none
! Local variables
!---------------------------------------------------------------------//
      integer :: dm
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_DEM")

! Determine the domain volume which is used to calculate the total
! number of particles and the number of particles in each phase.
! Values of DZ(1) or zlength are guaranteed at this point due to
! check_geometry_prereqs. If the user left both undefined and
! NO_K = .T., then they are set to ONE. If dz(1) is undefined but
! zlength is defined, then dz(1) is set to zlength (and vice versa).
! If both are defined they must be equal.
      IF(DIMN.EQ.2) THEN
         IF (DES_MMAX.EQ.1) THEN
! account for a possible offset index when using d_p0 and ro_s.
            DM = MMAX+1
! Warn the user if the domain depth is not equal to the particle
! diameter as it may cause problems for coupled simulations.
! The user should also be aware of this when interpreting
! volume/void fraction calculations (including bulk density).
            IF(.NOT.COMPARE(ZLENGTH,D_P0(DM))) THEN
               WRITE(ERR_MSG, 1000) iVal(d_p0(DM))
               CALL FLUSH_ERR_MSG
            ENDIF
         ELSE
! Let the user know basis of depth dimension for calculating number of
! particles. this will also be important when considering volume/void
! fraction calculations.
            WRITE(ERR_MSG, 1001)
            CALL FLUSH_ERR_MSG
         ENDIF
      ENDIF


 1000 FORMAT(' Message: ',&
      'WARNING: zlength or dz(1) is used to calculate the ',&
      'number of particles in the 2D simulation when ',&
      'GENER_PART_CONFIG is T and DIMN = 2.',/10X,'This depth ',&
      'does not equal D_P0 = ', A, '.')

 1001 FORMAT(' Message: ',&
      'WARNING: zlength or dz(1) is used to calculate the ',&
      'number of particles in the 2D simulation when ',&
      'GENER_PART_CONFIG is T and DIMN = 2.')


      IF (Gener_part_config.and.DES_CONTINUUM_HYBRID) THEN
         WRITE(ERR_MSG, 999)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 999  format('Error # 999: Gener_part_config set to', &
      ' true for DES_continuum hybrid', /, &
      ' This is not allowed, specify the initial particle', &
      ' configuration explicitly', /, &
      ' See MFIX readme', /,  &
      ' Please correct the data file.')

      CALL FINL_ERR_MSG

      END SUBROUTINE CHECK_IC_DEM
