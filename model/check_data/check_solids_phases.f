!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  SUBROUTINE: CHECK_SOLIDS_PHASES                                     !
!                                                                      !
!  Purpose: Driver routine for calls to solids phase checks.           !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_PHASES


! Global Variables:
!---------------------------------------------------------------------//
! Runtime flag specifying TFM solids
      use run, only: TFM_SOLIDS
! Runtime flag specifying DEM solids
      use run, only: DEM_SOLIDS
! Runtime flag specifying MPPIC solids
      use run, only: PIC_SOLIDS

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//

!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_PHASES")

! Impose the various model limitations.
      CALL CHECK_SOLIDS_MODEL_LIMITATIONS

! Checks common to all solids models.
      CALL CHECK_SOLIDS_COMMON_ALL

! Checks common to discrete solids phases (DEM, MPPIC).
      IF(DEM_SOLIDS .OR. PIC_SOLIDS) &
         CALL CHECK_SOLIDS_COMMON_DISCRETE

! Checks specific to the particular solids phase.
      IF(TFM_SOLIDS) CALL CHECK_SOLIDS_CONTINUUM
      IF(DEM_SOLIDS) CALL CHECK_SOLIDS_DEM
      IF(PIC_SOLIDS) CALL CHECK_SOLIDS_MPPIC

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_PHASES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  SUBROUTINE: CHECK_SOLIDS_MODEL_LIMITATIONS                          !
!                                                                      !
!  Purpose: Impose the limitations of the various solids models. These !
!  checks should be 'high level' in that they only ensure that models  !
!  are only used with phases that support them.                        !
!                                                                      !
!  Author: J.Musser                                   Date: 28-Feb-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_MODEL_LIMITATIONS

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: TFM solids present.
      use run, only: TFM_SOLIDS
! Runtime Flag: DEM solids present.
      use run, only: DEM_SOLIDS
! Runtime Flag: PIC solids present.
      use run, only: PIC_SOLIDS
! Runtime Flag: Invoke a cohesion model for DES simulation.
      use discretelement, only: USE_COHESION
! Runtime Flag: Sovle energy equations
      use run, only: ENERGY_EQ
! Runtime Flag: Sovle species equations
      use run, only: SPECIES_EQ

! Number of solid phases specified by the user/TFM model.
      use physprop, only: SMAX
! Number of discrete solids.
      use discretelement, only: DES_MMAX

! Global Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      implicit none


! Local Variables:
!---------------------------------------------------------------------//
! Total number of solids phases.
      INTEGER :: MMAX_TOT


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_MODEL_LIMITATIONS")


! Set up the total number of solids.
      MMAX_TOT = SMAX + DES_MMAX


! The cohesion model is only implemented for DEM simulations
      IF(USE_COHESION) THEN
         IF(TFM_SOLIDS .OR. PIC_SOLIDS) THEN
            WRITE(ERR_MSG, 2000)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2000 FORMAT('Error 2000: The solids cohesion model is only available',&
         ' for DEM',/' solids. Please correct the mfix.dat file.')
      ENDIF

! Place holder
      IF(ENERGY_EQ .AND. (TFM_SOLIDS .OR. PIC_SOLIDS)) THEN
!         WRITE(ERR_MSG, 2002)
!         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2002 FORMAT('Error 2002: The solids-solids conduction model is only', &
         ' available',/' for DEM only. Please correct the mfix.dat',   &
         ' file.')
      ENDIF


! This is only implemented for pure TFM or pure DEM simulations.
      IF(any(SPECIES_EQ(1:MMAX_TOT))) THEN
         IF(TFM_SOLIDS .AND. DEM_SOLIDS) THEN
            WRITE(ERR_MSG, 5000)
            CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
         ENDIF

 5000 FORMAT('Error 5000: Species equations are not available with',   &
         ' the hybrid',/'solids model. Please correct the mfix.dat',   &
         ' file.')
      ENDIF


      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_SOLIDS_MODEL_LIMITATIONS
