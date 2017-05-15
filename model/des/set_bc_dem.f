!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_DEM                                              !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Purpose: Check the data provided for the des mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the proveded information into the format       !
!  necessary for the dependent subrountines to function properly.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_DEM

      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE run
      USE mfix_pic
      use mpi_utility

      use bc

      use error_manager

      IMPLICIT NONE


      CALL INIT_ERR_MSG("SET_BC_DEM")

! The variable PARTICLES should already be set by this point if using
! gener_part_config option
      IF(PARTICLES == UNDEFINED_I) THEN
         PARTICLES = 0
      ENDIF

! If the system is started without any particles and an inlet is not
! specified, the run is likely aborted.
! Inlet/outlet for MPPIC are based off the regular mfix declarations,
! and so DEM_BCMI could still be zero.
      IF(PARTICLES == 0 .AND. DEM_BCMI == 0) THEN
         WRITE(ERR_MSG, 1202)
         CALL FLUSH_ERR_MSG
      ENDIF

 1202 FORMAT('WARNING 1202: The system is initiated with no particles',&
         ' and no',/'solids inlet was detected.')

      IF(DEM_BCMI > 0) CALL SET_BC_DEM_MI
      IF(DEM_BCMO > 0) CALL SET_BC_DEM_MO

! Set the flag that one or more DEM MI/MO exists.
      DEM_MIO = (DEM_BCMI /= 0 .OR. DEM_BCMO /= 0)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_DEM
