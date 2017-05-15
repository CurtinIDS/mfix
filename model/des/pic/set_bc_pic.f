!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_PIC                                              !
!  Author: R. Garg                                    Date: 11-Jun-14  !
!                                                                      !
!  Purpose: Check the data provided for the pic mass inflow boundary   !
!  condition and flag errors if the data is improper.  This module is  !
!  also used to convert the proveded information into the format       !
!  necessary for the dependent subrountines to function properly.      !
!  Rehack of set_bc_dem                                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_BC_PIC

! Modules
!---------------------------------------------------------------------//
      USE discretelement, only: particles
      use error_manager
      USE param1, only: undefined_i
      USE pic_bc, only: pic_bcmi, pic_bcmo, pic_mio
      IMPLICIT NONE
!......................................................................!

      CALL INIT_ERR_MSG("SET_BC_PIC")

! The variable PARTICLES should already be set by this point if using
! gener_part_config option
      IF(PARTICLES == UNDEFINED_I) THEN
         PARTICLES = 0
      ENDIF

! If the system is started without any particles and an inlet is not
! specified, the run is likely aborted.
! Inlet/outlet for MPPIC are based off the regular mfix declarations,
! and so PIC_BCMI could still be zero.
      IF(PARTICLES == 0 .AND. PIC_BCMI == 0) THEN
         WRITE(ERR_MSG, 1202)
         CALL FLUSH_ERR_MSG
      ENDIF

 1202 FORMAT('WARNING 1202: The system is initiated with no particles',&
         ' and no',/'solids inlet was detected.')

      IF(PIC_BCMI > 0) CALL SET_BC_PIC_MI
      IF(PIC_BCMO > 0) CALL SET_BC_PIC_MO

! Set the flag that one or more PIC MI/MO exists.
      PIC_MIO = (PIC_BCMI /= 0 .OR. PIC_BCMO /= 0)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE SET_BC_PIC
