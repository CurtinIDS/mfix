!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS2_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms are applied and the time step updated. The   !
!  The user may insert code in this routine or call user defined       !
!  subroutines.                                                        !
!                                                                      !
!  This routien is called from the time loop, but no indicies (fluid   !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR2_DES

      Use discretelement
      Use run
      Use usr

      IMPLICIT NONE

! Check if particle reached the peak of the bounce.
      IF(DES_VEL_new(1,2) < 0.0) THEN
         IF(yVELO > 0.0) THEN
            IF(BOUNCE_COUNT < MAX_BOUNCE) THEN
               BOUNCE_COUNT = BOUNCE_COUNT + 1
               MAX_HEIGHT(BOUNCE_COUNT) = max(yPOSO, DES_POS_NEW(1,2))
            ENDIF
         ENDIF
      ENDIF

      yVELO = DES_VEL_NEW(1,2)
      yPOSO = DES_POS_NEW(1,2)

      RETURN
      END SUBROUTINE USR2_DES
