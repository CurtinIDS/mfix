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

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use run
      Use usr

      use output, only: USR_DT

      IMPLICIT NONE

      DOUBLE PRECISION :: OVERSHOOT

      DOUBLE PRECISION, SAVE :: USR_TIME = 0.0


      IF(USR_TIME == 0.0) THEN
         USR_TIME = USR_DT(1)
         RETURN
      ENDIF


      OVERSHOOT = S_TIME + 0.1d0*DTSOLID

! Write special output, if needed.
      IF(OVERSHOOT >= USR_TIME) THEN
! Update the time to write the special output.
         USR_TIME = (INT((OVERSHOOT)/USR_DT(1))+1)*USR_DT(1)
! Write the output.
         CALL WRITE_GRAN_TEMP
      ENDIF

      RETURN
      END SUBROUTINE USR2_DES
