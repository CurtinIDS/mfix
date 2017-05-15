!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR0                                                   !
!  Author: J.Musser                                   Date: dd-mmm-yy  !
!  Purpose: This routine is called before the time loop starts and is  !
!           user-definable.  The user may insert code in this routine  !
!           or call appropriate user defined subroutines.  This        !
!           can be used for setting constants and checking errors in   !
!           data.  This routine is not called from an IJK loop, hence  !
!           all indices are undefined.                                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR0

      use discretelement
      use exit, only: mfix_exit
      use usr

      IMPLICIT NONE

! Initialize:
! Previoius position and velocity.
      yPOSO = DES_POS_NEW(1,2)
      yVELO = DES_VEL_NEW(1,2)
! Bounce counter
      BOUNCE_COUNT = 0
! Max bounce height
      MAX_HEIGHT(0) = yPOSO
! Drop height
      h0 = yPOSO

      return
      END SUBROUTINE USR0
