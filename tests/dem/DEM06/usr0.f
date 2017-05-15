!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!                                                                      C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR0

      use discretelement, only: DES_POS_NEW
      use discretelement, only: DES_VEL_NEW
      use discretelement, only: PIP
      use exit, only: mfix_exit
      use usr

      IMPLICIT NONE

      INTEGER :: NP

      IF(PIP /= 1) THEN
         write(*,"(3x, 'invalid setup for test case')")
         call mfix_exit(0)
      ENDIF

! Store the initial position and velocity
      RK4_POS(:) = DES_POS_NEW(1,:)
      RK4_VEL(:) = DES_VEL_NEW(1,:)

      return
      END SUBROUTINE USR0
