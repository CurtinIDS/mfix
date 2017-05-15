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

      use constant, only: PI
      use discretelement, only: PIP
      use discretelement, only: VEL => DES_VEL_NEW
      use exit, only: mfix_exit
      use usr

      IMPLICIT NONE

      INTEGER :: NP

      IF(PIP /= 93) THEN
         write(*,"(3x, 'invalid setup for test case')")
         call mfix_exit(0)
      ENDIF

! Store the collision angle and initial tangential velocity
      DO NP=1, 62
         INIT_VEL_T(NP) = sqrt(VEL(NP,1)**2 + VEL(NP,3)**2)
         INIT_ANGLE(NP) = abs(atan(INIT_VEL_T(NP)/VEL(NP,2)))*180.0/PI
      ENDDO

      return
      END SUBROUTINE USR0
