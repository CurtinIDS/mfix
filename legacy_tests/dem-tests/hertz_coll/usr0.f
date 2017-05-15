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

      use usr
      use discretelement, only: PIP
      use discretelement, only: VEL => DES_VEL_NEW
      use constant, only: PI

      IMPLICIT NONE

      DOUBLE PRECISION :: NP

! Store the collision angle and initial tangential velocity
      DO NP=1, PIP
         INIT_ANGLE(NP) = abs(atan(VEL(1,NP)/VEL(2,NP)))*180.0/PI
         INIT_VEL_T(NP) = VEL(1,NP)
      ENDDO



      return
      END SUBROUTINE USR0
