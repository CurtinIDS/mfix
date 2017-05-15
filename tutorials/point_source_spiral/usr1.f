!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR1                                                   C
!  Purpose: This routine is called from the time loop and is           C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting or checking errors in quantities   C
!           that vary with time.  This routine is not called from an   C
!           IJK loop, hence all indices are undefined.                 C               C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR1

      use ps
      use constant
      use physprop
      use run
      use usr

      IMPLICIT NONE

      DOUBLE PRECISION :: lRad

      DOUBLE PRECISION :: lU, lV, lW


      INTEGER :: M

! One reveloution over a three second period.
      lRad = (2.0d0*Pi*time)/3.0d0

! Calculate the normalized velocity components.
      lV = (0.12d0/0.15d0)
      lU = (0.09d0/0.15d0)*cos(lRad)
      lW = (0.09d0/0.15d0)*sin(lRad)

! Update the gas phase components.
      PS_V_g(1) = lV
      PS_U_g(1) = lU
      PS_W_g(1) = lW

! Update the solids phase components.
      do M=1,MMAX
         PS_V_s(1,M) = lV
         PS_U_s(1,M) = lU
         PS_W_s(1,M) = lW
      enddo

      RETURN
      END SUBROUTINE USR1
