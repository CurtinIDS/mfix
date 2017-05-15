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

      USE usr
      USE param
      USE param1
      USE parallel
      USE physprop
      USE geometry
      USE fldvar
      USE indices
      USE constant
      USE toleranc
      USE compar
      USE run
      USE turb
      USE sendrecv
      USE discretelement
      USE functions

      IMPLICIT NONE

!                      Indices
      INTEGER          I, J, K, IJK
      DOUBLE PRECISION XX, YY, ZZ, XM, YM, ZM
      DOUBLE PRECISION, PARAMETER :: T_per=0.25d0
      double precision nu, ftauf, wnot, kx, ky, kmag

      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         XX = XE(I)
         YY = YN(J)
         XM = XE(I) - 0.5d0*DX(I)
         YM = YN(J) - 0.5d0*DY(J)


         nu = mu_g0/ro_g0 ! kinematic viscosity
         ftauf = 1.0 ! flow time scale
         wnot = 1.0 ! initial vorticity
         kx = 200*pi; ky = 200*pi ! wavenumber of TG vortices
         kmag = 0.05

         u_g(ijk) = -cos(kx*xx)*sin(ky*ym)*kmag
         v_g(ijk) =  sin(kx*xm)*cos(ky*yy)*kmag

      END DO

      RETURN
      END SUBROUTINE USR1
