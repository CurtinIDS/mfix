!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR1                                                   C
!  Purpose: This routine is called from the time loop and is           C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting or checking errors in quantities   C
!           that vary with time.  This routine is not called from an   C
!           IJK loop, hence all indices are undefined.                 C         !                                                                      C
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
!
      SUBROUTINE USR1
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE toleranc
      USE constant
      USE compar
      USE funits
      USE fun_avg
      USE usr
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
      INTEGER IJK, M, I, IMJK, IJMK, IJKM
      DOUBLE PRECISION DIFF, EP_g2
      DOUBLE PRECISION Sc1o3, UGC, VGC, WGC, USCM, VSCM, WSCM, VREL, Re
      INCLUDE 'usrnlst.inc'
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!
!
!     Wall heat transfer coefficient in cal/cm^3.s.K
!
!
      DO M = 1, MMAX
!
!$omp  parallel do private(IJK)
         DO IJK = IJKSTART3, IJKEND3
            IF (FLUID_AT(IJK)) THEN
	      I = I_OF(IJK)
              IMJK  = IM_OF(IJK)
              IJMK  = JM_OF(IJK)
              IJKM  = KM_OF(IJK)
!
!   Calculate Sherwood number for solids phases (Gunn 1978)
!
              EP_g2 = EP_g(IJK) * EP_g(IJK)
              DIFF = 4.26 * ((T_g(IJK)/1800.)**1.75) * 1013000. / P_g(IJK)
              Sc1o3 = (MU_g(IJK)/(RO_g(IJK) * DIFF))**(1./3.)
              UGC = AVG_X_E(U_g(IMJK), U_g(IJK), I)
              VGC = AVG_Y_N(V_g(IJMK), V_g(IJK))
              WGC = AVG_Z_T(W_g(IJKM), W_g(IJK))
!
              USCM = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
              VSCM = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
              WSCM = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
!
              VREL = SQRT((UGC - USCM)**2 + (VGC-VSCM)**2   &
                                        + (WGC-WSCM)**2 )
              Re = EP_g(IJK) * D_p0(M) * VREL * RO_g(IJK) / MU_g(IJK)
              N_sh(IJK, M) = ( (7. - 10. * EP_g(IJK) + 5. * EP_g2) &
                             *(ONE + 0.7 * Re**0.2 * Sc1o3)      &
                            + (1.33 - 2.4*EP_g(IJK) + 1.2*EP_g2) &
                             * Re**0.7 * Sc1o3 )
            ENDIF
         END DO
      END DO
      RETURN
      END SUBROUTINE USR1
