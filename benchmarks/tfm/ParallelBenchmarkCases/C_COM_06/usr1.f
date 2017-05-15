!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR1                                                   C
!  Purpose: This routine is called from the time loop and is           C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting or checking errors in quantities   C
!           that vary with time.  This routine is not called from an   C
!           IJK loop, hence all indices are undefined.                 C
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

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits
      USE compar
      USE fun_avg

      USE usr
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Bounded Phase temperatues (K)
      DOUBLE PRECISION, parameter :: MAX_TEMP = 2.5d3

      INTEGER :: IJK

      DOUBLE PRECISION :: Diff
      DOUBLE PRECISION :: N_Re, N_Sc, lN_Sh
      DOUBLE PRECISION :: xTg   ! Gas

!-----------------------------------------------
      INCLUDE 'species.inc'

      N_Sh = ZERO
      DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN
! Calculte the bounded gas phase temperature (K)
            xTg  = min(MAX_TEMP, T_g(IJK))
! Diffusion coefficient of water vapor in air.
            Diff = 4.26d-4 * ((xTg/1.8d3)**1.75d0) / (P_g(IJK) / 101.325d3)
! Reynolds Number
            N_Re = cal_NRe(1)
! Schmidt Number
            N_Sc = cal_NSc(Diff)
! Sherwood Number (Ranz and Marshal, 1952)
            N_Sh(IJK, 1) = cal_NSh(N_Re, N_Sc)
         ENDIF
      ENDDO

      RETURN
      CONTAINS

!----------------------------------------------------------------------!
! Function: calc_NRe(M)                                                !
!                                                                      !
! Purpose: Calculate the Reynolds number.                              !
!                                                                      !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION cal_NRe(M)

      INTEGER, INTENT(in) ::  M

! Various fluid cell indicies
      INTEGER I, IMJK, IJMK, IJKM
! Gas Velocity - cell centered
      DOUBLE PRECISION UGC, VGC, WGC
! Solids Velocity - cell centered
      DOUBLE PRECISION USCM, VSCM, WSCM
! Relative velocity.
      DOUBLE PRECISION VREL

! Initialize fluid cell variables
      I =  I_OF(IJK)
      IMJK  = IM_OF(IJK)
      IJMK  = JM_OF(IJK)
      IJKM  = KM_OF(IJK)

! Calculate velocity components at i, j, k
! Gas
      UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
      VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
      WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
! Solids
      USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
      VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M))
      WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

! magnitude of gas-solids relative velocity
      VREL = SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)

! Reynods Number
      IF(MU_g(IJK) > ZERO) THEN
         cal_NRe = EP_g(IJK) * D_P(IJK,M) * VREL * RO_g(IJK) / MU_g(IJK)
      ELSE
         cal_NRe = LARGE_NUMBER
      ENDIF

      RETURN
      END FUNCTION cal_NRe

!----------------------------------------------------------------------!
! Function: cal_NSc                                                    !
!                                                                      !
! Purpose: Calculate the Schmidt Number.                               !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION cal_NSc(Diff_Coeff)
! Diffustion coefficient
      DOUBLE PRECISION, intent(IN) :: Diff_Coeff

! Schmidt Number
      cal_NSc = MU_g(IJK)/(RO_g(IJK)*Diff_Coeff)

      RETURN
      END FUNCTION cal_NSc


!----------------------------------------------------------------------!
! Function: cal_NSh                                                    !
!                                                                      !
! Purpose: Calculate the Sherwood Number.                              !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION cal_NSh(lRe, lSc)

! Reynolds and Schmidt Numbers.
      DOUBLE PRECISION, intent(IN) :: lRe, lSc

! Cube root of Schmidt Number.
      DOUBLE PRECISION :: cr_Sc
! Squre of gas phase volume fraction
      DOUBLE PRECISION :: sEPg

      cr_Sc = lSc**(1.0d0/3.0d0)
      sEPg  = EP_g(IJK)*EP_g(IJK)

! Sherwood Number: Ranz and Marshall, 1952
      cal_NSh = (7.0d0 - 10.0d0*EP_g(IJK) + 5.0d0*sEPg) +              &
         (ONE + 0.7d0*(lRe**0.2d0)*cr_Sc) +                            &
         (1.33d0 - 2.4d0*EP_g(IJK) +  1.2d0*sEPg)*(lRe**0.7d0)*cr_Sc

      RETURN

      END FUNCTION cal_NSh

      END SUBROUTINE USR1
