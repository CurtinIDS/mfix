!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_A_U_g(A_m, B_m, IER)                            C
!  Purpose: Handle the special case of the center coefficient in       C
!  U_g momentum eq. becoming zero.                                     C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date:  2-AUG-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
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
      SUBROUTINE ADJUST_A_U_G(A_M, B_M)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE run
      USE indices
      USE usr       !just to ensure that this module gets compiled early on
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          I, IP, IJK, IJKE, IMJK
!
!                      Phase index
      INTEGER          M
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Vector b_m
      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------

      M = 0
      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN
!
!!!!$omp parallel do private(I,IP,IJK,IJKE,IMJK)
      DO IJK = ijkstart3, ijkend3
         IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN
            A_M(IJK,east,M) = ZERO
            A_M(IJK,west,M) = ZERO
            A_M(IJK,north,M) = ZERO
            A_M(IJK,south,M) = ZERO
            A_M(IJK,top,M) = ZERO
            A_M(IJK,bottom,M) = ZERO
            A_M(IJK,0,M) = -ONE
            IF (B_M(IJK,M) < ZERO) THEN
               IJKE = EAST_OF(IJK)
               IP = IP1(I_OF(IJK))
               IF (ROP_G(IJKE)*AYZ_U(IJK) > SMALL_NUMBER) THEN
                  B_M(IJK,M) = SQRT((-B_M(IJK,M)/(ROP_G(IJKE)*AVG_X_E(ONE,ZERO,&
                     IP)*AYZ_U(IJK))))
               ELSE
                  B_M(IJK,M) = ZERO
               ENDIF
            ELSE IF (B_M(IJK,M) > ZERO) THEN
               I = I_OF(IJK)
               IMJK = IM_OF(IJK)
               IF (ROP_G(IJK)*AYZ_U(IMJK) > SMALL_NUMBER) THEN
                  B_M(IJK,M) = SQRT(B_M(IJK,M)/(ROP_G(IJK)*AVG_X_E(ZERO,ONE,I)*&
                     AYZ_U(IMJK)))
               ELSE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDIF
         ENDIF
      END DO

      RETURN
      END SUBROUTINE ADJUST_A_U_G

!// Comments on the modifications for DMP version implementation
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added sendrecv module for COMMunication
