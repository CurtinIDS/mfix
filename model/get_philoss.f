!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_PHILOSS(LOSS)                                      C
!  Purpose: Determine the total loss of sclar Phi through reactor wallsC
!                                                                      C
!  Author: M. Syamlal                                 Date: 9-MAY-97   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2, KMAX2, MMAX, ROP_s, DX, DY, DZ, C
!                        X, IMIN1, JMIN1. KMIN1                        C
!  Variables modified: I, J, K, M, IJK                                 C
!                                                                      C
!  Local variables:  None                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_PHILOSS(PHI,DIFF,BC_PHI_W,BC_HW_PHI,BC_C_PHI,LOSS)
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
      USE physprop
      USE geometry
      USE fldvar
      USE bc
      USE indices
      USE energy
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
!
!                      scalar diffusion coefficient
      DOUBLE PRECISION Diff(DIMENSION_3)
!
!                      scalar value at wall
      DOUBLE PRECISION BC_Phi_w(DIMENSION_BC)
!
!                      wall transfer coefficient
      DOUBLE PRECISION BC_hw_Phi(DIMENSION_BC)
!
!                      wall flux coefficient
      DOUBLE PRECISION BC_C_Phi(DIMENSION_BC)
!
!                      scalar loss thorugh wall
      DOUBLE PRECISION LOSS
!
!                      Indices
      INTEGER          L, I, J, K, IJK, IJKE, IJKW, IJKN, IJKS
      INTEGER          IJKT, IJKB
!
!                      Indices
      INTEGER          I1, I2, J1, J2, K1, K2
!-----------------------------------------------

      LOSS = ZERO

      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN
            IF (BC_TYPE_ENUM(L)==NO_SLIP_WALL .OR. BC_TYPE_ENUM(L)==FREE_SLIP_WALL&
                .OR. BC_TYPE_ENUM(L)==PAR_SLIP_WALL) THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF(.NOT.IS_ON_myPE_Owns(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
                           IJKE = EAST_OF(IJK)
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              LOSS = LOSS + DIFF(IJKE)*AYZ(IJK)*(PHI(EAST_OF(&
                                 IJKE))-PHI(IJKE))*ODX_E(IP1(I))
                           ELSE
                              LOSS = LOSS + DIFF(IJKE)*AYZ(IJK)*(BC_HW_PHI(L)*(&
                                 PHI(IJKE)-BC_PHI_W(L))-BC_C_PHI(L))
                           ENDIF
                        ELSE IF (FLUID_AT(WEST_OF(IJK))) THEN
                           IJKW = WEST_OF(IJK)
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              LOSS = LOSS - DIFF(IJKW)*AYZ(IJKW)*(PHI(IJKW)-PHI&
                                 (WEST_OF(IJKW)))*ODX_E(IM1(IM1(I)))
                           ELSE
                              LOSS = LOSS + DIFF(IJKW)*AYZ(IJKW)*(BC_HW_PHI(L)*&
                                 (PHI(IJKW)-BC_PHI_W(L))-BC_C_PHI(L))
                           ENDIF
                        ELSE IF (FLUID_AT(NORTH_OF(IJK))) THEN
                           IJKN = NORTH_OF(IJK)
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              LOSS = LOSS + DIFF(IJKN)*AXZ(IJK)*(PHI(NORTH_OF(&
                                 IJKN))-PHI(IJKN))*ODY_N(JP1(J))
                           ELSE
                              LOSS = LOSS + DIFF(IJKN)*AXZ(IJK)*(BC_HW_PHI(L)*(&
                                 PHI(IJKN)-BC_PHI_W(L))-BC_C_PHI(L))
                           ENDIF
                        ELSE IF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           IJKS = SOUTH_OF(IJK)
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              LOSS = LOSS + DIFF(IJKS)*AXZ(JM_OF(IJK))*(PHI(&
                                 IJKS)-PHI(SOUTH_OF(IJKS)))*ODY_N(JM1(JM1(J)))
                           ELSE
                              LOSS = LOSS + DIFF(IJKS)*AXZ(JM_OF(IJK))*(&
                                 BC_HW_PHI(L)*(PHI(IJKS)-BC_PHI_W(L))-BC_C_PHI(&
                                 L))
                           ENDIF
                        ELSE IF (FLUID_AT(TOP_OF(IJK))) THEN
                           IJKT = TOP_OF(IJK)
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              LOSS = LOSS + DIFF(IJKT)*AXY(IJK)*(PHI(TOP_OF(&
                                 IJKT))-PHI(IJKT))*OX(I)*ODZ_T(KP1(K))
                           ELSE
                              LOSS = LOSS + DIFF(IJKT)*AXY(IJK)*(BC_HW_PHI(L)*(&
                                 PHI(IJKT)-BC_PHI_W(L))-BC_C_PHI(L))
                           ENDIF
                        ELSE IF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           IJKB = BOTTOM_OF(IJK)
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              LOSS = LOSS + DIFF(IJKB)*AXY(KM_OF(IJK))*(PHI(&
                                 IJKB)-PHI(BOTTOM_OF(IJKB)))*OX(I)*ODZ_T(KM1(&
                                 KM1(K)))
                           ELSE
                              LOSS = LOSS + DIFF(IJKB)*AXY(KM_OF(IJK))*(&
                                 BC_HW_PHI(L)*(PHI(IJKB)-BC_PHI_W(L))-BC_C_PHI(&
                                 L))
                           ENDIF
                        ENDIF
                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF
      END DO
      RETURN
      END SUBROUTINE GET_PHILOSS

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
