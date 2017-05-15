!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_HLOSS(HLOSS)                                       C
!  Purpose: Determine the total heat loss from the reactor             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAR-95  C
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
      SUBROUTINE GET_HLOSS(HLOSS)
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
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Total heat loss from the reactor
      DOUBLE PRECISION HLOSS, HLOSSm
!
!                      Indices
      INTEGER          M, IER
!
!-----------------------------------------------
!
      CALL GET_PHILOSS (T_G, K_G, BC_TW_G, BC_HW_T_G, BC_C_T_G, HLOSSM, IER)
      HLOSS = HLOSSM
      DO M = 1, MMAX
         CALL GET_PHILOSS (T_S(1,M), K_S(1,M), BC_TW_S(1,M), BC_HW_T_S(1,M), &
            BC_C_T_S(1,M), HLOSSM, IER)
         HLOSS = HLOSS + HLOSSM
      END DO
      RETURN
      END SUBROUTINE GET_HLOSS
