!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INIT_RESID(IER)                                        C
!  Purpose: Initialize residuals                                       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 8-JUL-96   C
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
      SUBROUTINE INIT_RESID()
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
      USE residual
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      residual number
      INTEGER          L
!
!                      Phase index
      INTEGER          M
!
!-----------------------------------------------
!
      L = 1
      IF (NRESID > 0) THEN
         M = 0
         IF (MMAX + 1 > 0) THEN
            RESID(:NRESID,:MMAX) = ZERO
            M = MMAX + 1
         ENDIF
         L = NRESID + 1
      ENDIF
      RETURN
      END SUBROUTINE INIT_RESID
