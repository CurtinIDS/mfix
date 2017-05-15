!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADJUST_LEQ(RESID, LEQ_IT, LEQ_METHOD, LEQI, LEQM, IER) C
!  Purpose: Adjusts liner equation solver method and iterations        C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 23-MAY-97  C
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
      SUBROUTINE ADJUST_LEQ(RESID, LEQ_ITL, LEQ_METHODL, LEQI, LEQM)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE toleranc
      USE leqsol
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER LEQ_ITL, LEQ_METHODL, LEQI, LEQM
      DOUBLE PRECISION RESID
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!                      Linear equation solver parameters used when
!                      a particular equation set has converged
      INTEGER, PARAMETER :: LEQ_IT_CONV = 5
      INTEGER, PARAMETER :: LEQ_METHOD_CONV = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!  The adjustment is disabled, because it was adversely affecting species
!  conservation
!      IF (LEQ_ADJUST .AND. RESID<=TOL_RESID*0.1) THEN
!         LEQM = LEQ_METHOD_CONV
!         LEQI = MIN(LEQ_IT_CONV,LEQ_ITL)
!      ELSE
         LEQM = LEQ_METHODL
         LEQI = LEQ_ITL
!      ENDIF
!
      RETURN
      END SUBROUTINE ADJUST_LEQ
