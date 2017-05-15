!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UNDEF_2_0 (Var, IER)                                   C
!  Purpose: change undefined values to zero.  Otherwise linear equationC
!           solver does not work                                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-JUL-96  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE UNDEF_2_0(VARDUM)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: VARDUM

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER :: IJK
!-----------------------------------------------
!
!  Local variables
!
!
!      IJK = 1
!      IF (IJKEND3 > 0) THEN
!         WHERE (VAR(:IJKMAX2) == UNDEFINED) VAR(:IJKMAX2) = ZERO
!         WHERE (VARDUM(IJKSTART3:IJKEND3) == UNDEFINED) VARDUM(IJKSTART3:IJKEND3) = ZERO
         WHERE (VARDUM(1:IJKEND3) == UNDEFINED) VARDUM(1:IJKEND3) = ZERO
!         IJK = IJKMAX2 + 1
!      ENDIF
      RETURN
      END SUBROUTINE UNDEF_2_0

!// Comments on the modifications for DMP version implementation
!//PG Changed local variable name from VAR to VARDUM, due to conflict in PG
!// 120 Replaced the index for initialization :IJKMAX2 --> VARDUM(1:IJKEND3)
