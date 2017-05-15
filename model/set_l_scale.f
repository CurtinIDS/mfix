!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_L_scale                                            C
!  Purpose: Initialize length scale for turbulence model               C
!                                                                      C
!  Author: W. Sams                                    Date: 04-MAY-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified: L_scale0, L_scale                               C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_L_SCALE
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
      USE constant
      USE visc_g
      USE geometry
      USE indices
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IJK
!-----------------------------------------------
!
      IJK = 1

!     IF (IJKMAX2 > 0) THEN
         L_SCALE(IJKSTART3:IJKEND3) = L_SCALE0
!        IJK = IJKMAX3 + 1
!     ENDIF
      RETURN
      END SUBROUTINE SET_L_SCALE

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!// 120 Replaced the index for initialization : :IJKMAX2 --> L_SCALE(IJKSTART3:IJKEND3)
