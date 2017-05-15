!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: COPY_A_g(A_VEL, A_m, IER)                              C
!  Purpose: Copy A_VEL_g to A_m                                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
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
      SUBROUTINE COPY_A_G(A_VEL, A_M)
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
      USE geometry
      USE indices
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
!                      Stored coefficients
      DOUBLE PRECISION A_VEL(DIMENSION_3, -3:3)
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Indices
      INTEGER          IJK
!-----------------------------------------------

      IJK = 1
      IF (IJKMAX2 > 0) THEN
         A_M(:IJKMAX2,west,0) = A_VEL(:IJKMAX2,west)
         A_M(:IJKMAX2,east,0) = A_VEL(:IJKMAX2,east)
         A_M(:IJKMAX2,south,0) = A_VEL(:IJKMAX2,south)
         A_M(:IJKMAX2,north,0) = A_VEL(:IJKMAX2,north)
         A_M(:IJKMAX2,bottom,0) = A_VEL(:IJKMAX2,bottom)
         A_M(:IJKMAX2,top,0) = A_VEL(:IJKMAX2,top)
         IJK = IJKMAX2 + 1
      ENDIF
      RETURN
      END SUBROUTINE COPY_A_G
!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: COPY_A_s(A_VEL, A_m, M, IER)                           C
!  Purpose: Copy A_VEL_s to A_m                                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
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
      SUBROUTINE COPY_A_S(A_VEL, A_M, M)
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
      USE geometry
      USE indices
      USE physprop
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
!                      Stored coefficients
      DOUBLE PRECISION A_VEL(DIMENSION_3, -3:3, DIMENSION_M)
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Indices
      INTEGER          IJK, M
!-----------------------------------------------

      IJK = 1
      IF (IJKMAX2 > 0) THEN
         A_M(:IJKMAX2,west,M) = A_VEL(:IJKMAX2,west,M)
         A_M(:IJKMAX2,east,M) = A_VEL(:IJKMAX2,east,M)
         A_M(:IJKMAX2,south,M) = A_VEL(:IJKMAX2,south,M)
         A_M(:IJKMAX2,north,M) = A_VEL(:IJKMAX2,north,M)
         A_M(:IJKMAX2,bottom,M) = A_VEL(:IJKMAX2,bottom,M)
         A_M(:IJKMAX2,top,M) = A_VEL(:IJKMAX2,top,M)
         IJK = IJKMAX2 + 1
      ENDIF
      RETURN
      END SUBROUTINE COPY_A_S

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
