!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Init_Ab_m(A_m, b_m, IJKMAX2, M, IER)                   C                     C
!  Purpose:Initialiize the sparse matrix coefficients and the          C
!           source vector.                                             C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 16-MAY-96  C
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
      SUBROUTINE INIT_AB_M(A_M, B_M, IJKMAX2A, M)
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
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Phase index
      INTEGER          M
!
!                      Maximum dimension
      INTEGER          IJKMAX2A
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Source vector
      DOUBLE PRECISION b_m(DIMENSION_3, 0:DIMENSION_M)
!
!-----------------------------------------------
!
!      IJK = 1
      IF (IJKMAX2A > 0) THEN
!$omp parallel
!$omp sections
         A_M(:,bottom,M) = ZERO
!$omp section
         A_M(:,south,M) = ZERO
!$omp section
         A_M(:,west,M) = ZERO
!$omp section
         A_M(:,0,M) = -ONE
!$omp section
         A_M(:,east,M) = ZERO
!$omp section
         A_M(:,north,M) = ZERO
!$omp section
         A_M(:,top,M) = ZERO
!$omp section
         B_M(:,M) = ZERO
!$omp end sections
!$omp end parallel
      ENDIF
      RETURN
    END SUBROUTINE INIT_AB_M
