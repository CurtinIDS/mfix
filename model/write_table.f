!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_TABLE (LEGEND, ARRAY, DIST_MIN, LSTART, LEND)    C
!  Purpose: To write a table of DX, DY, DZ, and cell wall locations    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 09-JAN-92  C
!  Reviewer: S. Venkatesan                            Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: NROW, L, L1, L2, L3, DIST, ARRAY1, ARRAY3          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_TABLE(LEGEND, ARRAY, DIST_MIN, LSTART, LEND)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE funits
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Legend
      CHARACTER(LEN=*)    LEGEND(3)
!
!                      DX, DY, or DZ Array to be written


!                      Starting array index
      INTEGER          LSTART
!
!                      Ending array index
      INTEGER          LEND
!//EFD Nov/11, avoid use of (*)
!//      DOUBLE PRECISION ARRAY(*)
      DOUBLE PRECISION ARRAY((LSTART-1):(LEND+1))
!
!                      Starting value of distance
      DOUBLE PRECISION DIST_MIN
!
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!
!                      Number of columns in the table.  When this is changed
!                      remember to change the FORMAT statement also.
      INTEGER, PARAMETER :: NCOL = 5
!
!                      Some dimension large enough for I, J, and K.
      INTEGER, PARAMETER :: DIMENSION_1 = MAX(DIM_I, DIM_J, DIM_K)

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Indices
      INTEGER          ARRAY1(DIMENSION_1)
!
!                      Array3 to be written
      DOUBLE PRECISION ARRAY3(DIMENSION_1)
!
!                      Number of rows
      INTEGER          NROW
!
!                      Temporary storage for distance calculation
      DOUBLE PRECISION DIST
!
!                      Local array indices
      INTEGER          L, L1, L2, L3
!-----------------------------------------------
!
!
!  Fill arrays 1 and 3
!
      DIST = DIST_MIN
      DO L = LSTART, LEND
         ARRAY1(L) = L
         ARRAY3(L) = DIST
         IF (L < LEND) DIST = DIST + ARRAY(L+1)
      END DO
      NROW = (LEND - LSTART + 1)/NCOL
!
      L2 = LSTART - 1
      DO L = 1, NROW
         L1 = L2 + 1
         L2 = L1 + NCOL - 1
         WRITE (UNIT_OUT, 1010) LEGEND(1), (ARRAY1(L3),L3=L1,L2)
         WRITE (UNIT_OUT, 1020) LEGEND(2), (ARRAY(L3),L3=L1,L2)
         WRITE (UNIT_OUT, 1030) LEGEND(3), (ARRAY3(L3),L3=L1,L2)
      END DO
      IF (NROW*NCOL < LEND - LSTART + 1) THEN
         L1 = L2 + 1
         L2 = LEND
         WRITE (UNIT_OUT, 1010) LEGEND(1), (ARRAY1(L3),L3=L1,L2)
         WRITE (UNIT_OUT, 1020) LEGEND(2), (ARRAY(L3),L3=L1,L2)
         WRITE (UNIT_OUT, 1030) LEGEND(3), (ARRAY3(L3),L3=L1,L2)
      ENDIF
      RETURN
!
 1010 FORMAT(7X,A3,2X,5(4X,I3,5X,1X))
 1020 FORMAT(7X,A3,2X,5(G12.5,1X))
 1030 FORMAT(7X,A3,2X,5(G12.5,1X),/)
      END SUBROUTINE WRITE_TABLE
