!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY_KC (ARRAY, K)                                C
!  Purpose: print out a 2D (constant k-plane) array to standard output C
!           (character)                                                C
!                                                                      C
!  Author: P.Nicoletti                                Date: 02-DEC-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2                                  C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: NCOL, NTAB, LL1, LL2, LL3, IFORM1, IFORM2, IJK, IJ2
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_ARRAY_KC(ARRAY)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE fldvar
      USE physprop
      USE indices
      USE funits
      USE compar
      USE mpi_utility
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      2D array to print out
      CHARACTER(LEN=4) :: ARRAY(*)
!
! local variables
!
!                      A line of characters to print
      CHARACTER(LEN=132) :: LINE
!
!                      number of columns to print out across the page
      INTEGER          NCOL
!
!                      number of tables the 2D array must be split into
!                      for printing
      INTEGER          NTAB
!
!                      loop indices
      INTEGER          LL1, LL2, LL3, LL4
!
!                      start and end 'I' for current table
      INTEGER          IFORM1 , IFORM2
!
!                      start 'IJ' and end 'IJ' for a given 'J' to print out
      INTEGER          IJK , IJ2
!
!-----------------------------------------------

!
! NOTE:  IF NCOL IS CHANGED TO A NUMBER OTHER THAN 24, THEN THE "24"
!        IN FORMATS 5050 AND 5100 MUST BE CHANGED TO THAT NUMBER.
!
      NCOL = 24
      NTAB = IMAX2/NCOL + 1
      IF (MOD(IMAX2,NCOL) == 0) NTAB = NTAB - 1
!
      DO LL1 = 1, NTAB
         IFORM1 = 1 + NCOL*(LL1 - 1)
         IFORM2 = NCOL*LL1
         IFORM2 = MIN(IFORM2,IMAX2)
         WRITE (UNIT_OUT, 5050) (LL3,LL3=IFORM1,IFORM2)
         DO LL2 = JMAX2, 1, -1
            IJK = funijk_io(IFORM1,LL2,1)
            IJ2 = funijk_io(IFORM2,LL2,1)
!efd
            WRITE (LINE, 5100) LL2, (ARRAY(LL3),LL3=IJK,IJ2)

            LL4 = 12 + (IFORM2 - IFORM1 + 1)*5
            WRITE (UNIT_OUT, '(A)') LINE(1:LL4)
         END DO
      END DO


 5050 FORMAT(3X,'J',3X,'I=',3X,24(I3,2X))
 5100 FORMAT(1X,I3,8X,24(A4,1X))
      RETURN
      END SUBROUTINE OUT_ARRAY_KC
