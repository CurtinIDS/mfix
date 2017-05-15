!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY_K                                            C
!  Purpose: print out a 2D (constant k-plane) array to standard output C
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
!  Variables modified:                                                 C
!                                                                      C
!  Local variables: NCOL, NTAB, LL1, LL2, LL3, IFORM1, IFORM2, IJK, IJ2
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_ARRAY_K(ARRAY)
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
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
      DOUBLE PRECISION ARRAY(*)
!
!                      number of columns to print out across the page
      INTEGER          NCOL
!
!                      number of tables the 2D array must be split into
!                      for printing
      INTEGER          NTAB
!
!                      loop indices
      INTEGER          LL1, LL2, LL3
!
!                      start and end 'I' for current table
      INTEGER          IFORM1 , IFORM2
!
!                      start 'IJ' and end 'IJ' for a given 'J' to print out
      INTEGER          IJK , IJ2
!

!
!-----------------------------------------------
!
! NOTE:  IF NCOL IS CHANGED TO A NUMBER GREATER THAN 30, THEN THE "30"
!        IN FORMATS 5050 AND 5100 MUST BE CHANGED TO THAT NUMBER.
!
      NCOL = 10
!// Adjust for cyclic in x direction
    IF(CYCLIC_X) then
      NTAB = (IMAX2-1)/NCOL + 1
    ELSE
      NTAB = IMAX2/NCOL + 1
    ENDIF

!// Adjust for cyclic in x direction
    IF(CYCLIC_X) then
      IF (MOD(IMAX2-1,NCOL) == 0) NTAB = NTAB - 1
    ELSE
      IF (MOD(IMAX2,NCOL) == 0) NTAB = NTAB - 1
    ENDIF
!
      DO LL1 = 1, NTAB

!// Adjust for cyclic in x direction
       IF(CYCLIC_X.AND.LL1.eq.1) then
         IFORM1 = 2 + NCOL*(LL1 - 1)
       ELSE
         IFORM1 = 1 + NCOL*(LL1 - 1)
       ENDIF

         IFORM2 = NCOL*LL1

!// Adjust for cyclic in x direction
       IF(CYCLIC_X) then
         IFORM2 = MIN(IFORM2,IMAX2-1)
       ELSE
         IFORM2 = MIN(IFORM2,IMAX2)
       ENDIF

         WRITE (UNIT_OUT, 5050) (LL3,LL3=IFORM1,IFORM2)
         DO LL2 = JMAX2, 1, -1
            IJK = funijk_io(IFORM1,LL2,1)
            IJ2 = funijk_io(IFORM2,LL2,1)
!efd
!            WRITE (UNIT_OUT, 5100) LL2, (ARRAY(LL3),LL3=IJK,IJ2)
            WRITE (UNIT_OUT, 5100) LL2,  &
                     (ARRAY(funijk_io(LL3,LL2,1)),LL3=IFORM1,IFORM2)
         END DO
      END DO
 5050 FORMAT(3X,'J',3X,'I=',3X,10(I3,9X))
 5100 FORMAT(1X,I3,3X,10(1PE12.4))
      RETURN
      END SUBROUTINE OUT_ARRAY_K
