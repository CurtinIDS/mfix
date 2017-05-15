!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_BIN_512                                            C
!  Purpose: write out an array in chunks of 512 bytes    (DP WORDS)    C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables: NWORDS, DS, L, NSEG, NREM, LC, N1, N2              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_BIN_512(IUNIT, ARRAY, N, NEXT_REC)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE machine
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      array to write out
      DOUBLE PRECISION ARRAY(*)
!
!                      output unit number
      INTEGER          IUNIT
!
!                      number of elements in ARRAY
      INTEGER          N
!
!                      next record number in direct access output file
      INTEGER          NEXT_REC
!
! local variables
!
!                      number of words for 512 bytes
      INTEGER          NWORDS
!
!                      loop counter
      INTEGER          L
!
!                      number of full 512 byte segments need to write N
!                      double precision words
      INTEGER          NSEG
!
!                      number of double precision words in the partially
!                      filled last record
      INTEGER          NREM
!
!                      loop counter
      INTEGER          LC
!
!                      write out array elements N1 to N2
      INTEGER          N1 , N2
!
!-----------------------------------------------
!
      NWORDS = NWORDS_DP
      IF (N <= NWORDS) THEN
         WRITE (IUNIT, REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      ENDIF
!
      NSEG = N/NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
!
! write out the full 512 byte segments
!
      DO LC = 1, NSEG
         WRITE (IUNIT, REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
      END DO
      IF (NREM /= 0) THEN
         WRITE (IUNIT, REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      ENDIF
!
      RETURN
      END SUBROUTINE OUT_BIN_512
