!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: IN_BIN_512                                             C
!  Purpose: read in an array in chunks of 512 bytes    (DP WORDS)      C
!                                                                      C
!  Author: P. Nicoletti                               Date: 03-JAN-92  C
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
!  Local variables: NWORDS, L, NSEG, NREM, LC, N1, N2                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
MODULE IN_BINARY_512

CONTAINS

      SUBROUTINE IN_BIN_512(IUNIT, ARRAY, N, NEXT_REC)
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
!                      number of words for 512 bytes (nwords * 4 = 512)
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
!-----------------------------------------------
!
!
      NWORDS = NWORDS_DP
      IF (N <= NWORDS) THEN
         READ (IUNIT, REC=NEXT_REC) (ARRAY(L),L=1,N)
         NEXT_REC = NEXT_REC + 1
         RETURN
      ENDIF
!
      NSEG = N/NWORDS
      NREM = MOD(N,NWORDS)
      N1 = 1
      N2 = NWORDS
!
! read the full 512 byte segments
!
      DO LC = 1, NSEG
         READ (IUNIT, REC=NEXT_REC) (ARRAY(L),L=N1,N2)
         N1 = N1 + NWORDS
         N2 = N2 + NWORDS
         NEXT_REC = NEXT_REC + 1
      END DO
      IF (NREM /= 0) THEN
         READ (IUNIT, REC=NEXT_REC) (ARRAY(L),L=N1,N)
         NEXT_REC = NEXT_REC + 1
      ENDIF

      RETURN
      END SUBROUTINE IN_BIN_512

      subroutine convert_from_io_dp(arr_io,arr_internal,nn)

      use geometry
      use indices
      USE compar
      USE functions

      implicit none

      double precision   arr_io(*) , arr_internal(*)
      integer            nn,i,j,k,ijk,ijk_io

      do k = 1,kmax2
         do j = 1,jmax2
            do i = 1,imax2
               ijk = funijk_gl(i,j,k)
               ijk_io = funijk_io(i,j,k)
               arr_internal(ijk) = arr_io(ijk_io)
            end do
         end do
      end do

      return
    end subroutine convert_from_io_dp

      subroutine convert_to_io_dp(arr_internal,arr_io,nn)

      use geometry
      use indices
      USE compar
      USE functions

      implicit none

      double precision   arr_io(*) , arr_internal(*)
      integer            nn,i,j,k,ijk,ijk_io

      do k = 1,kmax2
         do j = 1,jmax2
            do i = 1,imax2
               ijk  = funijk_gl(i,j,k)
               ijk_io = funijk_io(i,j,k)
               arr_io(ijk_io) = arr_internal(ijk)
            end do
         end do
      end do
!
      return
    end subroutine convert_to_io_dp

END MODULE IN_BINARY_512

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
