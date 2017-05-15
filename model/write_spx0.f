!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_SPX0(L, unit_add)                                !
!  Author: P. Nicoletti                               Date: 13-DEC-91  !
!                                                                      !
!  Purpose: Write out the initial single percision data files.         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_SPX0(L, UNIT_ADD)

! Global Variables:
!---------------------------------------------------------------------//
! User defined run_name
      use run, only: RUN_NAME
! Calendar information for start of run.
      use machine, only: ID_MONTH, ID_DAY, ID_YEAR
! Time information for start of run.
      use machine, only: ID_HOUR, ID_MINUTE, ID_SECOND
! Base file unit for SPx files.
      use funits, only: UNIT_SPX
! Flag: Use distributed I/O
      USE cdist, only: bDIST_IO
! Rank of current process.
      use compar, only: myPE
! Rank assigned to serial I/O
      use compar, only: PE_IO

      IMPLICIT NONE

! Passed Variables:
!----------------------------------------------------------------------!
! Index of SPx file.
      INTEGER, INTENT(IN) :: L
! Offset for use in post_mfix
      INTEGER, INTENT(IN) :: UNIT_ADD

! Local Variables:
!----------------------------------------------------------------------!
! File version ID
      CHARACTER(LEN=512) :: VERSION
! UNIT_SPX + offset from post_mfix
      INTEGER :: USPX
! Generic SPx end characters.
      CHARACTER(len=15), PARAMETER :: EXT_END = '123456789ABCDEF'

!......................................................................!

! Serial I/O: only PE_IO writes to a single file [default]
! Distributed I/O: all ranks write their own output file.
      IF(myPE /= PE_IO .AND. .NOT.bDIST_IO) RETURN

! Construct the file version string.
      VERSION = 'SPx = 02.00'
      WRITE(VERSION(3:3),"(A1)") EXT_END(L:L)

! Calculate the SPx file unit
      USPX = UNIT_SPX + UNIT_ADD + L

! Write the SPx file header.
      WRITE(USPX, REC=1) VERSION
      WRITE(USPX, REC=2) RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR,          &
         ID_HOUR, ID_MINUTE, ID_SECOND

!  The first field contains the pointer to the next record.
!  The second field contains the number of records written each time step
!  (The 4 and -1 are overwritten in WRITE_SPX1)
      WRITE (USPX, REC=3) 4, -1

      IF(UNIT_ADD == 0) FLUSH(USPX)

      RETURN
      END SUBROUTINE WRITE_SPX0
