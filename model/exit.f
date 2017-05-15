#include "version.inc"
MODULE exit
   CONTAINS
      SUBROUTINE MFIX_EXIT(myID, normal_termination)

! File unit for .OUT file
      USE funits, only : UNIT_OUT
! File unit for .LOG files
      USE funits, only : UNIT_LOG

      use compar
      use mpi_utility

      implicit none

! Rank ID
      INTEGER, INTENT(IN) :: myID
! if present, normal termination (won't print error message)
      LOGICAL, INTENT(IN), OPTIONAL :: normal_termination
! Logical showing that a file unit is open.
      LOGICAL :: isOpen
! The value passed via the dummy argument or the process ID.
      INTEGER :: myID_l
! Process ID (myPE) converted to a character string.
      CHARACTER(len=64) :: myID_c

! Set the ID of the caller.
      myID_c=''; WRITE(myID_c,*) myID

! Write out that this routine was called.
      IF (.not. present(normal_termination)) THEN
         IF(myPE == PE_IO) WRITE(*,1000)
         IF(DMP_LOG) THEN
            INQUIRE(UNIT=UNIT_LOG,OPENED=isOpen)
            IF(isOPEN) WRITE(UNIT_LOG,1001) trim(adjustl(myID_c))
         ENDIF
      ENDIF

! Terminate MPI.
      CALL exitMPI(myID_l)

! Close any open files.
      CALL CLOSE_FILE(UNIT_OUT)
      CALL CLOSE_FILE(UNIT_LOG)

! Last gasp...
      IF(myPE == PE_IO) WRITE(*,1002)

! Hard Stop.
      ERROR_STOP 1

 1000 FORMAT(2/,1x,70('*'),/' Fatal error reported on one or more',    &
        ' processes. The .LOG file',/' may contain additional',        &
        ' information about the failure.',/1x,70('*'))

 1001 FORMAT(2/,1x,70('*'),/' Fatal error reported on PE ',  &
         A,'. The .LOG file may contain',/' additional ',     &
        'information about the failure.',/1x,70('*'))

 1002 FORMAT(2/,1x,'Program Terminated.',2/)

      END SUBROUTINE MFIX_EXIT

!``````````````````````````````````````````````````````````````````````!
! Subroutine: CLOSE_FILE                                               !
!                                                                      !
! Purpose: Close a file if it is open.                                 !
!......................................................................!
      SUBROUTINE CLOSE_FILE(UNIT_l)

! Global Variables.
!---------------------------------------------------------------------//
! NONE

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: UNIT_l

! Local Variables.
!---------------------------------------------------------------------//
! Retruned status of the specifed file unit
      INTEGER :: IOS
! Logical indicating if the file is open
      LOGICAL :: FOPEN

! If the file is open...
      INQUIRE(UNIT=UNIT_l, OPENED=FOPEN, IOSTAT=IOS )
! Close it.
      IF(FOPEN) CLOSE(UNIT_l)

      RETURN
      END SUBROUTINE CLOSE_FILE

END MODULE exit
