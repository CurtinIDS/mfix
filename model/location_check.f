!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: LOCATION_CHECK                                          !
!  Author: P. Nicoletti                               Date: 02-DEC-91  !
!                                                                      !
!  Purpose: Check calculated and given cell locations for consistency  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE LOCATION_CHECK(CELL_SPECIFIED, CELL_CALCULATED,       &
         COUNTER, MESSAGE)

! Global variables:
!---------------------------------------------------------------------//
! Runtime flags stating direction is not solved.
      USE geometry, only: NO_I, NO_J, NO_K
! Flag for reinitializing a run.
      use run, only: REINITIALIZING

! Module procedure for error message management.
      use error_manager

      IMPLICIT NONE

! Dummy arguments:
!---------------------------------------------------------------------//
! Cell index specified in the input file.
      INTEGER, INTENT(INOUT) :: CELL_SPECIFIED
! Cell index calculated for location coordinate.
      INTEGER, INTENT(IN) :: CELL_CALCULATED
! Index for BC, IC, or IS
      INTEGER, INTENT(IN) :: COUNTER
! Error message to print out
      CHARACTER(len=*) :: MESSAGE
!......................................................................!

! During a reinitializing, "old" data is in the I/J/K arrays. This check
! overwrites the old data with the current calculated values.
      IF(REINITIALIZING) CELL_SPECIFIED = CELL_CALCULATED

! Check that the cell_specified in the data input equals to the cell
! calculated.
      IF(CELL_SPECIFIED == CELL_CALCULATED) RETURN


      IF(NO_K .AND. (MESSAGE(6:6)=='b' .OR. MESSAGE(6:6)=='t')) RETURN
      IF(NO_J .AND. (MESSAGE(6:6)=='s' .OR. MESSAGE(6:6)=='n')) RETURN
      IF(NO_I .AND. (MESSAGE(6:6)=='w' .OR. MESSAGE(6:6)=='e')) RETURN

      CALL INIT_ERR_MSG('LOCATION_CHECK')

      WRITE(ERR_MSG, 1000) MESSAGE, COUNTER, CELL_SPECIFIED,           &
         CELL_CALCULATED
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1000 FORMAT('Error 1000: IC, BC, OR IS consistency error for: ',A,/,  &
      'IC/BC/IS No',5X,'= ',I6,/,'Cell specified',2x,'= ',I6,/,        &
      'Cell calculated',1x,'= ',I6)

      CALL FINL_ERR_MSG()

      RETURN
      END SUBROUTINE LOCATION_CHECK
