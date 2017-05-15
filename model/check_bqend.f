!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: CHECK_BATCH_QUEUE_END                                   !
!  Author: A.Gel                                      Date:            !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE CHECK_BATCH_QUEUE_END(pEXIT_SIGNAL)

      use compar, only: PE_IO
      use error_manager
      use machine, only: WALL_TIME
      use run, only: get_tunit
      use mpi_utility, only: BCAST
      use run, only: BATCH_WALLCLOCK
      use run, only: TERM_BUFFER
      use time_cpu, only: WALL_START

      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: pEXIT_SIGNAL

! Logical flags for hault cases.
      LOGICAL :: USER_HAULT, WALL_HAULT
! Elapsed wall time, and fancy formatted buffer/batch queue times.
      DOUBLE PRECISION :: WALL_STOP, FANCY_BUFF, FANCY_BATCH
! Time units for formatted output.
      CHARACTER(LEN=4) :: WT_UNIT, BF_UNIT, BC_UNIT

! Calculate the current elapsed wall time.
      WALL_STOP = WALL_TIME()
      WALL_STOP = WALL_STOP - WALL_START

! Set flags for wall time exceeded or user specified hault.
      WALL_HAULT = ((WALL_STOP+TERM_BUFFER) >= BATCH_WALLCLOCK)
      INQUIRE(file="MFIX.STOP", exist=USER_HAULT)

! Report that the max user wall time was reached and exit.
      IF(WALL_HAULT) THEN
         CALL GET_TUNIT(WALL_STOP,WT_UNIT)
         FANCY_BUFF = TERM_BUFFER
         CALL GET_TUNIT(FANCY_BUFF, BF_UNIT)
         FANCY_BATCH = BATCH_WALLCLOCK
         CALL GET_TUNIT(FANCY_BATCH, BC_UNIT)
         WRITE(ERR_MSG, 1100) WALL_STOP, WT_UNIT, FANCY_BUFF, BF_UNIT, &
            FANCY_BATCH, BC_UNIT
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1100 FORMAT(2/,15('='),' REQUESTED CPU TIME LIMIT REACHED ',('='),/   &
         'Batch Wall Time:',3X,F9.2,1X,A,/'Elapsed Wall Time: ',F9.2,  &
         1X,A,/'Term Buffer:',7X,F9.2,A,/15('='),' REQUESTED CPU ',    &
         'TIME LIMIT REACHED ',('='))

! Report that the hault signal was detected.
      IF(USER_HAULT) THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 1200 FORMAT(2/,19('='),' MFIX STOP SIGNAL DETECTED ',19('='),/'MFIX.',&
         'STOP file detected in run directory. Terminating MFIX.',/    &
         'Please DO NOT FORGET to erase the MFIX.STOP file before ',   &
         'restarting',/19('='),'MFIX STOP SIGNAL DETECTED ',19('='))

! This routine was restructured so all MPI ranks to the same action. As
! a result, broadcasting the BATCHQ flag may not be needed.
      pEXIT_SIGNAL = (WALL_HAULT .OR. USER_HAULT) .OR. pEXIT_SIGNAL
      call bcast (pEXIT_SIGNAL,PE_IO)

      END SUBROUTINE CHECK_BATCH_QUEUE_END
