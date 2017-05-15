!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: OPEN_FILES                                             !
!  Author: P. Nicoletti                               Date: 12-DEC-91  !
!                                                                      !
!  Purpose: open all the files for this run                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE OPEN_FILES(RUN_NAME, RUN_TYPE, N_SPX)

      USE machine
      USE funits
      USE compar
      USE cdist

      use error_manager

      IMPLICIT NONE

! Error index: 0 - no error, 1 could not open file
      INTEGER :: IER(0:numPEs-1)
! RUN_NAME (as specified in input file)
      CHARACTER(LEN=*) :: RUN_NAME
! Run_type (as specified in input file)
      CHARACTER(LEN=*) :: RUN_TYPE
! Number of single precision output files (param.inc)
      INTEGER :: N_SPX
! local variables
      CHARACTER(len=4) :: EXT
! run_name + extension
      CHARACTER(len=255) :: FILE_NAME
! Loop counter
      INTEGER :: LC
! index to first blank character in run_name
      INTEGER :: NB
      CHARACTER(len=35) :: EXT_END
      CHARACTER(len=10) :: CSTATUS
! Character error code.
      CHARACTER(len=32) :: CER
!-----------------------------------------------


! Initialize the error manager.
      CALL INIT_ERR_MSG("OPEN_FILES")

! Initialize the error flag array.
      IER = 0

! Initialize the generic SPx extension.
      EXT = '.SPx'

! Generic SPx end characters in order.
      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'

! Get the length of RUN_NAME. Note that len_trim would allow the
! name to still contain spaces. The following approach truncates
! RUN_NAME at the first blank character.
      NB = INDEX(RUN_NAME,' ')

! Only PE_IO opens the RUN_NAME.OUT file.
      IF(myPE == PE_IO) CALL OPEN_FILE (RUN_NAME, NB, UNIT_OUT, '.OUT',&
         FILE_NAME, 'UNKNOWN', 'SEQUENTIAL','FORMATTED',132, IER(myPE))

! Check if there was an error opening the file.
      IF(ERROR_OPENING(IER)) THEN
         WRITE(ERR_MSG,3000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! Open the RES and SPx files. By default, only PE_IO opens these files,
! but all ranks open a rank-specific copy for distributed IO runs.
      SELECT CASE (TRIM(RUN_TYPE))

! Open the RES and SPx files for a new run.
!......................................................................
      CASE ('NEW')

         IF(myPE==PE_IO .OR.  bDist_IO) THEN

! Open the RES file.
            CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, &
               'NEW', 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE))
! Report errors.
            IF (IER(myPE) == 100) THEN
               WRITE(ERR_MSG, 1000)'RES', 'NEW', trim(FILE_NAME)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ELSEIF(IER(myPE) /= 0) THEN
               CER=''; WRITE(CER,*)
               WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ENDIF

! Open the SPx files.
            DO LC = 1, N_SPX
               EXT(4:4) = ext_end(LC:LC)
               CALL OPEN_FILE(RUN_NAME, NB, UNIT_SPX+LC, EXT,FILE_NAME,&
                  'NEW', 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE))
! Report errors.
               IF (IER(myPE) == 100) THEN
                  WRITE(ERR_MSG, 1000)EXT(2:4), 'NEW', trim(FILE_NAME)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ELSEIF(IER(myPE) /= 0) THEN
                  CER=''; WRITE(CER,*)
                  WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ENDIF
            ENDDO
         ENDIF


! Open the RES and SPx files for a typical restart run.
!......................................................................
      CASE ('RESTART_1')

! Open the RES file.
         IF(myPE == PE_IO .or. bDist_IO) THEN
            CALL OPEN_FILE(RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME,  &
               'OLD', 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE))
! Report errors.
            IF (IER(myPE) == 101) THEN
               WRITE(ERR_MSG, 1001)'RES', 'RESTART_1',trim(FILE_NAME)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ELSEIF(IER(myPE) /= 0) THEN
               CER=''; WRITE(CER,*)
               WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ENDIF

! Open the SPx files.
            DO LC = 1, N_SPX
               EXT(4:4) = EXT_END(LC:LC)
               CALL OPEN_FILE (RUN_NAME,NB, UNIT_SPX+LC,EXT, FILE_NAME,&
                  'OLD', 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE))
! Report errors.
               IF (IER(myPE) == 101) THEN
                  WRITE(ERR_MSG, 1001) EXT(2:4), 'RESTART_1',         &
                     trim(FILE_NAME)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ELSEIF(IER(myPE) /= 0) THEN
                  CER=''; WRITE(CER,*)
                  WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ENDIF
            END DO
         ENDIF


! Open the RES and SPx files for a typical restart run.
!......................................................................
      CASE ('RESTART_2')
! Open the RES file.
         CSTATUS = 'OLD'
         IF(myPE == PE_IO .OR. bDist_IO) THEN
            IF(bStart_with_one_res) CSTATUS = 'UNKNOWN'
            CALL OPEN_FILE (RUN_NAME, NB, UNIT_RES, '.RES', FILE_NAME, &
               CSTATUS,'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE))
! Report errors.
            IF (IER(myPE) == 101) THEN
               WRITE(ERR_MSG, 1001)'RES', 'RESTART_2',trim(FILE_NAME)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ELSEIF(IER(myPE) /= 0) THEN
               CER=''; WRITE(CER,*)
               WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
               CALL FLUSH_ERR_MSG
               GO TO 100
            ENDIF

! Open the SPx files.
            DO LC = 1, N_SPX
               EXT(4:4) = EXT_END(LC:LC)
               CALL OPEN_FILE (RUN_NAME,NB,UNIT_SPX+LC, EXT, FILE_NAME,&
                  'NEW' , 'DIRECT', 'UNFORMATTED', OPEN_N1, IER(myPE))
! Report errors.
               IF (IER(myPE) == 100) THEN
                  WRITE(ERR_MSG, 1000)EXT(2:4), 'RESTART_2',          &
                     trim(FILE_NAME)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ELSEIF(IER(myPE) /= 0) THEN
                  CER=''; WRITE(CER,*)
                  WRITE(ERR_MSG, 2000) trim(FILE_NAME), trim(CER)
                  CALL FLUSH_ERR_MSG
                  GO TO 100
               ENDIF
            END DO
         ENDIF

      CASE DEFAULT
         WRITE(ERR_MSG, 3000)
         CALL FLUSH_ERR_MSG
         GO TO 100

      END SELECT

! If an error was detected, abort the run.
  100 IF(ERROR_OPENING(IER)) CALL MFIX_EXIT(myPE)

! Initialize the error manager.
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: ',A,' file detected but RUN_TYPE=',A/,       &
         'Cannot open file: ',A)

 1001 FORMAT('Error 1001: ',A,' file missing for RUN_TYPE=',A/,        &
         'Cannot open file: ',A)

 2000 FORMAT('Error 2000: Unknown error opening file ',A,/             &
         'Error code: ',A)

 3000 FORMAT('Error 3000: Unknown run type: ',A)


      CONTAINS


!``````````````````````````````````````````````````````````````````````!
! FUNCTION: ERROR_OPENING                                              !
! Purpose: Collect the error flags from all processes and sum them.    !
! RESULT: .TRUE.  :: Sum of IER over all processes is non-zero.        !
!         .FALSE. :: GLOBAL_ALL_SUM is zero.                           !
!                                                                      !
!......................................................................!
      LOGICAL FUNCTION ERROR_OPENING(IER_l)

! MPI Wrapper function.
      use mpi_utility, only: GLOBAL_ALL_SUM

! Array containing error flags from all ranks.
      INTEGER, INTENT(IN) :: IER_L(0:numPEs-1)
! Initialize error flags.
      ERROR_OPENING = .FALSE.
! Globally collect flags.
      CALL GLOBAL_ALL_SUM(IER)
! Report errors.
      IF(sum(IER_l) /= 0) ERROR_OPENING = .TRUE.

      RETURN
      END FUNCTION ERROR_OPENING

      END SUBROUTINE OPEN_FILES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: OPEN_PE_LOG                                            !
!  Author: P. Nicoletti                               Date: 12-DEC-91  !
!                                                                      !
!  Purpose: Every rank open a .LOG file for domain specific errors.    !
!  This routine should only be invoked before writing to the log and   !
!  exiting.                                                            !
!                                                                      !
!  This routine only opens files when the following are met:           !
!    (1) MFIX is run in DMP parallel (MPI)                             !
!    (2) ENABLE_DMP_LOG is not set in the mfix.dat file.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE OPEN_PE_LOG(IER)

! Global Variables:
!---------------------------------------------------------------------//
! File unit for LOG files.
      USE funits, only: UNIT_LOG
! User specifed run name
      USE run, only: RUN_NAME
! MPI Rank of current process.
      USE compar, only: myPE
! Total number of MPI ranks.
      USE compar, only: numPEs
! Flag: My rank reports errors.
      use funits, only: DMP_LOG
! Flag: The log had to be opened.
      use funits, only: LOG_WAS_CLOSED

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Error index.
      INTEGER, INTENT(inout) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Log file name.
      CHARACTER(len=255) :: LOGFILE
      CHARACTER(len=255) :: FILE_NAME
! Flag for LOG files that are already open.
      LOGICAL :: DO_NOTHING
! Index of first blank character in RUN_NAME
      INTEGER :: NB
!......................................................................!


! Enable output from this rank.
      DMP_LOG = .TRUE.

! Return to the caller if this rank is already connect to a log file.
      INQUIRE(UNIT=UNIT_LOG, OPENED=DO_NOTHING)
      IF(DO_NOTHING) RETURN

! Flag that the log had to be opened.
      LOG_WAS_CLOSED = .TRUE.

! Verify the length of user-provided name.
      LOGFILE = ''
      NB = INDEX(RUN_NAME,' ')

! Specify the .LOG file name based on MPI Rank extenion.
      IF(numPEs == 1) THEN
         WRITE(LOGFILE,"(A)")RUN_NAME(1:(NB-1))
      ELSEIF(numPEs <    10) THEN
         WRITE(LOGFILE,"(A,'_',I1.1)") RUN_NAME(1:(NB-1)), myPE
      ELSEIF(numPEs <   100) THEN
         WRITE(LOGFILE,"(A,'_',I2.2)") RUN_NAME(1:(NB-1)), myPE
      ELSEIF(numPEs <  1000) THEN
         WRITE(LOGFILE,"(A,'_',I3.3)") RUN_NAME(1:(NB-1)), myPE
      ELSEIF(numPEs < 10000) THEN
         WRITE(LOGFILE,"(A,'_',I4.4)") RUN_NAME(1:(NB-1)), myPE
      ELSE
         WRITE(LOGFILE,"(A,'_',I8.8)") RUN_NAME(1:(NB-1)), myPE
      ENDIF

! Open the .LOG file. From here forward, all routines should store
! error messages (at a minimum) in the .LOG file.
      NB = len_trim(LOGFILE)+1
      CALL OPEN_FILE(LOGFILE, NB, UNIT_LOG, '.LOG', FILE_NAME,         &
         'APPEND', 'SEQUENTIAL', 'FORMATTED', 132,  IER)

      RETURN
      END SUBROUTINE OPEN_PE_LOG


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CLOSE_PE_LOG                                           !
!  Author: P. Nicoletti                               Date: 12-DEC-91  !
!                                                                      !
!  Purpose: Every rank open a .LOG file for domain specific errors.    !
!  This routine should only be invoked before writing to the log and   !
!  exiting.                                                            !
!                                                                      !
!  This routine only opens files when the following are met:           !
!    (1) MFIX is run in DMP parallel (MPI)                             !
!    (2) ENABLE_DMP_LOG is not set in the mfix.dat file.               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CLOSE_PE_LOG

! Global Variables:
!---------------------------------------------------------------------//
! File unit for LOG files.
      USE funits, only: UNIT_LOG
! Flag: My rank reports errors.
      use funits, only: DMP_LOG
! Flag: The log had to be opened.
      use funits, only: LOG_WAS_CLOSED

      IMPLICIT NONE

!......................................................................!


! The log had to be opened for global error.
      IF(LOG_WAS_CLOSED) THEN
! Reset the flag.
         LOG_WAS_CLOSED = .FALSE.
! Disable output from this rank and close connection to LOG.
         DMP_LOG = .FALSE.
! Return to the caller if this rank is already connect to a log file.
         CLOSE(UNIT_LOG)
      ENDIF

      RETURN
      END SUBROUTINE CLOSE_PE_LOG
