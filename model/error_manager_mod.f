! -*- f90 -*-
!----------------------------------------------------------------------!
! Module: ERROR_MANAGER                                                !
!                                                                      !
! Purpose: Unify error message handeling.                              !
!                                                                      !
!----------------------------------------------------------------------!
      MODULE ERROR_MANAGER

      use, intrinsic :: ISO_C_BINDING
      use exit, only: mfix_exit

      implicit none

! Interface
!---------------------------------------------------------------------//
      interface iVal
         module procedure iVal_int
         module procedure iVal_dbl
         module procedure iVal_log
      end interface

! Maximum number of lines a message can have before a flush is needed.
      INTEGER, PARAMETER :: LINE_COUNT  = 32
! Maximum number of characters per line.
      INTEGER, PARAMETER :: LINE_LENGTH = 256

! Character string for storing the error message.
      CHARACTER(LEN=LINE_LENGTH), DIMENSION(LINE_COUNT) :: ERR_MSG

! Depth that the current call tree can go.
      INTEGER, PARAMETER, PRIVATE :: MAX_CALL_DEPTH = 16
! Current call depth.
      INTEGER, PRIVATE :: CALL_DEPTH

! The name of the calling routine. Set by calling: INIT_ERR_MSG
      CHARACTER(LEN=128), DIMENSION(MAX_CALL_DEPTH), PRIVATE :: CALLERS

! Flag for writing messages to the screen.
      LOGICAL, PRIVATE :: SCR_LOG

! Error Flag.
      INTEGER :: IER_EM

      contains

!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_ERROR_MANAGER                                       !
!                                                                      !
! Purpose: Initialize the error manager. This routine also opens the   !
! .LOG file(s) based on user input settings.                           !
!......................................................................!
      SUBROUTINE INIT_ERROR_MANAGER

! Global Variables:
!---------------------------------------------------------------------//
! Name given to current run.
      use run, only: RUN_NAME
! Flag: All ranks report errors.
      use output, only: ENABLE_DMP_LOG
! Flag: My rank reports errors.
      use funits, only: DMP_LOG
! Flag: Provide the full log.
      use output, only: FULL_LOG
! Rank ID of process
      use compar, only: myPE
! Rank ID for IO handeling
      use compar, only: PE_IO
! Number of ranks in parallel run.
      use compar, only: numPEs
! File unit for LOG messages.
      use funits, only: UNIT_LOG
! Undefined character string.
      use param1, only: UNDEFINED_C

! Global Routine Access:
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! Log file name.
      CHARACTER(len=255) :: LOGFILE
      CHARACTER(len=255) :: FILE_NAME
! First non-blank character in run_name.
      INTEGER :: NB
! Integer error flag
      INTEGER :: IER(0:numPEs-1)

! Initialize the error flags.
      IER = 0
      IER_EM = 0
! Initialize the call tree depth.
      CALL_DEPTH = 0
! Clear the error message storage container.
      ERR_MSG = ''
! Clear the caller routine information.
      CALLERS = ''

! This turns on error messaging from all processes.
      DMP_LOG = (myPE == PE_IO) .OR. ENABLE_DMP_LOG
! Flag for printing screen messages.
      SCR_LOG = (myPE == PE_IO) .AND. FULL_LOG

! Verify the length of user-provided name.
      LOGFILE = ''
      NB = INDEX(RUN_NAME,' ')
! RUN_NAME length too short.
      IF(RUN_NAME == UNDEFINED_C .OR. NB <= 1) THEN
         IF(myPE  == PE_IO) WRITE (*, 1000) 'short'
         CALL MFIX_EXIT(myPE)
! RUN_NAME length too long.
      ELSEIF(NB + 10 > LEN(LOGFILE)) THEN
         IF(myPE == PE_IO) WRITE (*, 1000) 'long'
         CALL MFIX_EXIT(myPE)
! RUN_NAME legnth just right.
      ELSE
! Specify the .LOG file name based on MPI Rank extenion.
         IF(numPEs == 1 .OR. .NOT.ENABLE_DMP_LOG) THEN
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
      ENDIF

! Open the .LOG file. From here forward, all routines should store
! error messages (at a minimum) in the .LOG file.
      IF(DMP_LOG) THEN
         NB = len_trim(LOGFILE)+1
         CALL OPEN_FILE(LOGFILE, NB, UNIT_LOG, '.LOG', FILE_NAME,      &
            'APPEND', 'SEQUENTIAL', 'FORMATTED', 132,  IER(myPE))
      ENDIF

! Verify that the .LOG file was successfully opened. Otherwise, flag the
! error and abort.
      CALL GLOBAL_ALL_SUM(IER)
      IF(sum(IER) /= 0) THEN
         IF(myPE == PE_IO) WRITE(*,1001) trim(FILE_NAME)
         CALL MFIX_EXIT(myPE)
      ENDIF

      RETURN

 1000 FORMAT(2/,1X,70('*')/' From: INIT_ERROR_MANAGER',/               &
         ' Error 1000: RUN_NAME too ',A,'. Please correct the',        &
         ' mfix.dat file.',/1x,70('*'),2/)

 1001 FORMAT(2/,1X,70('*')/' From: INIT_ERROR_MANAGER',/               &
         ' Error 1001: Failed to open log file: ',A,/' Aborting run.'/,&
         1x,70('*'),2/)

      END SUBROUTINE INIT_ERROR_MANAGER

!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_ERR_MSG                                             !
!                                                                      !
! Purpose: Initialize the error manager for the local routine. This    !
! call is needed to set the caller routines name for error messages.   !
!......................................................................!
      SUBROUTINE INIT_ERR_MSG(CALLER)

! Rank ID of process
      use compar, only: myPE
! Flag: My rank reports errors.
      use funits, only: DMP_LOG
! File unit for LOG messages.
      use funits, only: UNIT_LOG

      implicit none

      CHARACTER(LEN=*), intent(IN) :: CALLER

! Verify that the maximum call dept will not be exceeded.  If so, flag
! the error and exit.
      IF(CALL_DEPTH + 1 > MAX_CALL_DEPTH) THEN
         IF(SCR_LOG) WRITE(*,1000) CALL_DEPTH
         IF(DMP_LOG) WRITE(UNIT_LOG,1000) CALL_DEPTH
         CALL SHOW_CALL_TREE
         CALL MFIX_EXIT(myPE)
      ELSE
! Store the caller routines name.
         CALL_DEPTH = CALL_DEPTH + 1
         CALLERS(CALL_DEPTH) = trim(CALLER)
      ENDIF

! Clear out the error manager.
      ERR_MSG=''

      RETURN

 1000 FORMAT(/1X,70('*')/' From: ERROR_MANAGER --> INIT_ERR_MSG',/     &
         ' Error 1000: Invalid ERROR_MANAGER usage. The maximum call', &
         ' depth ',/' was exceeded. The calls to INIT_ERR_MSG should', &
         ' have corresponding',/' calls to FINL_ERR_MSG. The current', &
         ' CALL tree depth is: ',I4)

      END SUBROUTINE INIT_ERR_MSG

!``````````````````````````````````````````````````````````````````````!
! Subroutine: FINL_ERR_MSG                                             !
!                                                                      !
! Purpose: Finalize the error manager. The call is needed to clear out !
! old information and unset the lock.                                  !
!......................................................................!
      SUBROUTINE FINL_ERR_MSG

! Rank ID of process
      use compar, only: myPE
! Flag: My rank reports errors.
      use funits, only: DMP_LOG
! File unit for LOG messages.
      use funits, only: UNIT_LOG

      implicit none

! Single line.
      CHARACTER(LEN=LINE_LENGTH) :: LINE
! Line length with trailing space removed.
      INTEGER :: LENGTH
! Line Counter
      INTEGER :: LC
! Number of non-empty lines.
      INTEGER :: COUNT

! The current calling routine.
      CHARACTER(LEN=128) :: CALLER

! Verify that at the INIT routine was called.
      IF(CALL_DEPTH < 1) THEN
         IF(SCR_LOG) WRITE(*,1000)
         IF(DMP_LOG) WRITE(UNIT_LOG,1000)
         CALL MFIX_EXIT(myPE)
      ELSE
! Store the current caller, clear the array position, and decrement
! the counter.
         CALLER = CALLERS(CALL_DEPTH)
         CALLERS(CALL_DEPTH) = ''
         CALL_DEPTH = CALL_DEPTH - 1
      ENDIF

! Verify that the error message container is empty.
      COUNT = 0
      DO LC = 1, LINE_COUNT
         LINE = ERR_MSG(LC)
         LENGTH = len_trim(LINE)
         IF(0 < LENGTH .AND. LENGTH < 256 ) COUNT = COUNT + 1
      ENDDO

! If the error message container is not empty, report the error, dump
! the error message and abort MFIX.
      IF(COUNT /= 0) THEN
         IF(SCR_LOG) WRITE(*,1001) trim(CALLER)
         IF(DMP_LOG) WRITE(UNIT_LOG,1001) trim(CALLER)
! Write out the error message container contents.
         DO LC = 1, LINE_COUNT
            LINE = ERR_MSG(LC)
            LENGTH = len_trim(LINE)
            IF(0 < LENGTH .AND. LENGTH < 256 ) THEN
               IF(SCR_LOG) WRITE(*,1002)LC, LENGTH, trim(LINE)
               IF(DMP_LOG) WRITE(UNIT_LOG,1002)LC, LENGTH, trim(LINE)
            ENDIF
         ENDDO
         IF(SCR_LOG) WRITE(*,1003)
         IF(DMP_LOG) WRITE(UNIT_LOG, 1003)
         CALL MFIX_EXIT(myPE)
      ENDIF

! This shouldn't be needed, but it doesn't hurt.
      ERR_MSG = ''

      RETURN

 1000 FORMAT(/1X,70('*')/' From: ERROR_MANAGER --> FINL_ERR_MSG',/     &
         ' Error 1000: Ivalid ERROR_MANAGER usage. A call to FINL_ERR',&
         '_MSG was',/' made while the call tree is empty. This can',   &
         ' occur if a call to',/' FINL_ERR_MSG was made without a',    &
         ' corresponding call to INIT_ERR_MSG.',/' Aborting MFIX.'/    &
         1x,70('*'),2/)

 1001 FORMAT(/1X,70('*')/' From: ERROR_MANAGER --> FINL_ERR_MSG',/     &
         ' Error 1001: Error container ERR_MSG not empty.',/           &
         ' CALLERS: ',A,2/' Contents:')

 1002 FORMAT(' LC ',I2.2,': LEN: ',I3.3,1x,A)

 1003 FORMAT(/,1x,'Aborting MFIX.',1x,70('*'),2/)

      END SUBROUTINE FINL_ERR_MSG

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      SUBROUTINE FLUSH_ERR_MSG(DEBUG, HEADER, FOOTER, ABORT, LOG, &
         CALL_TREE)

! Rank ID of process
      use compar, only: myPE
! Flag: My rank reports errors.
      use funits, only: DMP_LOG
! File unit for LOG messages.
      use funits, only: UNIT_LOG
! Flag to reinitialize the code.
      use run, only: REINITIALIZING

! Dummy Arguments:
!---------------------------------------------------------------------//
! Debug flag.
      LOGICAL, INTENT(IN), OPTIONAL :: DEBUG
! Flag to suppress the message header.
      LOGICAL, INTENT(IN), OPTIONAL :: HEADER
! Flag to suppress the message footer.
      LOGICAL, INTENT(IN), OPTIONAL :: FOOTER
! Flag to abort execution by invoking MFIX_EXIT.
      LOGICAL, INTENT(IN), OPTIONAL :: ABORT
! Flag to force (or override) writing data to the log file.
      LOGICAL, INTENT(IN), OPTIONAL :: LOG
! Provide the call tree in error message.
      LOGICAL, INTENT(IN), OPTIONAL :: CALL_TREE

! Local Variables:
!---------------------------------------------------------------------//
! Single line.
      CHARACTER(LEN=LINE_LENGTH) :: LINE
! Line length with trailing space removed.
      INTEGER :: LENGTH
! Index of last line in the message.
      INTEGER :: LAST_LINE
! Line Counter
      INTEGER :: LC
! Local debug flag.
      LOGICAL :: D_FLAG
! Local flag to suppress writing the header.
      LOGICAL :: H_FLAG
! Local flag to suppress writing the footer.
      LOGICAL :: F_FLAG
! Local abort flag.
      LOGICAL :: A_FLAG
! Local call tree flag.
      LOGICAL :: CT_FLAG
! Local flag to store output to UNIT_LOG
      LOGICAL :: UNT_LOG

! The current calling routine.
      CHARACTER(LEN=128) :: CALLER

! Set the abort flag. Continue running by default.
      IF(PRESENT(ABORT))THEN
         A_FLAG = ABORT
      ELSE
         A_FLAG = .FALSE.
      ENDIF

! Set the local debug flag. Suppress debugging messages by default.
      IF(PRESENT(DEBUG)) THEN
         D_FLAG = DEBUG
      ELSE
         D_FLAG = .FALSE.
      ENDIF

! Set the header flag. Write the header by default.
      IF(PRESENT(HEADER)) THEN
         H_FLAG = HEADER
      ELSE
         H_FLAG = .TRUE.
      ENDIF

! Set the footer flag. Write the footer by default.
      IF(PRESENT(FOOTER))THEN
         F_FLAG = FOOTER
      ELSE
         F_FLAG = .TRUE.
      ENDIF

! Set the call tree flag. Suppress the call tree by default.
      IF(PRESENT(LOG)) THEN
         UNT_LOG = DMP_LOG .AND. LOG
      ELSE
         UNT_LOG = DMP_LOG
      ENDIF

! Set the call tree flag. Suppress the call tree by default.
      IF(PRESENT(CALL_TREE)) THEN
         CT_FLAG = CALL_TREE
      ELSE
         CT_FLAG = .FALSE.
      ENDIF

! Write out header infomration.
      IF(H_FLAG) THEN
! Set the current caller.
         CALLER = CALLERS(CALL_DEPTH)
         IF(D_FLAG) THEN
            IF(SCR_LOG) WRITE(*,2000) trim(CALLER)
            IF(UNT_LOG) WRITE(UNIT_LOG,2000) trim(CALLER)
         ELSE
            IF(SCR_LOG) WRITE(*,1000) trim(CALLER)
            IF(UNT_LOG) WRITE(UNIT_LOG,1000) trim(CALLER)
         ENDIF
      ENDIF

! Find the end of the message.
      LAST_LINE = 0
      DO LC = 1, LINE_COUNT
         LINE = ERR_MSG(LC)
         LENGTH = len_trim(LINE)
         IF(0 < LENGTH .AND. LENGTH < 256 ) LAST_LINE = LC
      ENDDO

! Write the message body.
      IF(D_FLAG)THEN
         DO LC = 1, LINE_COUNT
            LINE = ERR_MSG(LC)
            LENGTH = len_trim(LINE)
            IF(LENGTH == 0) THEN
               IF(SCR_LOG) WRITE(*,2001) LC, LENGTH, "EMPTY."
               IF(UNT_LOG) WRITE(UNIT_LOG,2001) LC, LENGTH, "EMPTY."
            ELSEIF(LENGTH >=  LINE_LENGTH)THEN
               IF(SCR_LOG) WRITE(*,2001) LC, LENGTH, "OVERFLOW."
               IF(UNT_LOG) WRITE(UNIT_LOG,2001) LC, LENGTH, "OVERFLOW."
            ELSE
               IF(SCR_LOG) WRITE(*,2001) LC, LENGTH, trim(LINE)
               IF(UNT_LOG) WRITE(UNIT_LOG,2001) LC, LENGTH, trim(LINE)
            ENDIF
         ENDDO
      ELSE
         DO LC = 1, LAST_LINE
            LINE = ERR_MSG(LC)
            LENGTH = len_trim(LINE)
            IF(0 < LENGTH .AND. LENGTH < 256 ) THEN
               IF(SCR_LOG) WRITE(*,1001) trim(LINE)
               IF(UNT_LOG) WRITE(UNIT_LOG,1001) trim(LINE)
            ELSE
               IF(SCR_LOG) WRITE(*,"('  ')")
               IF(UNT_LOG) WRITE(UNIT_LOG,"('  ')")
            ENDIF
         ENDDO
         IF(LAST_LINE == 0) THEN
            IF(SCR_LOG) WRITE(*,"('  ')")
            IF(UNT_LOG) WRITE(UNIT_LOG,"('  ')")
         ENDIF
      ENDIF

! Print footer.
      IF(F_FLAG) THEN
         IF(D_FLAG) THEN
            IF(SCR_LOG) WRITE(*, 2002)
            IF(UNT_LOG) WRITE(UNIT_LOG, 2002)
         ELSE
            IF(SCR_LOG) WRITE(*, 1002)
            IF(UNT_LOG) WRITE(UNIT_LOG, 1002)
         ENDIF
      ENDIF

! Clear the message array.
      ERR_MSG=''

! Abort the run if specified.
      IF(A_FLAG) THEN
         IF(REINITIALIZING)THEN
            IER_EM = 1
         ELSE
            IF(D_FLAG) WRITE(*,3000) myPE
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      RETURN

 1000 FORMAT(2/,1x,70('*'),/' From: ',A)
 1001 FORMAT(1x,A)
 1002 FORMAT(1x,70('*'))

 2000 FORMAT(2/,'--- HEADER ---> ',70('*'),/'--- HEADER ---> From: ',A)
 2001 FORMAT('LC ',I2.2,': LEN: ',I3.3,1x,A)
 2002 FORMAT('--- FOOTER --->',1x,70('*'))

 3000 FORMAT(2x,'Rank ',I5,' calling MFIX_EXIT from FLUSH_ERR_MSG.')

      END SUBROUTINE FLUSH_ERR_MSG


!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      SUBROUTINE SHOW_CALL_TREE(HEADER, FOOTER)

! Flag: My rank reports errors.
      use funits, only: DMP_LOG
! File unit for LOG messages.
      use funits, only: UNIT_LOG

! Dummy Arguments:
!---------------------------------------------------------------------//
! Flag to suppress the message header.
      LOGICAL, INTENT(IN), OPTIONAL :: HEADER
! Flag to suppress the message footer.
      LOGICAL, INTENT(IN), OPTIONAL :: FOOTER

! Local Variables:
!---------------------------------------------------------------------//
! Local flag to suppress writing the header.
      LOGICAL :: H_FLAG
! Local flag to suppress writing the footer.
      LOGICAL :: F_FLAG
! Generic loop counters.
      INTEGER ::  LC, SL

! Set the header flag. Write the header by default.
      H_FLAG = merge(HEADER, .TRUE., PRESENT(HEADER))
! Set the footer flag. Write the footer by default.
      F_FLAG = merge(FOOTER, .TRUE., PRESENT(FOOTER))

! Header
      IF(H_FLAG) THEN
         IF(SCR_LOG) WRITE(*,1000)
         IF(DMP_LOG) WRITE(UNIT_LOG,1000)
      ENDIF

! Call Tree
      DO LC=1,MAX_CALL_DEPTH
         DO SL=1,LC
            IF(SCR_LOG) WRITE(*,1001,ADVANCE='NO')
            IF(DMP_LOG) WRITE(UNIT_LOG,1001,ADVANCE='NO')
         ENDDO
         IF(SCR_LOG) WRITE(*,1002,ADVANCE='YES') CALLERS(LC)
         IF(DMP_LOG) WRITE(UNIT_LOG,1002,ADVANCE='YES') CALLERS(LC)
      ENDDO

! Footer.
      IF(F_FLAG) THEN
         IF(SCR_LOG) WRITE(*,1003)
         IF(DMP_LOG) WRITE(UNIT_LOG,1003)
      ENDIF

      RETURN

 1000 FORMAT(2/,1x,70('*'),' CALL TREE INFORMATION')
 1001 FORMAT(' ')
 1002 FORMAT('> ',A)
 1003 FORMAT(/1x,70('*'))

      END SUBROUTINE SHOW_CALL_TREE

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      CHARACTER(len=32) FUNCTION iVar(VAR, i1, i2, i3)

      CHARACTER(len=*), intent(in) :: VAR

      INTEGER,  intent(in) :: i1
      INTEGER, OPTIONAL, intent(in) :: i2
      INTEGER, OPTIONAL, intent(in) :: i3

      CHARACTER(len=16) :: iASc
      CHARACTER(len=64) :: tVAR

      iASc=''; WRITE(iASc,*)i1
      tVar=''; WRITE(tVar,"(A,'(',A)") &
         trim(adjustl(VAR)), trim(adjustl(iASc))

      IF(PRESENT(i2))THEN
         iASc=''; WRITE(iASc,*)i2
         WRITE(tVar,"(A,',',A)") trim(tVar), trim(adjustl(iASc))
      ENDIF

      IF(PRESENT(i3))THEN
         iASc=''; WRITE(iASc,*)i3
         WRITE(tVar,"(A,',',A)") trim(tVar), trim(adjustl(iASc))
      ENDIF

      WRITE(tVar,"(A,')')") trim(tVar)

      iVar = trim(adjustl(tVar))

      RETURN
      END FUNCTION iVar

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      CHARACTER(len=32) FUNCTION iVal_int(VAL)
      INTEGER, intent(in) :: VAL

      CHARACTER(len=32) :: iASc

      WRITE(iASc,*) VAL
      iVal_int = trim(adjustl(iASc))

      END FUNCTION iVal_int

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      CHARACTER(len=32) FUNCTION iVal_dbl(VAL)
      DOUBLE PRECISION, intent(in) :: VAL

      CHARACTER(len=32) :: dASc

      IF(abs(VAL) < 1.0d-2 .AND. abs(VAL) < 1.0d2) THEN
         WRITE(dASc,"(F18.4)") VAL
      ELSE
         WRITE(dASc,"(G18.4)") VAL
      ENDIF

      iVal_dbl = trim(adjustl(dASc))

      END FUNCTION iVal_dbl

!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!......................................................................!
      CHARACTER(len=32) FUNCTION iVal_log(VAL)
      LOGICAL, intent(in) :: VAL

      IF(VAL) THEN
         iVal_log = ".TRUE."
      ELSE
         iVal_log = ".FALSE."
      ENDIF

      RETURN
      END FUNCTION iVal_log

!``````````````````````````````````````````````````````````````````````!
! Function: Reports TRUE if one or more processes set an ABORT flag.   !
!......................................................................!
      LOGICAL FUNCTION REINIT_ERROR()

! Global Routine Access:
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM

      CALL GLOBAL_ALL_SUM(IER_EM)
      REINIT_ERROR = (IER_EM /= 0)
      IER_EM = 0
      RETURN
      END FUNCTION REINIT_ERROR

      END MODULE ERROR_MANAGER
