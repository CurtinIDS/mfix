!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_FROM_SPX (AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)  C
!                                                                      C
!  Purpose: OUTARR type output for a specified time from the SPX files C
!           USES N_SPX                                                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 15-MAR-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: TIME_FOR_OUT, FILE_NAME, PRINTED_MESS, MORE_DATA   C
!                   ERROR, L                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER,TIME_REAL)
!
!
      Use param
      Use param1
      Use fldvar
      Use geometry
      Use indices
      Use run
      Use funits
      Use post3d
      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
      REAL      TIME_FOUND
      LOGICAL   AT_EOF(*), READ_SPX(*)
      INTEGER   REC_POINTER(*)
      INTEGER   NSTEP_1
      REAL      TIME_REAL(*)
      LOGICAL   PRINTED_MESS(N_SPX)
      LOGICAL   ERROR
!
      INTEGER L
!
      IF (.NOT.DO_XFORMS) THEN
         WRITE (*,'(A)',ADVANCE='NO') 'Enter time to retrieve from Spx files > '
         READ  (*,*) TIME_FOR_OUT
         CALL GET_FILE_NAME(TEMP_FILE)
      END IF
!
      OPEN (UNIT=UNIT_OUT,FILE=TEMP_FILE,STATUS='UNKNOWN',convert='big_endian')
      CALL READ_RES0
!
      DO 10 L = 1,N_SPX
         READ_SPX(L) = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L) = .FALSE.
         PRINTED_MESS(L) = .FALSE.
10    CONTINUE
      DO 20 L = 1,N_SPX
         IF(.NOT.SPX_OPEN(L))GOTO 20
         READ_SPX(L) = .TRUE.
         IF(L .GT. 1)READ_SPX(L-1)=.FALSE.
         CALL SEEK_TIME(READ_SPX, TIME_FOR_OUT, REC_POINTER, TIME_FOUND)
         IF(TIME_FOUND .GE. ZERO) THEN
           WRITE (*,*) ' Found time in SPX file ', L
           PRINTED_MESS(L) = .TRUE.
         ENDIF
20    CONTINUE
      DO 30 L = 1,N_SPX
        IF(SPX_OPEN(L)) READ_SPX(L) = .TRUE.
30    CONTINUE
      CALL READ_SPX1(READ_SPX,REC_POINTER,AT_EOF, TIME_REAL,NSTEP_1)
!
      ERROR = .FALSE.
      DO L = 1,N_SPX
         IF (.NOT.PRINTED_MESS(L)) THEN
            WRITE(UNIT_OUT,*) ' Did not find time in SPX file ' , L
            ERROR = .TRUE.
         END IF
      END DO
!
      TIME  = TIME_FOR_OUT
      CALL WRITE_OUT1
      CLOSE (UNIT=UNIT_OUT)
!
      IF (ERROR) THEN
        WRITE (*,*)' WARNING: Some variables were not found !'
      END IF
!
      RETURN
      END
