!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RES_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)   C
!  Purpose: Create a RES file from the corresponding SPX files         C
!           NOTE : USES N_SPX                                          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 12-MAR-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IJKMAX2                                       C
!  Variables modified: IJK                                             C
!                                                                      C
!  Local variables: TIME_FOR_RES, L, PRINTED_MESS, ERROR               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE RES_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)
!
!
      Use param
      Use param1
      Use geometry
      Use indices
      Use run
      Use funits
      Use post3d
      Use physprop
      Use fldvar
      IMPLICIT NONE
!
      REAL     TIME_FOR_RES, TIME_FOUND
      LOGICAL  AT_EOF(*), READ_SPX(*)
      INTEGER  REC_POINTER(*), REC_POINTER_t(N_SPX)
      INTEGER  NSTEP_1 , ERROR_CODE
      CHARACTER(LEN=1) IANS
      REAL     TIME_REAL(*)
      LOGICAL  PRINTED_MESS(N_SPX)
      LOGICAL  ERROR
!
      INTEGER L, IJK
!
      ERROR = .FALSE.
!
      WRITE (*,'(A)',ADVANCE='NO') 'Enter time to retrieve from Spx files > '
      READ  (*,*) TIME_FOR_RES
!
      DO 10 L = 1,N_SPX
         READ_SPX(L) = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L) = .FALSE.
         PRINTED_MESS(L) = .FALSE.
10    CONTINUE
      DO 20 L = 1,N_SPX
         READ_SPX(L) = .TRUE.
         IF(L .GT. 1)READ_SPX(L-1)=.FALSE.
         CALL SEEK_TIME(READ_SPX, TIME_FOR_RES, REC_POINTER, TIME_FOUND)
         REC_POINTER_t(L) =  REC_POINTER(L)
         IF(TIME_FOUND .NE. TIME_FOR_RES) THEN
           WRITE(*,'(A, I2, A, G12.5)')&
             ' Did not find time in SPX file ' , L,&
            '.  Time found = ',TIME_FOUND
           ERROR = .TRUE.
         ELSE
           WRITE (*,'(A, I2)') ' Found time in SPX file ', L
           PRINTED_MESS(L) = .TRUE.
         ENDIF
20    CONTINUE
      DO 30 L = 1,N_SPX
         READ_SPX(L) = .TRUE.
30    CONTINUE
      CALL READ_SPX1(READ_SPX,REC_POINTER_t,AT_EOF, TIME_REAL,NSTEP_1)
!
      ERROR_CODE = 0
      IF (ERROR) THEN
         IF (MMAX.GT.1 .AND. .NOT.PRINTED_MESS(5)) THEN
            WRITE (*,*) ' '
!            WRITE (*,*) ' MMAX > 1 and did not find EP_s'
!            WRITE (*,*) ' can not create restart file'
!            RETURN
         END IF
         IF (.NOT.PRINTED_MESS(5) .AND. .NOT.PRINTED_MESS(1)) THEN
            WRITE (*,*) ' '
!            WRITE (*,*) ' Did not find EP_s or EP_g'
!            WRITE (*,*) ' can not create restart file'
!            RETURN
         END IF
         IF (.NOT.PRINTED_MESS(1) .AND. PRINTED_MESS(5)) THEN
            WRITE (*,*) ' '
            WRITE (*,*) ' Did not find  EP_g   will use EP_s'
            WRITE (*,*) ' '
         END IF
         IF (.NOT.PRINTED_MESS(5) .AND. PRINTED_MESS(1)) THEN
            WRITE (*,*) ' '
            WRITE (*,*) ' Did not find  EP_s   will use EP_g'
            ERROR_CODE = 1
            WRITE (*,*) ' '
         END IF
         WRITE (*,*) ' '
         WRITE (*,*) ' Create RESTART file anyway ? (Y for yes)'
         READ  (*,'(1A1)') IANS
         IF (IANS .NE. 'Y' .AND. IANS .NE. 'y') RETURN
         WRITE (*,'(A)',ADVANCE='NO') ' Time and DT ? '
         READ(*,*)TIME, DT

      END IF
!
!
!
      WRITE(*,*)
      WRITE(*,'(A)')' The old restart file will be over-written.'
      WRITE(*,'(A,G12.5,A)')' The records in SPx files after time = ',&
       TIME_FOR_RES, ' will be irrecoverably lost!'
      WRITE(*,'(A)',ADVANCE='NO')' Press Y to over write RES and SPX files '
      READ(*,'(1a1)') IANS
      IF (IANS .NE. 'Y' .AND. IANS .NE. 'y') RETURN
!
! CHANGE THE POINTER RECORD IN THE SPx files to correspond to TIME =
! TIME_FOR_RES.  Times greater than TIME_FOR_RES will be overwritten
! when MFIX is continued.  In SPX files with TIME_FOUND > TIME_FOR_RES
! time is reset as TIME_FOR_RES
!
      DO L = 1,N_SPX
        IF(ERROR .AND. PRINTED_MESS(L) ) THEN
          READ (UNIT_SPX+L,REC=REC_POINTER(L)) TIME_FOUND , NSTEP_1
        ENDIF
        IF(.NOT. PRINTED_MESS(L))THEN
          READ (UNIT_SPX+L,REC=REC_POINTER(L)) TIME_FOUND , NSTEP
          IF(TIME_FOUND .GT. TIME_FOR_RES) THEN
            WRITE (UNIT_SPX+L,REC=REC_POINTER(L))TIME_FOR_RES, NSTEP_1
          ENDIF
        ENDIF
        REC_POINTER(L) =  REC_POINTER(L) + NUM_REC(L)
        WRITE (UNIT_SPX+L,REC=3) REC_POINTER(L), NUM_REC(L)
      END DO
!
! set the time and cycle to the appropriate values in the RESTART file
!
      TIME  = TIME_FOR_RES
      NSTEP = NSTEP_1

      IF (ERROR_CODE.EQ.0) THEN
!
!       MAKE SURE THAT THE VOLUME FRACTIONS ADD UP TO 1.0 ... MIGHT NOT DUE TO
!       ROUND OFF ERROR IN GOING FROM DOUBLE TO SINGLE PRECISION
        DO IJK = 1,IJKMAX2
          CALL CALC_EP_g(IJK)
        END DO
      ELSE
         if(mmax .gt. 1.OR.ANY_SOLVE_ROs)Then
           WRITE(*,'(A)')' Cannot update SP5 file. Modify Post_mfix'
           return
         endif
         DO IJK = 1,IJKMAX2
            IF(EP_g(IJK) .NE. UNDEFINED) THEN
              ROP_s(IJK,1) = (ONE - EP_g(IJK)) * RO_s0(1)
            ELSE
              ROP_s(IJK,1) = ZERO
            ENDIF
         END DO
      END IF
!
! WRITE OUT THE RESTART FILE
!
      CALL WRITE_RES1
!
      WRITE(*,*)
      WRITE(*,'(A,A)',ADVANCE='NO')' RES and SPx files over written.',&
        ' Press any key to continue.'
      READ(*,'(1a1)') IANS
      RETURN
      END
