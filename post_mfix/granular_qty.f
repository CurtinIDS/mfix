!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_GRANULAR_QTY                                       C
!  Purpose: calculate MU_s/EP_s from the RES file                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-AUG-94  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX, IJKMAX2, EP_s, RUN_NAME, TIME           C
!  Variables modified: M, IJK                                          C
!                                                                      C
!  Local variables: FILE_NAME, ARRAY                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_GRANULAR_QTY

      USE compar
      USE constant
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE machine
      USE param
      USE param1
      USE physprop
      USE rdf
      USE run, only: any_solve_ros, nstep, run_name, time
      USE visc_s

      IMPLICIT NONE
      INCLUDE 'xforms.inc'

!  Local variables
!
      INTEGER      MM, M, IER, IJK
!
!                file version ID
      CHARACTER(LEN=512) :: VERSION
      INTEGER     REC_POINTER(N_SPX), L, NEXT_REC, NUM_REC, NSTEP_1
      LOGICAL     READ_SPX(N_SPX), AT_EOF(N_SPX)
      CHARACTER   ANSWER
      REAL        TIME_REAL(N_SPX), TIME_NOW
      REAL        TIME_IN_RES
      DOUBLE PRECISION K_1m, Tmp(DIMENSION_3)

      CALL READ_RES1
      TIME_IN_RES = TIME
!
      IF (DO_XFORMS) THEN
         C_E = E_PASS
      ELSE
         IF (C_E .EQ. UNDEFINED) THEN
            WRITE (*,'(A)',ADVANCE='NO')&
              ' Enter Coefficient of restitution value (e): '
            READ  (*,*) C_e
         END IF
      END IF
!
      IF(.NOT.DO_XFORMS)THEN
         CALC_GT = .FALSE.
         CALC_GV = .FALSE.
         WRITE (*,'(A)',ADVANCE='NO')&
           ' Do you need granular temperature ? (Y/N) '
         READ (*,'(A)') ANSWER
         IF(ANSWER .EQ. 'Y' .OR. ANSWER .EQ. 'y') CALC_GT = .TRUE.
      ENDIF
      IF (CALC_GT) THEN
        IF (.NOT.DO_XFORMS) CALL GET_FILE_NAME(TEMP_FILE2)
        OPEN (UNIT=70,FILE=TEMP_FILE2,STATUS='NEW',RECL=128,&
            ACCESS='DIRECT',FORM='UNFORMATTED',ERR=1000,convert='big_endian')
!
        VERSION = 'SP1 = 01.00'
        WRITE (70,REC=1) VERSION
        WRITE (70,REC=2)RUN_NAME,ID_MONTH,ID_DAY,ID_YEAR,&
                ID_HOUR,ID_MINUTE,ID_SECOND
        WRITE (70,REC=3) 4, -1
      ENDIF
!
      IF(.NOT.DO_XFORMS)THEN
         CALC_GV = .FALSE.
         WRITE (*,'(A)',ADVANCE='NO')&
           ' Do you need granular viscosity ? (Y/N) '
         READ (*,'(A)') ANSWER
         IF(ANSWER .EQ. 'Y' .OR. ANSWER .EQ. 'y') CALC_GV = .TRUE.
      ENDIF
      IF (CALC_GV) THEN
        IF (.NOT.DO_XFORMS) CALL GET_FILE_NAME(TEMP_FILE)
        OPEN (UNIT=71,FILE=TEMP_FILE,STATUS='NEW',RECL=128,&
            ACCESS='DIRECT',FORM='UNFORMATTED',ERR=1000,convert='big_endian')
!
        VERSION = 'SP1 = 01.00'
        WRITE (71,REC=1) VERSION
        WRITE (71,REC=2)RUN_NAME,ID_MONTH,ID_DAY,ID_YEAR,&
                ID_HOUR,ID_MINUTE,ID_SECOND
        WRITE (71,REC=3) 4, -1
      ENDIF
!
      IF(CALC_GT .OR. CALC_GV) THEN
        IF(MMAX .GT. 1) THEN
          IF (.NOT.DO_XFORMS) THEN
             WRITE(*, '(A)',ADVANCE='NO')' Enter solids phase number: '
             READ(*, *)M_USE
          END IF
        ELSE
          M_USE = 1
        ENDIF
      ELSE
        RETURN
      ENDIF
      MM = M_USE
!
10    IF (.NOT.DO_XFORMS) THEN
         WRITE (*,'(A)',ADVANCE='NO') ' Enter TIME_START and TIME_END: '
         READ (*, *)TIME_START, TIME_END
      END IF
!
      DO L = 1, N_SPX
         REC_POINTER(L) = 4
      END DO
!
      DO L = 1, N_SPX
         READ_SPX(L)    = .FALSE.
         AT_EOF(L)      = .FALSE.
      END DO
      READ_SPX(1)       = .TRUE.
      READ_SPX(4)       = .TRUE.
      IF(MMAX .GT. 1.OR.ANY_SOLVE_ROs) READ_SPX(5) = .TRUE.
      IF(ANY_SOLVE_ROs) READ_SPX(7) = .TRUE.  ! Will compute solids density
!
      IF(TIME_START .LT. TIME_IN_RES) THEN
        CALL SEEK_TIME(READ_SPX, TIME_START, REC_POINTER, TIME_NOW)
        IF(TIME_NOW .LT. ZERO) THEN
          WRITE(*,*)' Could not find record for TIME_START'
          GOTO 10
        ENDIF
      ENDIF
!
100   CONTINUE
      IF(TIME_START .LT. TIME_IN_RES) THEN
        CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,&
                          TIME_NOW, TIME_REAL, NSTEP_1)
      ELSE
        CALL READ_RES1
        TIME_NOW = TIME_IN_RES
      ENDIF
!
      IF (TIME_NOW .LT. ZERO) GOTO 1000
      IF (TIME_NOW .LT. TIME_START) GOTO 100
      IF(MMAX .EQ. 1.AND.(.NOT.ANY_SOLVE_ROs))  THEN
        DO IJK = 1, IJKMAX2
          ROP_s(IJK, 1) = RO_s0(1) * (1. - EP_g(IJK))
        END DO
      ENDIF
!
      DO M = 1, MMAX
        CALL CALC_MU_s(M, IER)
      ENDDO
!
      TIME = TIME_NOW
      NSTEP = NSTEP_1
      IF(CALC_GT) THEN
         DO IJK = 1, IJKMAX2
           IF( FLUID_AT(IJK) ) THEN
             K_1m = 2.D0 * (ONE + C_e) * RO_s(IJK,MM) * G_0(IJK, MM, MM)
             IF(K_1m .NE. 0.0 .AND. EP_s(IJK,MM) .NE. 0.0) THEN
               Tmp(IJK) = P_s(IJK,MM)/K_1m/EP_s(IJK, MM)**2
             ELSE
               Tmp(IJK) = 0.0
             ENDIF
           ELSE
             Tmp(IJK) = 0.0
           ENDIF
         END DO
         READ  (70,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (70,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL OUT_BIN_R(70,Tmp,IJKMAX2,NEXT_REC)
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(70,REC=3) NEXT_REC, NUM_REC
      ENDIF
!
      IF(CALC_GV) THEN
         DO IJK = 1, IJKMAX2
           IF( FLUID_AT(IJK) ) THEN
             Tmp(IJK) = MU_s(IJK, MM)
           ELSE
             Tmp(IJK) = 0.0
           ENDIF
         END DO
         READ  (71,REC=3) NEXT_REC, NUM_REC
         NUM_REC = NEXT_REC
         WRITE (71,REC=NEXT_REC) REAL(TIME) , NSTEP
         NEXT_REC = NEXT_REC + 1
         CALL OUT_BIN_R(71,MU_s(1,MM),IJKMAX2,NEXT_REC)
         NUM_REC = NEXT_REC - NUM_REC
         WRITE(71,REC=3) NEXT_REC, NUM_REC
      ENDIF
!
      IF(TIME_START .LT. TIME_IN_RES) GOTO 100
!
1000  CLOSE (UNIT=70)
      CLOSE (UNIT=71)
      RETURN
      END
