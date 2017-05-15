!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_CORR_TYPE_1                                       C
!  Purpose: Driver routine for correlation calculations : type # 1     C
!           avg(ep_g),sdv(ep_g),                                       C
!           avg(v_g) , sdv(v_g)                                        C
!           Uses N_SPX                                                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-AUG-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMIN1, IMAX1, JMIN1, JMAX1, KMIN1, KMAX1,NSUM C
!  Variables modified: STARTED, PLOT_TYPE, VAR_INDEX, LOC_X, LOC_Y     C
!                      LOC_Z, PLOT_TYPE2, I, J, K, IJK                 C
!                                                                      C
!  Local variables: TIME_1, TIME_2, FILE1_INDEX,FILE2_INDEX            C
!                   NX, NY, NZ, L, NT, FINISH                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CORR_TYPE_1
!
!
      Use param
      Use param1
      Use fldvar
      Use run
      Use geometry
      Use indices
      Use post3d
      Use physprop
      Use correl
      Use compar
      Use functions

      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
      REAL              TIME_REAL(N_SPX) , TIME_1 , TIME_2, TIME_FOUND
      REAL              TIME_NOW
      INTEGER           NSTEP_1
      INTEGER           NX , NY , NZ
      INTEGER           REC_POINTER(N_SPX) , L , NT
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX) , FINISH
      LOGICAL           INTER
      INTEGER           I, J, K, IJK
!
      INTER = .FALSE.
      IF (.NOT.DO_XFORMS) THEN
         WRITE (*,*) ' ENTER TIME_START,TIME_END'
         READ  (*,*) TIME_1 , TIME_2
      ELSE
         TIME_1 = TIME_START
         TIME_2 = TIME_END
      END IF
!
      DO L = 1,N_SPX
         READ_SPX(L)    = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
      READ_SPX(3) = .TRUE.    !  V_g
      READ_SPX(1) = .TRUE.    !  EP_g
      CALL SEEK_TIME(READ_SPX, TIME_1, REC_POINTER, TIME_FOUND)
      IF(TIME_FOUND .LT. ZERO) THEN
        WRITE(*,*)' Could not find record for TIME_START'
        RETURN
      ENDIF
      NT = 0
      FINISH = .FALSE.
      STARTED = 0
!
100   CONTINUE
      CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,&
                         TIME_NOW, TIME_REAL, NSTEP_1)
      IF (DO_XFORMS) THEN
         CALL CHECK_INTER(INTER)
         IF (INTER) RETURN
      END IF
      IF (TIME_NOW .LT. ZERO .OR. TIME_NOW .GT. TIME_2) GOTO 200
      IF (TIME_NOW .LT. TIME_1) GOTO 100
      NT = NT + 1
      CALL CALC_CORR_01 (FINISH,INTER)
      GOTO 100
!
200   IF (NT.EQ.0) THEN
         WRITE (*,*) ' No times found in common'
         RETURN
      END IF
!
      FINISH = .TRUE.
      CALL CALC_CORR_01(FINISH,INTER)
      IF (INTER) RETURN
!
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/,A)')' Average EP_g'
         CALL GET_FILE_NAME(TEMP_FILE4)
      END IF
      OPEN (UNIT=40,FILE=TEMP_FILE4,STATUS='UNKNOWN',CONVERT='BIG_ENDIAN')
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE (40,*)' Average EP_g'
      WRITE (40,'(1X, 2(A, 2X, G12.5))')&
        'Start time = ', TIME_1,'End time = ',TIME_2
      WRITE (40,'(1X, 3(2X, A, I4))')&
       'NZ = ', NZ, 'NY = ', NY, 'NX = ', NX
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               WRITE (40,*) AVG_EP_g(IJK)
            END DO
         END DO
      END DO
      CLOSE (UNIT=40)
!
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/,A)')' Standard Deviation of EP_g'
         CALL GET_FILE_NAME(TEMP_FILE2)
      END IF
      OPEN (UNIT=40,FILE=TEMP_FILE2,STATUS='UNKNOWN',CONVERT='BIG_ENDIAN')
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE(40,*)' Standard Deviation of EP_g'
      WRITE (40,'(1X, 2(A, 2X, G12.5))')&
        'Start time = ', TIME_1,'End time = ',TIME_2
      WRITE (40,'(1X, 3(2X, A, I4))')&
       'NZ = ', NZ, 'NY = ', NY, 'NX = ', NX
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               WRITE (40,*) SDV_EP_g(IJK)
            END DO
         END DO
      END DO
      CLOSE (UNIT=40)
!
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/,A)') ' Average V_g'
         CALL GET_FILE_NAME(TEMP_FILE3)
      END IF
      OPEN (UNIT=40,FILE=TEMP_FILE3,STATUS='UNKNOWN',CONVERT='BIG_ENDIAN')
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE(40,*)' Average V_g'
      WRITE (40,'(1X, 2(A, 2X, G12.5))')&
        'Start time = ', TIME_1,'End time = ',TIME_2
      WRITE (40,'(1X, 3(2X, A, I4))')&
       'NZ = ', NZ, 'NY = ', NY, 'NX = ', NX
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               WRITE (40,*) AVG_V_g(IJK)
            END DO
         END DO
      END DO
      CLOSE (UNIT=40)
!
      IF (.NOT.DO_XFORMS) THEN
         WRITE(*,'(/,A)') ' Standard deviation of V_g'
         CALL GET_FILE_NAME(TEMP_FILE)
      END IF
      OPEN (UNIT=40,FILE=TEMP_FILE,STATUS='UNKNOWN',CONVERT='BIG_ENDIAN')
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE(40,*)' Standard Deviation of V_g'
      WRITE (40,'(1X, 2(A, 2X, G12.5))')&
        'Start time = ', TIME_1,'End time = ',TIME_2
      WRITE (40,'(1X, 3(2X, A, I4))')&
       'NZ = ', NZ, 'NY = ', NY, 'NX = ', NX
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               WRITE (40,*) SDV_V_g(IJK)
            END DO
         END DO
      END DO
      CLOSE (UNIT=40)
!
      RETURN
      END
