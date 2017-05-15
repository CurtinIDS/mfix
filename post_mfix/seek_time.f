!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SEEK_TIME                                              C
!  Purpose: Subroutine for reading SPX files for a given time          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 03-NOV-93  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SEEK_TIME(READ_SPX, TIME_NEEDED, REC_POINTER, &
                           TIME_FOUND)
      USE param
      USE param1
      USE post3d
      USE funits
      USE toleranc

      IMPLICIT NONE

      INTEGER   L, REC_POINTER(*), NEXT_REC, NSTEP, UPorDOWN
      LOGICAL   READ_SPX(*), NOT_FOUND
      REAL      TIME_NEEDED, TIME_FOUND, TIME_REAL(N_SPX)
!
      TIME_FOUND = TIME_NEEDED
      NOT_FOUND  = .FALSE.
!
      DO 100 L = 1, N_SPX
        IF(READ_SPX(L)) THEN
          IF(.NOT.SPX_OPEN(L)) THEN
            WRITE(*,'(A,I2,A)')' SP',L,' file is not open'
            TIME_FOUND = -ONE
            RETURN
          ENDIF
!
          UPorDOWN = 0
          NEXT_REC = REC_POINTER(L)
          IF(NEXT_REC .LT. 4) NEXT_REC = 4
          IF(NEXT_REC .GT. (LAST_REC(L) - NUM_REC(L)))&
             NEXT_REC = LAST_REC(L) - NUM_REC(L)
10        CONTINUE
          READ (UNIT_SPX+L,REC=NEXT_REC) TIME_REAL(L),NSTEP
          IF(.NOT. COMPARE(DBLE(TIME_REAL(L)), DBLE(TIME_FOUND))) THEN
            IF(UPorDOWN .EQ. 0)THEN
              IF(TIME_REAL(L) .GT. TIME_FOUND) THEN
                UPorDOWN = -1
              ELSE
                UPorDOWN = 1
              ENDIF
            ENDIF
            IF(UPorDOWN .EQ. 1)THEN
              IF(TIME_REAL(L) .GT. TIME_FOUND) THEN
                IF(NOT_FOUND) THEN
                  TIME_FOUND = -ONE
                  RETURN
                ENDIF
                TIME_FOUND = TIME_REAL(L)
                REC_POINTER(L) = NEXT_REC
                NOT_FOUND = .TRUE.
              ELSE
                NEXT_REC = NEXT_REC + NUM_REC(L)
                IF(NEXT_REC .LT. LAST_REC(L))GOTO 10
                TIME_FOUND = -ONE
                RETURN
              ENDIF
            ELSEIF(UPorDOWN .EQ.-1)THEN
              IF(TIME_REAL(L) .LT. TIME_FOUND) THEN
                IF(NOT_FOUND) THEN
                  TIME_FOUND = -ONE
                  RETURN
                ENDIF
                TIME_FOUND = TIME_REAL(L)
                REC_POINTER(L) = NEXT_REC
                NOT_FOUND = .TRUE.
              ELSE
                NEXT_REC = NEXT_REC - NUM_REC(L)
                IF(NEXT_REC .GE. 4)GOTO 10
                TIME_FOUND = -ONE
                RETURN
              ENDIF
            ENDIF
          ELSE
           REC_POINTER(L) = NEXT_REC
          ENDIF
        ENDIF
100   CONTINUE
!
      RETURN
      END
