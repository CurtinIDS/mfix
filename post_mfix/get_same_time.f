!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,          C
!                              TIME_NOW, TIME_REAL, NSTEP_1)           C
!                                                                      C
!  Purpose: Synchronize the files enabled by READ_SPX                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 27-OCT-93  C
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
!  Local variables: TDIFF                                              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_SAME_TIME (READ_SPX, REC_POINTER,&
                                AT_EOF, TIME_NOW, TIME_REAL, NSTEP_1)

      USE param
      USE param1
      USE toleranc

      IMPLICIT NONE

      INTEGER   REC_POINTER(*) , NSTEP_1
      INTEGER   L
      LOGICAL   READ_SPX(*) , READ_SPX_STORE(N_SPX), AT_EOF(*)
      REAL      TIME_REAL(*), TIME_NOW
!
      TIME_NOW = -ONE
!
!  Store READ_SPX array
!
      DO 10 L = 1, N_SPX
        READ_SPX_STORE(L) = READ_SPX(L)
10    CONTINUE
!
100   CONTINUE
      CALL READ_SPX1(READ_SPX,REC_POINTER,AT_EOF, TIME_REAL, NSTEP_1)
!
!  Compare times
!
      DO 120 L = 1, N_SPX
        IF(READ_SPX(L)) THEN
          IF(AT_EOF(L)) THEN
            TIME_NOW = -ONE
            RETURN
          ENDIF
          TIME_NOW = MAX(TIME_NOW, TIME_REAL(L))
        ENDIF
120   CONTINUE
!
      DO 140 L = 1, N_SPX
        IF(READ_SPX(L)) THEN
          IF( COMPARE( DBLE(TIME_REAL(L)), DBLE(TIME_NOW))) THEN
            READ_SPX(L) = .FALSE.
          ENDIF
        ENDIF
140   CONTINUE
!
      DO 160 L = 1, N_SPX
        IF(READ_SPX(L)) GOTO 100
160   CONTINUE
!
!  Restore READ_SPX array
!
      DO 500 L = 1, N_SPX
        READ_SPX(L) = READ_SPX_STORE(L)
500   CONTINUE
!
      RETURN
!
      END
