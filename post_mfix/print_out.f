!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PRINT_OUT                                              C
!  Purpose: Create ASCII outputs of all arrays                         C
!                                                                      C
!  Author: P. Nicoletti                               Date: dd-mmm-yy  C
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
      SUBROUTINE PRINT_OUT
!
      Use param
      Use param1
      Use post3d
      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
      REAL              TIME_REAL(N_SPX)
      INTEGER           REC_POINTER(N_SPX)
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX)
!
      IF (DO_XFORMS) THEN
         SELECTION = XCODE
         GOTO 20
      END IF
!
10    WRITE (*,*) &
         '  0   - Return to main menu'
      WRITE (*,*) &
         '  1   - Print out variables from RES file'
      WRITE (*,*) &
         '  2   - Print out variables from SPX files'
!
      CALL GET_SELECTION (SELECTION)
!
 20   CONTINUE
      IF(SELECTION .EQ. 0) THEN
        RETURN
      ELSEIF(SELECTION .EQ. 1) THEN
         CALL OUT_FROM_RES(TEMP_FILE)
      ELSEIF(SELECTION .EQ. 2)THEN
        CALL OUT_FROM_SPX(AT_EOF,READ_SPX,REC_POINTER, TIME_REAL)
      ENDIF
      IF (DO_XFORMS) RETURN
      GOTO 10
!
      END
