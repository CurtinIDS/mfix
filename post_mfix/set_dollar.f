!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_DOLLAR (LINE,LINE_LEN)                             C
!  Purpose: Append a $ to the end of a string (for NCAR graphics)      C
!                                                                      C
!  Author: P. Nicoletti                               Date: 12-APR-92  C
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
!  Local variables: L, LAST_CHAR                                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_DOLLAR(LINE,LINE_LEN)
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) LINE
      INTEGER       LINE_LEN, L, LAST_CHAR
!
      LAST_CHAR = LINE_LEN
      DO L = 1,LINE_LEN-1
         IF (LINE(L:L).NE.' ') LAST_CHAR = L
      END DO
!
      IF (LAST_CHAR.LT.LINE_LEN) THEN
         LINE(LAST_CHAR+1:LAST_CHAR+1) = '$'
      ELSE
         LINE(LAST_CHAR:LAST_CHAR) = '$'
      END IF
!
      RETURN
      END
