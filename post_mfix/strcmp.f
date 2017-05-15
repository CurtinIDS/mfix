!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STRCMP                                                 C
!  Purpose: Compare two strings                                        C
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
      LOGICAL FUNCTION STRCMP(STRING1, STRING2)
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: STRING1, STRING2
      INTEGER LEN1, LEN2, L
!
      STRCMP = .FALSE.
      LEN1 = LEN(STRING1)
      LEN2 = LEN(STRING2)
      IF(LEN1 .NE. LEN2) RETURN
      DO 10 L = 1, LEN1
        IF(STRING1(L:L) .NE. STRING2(L:L))RETURN
10    CONTINUE
      STRCMP = .TRUE.
      RETURN
      END
