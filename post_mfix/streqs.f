!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STREQS(STRING1, STRING2)                               C
!  Purpose: Transfer string2 to string1                                C
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
      SUBROUTINE STREQS(STRING1, STRING2)
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) :: STRING1, STRING2
      INTEGER LEN1, L
!
      LEN1 = MIN(LEN(STRING1), LEN(STRING2))
      DO 10 L = 1, LEN1
        STRING1(L:L) = STRING2(L:L)
10    CONTINUE
      RETURN
      END
