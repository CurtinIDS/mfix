!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ANY_MORE_DATA (READ_SPX,AT_EOF)                        C
!  Purpose: Determine whether all data has been read from the          C
!           requested/needed SPX files ... uses N_SPX                  C
!                                                                      C
!  Author: P. Nicoletti                               Date: 20-MAR-92  C
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
!  Local variables: L                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      LOGICAL FUNCTION ANY_MORE_DATA(READ_SPX,AT_EOF)
!
!
        Use param
        Use param1
      IMPLICIT NONE
!
!     passed arguments
!
      LOGICAL READ_SPX(*) , AT_EOF(*)
!
!     local variables
!
      INTEGER L
!
      ANY_MORE_DATA = .FALSE.
      DO 100 L = 1,N_SPX
         IF (READ_SPX(L) .AND. .NOT.AT_EOF(L)) ANY_MORE_DATA = .TRUE.
100   CONTINUE
!
      RETURN
      END
