!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_FILE_NAME(FILE_NAME)                               C
!  Purpose: GET A FILENAME FOR OUTPUT DATA                             C
!                                                                      C
!  Author: P. Nicoletti                               Date: 27-FEB-92  C
!  Reviewer:                                                           C
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
      SUBROUTINE GET_FILE_NAME (FILE_NAME)
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*) FILE_NAME
      LOGICAL       FILE_EXIST
      CHARACTER   ANSWER
!
100   WRITE (*,'(A)',ADVANCE='NO') 'Enter Filename for this data > '
      READ  (*,'(A)') FILE_NAME
      INQUIRE (FILE=FILE_NAME,EXIST=FILE_EXIST)
      IF (FILE_EXIST) THEN
         WRITE (*,*) ' '
         WRITE (*,'(A)',ADVANCE='NO')&
             'File already exists.  Overwrite ? (Y/N) > '
         READ (*,'(A)') ANSWER
         IF(ANSWER .EQ. 'Y' .OR. ANSWER .EQ. 'y') RETURN
         GOTO 100
      END IF
!
      RETURN
      END
