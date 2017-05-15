!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_FROM_RES                                           C
!  Purpose: OUTARR type output from the RES file                       C
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
!  Variables referenced: TIME                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: FILE_NAME                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE OUT_FROM_RES(FILE_NAME)
!
!
      Use param
      Use param1
      Use run
      Use funits
      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
      CHARACTER FILE_NAME*(*)
!
      IF (.NOT.DO_XFORMS) CALL GET_FILE_NAME(FILE_NAME)
      OPEN (UNIT=UNIT_OUT,FILE=FILE_NAME,STATUS='UNKNOWN',convert='big_endian')
      CALL READ_RES0
      CALL READ_RES1
      WRITE (*,*) ' time in RES file = ' , TIME
      CALL WRITE_OUT1
      CLOSE (UNIT=UNIT_OUT)
!
      RETURN
      END
