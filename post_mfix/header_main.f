!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: HEADER_MAIN                                            C
!  Purpose: Write out the main selection menu and get user selection   C
!                                                                      C
!  Author: P. Nicoletti                               Date: 05-FEB-92  C
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
!  Variables modified: SELECTION                                       C
!                                                                      C
!  Local variables: NERROR                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE HEADER_MAIN
!
      Use param
      Use param1
      Use post3d
      IMPLICIT NONE
!
!             NUMBER OF INCORRECT INPUT SELECTIONS MADE
      INTEGER NERROR
!
      NERROR = -1
!
10    NERROR = NERROR + 1
      IF (NERROR.GT.10) THEN
         WRITE (*,*) ' HEADER_MAIN : TOO MANY INCORRECT INPUTS'
         STOP
      END IF
      WRITE (*,*)&
        ' *************************************************'
      WRITE (*,*)&
        '  0   - Exit POST_MFIX'
      WRITE (*,*)&
        '  1   - Examine/print data'
      WRITE (*,*)&
        '  2   - Write .RES from data in .SPx files'
      WRITE (*,*)&
        '  3   - Write .RES for a new grid, using old data'
      WRITE (*,*)&
        '  4   - Calculate miscellaneous quantities'
      WRITE (*,*)&
        '  5   - Print out variables'
      WRITE (*,*)&
        '  6   - Call user defined subroutine USR_POST'
      WRITE (*,*)&
        '  7   - Write a new SPx file with selected records'
      WRITE (*,*)&
        '  8   - Write new SPx files with time averaged data'
      WRITE (*,*)&
        '  9   - Perform ORNL calculations'
      WRITE (*,*)&
        ' 10   - run scavenger code'
       WRITE (*,*)&
        ' *************************************************'
!
      CALL GET_SELECTION (SELECTION)
!
      RETURN
      END
