!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_SELECTION (SELECTION)                              C
!  Purpose: GET A SELECTION INTEGER FROM THE USER                      C
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
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_SELECTION (SELECTION)
!
      IMPLICIT NONE
!
      INTEGER SELECTION
      character(LEN=2) :: sel
!
      WRITE (*,*) ' '
      WRITE (*,'(A)',ADVANCE='NO') ' Enter menu selection > '
!
!      READ (*,'(I)',ERR=10) SELECTION
      READ (*,'(A2)',ERR=10) sel
      if (sel(2:2) .eq. ' ') then
         sel(2:2) = sel(1:1)
         sel(1:1) = '0'
      end if
      read(sel,'(i2)') selection
      RETURN
10    SELECTION = -1000000
      RETURN
      END
