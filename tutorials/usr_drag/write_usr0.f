!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR0                                             C
!  Purpose: Write initial part of user-defined output                  C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
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
      SUBROUTINE WRITE_USR0

      use discretelement, only: DES_MMAX

      IMPLICIT NONE

      INTEGER :: M
      INTEGER :: lUNIT
      CHARACTER(len=32) :: FNAME

      DO M=1, DES_MMAX
         FNAME=''; WRITE(FNAME,"('POST_PAR_',I1'.dat')")M
         lUNIT = 750 + M
         open(UNIT=lUNIT, FILE=trim(adjustl(FNAME)), STATUS='NEW')
      ENDDO

      RETURN
      END SUBROUTINE WRITE_USR0
