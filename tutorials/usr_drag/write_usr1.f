!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR1 (L)                                         C
!  Purpose: Write user-defined output                                  C
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
      SUBROUTINE WRITE_USR1(L)

      use run, only: TIME

      use discretelement, only: DES_MMAX
      use discretelement, only: DES_VEL_NEW
      use discretelement, only: DES_POS_NEW
      use discretelement, only: PIJK, PIP
      use functions, only: IS_NONEXISTENT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L

      INTEGER :: lUNIT
      INTEGER :: NP, M

      IF( L /= 1) RETURN

      DO NP=1, PIP
         IF(is_nonexistent(NP)) CYCLE
         lUNIT = 750 + PIJK(NP,5)
         WRITE(lUNIT,1100) TIME, DES_POS_NEW(NP,2),                    &
            sqrt(DOT_PRODUCT(DES_VEL_NEW(NP,:),DES_VEL_NEW(NP,:)))
      ENDDO

 1100 FORMAT(3x,F6.2,3x,F7.2,3x,F5.2)

      RETURN
      END SUBROUTINE WRITE_USR1
