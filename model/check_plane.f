!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_PLANE                                            C
!  Purpose: make sure the flow boundary condition or internal surface  C
!           is a plane                                                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 10-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_PLANE(X_CONSTANT, Y_CONSTANT, Z_CONSTANT, BC, NAME)

!-----------------------------------------------
! Modules
!-----------------------------------------------
         USE compar
         USE exit, only: mfix_exit
         USE funits

      IMPLICIT NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! surface indicators
      LOGICAL, INTENT(IN) :: X_CONSTANT,Y_CONSTANT,Z_CONSTANT
! boundary condition or internal surface index
      INTEGER, INTENT(IN) ::  BC
! BC or IS
      CHARACTER(LEN=2), INTENT(IN) :: NAME
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! number of directions that are not constant (must equal 2)
      INTEGER :: N
!-----------------------------------------------


! number of directions that are not constant (must equal 2)
      N = 3
      IF (X_CONSTANT) N = N - 1
      IF (Y_CONSTANT) N = N - 1
      IF (Z_CONSTANT) N = N - 1

      IF (N /= 2) THEN
         IF(DMP_LOG)WRITE (UNIT_LOG, 1000) NAME, BC
         call mfix_exit(myPE)
      ENDIF

      RETURN

 1000 FORMAT(/70('*')//' From: CHECK_PLANE',/'Message: ',A,' No ',I3,&
         ' is not a plane',/70('*')/)
      END SUBROUTINE CHECK_PLANE

