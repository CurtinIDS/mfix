!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EVAL_USR_FCT                                           C
!  Purpose: Evaluates a user-defined function f_usr(x1,x2,x3)          C
!           where (x1,x2,x3) are the Cartesian coordinates of any      C
!           point in the computational domain                          C
! Regions where f_usr < 0 are part of the computational domain.        C
! Regions where f_usr > 0 are excluded from the computational domain.  C
! The integer Q is used to assign a boundary condition to the cut cell C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE EVAL_USR_FCT(x1,x2,x3,Q,f_usr,CLIP_FLAG)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE fldvar
      USE quadric
      USE cutcell

      IMPLICIT NONE


!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
      DOUBLE PRECISION :: x1,x2,x3
      DOUBLE PRECISION :: f_usr
      INTEGER :: Q
      LOGICAL :: CLIP_FLAG
!
!
!  Include files defining statement functions here
!
      IF(N_USR_DEF < 1) RETURN
!
!  Insert user-defined code here
!  The code must define a value for f_usr and Q
!

      CLIP_FLAG = .TRUE.
      RETURN


      END SUBROUTINE EVAL_USR_FCT

