!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EVAL_USR_FCT                                           C
!  Purpose: Evaluates a user-defined function f_usr(x1,x2,x3)          C
!           where (x1,x2,x3) are the Cartesian coordinates of any      C
!           point in the computational domain                          C
! Regions where f_usr < 0 are part of the computational domain.        C
! Regions where f_usr > 0 are excluded from the computational domain.  C
! The integer BC_ID is used to assign a boundary condition to the      C
! cut cell                                                             C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE EVAL_USR_FCT(x1,x2,x3,BCID,f_usr,CLIP_FLAG)

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
      DOUBLE PRECISION, PARAMETER :: eps = 1.0D-9
      INTEGER :: BCID
      LOGICAL :: CLIP_FLAG

      DOUBLE PRECISION :: TOL_UDF
      LOGICAL :: test1,test2
      DOUBLE PRECISION :: Amplitude, lambda,y_boundary
!
!
!  Include files defining statement functions here
!
      IF(N_USR_DEF < 1) RETURN
!
!  Insert user-defined code here
!  The code must define a value for f_usr and BC_ID
!
      TOL_UDF = 1.0D-9

      Amplitude = 0.0025
      lambda = 0.05


      IF(5.0*lambda<=x1.and.x1<=15.0*lambda) THEN
         y_boundary = Amplitude * DSIN(2.0*PI*(x1-0.25*lambda)/lambda) +   Amplitude
      ELSE
         y_boundary = ZERO
      ENDIF

      test1 = (x2 > y_boundary + TOL_UDF)
      test2 = (x2 < y_boundary - TOL_UDF)

      IF(test1) then
         f_usr = -ONE
      ELSEIF(test2) then
         f_usr = ONE
      ELSE
         f_usr = ZERO
      ENDIF

      BCID = 12

      CLIP_FLAG = .TRUE.
      RETURN


      END SUBROUTINE EVAL_USR_FCT

