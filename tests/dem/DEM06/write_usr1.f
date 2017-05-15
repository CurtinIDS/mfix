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

      use usr, only: UPDATE_RK4_SOL
      use run, only: TIME

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L


      SELECT CASE(L)
      CASE(1)
         CALL UPDATE_RK4_SOL(TIME)
         CALL WRITE_DES_OUT(TIME)
      END SELECT

      RETURN
      END SUBROUTINE WRITE_USR1


!......................................................................!
!  Subroutine: WRITE_DES_Out                                           !
!                                                                      !
!  Purpose: Calculate the position and velocity (1D Y-axis) of a free  !
!  falling particle. Compare the results to the MFIX-DEM solultion.    !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!......................................................................!
      SUBROUTINE WRITE_DES_Out(lTime)

      Use discretelement, only: DES_POS_NEW, DES_VEL_NEW
      Use discretelement, only: DES_USR_VAR
      Use usr, only: RK4_POS, RK4_VEL

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: lTime


! Local variables
!---------------------------------------------------------------------//
! file unit for heat transfer data
      INTEGER, PARAMETER :: lUNIT = 2030


! Open the file.
      OPEN(UNIT=lUNIT, FILE='POST_POS.dat',                            &
         POSITION="APPEND", STATUS='OLD')
! Write the results to file.
      WRITE(lUNIT,1000) lTime, RK4_POS(2), DES_POS_NEW(1,2),           &
         ABS_ERR(RK4_POS(2), DES_POS_NEW(1,2))
! Close the output file.
      CLOSE(lUNIT)


! Open the file.
      OPEN(UNIT=lUNIT, FILE='POST_VEL.dat',                            &
         POSITION="APPEND", STATUS='OLD')
! Write the results to file.
      WRITE(lUNIT,1000) lTime, RK4_VEL(2), DES_VEL_NEW(1,2),           &
         ABS_ERR(RK4_VEL(2), DES_VEL_NEW(1,2))
! Close the output file.
      CLOSE(lUNIT)

      RETURN

 1000 FORMAT(3x,F15.8,5X,F15.8,2(3x,F15.8))

      CONTAINS


!......................................................................!
!                                                                      !
!  Subroutine name: ABS_ERR                                            !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate either the absolute percent relative error or    !
!  the absolute error.                                                 !
!                                                                      !
!......................................................................!
      DOUBLE PRECISION FUNCTION ABS_ERR(EXT, NUM)

      use param1, only: SMALL_NUMBER

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: EXT, NUM

      IF(ABS(EXT) > SMALL_NUMBER) THEN
         ABS_ERR = ABS((EXT - NUM)/EXT)*1.0d2
      ELSE
         ABS_ERR = ABS(EXT - NUM)*1.0d2
      ENDIF
      END FUNCTION ABS_ERR

      END SUBROUTINE WRITE_DES_Out
