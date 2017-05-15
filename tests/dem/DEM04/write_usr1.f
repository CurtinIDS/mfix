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

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: L

      SELECT CASE(L)
      CASE(1); CALL WRITE_DES_OUT(TIME)
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

      Use discretelement
      Use run
      Use usr
      use toleranc, only: COMPARE

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: lTime


! Local variables
!---------------------------------------------------------------------//
! file unit for heat transfer data
      INTEGER, PARAMETER :: lUNIT = 2030
! Slip velocity at contact, error, non-dimensional values.
      DOUBLE PRECISION :: SLIP, ERR, ANL_ND, DEM_ND
! Flag that rolling friction already ended.
      LOGICAL, SAVE :: ROLLFRIC_END = .FALSE.

! Return: Rolling friction already ended.
      IF(ROLLFRIC_END) RETURN

! Calculate the slip velocity.
      SLIP = DES_VEL_NEW(1,1) + OMEGA_NEW(1,3)*DES_RADIUS(1)

! Check for a sign flip or a small difference.
      IF(COMPARE(abs(SLIP),1.0d-6) .OR. SLIP < ZERO) THEN
         ROLLFRIC_END = .TRUE.

! Open the files.
         OPEN(UNIT=lUnit,FILE='POST_TIME.dat', &
            POSITION="APPEND",STATUS='OLD')

! Calculate the non-dimensional end slip times
         ANL_ND = 2.0d0/7.0d0
         DEM_ND = abs(MEW*grav(2)/u0) * lTime

         Err = (abs(ANL_ND-DEM_ND)/abs(ANL_ND) )*100.

! Write the results to a file.
         WRITE(lUNIT,1000) MEW_W, ANL_ND, DEM_ND, Err
         CLOSE(lUNIT)

! Open the files.
         OPEN(UNIT=lUnit,FILE='POST_TVEL.dat', &
            POSITION="APPEND",STATUS='OLD')

! Calculate the non-dimensional translational velocity.
         ANL_ND = 5.0d0/7.0d0
         DEM_ND = abs(DES_VEL_NEW(1,1)/u0)

         Err = (abs(ANL_ND-DEM_ND)/abs(ANL_ND) )*100.

! Write the results to a file.
         WRITE(lUNIT,1000) MEW_W, ANL_ND, DEM_ND, Err
         CLOSE(lUNIT)

! Open the files.
         OPEN(UNIT=lUnit,FILE='POST_AVEL.dat', &
            POSITION="APPEND",STATUS='OLD')

! Calculate the non-dimensional angular velocity.
         ANL_ND = 5.0d0/7.0d0
         DEM_ND = abs(OMEGA_NEW(1,3)*DES_RADIUS(1)/u0)

         Err = (abs(ANL_ND-DEM_ND)/abs(ANL_ND) )*100.

! Write the results to a file.
         WRITE(lUNIT,1000) MEW_W, ANL_ND, DEM_ND, Err
         CLOSE(lUNIT)

      ENDIF

      RETURN

 1000 FORMAT(3x,F15.8,5X,F15.8,2(3x,F15.8))

      END SUBROUTINE WRITE_DES_Out
