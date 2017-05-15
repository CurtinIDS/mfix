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

      use discretelement, only: DES_USR_VAR
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
      use compar

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: lTime

! Local variables
!---------------------------------------------------------------------//
! file name
      CHARACTER*64 :: FNAME1, FNAME2
! logical used for testing is the data file already exists
      LOGICAL :: F_EXISTS1, F_EXISTS2
! file unit for heat transfer data
      INTEGER, PARAMETER :: uPos = 2030
      INTEGER, PARAMETER :: uVel = 2031

! Analytic position and velocity
      double precision :: lPos_Y, lVel_Y
! Absolute relative error between MFIX solution and analytic solution.
      double precision :: Pos_aErr, Pos_rErr
      double precision :: Vel_aErr, Vel_rErr

      double precision :: lGrav
      double precision :: lRad
      integer :: lStage

      IF(myPE /= PE_IO) RETURN

! Open the files.
      OPEN(UNIT=uPOS,FILE='POST_POS.dat',POSITION="APPEND",STATUS='OLD')
      OPEN(UNIT=uVEL,FILE='POST_VEL.dat',POSITION="APPEND",STATUS='OLD')

! Set local variables.
      lGrav = -grav(2)
      lRad  = des_radius(1)
      lStage = 0

! Calculate the position and velocity of the particle

! Stage 1: Free fall
      if(lTime < time_c) then
         lStage = 1
         lPos_Y = y_s1(h0, lGrav, lTime)
         lVel_Y = dydt_s1(lGrav, lTime)
! Stage 2: Contact
      elseif( lTime < time_r) then
         lStage = 2
         lPos_Y = y_s2(h0, lRad, b_r, w0_r, lGrav, lTime)
         lVel_Y = dydt_s2(h0, lRad, b_r, w0_r, lGrav, lTime)
! Stage 3: Rebound
      else
         lStage = 3
         lPos_Y = y_s3(lRad, lGrav, lTime)
         lVel_Y = dydt_s3(lGrav, lTime)
      endif

! Calculate the absolute and absolute relative errors.
       Pos_aErr = abs(lPos_Y - DES_POS_new(1,2))
       if(lPos_Y > 0.0d0) then
          Pos_rErr = Pos_aErr/lPos_Y*100
       else
          Pos_rErr = Pos_aErr*100
       endif

       Vel_aErr = abs(lVel_Y - DES_VEL_new(1,2))
       if(lVel_Y > 0.0d0) then
          Vel_rErr = Vel_aErr/lVel_Y*100
       else
          Vel_rErr = Vel_aErr*100
       endif

! Write the results to a file.
      WRITE(uPos,"(F15.8,5x,I1,5X,F15.8,3(3x,F15.8))") lTime, &
         lStage, lPos_Y,DES_POS_new(1,2),Pos_rErr, Pos_aErr
      CLOSE(uPos)


      WRITE(uVel,"(F15.8,5x,I1,5X,F15.8,3(3x,F15.8))")lTime, &
         lStage, lVel_Y, DES_VEL_new(1,2),Vel_rErr, Vel_aErr
      CLOSE(uVel)


      RETURN
      END SUBROUTINE WRITE_DES_Out
