!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Subroutine: CFFCTOWALL
!  Purpose: Calculate the total force and torque on a particle
!
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04
!  Reviewer: Rahul Garg                               Date: 02-Aug-07!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE CFFCTOWALL(L, NORM, DIST_LI, FN, FT)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! particle index
      INTEGER, INTENT(IN) :: L
! distance between particle center and wall
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI
! unit normal vector along the line of contact pointing from
! particle L to wall
      DOUBLE PRECISION, INTENT(IN) :: NORM(3)
! normal and tangential force
      DOUBLE PRECISION, INTENT(IN) :: FN(3), FT(3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! local variable for calculating torque on particle
      DOUBLE PRECISION :: CROSSP(3)
! distance from the contact point to the particle center
      DOUBLE PRECISION DIST_CL
!------------------------------------------------

! total contact force
      FC(L,:) = FC(L,:) + FN(:) + FT(:)

! calculate the distance from the particle center to the wall
      DIST_CL = DIST_LI - DES_RADIUS(L)

! total torque
      CROSSP = DES_CROSSPRDCT(NORM, FT)
      TOW(L,:) = TOW(L,:) + DIST_CL*CROSSP(:)

      RETURN
      END SUBROUTINE CFFCTOWALL


