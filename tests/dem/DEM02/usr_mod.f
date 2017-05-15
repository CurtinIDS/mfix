      MODULE usr

! a dummy variable listed in usrnlst.inc
      DOUBLE PRECISION DUMMY_DP

! Maximum number of bounces.
      integer, parameter :: MAX_BOUNCE = 50
! The total number of bounces
      integer :: BOUNCE_COUNT
! Record of the maximum bounce heights
      double precision :: MAX_HEIGHT(0:MAX_BOUNCE)

! Initial drop height
      double precision :: h0

!

! The previous position and velocity.
      double precision :: yPOSO
      double precision :: yVELO


      contains


!......................................................................!
!  Function name: y_s1                                                 !
!                                                                      !
!  Purpose: Calculate the maximum heat of the particle center attained !
!  after k collisions based on the hard sphere model.                  !
!                                                                      !
!  Author: J.Musser                                   Date:  Apr-15    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 21, Equation (65).                                             !
!......................................................................!
      double precision function MAX_HEIGHT_HS(k)

      use discretelement, only: DES_RADIUS
      use discretelement, only: DES_EN_WALL_INPUT

      implicit none

! K'th bounce
      integer, intent(in) :: k

      if(k == 0) then
         MAX_HEIGHT_HS = h0
      else
         MAX_HEIGHT_HS = DES_RADIUS(1) + &
            (h0 - DES_RADIUS(1)) * (DES_EN_WALL_INPUT(1)**(2*k))
      endif

      return 
      end function MAX_HEIGHT_HS


      END MODULE usr
