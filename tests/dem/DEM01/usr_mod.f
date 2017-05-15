      MODULE usr

        Use param
        Use param1


! a dummy variable listed in usrnlst.inc
      DOUBLE PRECISION DUMMY_DP


! Initial Starting height of particle
      double precision :: h0

! start time of contact stage (sec)
      double precision :: time_c
! velocity at start of contact stage (cm/sec)
      double precision :: vel_c

! start time of rebound stage
      double precision :: time_r
! velocity at start of rebound stage (cm/sec)
      double precision :: vel_r

! Parameters from particle collision properties.
      double precision :: b_r, w0_r

      contains

!......................................................................!
!  Function name: y_s1                                                 !
!                                                                      !
!  Purpose: Calculate the position (1D in Y-axis direction) during     !
!  stage 1, free fall.                                                 !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 18, Equation (55).                                             !
!......................................................................!
      double precision function y_s1(h0, g, t1)

      double precision, intent(in) :: h0  ! initial particle height
      double precision, intent(in) :: g   ! gravity (down is positive)
      double precision, intent(in) :: t1  ! simulation time

      double precision :: lt              ! local adjusted time

      lt = t1
      y_s1 = h0 - 0.5d0*g*(t1*t1)

      end function y_s1


!......................................................................!
!  Function name: dydt_s1                                              !
!                                                                      !
!  Purpose: Calculate the velocity (1D in Y-axis direction) during     !
!  stage 1, free fall.                                                 !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 18, Equation (54)                                              !
!......................................................................!
      double precision function dydt_s1(g, t1)

      double precision, intent(in) :: g  ! gravity (down is +)
      double precision, intent(in) :: t1 ! simulation time

      double precision :: lt             ! local adjusted time

      lt = t1
      dydt_s1 = -g*lt

      end function dydt_s1


!......................................................................!
!  Function name: y_s2                                                 !
!                                                                      !
!  Purpose: Calculate the position (1D in Y-axis direction) during     !
!  stage 2, contact.                                                   !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf          !
!                                                                      !
!  page 19, Equation (60)                                              !
!......................................................................!
      double precision function y_s2(h0, rp, b, w0, g, t1)

      double precision, intent(in) :: h0  ! initial particle height
      double precision, intent(in) :: rp  ! particle radius
      double precision, intent(in) :: b   ! physical parameter (p 18)
      double precision, intent(in) :: w0  ! physical parameter (p 18)
      double precision, intent(in) :: g   ! gravity (down is +)
      double precision, intent(in) :: t1  ! time

      double precision :: c0, c1, c2      ! coefficients
      double precision :: lt              ! local adjusted time

      lt = t1 - time_c
      c0 = dsqrt(1.0d0 - b*b)*w0
      c1 = g/(w0*w0)
      c2 = ((b*g)/w0 - dsqrt(2.0d0*g*(h0-rp)))/ (w0*dsqrt(1.0d0 - b*b))

      y_s2 = (c1*cos(c0*lt) + c2*sin(c0*lt))*exp(-b*w0*lt)             &
         + rp - g/(w0*w0)

      end function y_s2


!......................................................................!
!  Function name: dydt_s2                                              !
!                                                                      !
!  Purpose: Calculate the velocity (1D in Y-axis direction) during     !
!  stage 2, contact.                                                   !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf          !
!                                                                      !
!  page 19, Equation (61)                                              !
!......................................................................!
      double precision function dydt_s2(h0, rp, b, w0, g, t1)

      double precision, intent(in) :: h0  ! initial particle height
      double precision, intent(in) :: rp  ! particle radius
      double precision, intent(in) :: b   ! physical parameter (p 18)
      double precision, intent(in) :: w0  ! physical parameter (p 18)
      double precision, intent(in) :: g   ! gravity (down is +)
      double precision, intent(in) :: t1  ! time

      double precision :: c0, c1, c2      ! coefficients
      double precision :: lt              ! local adjusted time

      lt = t1 - time_c
      c0 = w0*dsqrt(1.0d0 - b*b)
      c1 = -dsqrt(2.0d0*g*(h0 - rp))
      c2 = (b*w0*dsqrt(2.0d0*g*(h0 - rp)) - g)/(w0*dsqrt(1.0d0 - b*b))

      dydt_s2 = (c1*cos(c0*lt) + c2*sin(c0*lt)) * exp(-b*w0*lt)

      end function dydt_s2


!......................................................................!
!  Function name: y_s3                                                 !
!                                                                      !
!  Purpose: Calculate the position (1D in Y-axis direction) during     !
!  stage 3, rebound.                                                   !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf          !
!                                                                      !
!  page 19, Equation (63)                                              !
!......................................................................!
      double precision function y_s3(rp, g, t1)

      double precision, intent(in) :: rp  ! particle radius
      double precision, intent(in) :: g   ! gravity (down is +)
      double precision, intent(in) :: t1  ! time

      double precision :: lt              ! local adjusted time

      lt = t1 - time_r

      y_s3 = rp + vel_r*lt - 0.5d0*g*(lt*lt)

      end function y_s3


!......................................................................!
!  Function name: dydt_s3                                              !
!                                                                      !
!  Purpose: Calculate the velocity (1D in Y-axis direction) during     !
!  stage 3, rebound.                                                   !
!                                                                      !
!  Author: J.Musser                                   Date:  Jan-13    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf          !
!                                                                      !
!  page 19, Equation (62)                                              !
!......................................................................!
      double precision function dydt_s3(g, t1)

      double precision, intent(in) :: g   ! gravity (down is +)
      double precision, intent(in) :: t1  ! time

      double precision :: lt              ! local adjusted time

      lt = t1 - time_r
      dydt_s3 = vel_r - g*lt


      end function dydt_s3


      END MODULE usr
