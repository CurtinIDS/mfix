      MODULE usr

      Use param
      Use param1


! a dummy variable listed in usrnlst.inc
      DOUBLE PRECISION DUMMY_DP

! Position and velocity for direct integration
      DOUBLE PRECISION :: gX1, gX2
      DOUBLE PRECISION :: gY1, gY2

! Body forces. (Gravity)
      double precision :: F1b
      double precision :: F2b

      contains

!......................................................................!
!                                                                      !
!  Function name: F1kw                                                 !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate particle 1-wall spring force.                    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      double precision function F1kw(y1)

      use discretelement, only: KN_W
      use discretelement, only: PMASS, DES_RADIUS

      implicit none

      double precision, intent(in) :: y1  ! particle height.

      F1kw = -(KN_W/PMASS(1))*(y1 - DES_RADIUS(1))

      return
      end function F1kw


!......................................................................!
!                                                                      !
!  Function name: F1dw                                                 !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate particle 1-wall damping force.                   !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      double precision function F1dw(x1)

      use discretelement, only: DES_ETAN_WALL
      use discretelement, only: PMASS

      implicit none

      double precision, intent(in) :: x1 ! particle velocity

      F1dw = -(DES_ETAN_WALL(1)/PMASS(1)) * x1

      return
      end function F1dw


!......................................................................!
!                                                                      !
!  Function name: F12k                                                 !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate particle-particle spring force.                  !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      double precision function F12k(y1, y2)

      use discretelement, only: KN
      use discretelement, only: PMASS, DES_RADIUS

      implicit none

      double precision, intent(in) :: y1, y2 ! particle height

      F12k = -(KN/PMASS(1))*(2.0d0*DES_RADIUS(1) - (y2 - y1))

      return
      end function F12k



!......................................................................!
!                                                                      !
!  Function name: F12d                                                 !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate particle 1 - particle 2 damping force.           !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      double precision function F12d(x1, x2)

      use discretelement, only: KN
      use discretelement, only: PMASS, DES_ETAN

      implicit none

      double precision, intent(in) :: x1, x2 ! particle height

      F12d = -(DES_ETAN(1,2)/PMASS(1))*(x1 - x2)

      return
      end function F12d


!......................................................................!
!                                                                      !
!  Function name: F2kw                                                 !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate particle 2-wall spring force.                    !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      double precision function F2kw(y2)

      use discretelement, only: KN_W
      use discretelement, only: PMASS, DES_RADIUS
      use geometry, only: YLENGTH

      implicit none

      double precision, intent(in) :: y2  ! particle height.

      F2kw = -(KN_W/PMASS(2))*(DES_RADIUS(2) - (YLENGTH - y2))

      return
      end function F2kw



!......................................................................!
!                                                                      !
!  Function name: F2dw                                                 !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate particle 2-wall damping force.                   !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      double precision function F2dw(x2)

      use discretelement, only: DES_ETAN_WALL
      use discretelement, only: PMASS

      implicit none

      double precision, intent(in) :: x2 ! particle velocity

      F2dw = -(DES_ETAN_WALL(2)/PMASS(2)) * x2

      return
      end function F2dw



!......................................................................!
!                                                                      !
!  Function name: F21k                                                 !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate particle 2 - particle 1 spring force.            !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      double precision function F21k(y1, y2)

      use discretelement, only: PMASS

      implicit none

      double precision, intent(in) :: y1, y2 ! particle height

      F21k = -(PMASS(1)/PMASS(2))* F12k(y1, y2)

      return
      end function F21k



!......................................................................!
!                                                                      !
!  Function name: F21d                                                 !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate particle 1 - particle 2 damping force.           !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      double precision function F21d(x1, x2)

      use discretelement, only: PMASS

      implicit none

      double precision, intent(in) :: x1, x2 ! particle height

      F21d = -(PMASS(1)/PMASS(2)) * F12d(x1,x2)

      return
      end function F21d



!......................................................................!
!                                                                      !
!  Subroutine name: F21d                                               !
!  Author: J.Musser                                   Date:  May-2014  !
!                                                                      !
!  Purpose: Calculate particle 1 - particle 2 damping force.           !
!                                                                      !
!  Ref: R. Garg, J. Galvin, T. Li, and S. Pannala, Documentation of    !
!  open-source MFIX-DEM software for gas-solids flows," from URL:      !
!  https://mfix.netl.doe.gov/documentation/dem_doc_2012-1.pdf,         !
!  page 22, Equation (67).                                             !
!......................................................................!
      SUBROUTINE RK4_V4(DT, y1_a, x1_a, y2_a, x2_a)

      implicit none

      double precision, intent(in) :: DT

      double precision, intent(inout) :: y1_a, x1_a
      double precision, intent(inout) :: y2_a, x2_a

      double precision :: y1, x1
      double precision :: y2, x2

      double precision :: x1_k1, x1_k2, x1_k3, x1_k4, x2_k1, x2_k2, x2_k3, x2_k4
      double precision :: y1_k1, y1_k2, y1_k3, y1_k4, y2_k1, y2_k2, y2_k3, y2_k4

      y1 = y1_a
      y2 = y2_a
      x1 = x1_a
      x2 = x2_a

      y1_K1 = DT * x1
      y2_K1 = DT * x2
      x1_K1 = DT * (F1b+F1kw(y1)+F1dw(x1)+F12k(y1,y2)+F12d(x1,x2))
      x2_K1 = DT * (F2b+F2kw(y2)+F2dw(x2)+F21k(y1,y2)+F21d(x1,x2))


      y1 = y1_a + y1_K1/2.0d0
      y2 = y2_a + y2_K1/2.0d0
      x1 = x1_a + x1_K1/2.0d0
      x2 = x2_a + x2_K1/2.0d0

      y1_K2 = DT * x1
      y2_K2 = DT * x2
      x1_K2 = DT * (F1b+F1kw(y1)+F1dw(x1)+F12k(y1,y2)+F12d(x1,x2))
      x2_K2 = DT * (F2b+F2kw(y2)+F2dw(x2)+F21k(y1,y2)+F21d(x1,x2))


      y1 = y1_a + y1_K2/2.0d0
      y2 = y2_a + y2_K2/2.0d0
      x1 = x1_a + x1_K2/2.0d0
      x2 = x2_a + x2_K2/2.0d0

      y1_K3 = DT * x1
      y2_K3 = DT * x2
      x1_K3 = DT * (F1b+F1kw(y1)+F1dw(x1)+F12k(y1,y2)+F12d(x1,x2))
      x2_K3 = DT * (F2b+F2kw(y2)+F2dw(x2)+F21k(y1,y2)+F21d(x1,x2))


      y1 = y1_a + y1_K3
      y2 = y2_a + y2_K3
      x1 = x1_a + x1_K3
      x2 = x2_a + x2_K3

      y1_K4 = DT * x1
      y2_K4 = DT * x2
      x1_K4 = DT * (F1b+F1kw(y1)+F1dw(x1)+F12k(y1,y2)+F12d(x1,x2))
      x2_K4 = DT * (F2b+F2kw(y2)+F2dw(x2)+F21k(y1,y2)+F21d(x1,x2))

      y1_a = y1_a + (y1_K1 + 2.0d0*y1_K2 + 2.0d0*y1_K3 + y1_K4)/6.0d0
      y2_a = y2_a + (y2_K1 + 2.0d0*y2_K2 + 2.0d0*y2_K3 + y2_K4)/6.0d0
      x1_a = x1_a + (x1_K1 + 2.0d0*x1_K2 + 2.0d0*x1_K3 + x1_K4)/6.0d0
      x2_a = x2_a + (x2_K1 + 2.0d0*x2_K2 + 2.0d0*x2_K3 + x2_K4)/6.0d0

      RETURN
      END SUBROUTINE RK4_V4

      END MODULE usr
