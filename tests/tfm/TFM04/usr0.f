!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR0                                                   C
!  Purpose: This routine is called before the time loop starts and is  C
!           user-definable.  The user may insert code in this routine  C
!           or call appropriate user defined subroutines.  This        C
!           can be used for setting constants and checking errors in   C
!           data.  This routine is not called from an IJK loop, hence  C
!           all indices are undefined.                                 C
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
!
      SUBROUTINE USR0
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
      Use usr, only     : xtr, ytr, ztr
      Use compar, only  : ijkstart3, ijkend3
      Use indices, only : i_of, j_of
      Use geometry, only: dx, dy, dz      
      Use param, only   : dimension_3
      Use param1, only  : half
      use fldvar, only  : p_g, u_g, v_g
      use usr, only     : p_g_ex, u_g_ex, v_g_ex
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!

! looping indices
        integer :: ijk, i, j, k

        integer :: ii, jj, kk

! temporary variables
        double precision  :: xt, yt      

! for function call
        double precision  :: gresho_ic
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!


! allocate variables defined in usr_mod.f
        allocate(xtr(dimension_3))
        allocate(ytr(dimension_3))
        allocate(ztr(dimension_3))
        allocate(p_g_ex(dimension_3))
        allocate(u_g_ex(dimension_3))
        allocate(v_g_ex(dimension_3))

       
        call generate_grid_locations

! Set initial conditions for the Gresho problem
      DO ijk = ijkstart3, ijkend3

        i = i_of(IJK)
        j = j_of(IJK)

! u_g         
        xt = xtr(ijk)
        yt = ytr(ijk) - dy(j)*half
        u_g(ijk) = gresho_ic(xt,yt,'u_g')
        u_g_ex(ijk) = u_g(ijk)

! v_g
        xt = xtr(ijk) - dx(i)*half
        yt = ytr(ijk)
        v_g(ijk) = gresho_ic(xt,yt,'v_g')
        v_g_ex(ijk) = v_g(ijk)

! p_g
        xt = xtr(ijk) - dx(i)*half
        yt = ytr(ijk) - dy(j)*half
        p_g(ijk) = gresho_ic(xt,yt,'p_g')
        p_g_ex(ijk) = p_g(ijk)
         
      END DO

      RETURN
      END SUBROUTINE USR0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: generate_grid_locations                                !
!  Purpose: Generate grid node location arrays from dx, dy dz          !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jun 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE generate_grid_locations
      Use compar, only    : ijkstart3, ijkend3
      Use indices, only   : i_of, j_of, k_of
      Use geometry, only  : dx, dy, dz
      Use usr, only       : xtr, ytr, ztr
      Use param1, only    : zero
      IMPLICIT NONE

! looping indices
        integer :: ijk, i, j, k

        integer :: ii, jj, kk

! temporary variables
        double precision  :: xt, yt, zt


! generate grid locations for exact solution calculation
        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          xt = zero - dx(1)
          yt = zero - dy(1)
          zt = zero - dz(1)

          do ii = 1, i
            xt = xt + dx(ii)
          end do

          do jj = 1, j
            yt = yt + dy(jj)
          end do

          do kk = 1, k
            zt = zt + dz(kk)
          end do

          xtr(ijk) = xt
          ytr(ijk) = yt
          ztr(ijk) = zt
        end do


      END SUBROUTINE generate_grid_locations


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: gresho_ic                                              !
!  Purpose: Function to return gresho problem inital conditions.       !
!                                                                      !
!  Reference:                                                          ! 
!   [1] Liska, R. & Wendroff, B. (2003). Comparison of Several         !
!   Difference Schemes on 1D and 2D Test Problems for the              !
!   Euler Equations.                                                   !
!   SIAM J. Sci. Comput., 25, 995--1017.                               !
!   doi: 10.1137/s1064827502402120                                     !
!                                                                      ! 
!  Author: Aniruddha Choudhary                        Date: Jun 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      DOUBLE PRECISION FUNCTION gresho_ic(xt,yt,ch)
      IMPLICIT NONE

! input x      
      double precision, intent(in)  :: xt
! input y      
      double precision, intent(in)  :: yt
! variable for which return value is desired      
      character(len=*), intent(in)  :: ch

!! local variables
! vortex center coordinates
      double precision  :: x0, y0

! radius at (xt,yt)
      double precision  :: r

! u_phi as a function of r at (xt,y)      
      double precision  :: uphir

! pressure as a function of r at (xt,yt)
      double precision  :: pr

! temporary return variables
      double precision  :: p_g_return
      double precision  :: u_g_return     
      double precision  :: v_g_return


! center of vortex        
        x0 = 0.5d0
        y0 = 0.5d0

! radius at (xt,yt)       
        r = sqrt((xt-x0)**2 + (yt-y0)**2)

! velocity and pressure distribution based upon 
! Reference: Liska and Wendroff (2003)        
        if((r.ge.0.0d0).and.(r.lt.0.2d0)) then
          uphir = 5.0d0*r
          pr = 5.0d0 + 25.0d0/2.0d0*r**2
        elseif((r.ge.0.2d0).and.(r.lt.0.4d0)) then
          uphir = 2.0d0 - 5.0d0*r
          pr = 9.0d0 - 4.0d0*log(0.2d0) + 25.0d0/2.0d0*r**2 &
              - 20.0d0*r + 4.0d0*log(r)
        elseif((r.ge.0.4d0)) then
          uphir = 0.0d0
          pr = 3.0d0 + 4.0d0*log(2.0d0)
        else
          write(*,*) "Check IC setup in usr0.f"
          ERROR STOP
        end if

! set return values        
        if(.not.(r.eq.0.0d0)) then
          u_g_return = -uphir*((yt-y0)/r)
          v_g_return = uphir*((xt-x0)/r)
        else
          u_g_return = 0.0d0
          v_g_return = 0.0d0
        endif

        p_g_return = pr

! return desired value
        select case(trim(ch))
        case('u_g')
          gresho_ic = u_g_return
        case('v_g')
          gresho_ic = v_g_return
        case('p_g')
          gresho_ic = p_g_return
        case default
          write(*,*) "Check gresho_ic in usr0.f"
          ERROR STOP
        end select

        RETURN

      END FUNCTION gresho_ic
