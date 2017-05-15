!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: MMS_SRC                                                !
!  Purpose: Global storage containers for MMS variables.               !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: 17-Oct-11  !
!  email: anirudd@vt.edu                                               !
!                                                                      !
!  Reviewer: J.Musser                                 Date: 04-Dec-13  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE mms

! By default the MMS functions are unavailable.
      LOGICAL :: USE_MMS = .TRUE.

!! Method of Manufactured Solutions (MMS) and Tecplot variables :

! Gas volume fraction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_Ep_g

! Gas pressure
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_P_g

! Gas velocity components
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_g

! Gas temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_g

! Solid bulk density
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_ROP_s

! Solids velocity components
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_s
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_s

! Solids temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_s

! Granular temperature
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_Theta_m

! Gas continuity MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_ROP_g_Src

! Gas Momentum source terms
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_g_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_g_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_g_Src

! Gas energy equation MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_g_Src

! Solid continuity MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_ROP_s_Src

! Solid momentum MMS source terms
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_U_s_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_V_s_Src
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_W_s_Src

! Solid energy equation MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_T_s_Src

! Granular energy MMS source term
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MMS_Theta_m_Src

! Temporary variable for pressure shifting while plotting and
! discretization error norm calculation
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  P_G_Sh

! Index for pressure shifting
      INTEGER :: IJK_Sh

!  x, y, z locations of top-right corner of a cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  xtr
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ytr
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ztr



      contains



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_MMS_VARS                                      !
!  Purpose: Allocate memory for allocatable variables defined inside   !
!  MMS module.                                                         !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Feb 2015   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_MMS_VARS
      use param, only   : DIMENSION_3

      IMPLICIT NONE

! Note: both fluid and solid phase variables can be of size DIMENSION_3
! since MMS tests do not work for more than one solid phase.

        allocate(MMS_Ep_g(DIMENSION_3))
        allocate(MMS_P_g(DIMENSION_3))
        allocate(MMS_U_g(DIMENSION_3))
        allocate(MMS_V_g(DIMENSION_3))
        allocate(MMS_W_g(DIMENSION_3))
        allocate(MMS_T_g(DIMENSION_3))

        allocate(MMS_ROP_s(DIMENSION_3))
        allocate(MMS_U_s(DIMENSION_3))
        allocate(MMS_V_s(DIMENSION_3))
        allocate(MMS_W_s(DIMENSION_3))
        allocate(MMS_T_s(DIMENSION_3))
        allocate(MMS_Theta_m(DIMENSION_3))

        allocate(MMS_ROP_g_Src(DIMENSION_3))
        allocate(MMS_U_g_Src(DIMENSION_3))
        allocate(MMS_V_g_Src(DIMENSION_3))
        allocate(MMS_W_g_Src(DIMENSION_3))
        allocate(MMS_T_g_Src(DIMENSION_3))

        allocate(MMS_ROP_s_Src(DIMENSION_3))
        allocate(MMS_U_s_Src(DIMENSION_3))
        allocate(MMS_V_s_Src(DIMENSION_3))
        allocate(MMS_W_s_Src(DIMENSION_3))
        allocate(MMS_T_s_Src(DIMENSION_3))
        allocate(MMS_Theta_m_Src(DIMENSION_3))

        allocate(P_g_Sh(DIMENSION_3))

        allocate(xtr(DIMENSION_3))
        allocate(ytr(DIMENSION_3))
        allocate(ztr(DIMENSION_3))


      RETURN
      END SUBROUTINE ALLOCATE_MMS_VARS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DEALLOCATE_MMS_VARS                                    !
!  Purpose: Deallocate memory for allocatable variables defined inside !
!  MMS module.                                                         !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Feb 2015   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DEALLOCATE_MMS_VARS

      IMPLICIT NONE


        deallocate(MMS_Ep_g)
        deallocate(MMS_P_g)
        deallocate(MMS_U_g)
        deallocate(MMS_V_g)
        deallocate(MMS_W_g)
        deallocate(MMS_T_g)

        deallocate(MMS_ROP_s)
        deallocate(MMS_U_s)
        deallocate(MMS_V_s)
        deallocate(MMS_W_s)
        deallocate(MMS_T_s)
        deallocate(MMS_Theta_m)

        deallocate(MMS_ROP_g_Src)
        deallocate(MMS_U_g_Src)
        deallocate(MMS_V_g_Src)
        deallocate(MMS_W_g_Src)
        deallocate(MMS_T_g_Src)

        deallocate(MMS_ROP_s_Src)
        deallocate(MMS_U_s_Src)
        deallocate(MMS_V_s_Src)
        deallocate(MMS_W_s_Src)
        deallocate(MMS_T_s_Src)
        deallocate(MMS_Theta_m_Src)

        deallocate(P_g_Sh)

        deallocate(xtr)
        deallocate(ytr)
        deallocate(ztr)


      RETURN
      END SUBROUTINE DEALLOCATE_MMS_VARS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name:CALCULATE_MMS                                           !
!  Purpose: Generic stub/placeholder for MMS analytical solutions.     !
!                                                                      !
!  Author: J.Musser                                   Date: 04-Dec-13  !
!                                                                      !
! Revision: 1                                                          !
! Purpose: Add case specific code                                      !
! Author: Aniruddha Choudhary                        Date: 04-Feb-15   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALCULATE_MMS
      use compar, only      : ijkstart3, ijkend3
      use compar, only      : myPE, PE_IO
      use fldvar, only      : p_g
      use functions, only   : IS_ON_myPE_owns
      use functions, only   : funijk_gl
      use geometry, only    : dx, dy, dz
      use geometry, only    : imax1, jmax1, kmin1
      use indices, only     : i_of, j_of, k_of
      use param1, only      : zero, half
      IMPLICIT NONE

! indices
      integer               :: ijk, i, j, k, ii, jj, kk

! temporary location variables
      double precision      :: xt, yt, zt


      if(myPE == PE_IO) write(*,"(3x, 'Calculating MMS')")

! allocate mms variables here
      call allocate_mms_vars

! set reference point for shifting pressure
      ijk_sh = funijk_gl( imax1/2+1, jmax1/2+1, kmin1)  ! for 2D

! generate grid locations for plotting and mms calculations
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

! set MMS analytical solutions
      do ijk = ijkstart3, ijkend3

         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)

! scalar variables
         xt = xtr(ijk) - dx(i)*half
         yt = ytr(ijk) - dy(j)*half
         zt = ztr(ijk) - dz(k)*half
         mms_p_g(ijk) = mms_function(xt, yt, zt, 1)
         mms_t_g(ijk) = mms_function(xt, yt, zt, 8)
         mms_t_s(ijk) = mms_function(xt, yt, zt, 9)
         mms_ep_g(ijk) = mms_function(xt, yt, zt, 10)
         mms_rop_s(ijk) = mms_function(xt, yt, zt, 11)
         mms_theta_m(ijk) = mms_function(xt, yt, zt, 12)

! vector variables - (x)
         xt = xtr(ijk)
         yt = ytr(ijk) - dy(j)*half
         zt = ztr(ijk) - dz(k)*half
         mms_u_g(ijk) = mms_function(xt, yt, zt, 2)
         mms_u_s(ijk) = mms_function(xt, yt, zt, 5)

! vector variables - (y)
         xt = xtr(ijk) - dx(i)*half
         yt = ytr(ijk)
         zt = ztr(ijk) - dz(k)*half
         mms_v_g(ijk) = mms_function(xt, yt, zt, 3)
         mms_v_s(ijk) = mms_function(xt, yt, zt, 6)

! vector variables - (z)
         xt = xtr(ijk) - dx(i)*half
         yt = ytr(ijk) - dy(j)*half
         zt = ztr(ijk)
         mms_w_g(ijk) = mms_function(xt, yt, zt, 4)
         mms_w_s(ijk) = mms_function(xt, yt, zt, 7)

      enddo


      RETURN
      END SUBROUTINE CALCULATE_MMS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name:CALCULATE_MMS_SOURCE                                    !
!  Purpose: Generic stub/placeholder for MMS source terms.             !
!                                                                      !
!  Author: J.Musser                                   Date: 04-Dec-13  !
!                                                                      !
! Revision: 1                                                          !
! Purpose: Add case specific code                                      !
! Author: Aniruddha Choudhary                         Date: 04-Feb-15  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALCULATE_MMS_SOURCE
      use compar, only      : ijkstart3, ijkend3
      use geometry, only    : dx, dy, dz
      use indices, only     : i_of, j_of, k_of
      use param1, only      : half
      IMPLICIT NONE

! indices
      integer               :: ijk, i, j, k

! temporary location variables
      double precision      :: xt, yt, zt


      do ijk = ijkstart3, ijkend3

        i = i_of(ijk)
        j = j_of(ijk)
        k = k_of(ijk)

! scalar variables
        xt = xtr(ijk) - dx(i)*half
        yt = ytr(ijk) - dy(j)*half
        zt = ztr(ijk) - dz(k)*half
        mms_t_g_src(ijk) = mms_source(xt, yt, zt, 8)
        mms_t_s_src(ijk) = mms_source(xt, yt, zt, 9)
        mms_rop_g_src(ijk) = mms_source(xt, yt, zt, 10)
        mms_rop_s_src(ijk) = mms_source(xt, yt, zt, 11)
        mms_theta_m_src(ijk) = mms_source(xt, yt, zt, 12)

! vector variable - (x)
        xt = xtr(ijk)
        yt = ytr(ijk) - dy(j)*half
        zt = ztr(ijk) - dz(k)*half
        mms_u_g_src(ijk) = mms_source(xt, yt, zt, 2)
        mms_u_s_src(ijk) = mms_source(xt, yt, zt, 5)


! vector variable - (y)
        xt = xtr(ijk) - dx(i)*half
        yt = ytr(ijk)
        zt = ztr(ijk) - dz(k)*half
        mms_v_g_src(ijk) = mms_source(xt, yt, zt, 3)
        mms_v_s_src(ijk) = mms_source(xt, yt, zt, 6)

! vector variable - (z)
        xt = xtr(ijk) - dx(i)*half
        yt = ytr(ijk) - dy(j)*half
        zt = ztr(ijk)
        mms_w_g_src(ijk) = mms_source(xt, yt, zt, 4)
        mms_w_s_src(ijk) = mms_source(xt, yt, zt, 7)

      enddo


      RETURN
      End SUBROUTINE CALCULATE_MMS_SOURCE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Function: MMS_Function                                               !
! Purpose: Add analytical solutions here                               !
!                                                                      !
! Author: Aniruddha Choudhary                        Date: 17-Oct-11   !
! email: anirudd@vt.edu                                                !
!                                                                      !
! Reviewer: J.Musser                                 Date: 04-Dec-13   !
!                                                                      !
! Revision: 1                                                          !
! Purpose: Cleanup                                                     !
! Author: Aniruddha Choudhary                        Date: 04-Feb-15   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION MMS_Function(xt, yt, zt, idx)
      use constant, only    : pi
      use param1, only      : zero
      IMPLICIT NONE

! temporary coordinates
      double precision      :: xt, yt, zt

! index for the variable
! index: 1=Pg, 2=Ug, 3=Vg, 4=Wg 5=Us 6=Vs 7=Ws 8=Tg 9=Ts
!        10=Ep_g 11=ROP_s 12=Theta_m
      integer               :: idx

! local variables: coefficients in manufactured solutions
      double precision      :: pg0=100.0d0
      double precision      :: ug0=5.0d0


      select case(idx)
      case(1)
      !pg!
        mms_function = pg0*cos(2.0d0*pi*(xt + yt))
      case(2)
      !ug!
        mms_function = ug0*sin(2.0d0*pi*(xt + yt))**2
      case(3)
      !vg!
        mms_function = ug0*cos(2.0d0*pi*(xt + yt))**2
      case(4)
      !wg!
        mms_function = zero
      case(5)
      !us!
        mms_function = zero
      case(6)
      !vs!
        mms_function = zero
      case(7)
      !ws!
        mms_function = zero
      case(8)
      !tg!
        mms_function = zero
      case(9)
      !ts!
        mms_function = zero
      case(10)
      !ep_g!
        mms_function = 1.0d0
      case(11)
      !rop_s!
        mms_function = zero
      case(12)
      !theta_m!
        mms_function = zero
      end select


      END FUNCTION MMS_Function


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Function: MMS_Source                                                 !
! Purpose: Add MMS source terms here                                   !
!                                                                      !
! Author: Aniruddha Choudhary                        Date: 17-Oct-11   !
! email: anirudd@vt.edu                                                !
!                                                                      !
! Reviewer: J.Musser                                 Date: 04-Dec-13   !
!                                                                      !
! Revision: 1                                                          !
! Purpose: Cleanup                                                     !
! Author: Aniruddha Choudhary                        Date: 04-Feb-15   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      DOUBLE PRECISION FUNCTION MMS_Source(xt, yt, zt, idx)
      use constant, only    : pi
      use physprop, only    : MU_g0
      use param1, only      : zero
      IMPLICIT NONE

! temporary coordinates
      double precision      :: xt, yt, zt

! index for the variable
! index: 1=Pg, 2=Ug, 3=Vg, 4=Wg 5=Us 6=Vs 7=Ws 8=Tg 9=Ts
!        10=Ep_g 11=ROP_s 12=Theta_m
      integer               :: idx

! local variables: coefficients in manufactured solutions
      double precision      :: pg0=100.0d0
      double precision      :: ug0=5.0d0

! local variables: physical constants
      double precision      :: mug


      mug = MU_g0

      select case(idx)
      case(1)
      !pgsrc -> no source terms for pressure correction equation !
        write(*,*) "wrong index in mms_source. stop."
        stop
      case(2)
      !ugsrc!
        mms_source = -16*mug*pi**2*ug0*cos(2*pi*(xt + yt))**2 - &
        2*pg0*pi*sin(2*pi*(xt + yt)) + &
        4*pi*ug0**2*cos(2*pi*(xt + yt))**3*sin(2*pi*(xt + yt)) + &
        16*mug*pi**2*ug0*sin(2*pi*(xt + yt))**2 + &
        4*pi*ug0**2*cos(2*pi*(xt + yt))*sin(2*pi*(xt + yt))**3
      case(3)
      !vgsrc!
        mms_source = 16*mug*pi**2*ug0*cos(2*pi*(xt + yt))**2 - &
        2*pg0*pi*sin(2*pi*(xt + yt)) - &
        4*pi*ug0**2*cos(2*pi*(xt + yt))**3*sin(2*pi*(xt + yt)) - &
        16*mug*pi**2*ug0*sin(2*pi*(xt + yt))**2 - &
        4*pi*ug0**2*cos(2*pi*(xt + yt))*sin(2*pi*(xt + yt))**3
      case(4)
      !wgsrc!
        mms_source = zero
      case(5)
      !ussrc!
        mms_source = zero
      case(6)
      !vssrc!
        mms_source = zero
      case(7)
      !wssrc!
        mms_source = zero
      case(8)
      !tgsrc!
        mms_source = zero
      case(9)
      !tssrc!
        mms_source = zero
      case(10)
      !ropgsrc!
        mms_source = zero
      case(11)
      !ropssrc!
        mms_source = zero
      case(12)
      !thssrc!
        mms_source = zero
      end select


      END FUNCTION MMS_Source



      END MODULE mms
