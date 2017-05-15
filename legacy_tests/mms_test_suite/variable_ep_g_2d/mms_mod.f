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
      use compar, only      : myPE, PE_IO
      use compar, only      : ijkstart3, ijkend3
      use functions, only   : funijk_gl
      use indices, only     : i_of, j_of, k_of
      use functions, only   : IS_ON_myPE_owns
      use geometry, only    : imax1, jmax1, kmin1
      use geometry, only    : dx, dy, dz
      use param1, only      : zero, half
      use fldvar, only      : p_g
      IMPLICIT NONE

! indices
      integer               :: ijk, i, j, k, ii, jj, kk

! temporary location variables
      double precision      :: xt, yt, zt


      if(myPE == PE_IO) write(*,"(3x, 'Calculating MMS')")

! allocate mms variables here
!      call allocate_mms_vars

! set reference point for shifting pressure
      ijk_sh = funijk_gl( imax1/2+1, jmax1/2+1, kmin1)  ! for2D

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
      use indices, only     : i_of, j_of, k_of
      use geometry, only    : dx, dy, dz
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
      use physprop, only    : ro_s0
      IMPLICIT NONE

! temporary coordinates
      double precision      :: xt, yt, zt

! index for the variable
! index: 1=Pg, 2=Ug, 3=Vg, 4=Wg 5=Us 6=Vs 7=Ws 8=Tg 9=Ts
!        10=Ep_g 11=ROP_s 12=Theta_m
      integer               :: idx

! local variables: coefficients in manufactured solutions
      double precision      :: pg0=100.0d0, &
                pgx=20.0d0, pgy=-50.0d0, pgz=20.0d0, &
                pgxy=-25.0d0, pgyz=-10.0d0, pgzx=10.0d0, &
                apgx=0.4d0, apgy=0.45d0, apgz=0.85d0, &
                apgxy=0.75d0, apgyz=0.7d0, apgzx=0.8d0
      double precision      :: ug0=7.0d0, &
                ugx=3.0d0, ugy=-4.0d0, ugz=-3.0d0, &
                ugxy=2.0d0, ugyz=1.5d0, ugzx=-2.0d0, &
                augx=0.5d0, augy=0.85d0, augz=0.4d0, &
                augxy=0.6d0, augyz=0.8d0, augzx=0.9d0
      double precision      :: vg0=9.0d0, &
                vgx=-5.0d0, vgy=4.0d0, vgz=5.0d0, &
                vgxy=-3.0d0, vgyz=2.5d0, vgzx=3.5d0, &
                avgx=0.8d0, avgy=0.8d0, avgz=0.5d0, &
                avgxy=0.9d0, avgyz=0.4d0, avgzx=0.6d0
      double precision      :: wg0=8.0d0, &
                wgx=-4.0d0, wgy=3.5d0, wgz=4.2d0, &
                wgxy=-2.2d0, wgyz=2.1d0, wgzx=2.5d0, &
                awgx=0.85d0, awgy=0.9d0, awgz=0.5d0, &
                awgxy=0.4d0, awgyz=0.8d0, awgzx=0.75d0
      double precision      :: us0=5.0d0
      double precision      :: vs0=5.0d0
      double precision      :: ws0=5.0d0
      double precision      :: Tg0=350.0d0, &
                Tgx=10.0d0, Tgy=-30.0d0, Tgz=20.0d0, &
                Tgxy=-12.0d0, Tgyz=10.0d0, Tgzx=8.0d0, &
                aTgx=0.75d0, aTgy=1.25d0, aTgz=0.8d0, &
                aTgxy=0.65d0, aTgyz=0.5d0, aTgzx=0.6d0
      double precision      :: Ts0=300.0d0, &
                Tsx=15.0d0, Tsy=-20.0d0, Tsz=15.0d0, &
                Tsxy=-10.0d0, Tsyz=12.0d0, Tszx=10.0d0, &
                aTsx=0.5d0, aTsy=0.9d0, aTsz=0.8d0, &
                aTsxy=0.5d0, aTsyz=0.65d0, aTszx=0.4d0
      double precision      :: es0=0.3d0, &
                esx=0.06d0, esy=-0.1d0, esz=0.06d0, &
                esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
                aesx=0.4d0, aesy=0.5d0, aesz=0.5d0, &
                aesxy=0.4d0, aesyz=0.4d0, aeszx=0.4d0
!      double precision      :: es0=0.3d0, &
!                esx=0.0d0, esy=0.0d0, esz=0.0d0, &
!                esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
!                aesx=0.4d0, aesy=0.5d0, aesz=0.5d0, &
!                aesxy=0.4d0, aesyz=0.4d0, aeszx=0.4d0
      double precision      :: Ths0=100.0d0, &
                Thsx=5.0d0, Thsy=-10.0d0, Thsz=12.0d0, &
                Thsxy=-8.0d0, Thsyz=10.0d0, Thszx=7.0d0, &
                aThsx=0.8d0, aThsy=1.25d0, aThsz=0.7d0, &
                aThsxy=0.5d0, aThsyz=0.6d0, aThszx=0.7d0
      double precision      :: ros


      ros = ro_s0(1)

      select case(idx)
      case(1)
      !pg!
        mms_function = xt**2 + 10.0
      case(2)
      !ug!
        mms_function = (3.0*yt**2 + 2.0)/(-0.1*xt**2 + 0.7)
      case(3)
      !vg!
        mms_function = (3.0*xt**2 + 2.0)/(-0.1*xt**2 + 0.7)
      case(4)
      !wg!
        mms_function = zero
      case(5)
      !us!
        mms_function = (2.0*yt**2 + 1.0)/(0.1*xt**2 + 0.3)
      case(6)
      !vs!
        mms_function = (2.0*xt**2 + 1.0)/(0.1*xt**2 + 0.3)
      case(7)
      !ws!
        mms_function = zero
      case(8)
      !tg!
        mms_function = 20.0*xt*yt + 300.0
      case(9)
      !ts!
        mms_function = 30.0*xt*yt + 350.0
      case(10)
      !ep_g!
        mms_function = 1.0d0 - ( 0.1*xt**2 + 0.3 )
      case(11)
      !rop_s!
        mms_function = ros*(0.1*xt**2 + 0.3)
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
      use constant, only    : pi, gas_const
      use physprop, only    : ro_s0, mu_g0, ro_g0, MW_avg, &
                                     mu_s0, C_pg0, K_g0, C_ps0, K_s0
      use param1, only      : zero
      IMPLICIT NONE

! temporary coordinates
      double precision      :: xt, yt, zt

! index for the variable
! index: 1=Pg, 2=Ug, 3=Vg, 4=Wg 5=Us 6=Vs 7=Ws 8=Tg 9=Ts
!        10=Ep_g 11=ROP_s 12=Theta_m
      integer               :: idx

! local variables: coefficients in manufactured solutions
      double precision      :: pg0=100.0d0, &
                pgx=20.0d0, pgy=-50.0d0, pgz=20.0d0, &
                pgxy=-25.0d0, pgyz=-10.0d0, pgzx=10.0d0, &
                apgx=0.4d0, apgy=0.45d0, apgz=0.85d0, &
                apgxy=0.75d0, apgyz=0.7d0, apgzx=0.8d0
      double precision      :: ug0=7.0d0, &
                ugx=3.0d0, ugy=-4.0d0, ugz=-3.0d0, &
                ugxy=2.0d0, ugyz=1.5d0, ugzx=-2.0d0, &
                augx=0.5d0, augy=0.85d0, augz=0.4d0, &
                augxy=0.6d0, augyz=0.8d0, augzx=0.9d0
      double precision      :: vg0=9.0d0, &
                vgx=-5.0d0, vgy=4.0d0, vgz=5.0d0, &
                vgxy=-3.0d0, vgyz=2.5d0, vgzx=3.5d0, &
                avgx=0.8d0, avgy=0.8d0, avgz=0.5d0, &
                avgxy=0.9d0, avgyz=0.4d0, avgzx=0.6d0
      double precision      :: wg0=8.0d0, &
                wgx=-4.0d0, wgy=3.5d0, wgz=4.2d0, &
                wgxy=-2.2d0, wgyz=2.1d0, wgzx=2.5d0, &
                awgx=0.85d0, awgy=0.9d0, awgz=0.5d0, &
                awgxy=0.4d0, awgyz=0.8d0, awgzx=0.75d0
      double precision      :: us0=5.0d0
      double precision      :: vs0=5.0d0
      double precision      :: ws0=5.0d0
      double precision      :: Tg0=350.0d0, &
                Tgx=10.0d0, Tgy=-30.0d0, Tgz=20.0d0, &
                Tgxy=-12.0d0, Tgyz=10.0d0, Tgzx=8.0d0, &
                aTgx=0.75d0, aTgy=1.25d0, aTgz=0.8d0, &
                aTgxy=0.65d0, aTgyz=0.5d0, aTgzx=0.6d0
      double precision      :: Ts0=300.0d0, &
                Tsx=15.0d0, Tsy=-20.0d0, Tsz=15.0d0, &
                Tsxy=-10.0d0, Tsyz=12.0d0, Tszx=10.0d0, &
                aTsx=0.5d0, aTsy=0.9d0, aTsz=0.8d0, &
                aTsxy=0.5d0, aTsyz=0.65d0, aTszx=0.4d0
      double precision      :: es0=0.3d0, &
                esx=0.06d0, esy=-0.1d0, esz=0.06d0, &
                esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
                aesx=0.4d0, aesy=0.5d0, aesz=0.5d0, &
                aesxy=0.4d0, aesyz=0.4d0, aeszx=0.4d0
!      double precision      :: es0=0.3d0, &
!                esx=0.0d0, esy=0.0d0, esz=0.0d0, &
!                esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
!                aesx=0.4d0, aesy=0.5d0, aesz=0.5d0, &
!                aesxy=0.4d0, aesyz=0.4d0, aeszx=0.4d0
      double precision      :: Ths0=100.0d0, &
                Thsx=5.0d0, Thsy=-10.0d0, Thsz=12.0d0, &
                Thsxy=-8.0d0, Thsyz=10.0d0, Thszx=7.0d0, &
                aThsx=0.8d0, aThsy=1.25d0, aThsz=0.7d0, &
                aThsxy=0.5d0, aThsyz=0.6d0, aThszx=0.7d0

! local variables within source functions
      double precision      :: ros, mug, rog, MW, Rg, mus, Cpg, kg, &
                                Cps, ks


      mug   = MU_g0
      rog   = RO_g0
      Cpg   = C_pg0
      kg    = K_g0
      MW    = MW_AVG
      Rg    = Gas_Const
      mus   = MU_s0(1)
      ros   = ro_s0(1)
      Cps   = C_ps0(1)
      ks    = K_s0(1)

      select case(idx)
      case(1)
      !pgsrc -> no source terms for pressure correction equation !
        write(*,*) "wrong index in mms_source. stop."
        stop
      case(2)
      !ugsrc!
        mms_source = -0.0266666666666667*mug*xt**2*(3.0*yt**2 + 2.0)/(-0.1*xt**2 + 0.7)**3&
         - mug*(-6.0/(0.1*xt**2 - 0.7) + (3.0*yt**2 + 2.0)*(-0.08*xt**2/(0.1*xt**2 -&
         0.7) + 0.2)/(0.1*xt**2 - 0.7)**2) - 0.0666666666666667*mug*(3.0*yt**2 +&
         2.0)/(-0.1*xt**2 + 0.7)**2 + rog*(0.2*xt*(3.0*yt**2 + 2.0)**2/(-0.1*xt**2 +&
         0.7)**2 + 6.0*yt*(3.0*xt**2 + 2.0)/(-0.1*xt**2 + 0.7)) + 2*xt*(-0.1*xt**2 +&
         0.7)
      case(3)
      !vgsrc!
        mms_source = -0.4*mug*xt*yt/(-0.1*xt**2 + 0.7)**2 - mug*(2.4*xt**2/(0.1*xt**2 -&
         0.7) - 0.08*xt**2*(3.0*xt**2 + 2.0)/(0.1*xt**2 - 0.7)**2 - 6.0 + (0.6*xt**2 +&
         0.4)/(0.1*xt**2 - 0.7))/(0.1*xt**2 - 0.7) + rog*(6.0*xt*(3.0*yt**2 +&
         2.0)/(-0.1*xt**2 + 0.7) + 0.2*xt*(3.0*xt**2 + 2.0)*(3.0*yt**2+&
         2.0)/(-0.1*xt**2 + 0.7)**2)
      case(4)
      !wgsrc!
        mms_source = zero
      case(5)
      !ussrc!
        mms_source = -0.0266666666666667*mus*xt**2*(2.0*yt**2 + 1.0)/(0.1*xt**2 + 0.3)**3 -&
         mus*(4.0/(0.1*xt**2 + 0.3) + (2.0*yt**2 + 1.0)*(0.08*xt**2/(0.1*xt**2 + 0.3) -&
         0.2)/(0.1*xt**2 + 0.3)**2) + 0.0666666666666667*mus*(2.0*yt**2 +&
         1.0)/(0.1*xt**2 + 0.3)**2 + ros*(-0.2*xt*(2.0*yt**2 + 1.0)**2/(0.1*xt**2 +&
         0.3)**2 + 4.0*yt*(2.0*xt**2 + 1.0)/(0.1*xt**2 + 0.3)) + 2*xt*(0.1*xt**2 + 0.3)
      case(6)
      !vssrc!
        mms_source = 0.266666666666667*mus*xt*yt/(0.1*xt**2 + 0.3)**2 -&
         mus*(-1.6*xt**2/(0.1*xt**2 + 0.3) + 0.08*xt**2*(2.0*xt**2 + 1.0)/(0.1*xt**2 +&
         0.3)**2 + (-0.4*xt**2 - 0.2)/(0.1*xt**2 + 0.3) + 4.0)/(0.1*xt**2 + 0.3) +&
         ros*(4.0*xt*(2.0*yt**2 + 1.0)/(0.1*xt**2 + 0.3) - 0.2*xt*(2.0*xt**2 +&
         1.0)*(2.0*yt**2 + 1.0)/(0.1*xt**2 + 0.3)**2)
      case(7)
      !wssrc!
        mms_source = zero
      case(8)
      !tgsrc!
        mms_source = Cpg*rog*(-0.1*xt**2 + 0.7)*(20.0*xt*(3.0*xt**2 + 2.0)/(-0.1*xt**2 +&
         0.7) + 20.0*yt*(3.0*yt**2 + 2.0)/(-0.1*xt**2 + 0.7))
      case(9)
      !tssrc!
        mms_source = Cps*ros*(0.1*xt**2 + 0.3)*(30.0*xt*(2.0*xt**2 + 1.0)/(0.1*xt**2 + 0.3)&
         + 30.0*yt*(2.0*yt**2 + 1.0)/(0.1*xt**2 + 0.3))
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
