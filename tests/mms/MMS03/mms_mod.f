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
      use geometry, only    : imax1, jmax1, kmax1
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
!      call allocate_mms_vars !! FLAGMMS

! set reference point for shifting pressure
      ijk_sh = funijk_gl( imax1/2+1, jmax1/2+1, kmax1/2+1)  ! for3D

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
        mms_function = pg0 + pgx*cos(apgx*pi*xt) + pgxy*cos(apgxy*pi*xt*yt) +&
         pgy*cos(apgy*pi*yt) + pgyz*sin(apgyz*pi*yt*zt) + pgz*sin(apgz*pi*zt) +&
         pgzx*cos(apgzx*pi*xt*zt)
      case(2)
      !ug!
        mms_function = (-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) + avgz*pi*vgz*sin(avgz*pi*zt) +&
         avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) + awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) +&
         awgy*pi*wgy*cos(awgy*pi*yt) + awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)
      case(3)
      !vg!
        mms_function = (augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) - augz*pi*ugz*sin(augz*pi*zt) -&
         augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) + awgx*pi*wgx*sin(awgx*pi*xt) -&
         awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)
      case(4)
      !wg!
        mms_function = (augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) + augy*pi*ugy*sin(augy*pi*yt) -&
         augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) + avgx*pi*vgx*cos(avgx*pi*xt) -&
         avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)
      case(5)
      !us!
        mms_function = us0*sin(0.5*pi*(xt + yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))
      case(6)
      !vs!
        mms_function = vs0*cos(0.5*pi*(xt + yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))
      case(7)
      !ws!
        mms_function = ws0/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))
      case(8)
      !tg!
        mms_function = tg0 + tgx*cos(atgx*pi*xt) + tgxy*cos(atgxy*pi*xt*yt) +&
         tgy*cos(atgy*pi*yt) + tgyz*sin(atgyz*pi*yt*zt) + tgz*sin(atgz*pi*zt) +&
         tgzx*cos(atgzx*pi*xt*zt)
      case(9)
      !ts!
        mms_function = ts0 + tsx*cos(atsx*pi*xt) + tsxy*cos(atsxy*pi*xt*yt) +&
         tsy*cos(atsy*pi*yt) + tsyz*sin(atsyz*pi*yt*zt) + tsz*sin(atsz*pi*zt) +&
         tszx*cos(atszx*pi*xt*zt)
      case(10)
      !ep_g!
        mms_function = 1.0d0 - (es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + &
         esz*sin(aesz*pi*zt))
      case(11)
      !rop_s!
        mms_function = ros*(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))
      case(12)
      !theta_m!
        mms_function = ths0 + thsx*cos(athsx*pi*xt) + thsxy*cos(athsxy*pi*xt*yt) +&
         thsy*cos(athsy*pi*yt) + thsyz*sin(athsyz*pi*yt*zt) + thsz*sin(athsz*pi*zt) +&
         thszx*cos(athszx*pi*xt*zt)
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
        mms_source = -mug*(pi**2*(-2*aesx**2*esx**2*pi*(-avgyz*vgyz*yt*cos(avgyz*pi*yt*zt)&
         + avgz*vgz*sin(avgz*pi*zt) + avgzx*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*wgy*cos(awgy*pi*yt) +&
         awgyz*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)**2 -&
         aesx**2*esx*pi*(-avgyz*vgyz*yt*cos(avgyz*pi*yt*zt) + avgz*vgz*sin(avgz*pi*zt)&
         + avgzx*vgzx*xt*sin(avgzx*pi*xt*zt) + awgxy*wgxy*xt*cos(awgxy*pi*xt*yt) +&
         awgy*wgy*cos(awgy*pi*yt) +&
         awgyz*wgyz*zt*cos(awgyz*pi*yt*zt))*cos(aesx*pi*xt)/(es0 + esx*cos(aesx*pi*xt)&
         + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) -&
         2*aesx*esx*(avgzx**2*pi*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*vgzx*sin(avgzx*pi*xt*zt) - awgxy**2*pi*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*wgxy*cos(awgxy*pi*xt*yt))*sin(aesx*pi*xt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         avgzx**3*pi*vgzx*xt*zt**2*sin(avgzx*pi*xt*zt) -&
         2*avgzx**2*vgzx*zt*cos(avgzx*pi*xt*zt) +&
         awgxy**3*pi*wgxy*xt*yt**2*cos(awgxy*pi*xt*yt) +&
         2*awgxy**2*wgxy*yt*sin(awgxy*pi*xt*yt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         pi**2*(-2*aesy**2*esy**2*pi*(-avgyz*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*vgz*sin(avgz*pi*zt) + avgzx*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*wgy*cos(awgy*pi*yt) +&
         awgyz*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesy*pi*yt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)**2 -&
         aesy**2*esy*pi*(-avgyz*vgyz*yt*cos(avgyz*pi*yt*zt) + avgz*vgz*sin(avgz*pi*zt)&
         + avgzx*vgzx*xt*sin(avgzx*pi*xt*zt) + awgxy*wgxy*xt*cos(awgxy*pi*xt*yt) +&
         awgy*wgy*cos(awgy*pi*yt) +&
         awgyz*wgyz*zt*cos(awgyz*pi*yt*zt))*cos(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt)&
         + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         2*aesy*esy*(-avgyz**2*pi*vgyz*yt*zt*sin(avgyz*pi*yt*zt) +&
         avgyz*vgyz*cos(avgyz*pi*yt*zt) + awgxy**2*pi*wgxy*xt**2*sin(awgxy*pi*xt*yt) +&
         awgy**2*pi*wgy*sin(awgy*pi*yt) +&
         awgyz**2*pi*wgyz*zt**2*sin(awgyz*pi*yt*zt))*sin(aesy*pi*yt)/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) -&
         avgyz**3*pi*vgyz*yt*zt**2*cos(avgyz*pi*yt*zt) -&
         2*avgyz**2*vgyz*zt*sin(avgyz*pi*yt*zt) +&
         awgxy**3*pi*wgxy*xt**3*cos(awgxy*pi*xt*yt) + awgy**3*pi*wgy*cos(awgy*pi*yt) +&
         awgyz**3*pi*wgyz*zt**3*cos(awgyz*pi*yt*zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         pi**2*(-2*aesz**2*esz**2*pi*(-avgyz*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*vgz*sin(avgz*pi*zt) + avgzx*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*wgy*cos(awgy*pi*yt) +&
         awgyz*wgyz*zt*cos(awgyz*pi*yt*zt))*cos(aesz*pi*zt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)**2 -&
         aesz**2*esz*pi*(-avgyz*vgyz*yt*cos(avgyz*pi*yt*zt) + avgz*vgz*sin(avgz*pi*zt)&
         + avgzx*vgzx*xt*sin(avgzx*pi*xt*zt) + awgxy*wgxy*xt*cos(awgxy*pi*xt*yt) +&
         awgy*wgy*cos(awgy*pi*yt) +&
         awgyz*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt)&
         + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         2*aesz*esz*(avgyz**2*pi*vgyz*yt**2*sin(avgyz*pi*yt*zt) +&
         avgz**2*pi*vgz*cos(avgz*pi*zt) + avgzx**2*pi*vgzx*xt**2*cos(avgzx*pi*xt*zt) -&
         awgyz**2*pi*wgyz*yt*zt*sin(awgyz*pi*yt*zt) +&
         awgyz*wgyz*cos(awgyz*pi*yt*zt))*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) -&
         avgyz**3*pi*vgyz*yt**3*cos(avgyz*pi*yt*zt) + avgz**3*pi*vgz*sin(avgz*pi*zt) +&
         avgzx**3*pi*vgzx*xt**3*sin(avgzx*pi*xt*zt) +&
         awgyz**3*pi*wgyz*yt**2*zt*cos(awgyz*pi*yt*zt) +&
         2*awgyz**2*wgyz*yt*sin(awgyz*pi*yt*zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)) -&
         mug*(2*aesx**2*esx**2*pi**2*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)**2/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesx**2*esx*pi**2*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*cos(aesx*pi*xt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         2*aesx*esx*pi*(avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt) -&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*sin(aesx*pi*xt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (-avgzx**3*pi**3*vgzx*xt*zt**2*sin(avgzx*pi*xt*zt) +&
         2*avgzx**2*pi**2*vgzx*zt*cos(avgzx*pi*xt*zt) -&
         awgxy**3*pi**3*wgxy*xt*yt**2*cos(awgxy*pi*xt*yt) -&
         2*awgxy**2*pi**2*wgxy*yt*sin(awgxy*pi*xt*yt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) -&
         mug*(2*aesx*aesy*esx*esy*pi**2*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesx*pi*xt)*sin(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesx*esx*pi*(-augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) +&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) -&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*sin(aesx*pi*xt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         aesy*esy*pi*(-augzx**2*pi**2*ugzx*xt*zt*cos(augzx*pi*xt*zt) -&
         augzx*pi*ugzx*sin(augzx*pi*xt*zt) + awgx**2*pi**2*wgx*cos(awgx*pi*xt) +&
         awgxy**2*pi**2*wgxy*yt**2*sin(awgxy*pi*xt*yt) +&
         awgzx**2*pi**2*wgzx*zt**2*cos(awgzx*pi*xt*zt))*sin(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (awgxy**3*pi**3*wgxy*xt*yt**2*cos(awgxy*pi*xt*yt) +&
         2*awgxy**2*pi**2*wgxy*yt*sin(awgxy*pi*xt*yt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) -&
         mug*(-2*aesx*aesz*esx*esz*pi**2*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*sin(aesx*pi*xt)*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesx*esx*pi*(augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) -&
         avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) -&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt))*sin(aesx*pi*xt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(augxy**2*pi**2*ugxy*xt*yt*cos(augxy*pi*xt*yt) +&
         augxy*pi*ugxy*sin(augxy*pi*xt*yt) - avgx**2*pi**2*vgx*sin(avgx*pi*xt) -&
         avgxy**2*pi**2*vgxy*yt**2*cos(avgxy*pi*xt*yt) -&
         avgzx**2*pi**2*vgzx*zt**2*cos(avgzx*pi*xt*zt))*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (avgzx**3*pi**3*vgzx*xt*zt**2*sin(avgzx*pi*xt*zt) -&
         2*avgzx**2*pi**2*vgzx*zt*cos(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) +&
         2*mug*(2*aesx**2*esx**2*pi**2*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)**2/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesx**2*esx*pi**2*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*cos(aesx*pi*xt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         2*aesx*aesy*esx*esy*pi**2*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesx*pi*xt)*sin(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         2*aesx*aesz*esx*esz*pi**2*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*sin(aesx*pi*xt)*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesx*esx*pi*(-augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) +&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) -&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*sin(aesx*pi*xt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         aesx*esx*pi*(augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) -&
         avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) -&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt))*sin(aesx*pi*xt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         2*aesx*esx*pi*(avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt) -&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*sin(aesx*pi*xt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         aesy*esy*pi*(-augzx**2*pi**2*ugzx*xt*zt*cos(augzx*pi*xt*zt) -&
         augzx*pi*ugzx*sin(augzx*pi*xt*zt) + awgx**2*pi**2*wgx*cos(awgx*pi*xt) +&
         awgxy**2*pi**2*wgxy*yt**2*sin(awgxy*pi*xt*yt) +&
         awgzx**2*pi**2*wgzx*zt**2*cos(awgzx*pi*xt*zt))*sin(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(augxy**2*pi**2*ugxy*xt*yt*cos(augxy*pi*xt*yt) +&
         augxy*pi*ugxy*sin(augxy*pi*xt*yt) - avgx**2*pi**2*vgx*sin(avgx*pi*xt) -&
         avgxy**2*pi**2*vgxy*yt**2*cos(avgxy*pi*xt*yt) -&
         avgzx**2*pi**2*vgzx*zt**2*cos(avgzx*pi*xt*zt))*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (avgzx**3*pi**3*vgzx*xt*zt**2*sin(avgzx*pi*xt*zt) -&
         2*avgzx**2*pi**2*vgzx*zt*cos(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (awgxy**3*pi**3*wgxy*xt*yt**2*cos(awgxy*pi*xt*yt) +&
         2*awgxy**2*pi**2*wgxy*yt*sin(awgxy*pi*xt*yt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (-avgzx**3*pi**3*vgzx*xt*zt**2*sin(avgzx*pi*xt*zt) +&
         2*avgzx**2*pi**2*vgzx*zt*cos(avgzx*pi*xt*zt) -&
         awgxy**3*pi**3*wgxy*xt*yt**2*cos(awgxy*pi*xt*yt) -&
         2*awgxy**2*pi**2*wgxy*yt*sin(awgxy*pi*xt*yt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0))/3 +&
         rog*(-aesx*esx*pi*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))**2*sin(aesx*pi*xt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         aesy*esy*pi*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (-augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) +&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) -&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) -&
         avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) -&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (2*avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         2*avgzx*pi*vgzx*sin(avgzx*pi*xt*zt) -&
         2*awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         2*awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (avgyz**2*pi**2*vgyz*yt**2*sin(avgyz*pi*yt*zt) +&
         avgz**2*pi**2*vgz*cos(avgz*pi*zt) +&
         avgzx**2*pi**2*vgzx*xt**2*cos(avgzx*pi*xt*zt) -&
         awgyz**2*pi**2*wgyz*yt*zt*sin(awgyz*pi*yt*zt) +&
         awgyz*pi*wgyz*cos(awgyz*pi*yt*zt))*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (avgyz**2*pi**2*vgyz*yt*zt*sin(avgyz*pi*yt*zt) -&
         avgyz*pi*vgyz*cos(avgyz*pi*yt*zt) -&
         awgxy**2*pi**2*wgxy*xt**2*sin(awgxy*pi*xt*yt) -&
         awgy**2*pi**2*wgy*sin(awgy*pi*yt) -&
         awgyz**2*pi**2*wgyz*zt**2*sin(awgyz*pi*yt*zt))*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt)&
         - augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) +&
         (-apgx*pgx*pi*sin(apgx*pi*xt) - apgxy*pgxy*pi*yt*sin(apgxy*pi*xt*yt) -&
         apgzx*pgzx*pi*zt*sin(apgzx*pi*xt*zt))*(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)
      case(3)
      !vgsrc!
        mms_source = -mug*(pi**2*(-2*aesy**2*esy**2*pi*(augyz*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*ugz*sin(augz*pi*zt) - augzx*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*wgx*sin(awgx*pi*xt) - awgxy*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesy*pi*yt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)**2 -&
         aesy**2*esy*pi*(augyz*ugyz*yt*cos(augyz*pi*yt*zt) - augz*ugz*sin(augz*pi*zt) -&
         augzx*ugzx*xt*sin(augzx*pi*xt*zt) + awgx*wgx*sin(awgx*pi*xt) -&
         awgxy*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*wgzx*zt*sin(awgzx*pi*xt*zt))*cos(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt)&
         + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         2*aesy*esy*(augyz**2*pi*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         augyz*ugyz*cos(augyz*pi*yt*zt) - awgxy**2*pi*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*wgxy*cos(awgxy*pi*xt*yt))*sin(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         augyz**3*pi*ugyz*yt*zt**2*cos(augyz*pi*yt*zt) +&
         2*augyz**2*ugyz*zt*sin(augyz*pi*yt*zt) -&
         awgxy**3*pi*wgxy*xt**2*yt*cos(awgxy*pi*xt*yt) -&
         2*awgxy**2*wgxy*xt*sin(awgxy*pi*xt*yt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         pi**2*(-2*aesx**2*esx**2*pi*(augyz*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*ugz*sin(augz*pi*zt) - augzx*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*wgx*sin(awgx*pi*xt) - awgxy*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesx*pi*xt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)**2 -&
         aesx**2*esx*pi*(augyz*ugyz*yt*cos(augyz*pi*yt*zt) - augz*ugz*sin(augz*pi*zt) -&
         augzx*ugzx*xt*sin(augzx*pi*xt*zt) + awgx*wgx*sin(awgx*pi*xt) -&
         awgxy*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*wgzx*zt*sin(awgzx*pi*xt*zt))*cos(aesx*pi*xt)/(es0 + esx*cos(aesx*pi*xt)&
         + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) -&
         2*aesx*esx*(-augzx**2*pi*ugzx*xt*zt*cos(augzx*pi*xt*zt) -&
         augzx*ugzx*sin(augzx*pi*xt*zt) + awgx**2*pi*wgx*cos(awgx*pi*xt) +&
         awgxy**2*pi*wgxy*yt**2*sin(awgxy*pi*xt*yt) +&
         awgzx**2*pi*wgzx*zt**2*cos(awgzx*pi*xt*zt))*sin(aesx*pi*xt)/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) -&
         augzx**3*pi*ugzx*xt*zt**2*sin(augzx*pi*xt*zt) +&
         2*augzx**2*ugzx*zt*cos(augzx*pi*xt*zt) + awgx**3*pi*wgx*sin(awgx*pi*xt) -&
         awgxy**3*pi*wgxy*yt**3*cos(awgxy*pi*xt*yt) +&
         awgzx**3*pi*wgzx*zt**3*sin(awgzx*pi*xt*zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) -&
         pi**2*(2*aesz**2*esz**2*pi*(augyz*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*ugz*sin(augz*pi*zt) - augzx*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*wgx*sin(awgx*pi*xt) - awgxy*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*wgzx*zt*sin(awgzx*pi*xt*zt))*cos(aesz*pi*zt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)**2 +&
         aesz**2*esz*pi*(augyz*ugyz*yt*cos(augyz*pi*yt*zt) - augz*ugz*sin(augz*pi*zt) -&
         augzx*ugzx*xt*sin(augzx*pi*xt*zt) + awgx*wgx*sin(awgx*pi*xt) -&
         awgxy*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt)&
         + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         2*aesz*esz*(augyz**2*pi*ugyz*yt**2*sin(augyz*pi*yt*zt) +&
         augz**2*pi*ugz*cos(augz*pi*zt) + augzx**2*pi*ugzx*xt**2*cos(augzx*pi*xt*zt) -&
         awgzx**2*pi*wgzx*xt*zt*cos(awgzx*pi*xt*zt) -&
         awgzx*wgzx*sin(awgzx*pi*xt*zt))*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) -&
         augyz**3*pi*ugyz*yt**3*cos(augyz*pi*yt*zt) + augz**3*pi*ugz*sin(augz*pi*zt) +&
         augzx**3*pi*ugzx*xt**3*sin(augzx*pi*xt*zt) -&
         awgzx**3*pi*wgzx*xt**2*zt*sin(awgzx*pi*xt*zt) +&
         2*awgzx**2*wgzx*xt*cos(awgzx*pi*xt*zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)) -&
         mug*(2*aesy**2*esy**2*pi**2*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesy*pi*yt)**2/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesy**2*esy*pi**2*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*cos(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         2*aesy*esy*pi*(-augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) +&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) -&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*sin(aesy*pi*yt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (-augyz**3*pi**3*ugyz*yt*zt**2*cos(augyz*pi*yt*zt) -&
         2*augyz**2*pi**2*ugyz*zt*sin(augyz*pi*yt*zt) +&
         awgxy**3*pi**3*wgxy*xt**2*yt*cos(awgxy*pi*xt*yt) +&
         2*awgxy**2*pi**2*wgxy*xt*sin(awgxy*pi*xt*yt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) -&
         mug*(2*aesx*aesy*esx*esy*pi**2*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)*sin(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesx*esx*pi*(avgyz**2*pi**2*vgyz*yt*zt*sin(avgyz*pi*yt*zt) -&
         avgyz*pi*vgyz*cos(avgyz*pi*yt*zt) -&
         awgxy**2*pi**2*wgxy*xt**2*sin(awgxy*pi*xt*yt) -&
         awgy**2*pi**2*wgy*sin(awgy*pi*yt) -&
         awgyz**2*pi**2*wgyz*zt**2*sin(awgyz*pi*yt*zt))*sin(aesx*pi*xt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         aesy*esy*pi*(avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt) -&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*sin(aesy*pi*yt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (-awgxy**3*pi**3*wgxy*xt**2*yt*cos(awgxy*pi*xt*yt) -&
         2*awgxy**2*pi**2*wgxy*xt*sin(awgxy*pi*xt*yt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) -&
         mug*(-2*aesy*aesz*esy*esz*pi**2*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*sin(aesy*pi*yt)*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesy*esy*pi*(augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) -&
         avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) -&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt))*sin(aesy*pi*yt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(augxy**2*pi**2*ugxy*xt**2*cos(augxy*pi*xt*yt) +&
         augy**2*pi**2*ugy*cos(augy*pi*yt) +&
         augyz**2*pi**2*ugyz*zt**2*sin(augyz*pi*yt*zt) -&
         avgxy**2*pi**2*vgxy*xt*yt*cos(avgxy*pi*xt*yt) -&
         avgxy*pi*vgxy*sin(avgxy*pi*xt*yt))*cos(aesz*pi*zt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (augyz**3*pi**3*ugyz*yt*zt**2*cos(augyz*pi*yt*zt) +&
         2*augyz**2*pi**2*ugyz*zt*sin(augyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) +&
         2*mug*(2*aesx*aesy*esx*esy*pi**2*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)*sin(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesx*esx*pi*(avgyz**2*pi**2*vgyz*yt*zt*sin(avgyz*pi*yt*zt) -&
         avgyz*pi*vgyz*cos(avgyz*pi*yt*zt) -&
         awgxy**2*pi**2*wgxy*xt**2*sin(awgxy*pi*xt*yt) -&
         awgy**2*pi**2*wgy*sin(awgy*pi*yt) -&
         awgyz**2*pi**2*wgyz*zt**2*sin(awgyz*pi*yt*zt))*sin(aesx*pi*xt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         2*aesy**2*esy**2*pi**2*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesy*pi*yt)**2/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesy**2*esy*pi**2*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*cos(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         2*aesy*aesz*esy*esz*pi**2*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*sin(aesy*pi*yt)*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         2*aesy*esy*pi*(-augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) +&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) -&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*sin(aesy*pi*yt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         aesy*esy*pi*(augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) -&
         avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) -&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt))*sin(aesy*pi*yt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         aesy*esy*pi*(avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt) -&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*sin(aesy*pi*yt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(augxy**2*pi**2*ugxy*xt**2*cos(augxy*pi*xt*yt) +&
         augy**2*pi**2*ugy*cos(augy*pi*yt) +&
         augyz**2*pi**2*ugyz*zt**2*sin(augyz*pi*yt*zt) -&
         avgxy**2*pi**2*vgxy*xt*yt*cos(avgxy*pi*xt*yt) -&
         avgxy*pi*vgxy*sin(avgxy*pi*xt*yt))*cos(aesz*pi*zt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (augyz**3*pi**3*ugyz*yt*zt**2*cos(augyz*pi*yt*zt) +&
         2*augyz**2*pi**2*ugyz*zt*sin(augyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (-awgxy**3*pi**3*wgxy*xt**2*yt*cos(awgxy*pi*xt*yt) -&
         2*awgxy**2*pi**2*wgxy*xt*sin(awgxy*pi*xt*yt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (-augyz**3*pi**3*ugyz*yt*zt**2*cos(augyz*pi*yt*zt) -&
         2*augyz**2*pi**2*ugyz*zt*sin(augyz*pi*yt*zt) +&
         awgxy**3*pi**3*wgxy*xt**2*yt*cos(awgxy*pi*xt*yt) +&
         2*awgxy**2*pi**2*wgxy*xt*sin(awgxy*pi*xt*yt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0))/3 +&
         rog*(-aesx*esx*pi*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         aesy*esy*pi*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))**2*sin(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (-2*augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         2*augyz*pi*ugyz*cos(augyz*pi*yt*zt) +&
         2*awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) -&
         2*awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) -&
         avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) -&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt))*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt) -&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (-augyz**2*pi**2*ugyz*yt**2*sin(augyz*pi*yt*zt) -&
         augz**2*pi**2*ugz*cos(augz*pi*zt) -&
         augzx**2*pi**2*ugzx*xt**2*cos(augzx*pi*xt*zt) +&
         awgzx**2*pi**2*wgzx*xt*zt*cos(awgzx*pi*xt*zt) +&
         awgzx*pi*wgzx*sin(awgzx*pi*xt*zt))*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (-augzx**2*pi**2*ugzx*xt*zt*cos(augzx*pi*xt*zt) -&
         augzx*pi*ugzx*sin(augzx*pi*xt*zt) + awgx**2*pi**2*wgx*cos(awgx*pi*xt) +&
         awgxy**2*pi**2*wgxy*yt**2*sin(awgxy*pi*xt*yt) +&
         awgzx**2*pi**2*wgzx*zt**2*cos(awgzx*pi*xt*zt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt)&
         + avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) +&
         (-apgxy*pgxy*pi*xt*sin(apgxy*pi*xt*yt) - apgy*pgy*pi*sin(apgy*pi*yt) +&
         apgyz*pgyz*pi*zt*cos(apgyz*pi*yt*zt))*(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)
      case(4)
      !wgsrc!
        mms_source = -mug*(-pi**2*(2*aesz**2*esz**2*pi*(augxy*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*ugy*sin(augy*pi*yt) - augyz*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*vgx*cos(avgx*pi*xt) - avgxy*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*vgzx*zt*sin(avgzx*pi*xt*zt))*cos(aesz*pi*zt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)**2 +&
         aesz**2*esz*pi*(augxy*ugxy*xt*sin(augxy*pi*xt*yt) + augy*ugy*sin(augy*pi*yt) -&
         augyz*ugyz*zt*cos(augyz*pi*yt*zt) + avgx*vgx*cos(avgx*pi*xt) -&
         avgxy*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*vgzx*zt*sin(avgzx*pi*xt*zt))*sin(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt)&
         + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         2*aesz*esz*(-augyz**2*pi*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         augyz*ugyz*cos(augyz*pi*yt*zt) + avgzx**2*pi*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*vgzx*sin(avgzx*pi*xt*zt))*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         augyz**3*pi*ugyz*yt**2*zt*cos(augyz*pi*yt*zt) +&
         2*augyz**2*ugyz*yt*sin(augyz*pi*yt*zt) +&
         avgzx**3*pi*vgzx*xt**2*zt*sin(avgzx*pi*xt*zt) -&
         2*avgzx**2*vgzx*xt*cos(avgzx*pi*xt*zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         pi**2*(-2*aesx**2*esx**2*pi*(augxy*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*ugy*sin(augy*pi*yt) - augyz*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*vgx*cos(avgx*pi*xt) - avgxy*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*vgzx*zt*sin(avgzx*pi*xt*zt))*sin(aesx*pi*xt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)**2 -&
         aesx**2*esx*pi*(augxy*ugxy*xt*sin(augxy*pi*xt*yt) + augy*ugy*sin(augy*pi*yt) -&
         augyz*ugyz*zt*cos(augyz*pi*yt*zt) + avgx*vgx*cos(avgx*pi*xt) -&
         avgxy*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*vgzx*zt*sin(avgzx*pi*xt*zt))*cos(aesx*pi*xt)/(es0 + esx*cos(aesx*pi*xt)&
         + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         2*aesx*esx*(-augxy**2*pi*ugxy*xt*yt*cos(augxy*pi*xt*yt) -&
         augxy*ugxy*sin(augxy*pi*xt*yt) + avgx**2*pi*vgx*sin(avgx*pi*xt) +&
         avgxy**2*pi*vgxy*yt**2*cos(avgxy*pi*xt*yt) +&
         avgzx**2*pi*vgzx*zt**2*cos(avgzx*pi*xt*zt))*sin(aesx*pi*xt)/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         augxy**3*pi*ugxy*xt*yt**2*sin(augxy*pi*xt*yt) -&
         2*augxy**2*ugxy*yt*cos(augxy*pi*xt*yt) + avgx**3*pi*vgx*cos(avgx*pi*xt) -&
         avgxy**3*pi*vgxy*yt**3*sin(avgxy*pi*xt*yt) -&
         avgzx**3*pi*vgzx*zt**3*sin(avgzx*pi*xt*zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         pi**2*(-2*aesy**2*esy**2*pi*(augxy*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*ugy*sin(augy*pi*yt) - augyz*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*vgx*cos(avgx*pi*xt) - avgxy*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*vgzx*zt*sin(avgzx*pi*xt*zt))*sin(aesy*pi*yt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)**2 -&
         aesy**2*esy*pi*(augxy*ugxy*xt*sin(augxy*pi*xt*yt) + augy*ugy*sin(augy*pi*yt) -&
         augyz*ugyz*zt*cos(augyz*pi*yt*zt) + avgx*vgx*cos(avgx*pi*xt) -&
         avgxy*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*vgzx*zt*sin(avgzx*pi*xt*zt))*cos(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt)&
         + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) -&
         2*aesy*esy*(augxy**2*pi*ugxy*xt**2*cos(augxy*pi*xt*yt) +&
         augy**2*pi*ugy*cos(augy*pi*yt) + augyz**2*pi*ugyz*zt**2*sin(augyz*pi*yt*zt) -&
         avgxy**2*pi*vgxy*xt*yt*cos(avgxy*pi*xt*yt) -&
         avgxy*vgxy*sin(avgxy*pi*xt*yt))*sin(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0) +&
         augxy**3*pi*ugxy*xt**3*sin(augxy*pi*xt*yt) + augy**3*pi*ugy*sin(augy*pi*yt) -&
         augyz**3*pi*ugyz*zt**3*cos(augyz*pi*yt*zt) -&
         avgxy**3*pi*vgxy*xt**2*yt*sin(avgxy*pi*xt*yt) +&
         2*avgxy**2*vgxy*xt*cos(avgxy*pi*xt*yt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) - 1.0)) -&
         mug*(2*aesz**2*esz**2*pi**2*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*cos(aesz*pi*zt)**2/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesz**2*esz*pi**2*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*sin(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         2*aesz*esz*pi*(augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) -&
         avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) -&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt))*cos(aesz*pi*zt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (augyz**3*pi**3*ugyz*yt**2*zt*cos(augyz*pi*yt*zt) +&
         2*augyz**2*pi**2*ugyz*yt*sin(augyz*pi*yt*zt) +&
         avgzx**3*pi**3*vgzx*xt**2*zt*sin(avgzx*pi*xt*zt) -&
         2*avgzx**2*pi**2*vgzx*xt*cos(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) -&
         mug*(-2*aesx*aesz*esx*esz*pi**2*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesx*esx*pi*(avgyz**2*pi**2*vgyz*yt**2*sin(avgyz*pi*yt*zt) +&
         avgz**2*pi**2*vgz*cos(avgz*pi*zt) +&
         avgzx**2*pi**2*vgzx*xt**2*cos(avgzx*pi*xt*zt) -&
         awgyz**2*pi**2*wgyz*yt*zt*sin(awgyz*pi*yt*zt) +&
         awgyz*pi*wgyz*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt) -&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*cos(aesz*pi*zt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (-avgzx**3*pi**3*vgzx*xt**2*zt*sin(avgzx*pi*xt*zt) +&
         2*avgzx**2*pi**2*vgzx*xt*cos(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) -&
         mug*(-2*aesy*aesz*esy*esz*pi**2*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesy*pi*yt)*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesy*esy*pi*(-augyz**2*pi**2*ugyz*yt**2*sin(augyz*pi*yt*zt) -&
         augz**2*pi**2*ugz*cos(augz*pi*zt) -&
         augzx**2*pi**2*ugzx*xt**2*cos(augzx*pi*xt*zt) +&
         awgzx**2*pi**2*wgzx*xt*zt*cos(awgzx*pi*xt*zt) +&
         awgzx*pi*wgzx*sin(awgzx*pi*xt*zt))*sin(aesy*pi*yt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(-augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) +&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) -&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*cos(aesz*pi*zt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (-augyz**3*pi**3*ugyz*yt**2*zt*cos(augyz*pi*yt*zt) -&
         2*augyz**2*pi**2*ugyz*yt*sin(augyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) +&
         2*mug*(-2*aesx*aesz*esx*esz*pi**2*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesx*esx*pi*(avgyz**2*pi**2*vgyz*yt**2*sin(avgyz*pi*yt*zt) +&
         avgz**2*pi**2*vgz*cos(avgz*pi*zt) +&
         avgzx**2*pi**2*vgzx*xt**2*cos(avgzx*pi*xt*zt) -&
         awgyz**2*pi**2*wgyz*yt*zt*sin(awgyz*pi*yt*zt) +&
         awgyz*pi*wgyz*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         2*aesy*aesz*esy*esz*pi**2*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesy*pi*yt)*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesy*esy*pi*(-augyz**2*pi**2*ugyz*yt**2*sin(augyz*pi*yt*zt) -&
         augz**2*pi**2*ugz*cos(augz*pi*zt) -&
         augzx**2*pi**2*ugzx*xt**2*cos(augzx*pi*xt*zt) +&
         awgzx**2*pi**2*wgzx*xt*zt*cos(awgzx*pi*xt*zt) +&
         awgzx*pi*wgzx*sin(awgzx*pi*xt*zt))*sin(aesy*pi*yt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         2*aesz**2*esz**2*pi**2*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*cos(aesz*pi*zt)**2/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**3 -&
         aesz**2*esz*pi**2*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*sin(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(-augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) +&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) -&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*cos(aesz*pi*zt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         2*aesz*esz*pi*(augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) -&
         avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) -&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt))*cos(aesz*pi*zt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt) -&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*cos(aesz*pi*zt)/(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (-augyz**3*pi**3*ugyz*yt**2*zt*cos(augyz*pi*yt*zt) -&
         2*augyz**2*pi**2*ugyz*yt*sin(augyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (-avgzx**3*pi**3*vgzx*xt**2*zt*sin(avgzx*pi*xt*zt) +&
         2*avgzx**2*pi**2*vgzx*xt*cos(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (augyz**3*pi**3*ugyz*yt**2*zt*cos(augyz*pi*yt*zt) +&
         2*augyz**2*pi**2*ugyz*yt*sin(augyz*pi*yt*zt) +&
         avgzx**3*pi**3*vgzx*xt**2*zt*sin(avgzx*pi*xt*zt) -&
         2*avgzx**2*pi**2*vgzx*xt*cos(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0))/3 +&
         rog*(-aesx*esx*pi*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))*sin(aesx*pi*xt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 -&
         aesy*esy*pi*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))*sin(aesy*pi*yt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         aesz*esz*pi*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))**2*cos(aesz*pi*zt)/(-es0 -&
         esx*cos(aesx*pi*xt) - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)**2 +&
         (-augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) +&
         augyz*pi*ugyz*cos(augyz*pi*yt*zt) +&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) -&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (2*augyz**2*pi**2*ugyz*yt*zt*sin(augyz*pi*yt*zt) -&
         2*augyz*pi*ugyz*cos(augyz*pi*yt*zt) -&
         2*avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) -&
         2*avgzx*pi*vgzx*sin(avgzx*pi*xt*zt))*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (avgzx**2*pi**2*vgzx*xt*zt*cos(avgzx*pi*xt*zt) +&
         avgzx*pi*vgzx*sin(avgzx*pi*xt*zt) -&
         awgxy**2*pi**2*wgxy*xt*yt*sin(awgxy*pi*xt*yt) +&
         awgxy*pi*wgxy*cos(awgxy*pi*xt*yt))*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (augxy**2*pi**2*ugxy*xt**2*cos(augxy*pi*xt*yt) +&
         augy**2*pi**2*ugy*cos(augy*pi*yt) +&
         augyz**2*pi**2*ugyz*zt**2*sin(augyz*pi*yt*zt) -&
         avgxy**2*pi**2*vgxy*xt*yt*cos(avgxy*pi*xt*yt) -&
         avgxy*pi*vgxy*sin(avgxy*pi*xt*yt))*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (augxy**2*pi**2*ugxy*xt*yt*cos(augxy*pi*xt*yt) +&
         augxy*pi*ugxy*sin(augxy*pi*xt*yt) - avgx**2*pi**2*vgx*sin(avgx*pi*xt) -&
         avgxy**2*pi**2*vgxy*yt**2*cos(avgxy*pi*xt*yt) -&
         avgzx**2*pi**2*vgzx*zt**2*cos(avgzx*pi*xt*zt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt)&
         + avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)) +&
         (apgyz*pgyz*pi*yt*cos(apgyz*pi*yt*zt) + apgz*pgz*pi*cos(apgz*pi*zt) -&
         apgzx*pgzx*pi*xt*sin(apgzx*pi*xt*zt))*(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0)
      case(5)
      !ussrc!
        mms_source = 2*aesx*aesz*esx*esz*mus*pi**2*ws0*sin(aesx*pi*xt)*cos(aesz*pi*zt)/(es0&
         + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 -&
         mus*(pi**2*us0*(2*aesx**2*esx**2*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesx*pi*xt)**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + aesx**2*esx*sin(0.5*pi*(xt + yt +&
         zt))**2*cos(aesx*pi*xt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) + 2.0*aesx*esx*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) - 0.5*sin(0.5*pi*(xt + yt + zt))**2&
         + 0.5*cos(0.5*pi*(xt + yt + zt))**2)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         pi**2*us0*(2*aesy**2*esy**2*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesy*pi*yt)**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + aesy**2*esy*sin(0.5*pi*(xt + yt +&
         zt))**2*cos(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) + 2.0*aesy*esy*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) - 0.5*sin(0.5*pi*(xt + yt + zt))**2&
         + 0.5*cos(0.5*pi*(xt + yt + zt))**2)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         pi**2*us0*(2*aesz**2*esz**2*sin(0.5*pi*(xt + yt +&
         zt))**2*cos(aesz*pi*zt)**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + aesz**2*esz*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) - 2.0*aesz*esz*sin(0.5*pi*(xt + yt + zt))*cos(0.5*pi*(xt&
         + yt + zt))*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) - 0.5*sin(0.5*pi*(xt + yt + zt))**2 + 0.5*cos(0.5*pi*(xt&
         + yt + zt))**2)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))) - mus*(2*aesx**2*esx**2*pi**2*us0*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesx*pi*xt)**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**3 + aesx**2*esx*pi**2*us0*sin(0.5*pi*(xt + yt +&
         zt))**2*cos(aesx*pi*xt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + 2.0*aesx*esx*pi**2*us0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 - 0.5*pi**2*us0*sin(0.5*pi*(xt +&
         yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) + 0.5*pi**2*us0*cos(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))) -&
         mus*(2*aesx*aesy*esx*esy*pi**2*vs0*sin(aesx*pi*xt)*sin(aesy*pi*yt)*cos(0.5*pi*(xt&
         + yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**3 - 1.0*aesx*esx*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 -&
         1.0*aesy*esy*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 + 0.5*pi**2*vs0*sin(0.5*pi*(xt +&
         yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) - 0.5*pi**2*vs0*cos(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))) +&
         2*mus*(2*aesx**2*esx**2*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesx*pi*xt)**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**3 + aesx**2*esx*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))**2*cos(aesx*pi*xt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 +&
         2*aesx*aesy*esx*esy*pi**2*vs0*sin(aesx*pi*xt)*sin(aesy*pi*yt)*cos(0.5*pi*(xt +&
         yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**3 -&
         2*aesx*aesz*esx*esz*pi**2*ws0*sin(aesx*pi*xt)*cos(aesz*pi*zt)/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 +&
         1.0*aesx*esx*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 -&
         1.0*aesy*esy*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2)/3 +&
         ros*(aesx*esx*pi*us0**2*sin(0.5*pi*(xt + yt + zt))**4*sin(aesx*pi*xt)/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 +&
         aesy*esy*pi*us0*vs0*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 -&
         aesz*esz*pi*us0*ws0*sin(0.5*pi*(xt + yt + zt))**2*cos(aesz*pi*zt)/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 +&
         2.0*pi*us0**2*sin(0.5*pi*(xt + yt + zt))**3*cos(0.5*pi*(xt + yt + zt))/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) -&
         1.0*pi*us0*vs0*sin(0.5*pi*(xt + yt + zt))**3*cos(0.5*pi*(xt + yt + zt))/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         1.0*pi*us0*vs0*sin(0.5*pi*(xt + yt + zt))*cos(0.5*pi*(xt + yt + zt))**3/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         1.0*pi*us0*ws0*sin(0.5*pi*(xt + yt + zt))*cos(0.5*pi*(xt + yt + zt))/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))) +&
         (-apgx*pgx*pi*sin(apgx*pi*xt) - apgxy*pgxy*pi*yt*sin(apgxy*pi*xt*yt) -&
         apgzx*pgzx*pi*zt*sin(apgzx*pi*xt*zt))*(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))
      case(6)
      !vssrc!
        mms_source = 2*aesy*aesz*esy*esz*mus*pi**2*ws0*sin(aesy*pi*yt)*cos(aesz*pi*zt)/(es0&
         + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 -&
         mus*(pi**2*vs0*(2*aesx**2*esx**2*sin(aesx*pi*xt)**2*cos(0.5*pi*(xt + yt +&
         zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + aesx**2*esx*cos(0.5*pi*(xt + yt +&
         zt))**2*cos(aesx*pi*xt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) - 2.0*aesx*esx*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) + 0.5*sin(0.5*pi*(xt + yt + zt))**2&
         - 0.5*cos(0.5*pi*(xt + yt + zt))**2)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         pi**2*vs0*(2*aesy**2*esy**2*sin(aesy*pi*yt)**2*cos(0.5*pi*(xt + yt +&
         zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + aesy**2*esy*cos(0.5*pi*(xt + yt +&
         zt))**2*cos(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) - 2.0*aesy*esy*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) + 0.5*sin(0.5*pi*(xt + yt + zt))**2&
         - 0.5*cos(0.5*pi*(xt + yt + zt))**2)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         pi**2*vs0*(2*aesz**2*esz**2*cos(0.5*pi*(xt + yt +&
         zt))**2*cos(aesz*pi*zt)**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + aesz**2*esz*sin(aesz*pi*zt)*cos(0.5*pi*(xt + yt +&
         zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) + 2.0*aesz*esz*sin(0.5*pi*(xt + yt + zt))*cos(0.5*pi*(xt&
         + yt + zt))*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) + 0.5*sin(0.5*pi*(xt + yt + zt))**2 - 0.5*cos(0.5*pi*(xt&
         + yt + zt))**2)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))) -&
         mus*(2*aesy**2*esy**2*pi**2*vs0*sin(aesy*pi*yt)**2*cos(0.5*pi*(xt + yt +&
         zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**3 + aesy**2*esy*pi**2*vs0*cos(0.5*pi*(xt + yt +&
         zt))**2*cos(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 - 2.0*aesy*esy*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 + 0.5*pi**2*vs0*sin(0.5*pi*(xt +&
         yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) - 0.5*pi**2*vs0*cos(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))) -&
         mus*(2*aesx*aesy*esx*esy*pi**2*us0*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesx*pi*xt)*sin(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 +&
         1.0*aesx*esx*pi**2*us0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 +&
         1.0*aesy*esy*pi**2*us0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 - 0.5*pi**2*us0*sin(0.5*pi*(xt +&
         yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt)) + 0.5*pi**2*us0*cos(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))) +&
         2*mus*(2*aesx*aesy*esx*esy*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesx*pi*xt)*sin(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 +&
         1.0*aesx*esx*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 +&
         2*aesy**2*esy**2*pi**2*vs0*sin(aesy*pi*yt)**2*cos(0.5*pi*(xt + yt +&
         zt))**2/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**3 + aesy**2*esy*pi**2*vs0*cos(0.5*pi*(xt + yt +&
         zt))**2*cos(aesy*pi*yt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 -&
         2*aesy*aesz*esy*esz*pi**2*ws0*sin(aesy*pi*yt)*cos(aesz*pi*zt)/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 -&
         1.0*aesy*esy*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2)/3 +&
         ros*(aesx*esx*pi*us0*vs0*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 +&
         aesy*esy*pi*vs0**2*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))**4/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 -&
         aesz*esz*pi*vs0*ws0*cos(0.5*pi*(xt + yt + zt))**2*cos(aesz*pi*zt)/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 -&
         1.0*pi*us0*vs0*sin(0.5*pi*(xt + yt + zt))**3*cos(0.5*pi*(xt + yt + zt))/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         1.0*pi*us0*vs0*sin(0.5*pi*(xt + yt + zt))*cos(0.5*pi*(xt + yt + zt))**3/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) -&
         2.0*pi*vs0**2*sin(0.5*pi*(xt + yt + zt))*cos(0.5*pi*(xt + yt + zt))**3/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) -&
         1.0*pi*vs0*ws0*sin(0.5*pi*(xt + yt + zt))*cos(0.5*pi*(xt + yt + zt))/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))) +&
         (-apgxy*pgxy*pi*xt*sin(apgxy*pi*xt*yt) - apgy*pgy*pi*sin(apgy*pi*yt) +&
         apgyz*pgyz*pi*zt*cos(apgyz*pi*yt*zt))*(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))
      case(7)
      !wssrc!
        mms_source = -2*aesz**2*esz**2*mus*pi**2*ws0*cos(aesz*pi*zt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 -&
         aesz**2*esz*mus*pi**2*ws0*sin(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 -&
         mus*(aesx**2*esx*pi**2*ws0*(2*esx*sin(aesx*pi*xt)**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         cos(aesx*pi*xt))/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + aesy**2*esy*pi**2*ws0*(2*esy*sin(aesy*pi*yt)**2/(es0&
         + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         cos(aesy*pi*yt))/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + aesz**2*esz*pi**2*ws0*(2*esz*cos(aesz*pi*zt)**2/(es0&
         + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         sin(aesz*pi*zt))/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2) - mus*(-2*aesx*aesz*esx*esz*pi**2*us0*sin(0.5*pi*(xt&
         + yt + zt))**2*sin(aesx*pi*xt)*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 +&
         1.0*aesx*esx*pi**2*us0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 -&
         1.0*aesz*esz*pi**2*us0*sin(0.5*pi*(xt + yt + zt))*cos(0.5*pi*(xt + yt +&
         zt))*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 - 0.5*pi**2*us0*sin(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         0.5*pi**2*us0*cos(0.5*pi*(xt + yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))) -&
         mus*(-2*aesy*aesz*esy*esz*pi**2*vs0*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt +&
         zt))**2*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**3 - 1.0*aesy*esy*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 +&
         1.0*aesz*esz*pi**2*vs0*sin(0.5*pi*(xt + yt + zt))*cos(0.5*pi*(xt + yt +&
         zt))*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**2 + 0.5*pi**2*vs0*sin(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) -&
         0.5*pi**2*vs0*cos(0.5*pi*(xt + yt + zt))**2/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))) +&
         2*mus*(-2*aesx*aesz*esx*esz*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))**2*sin(aesx*pi*xt)*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 +&
         1.0*aesx*esx*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesx*pi*xt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 -&
         2*aesy*aesz*esy*esz*pi**2*vs0*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt +&
         zt))**2*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) +&
         esz*sin(aesz*pi*zt))**3 - 1.0*aesy*esy*pi**2*vs0*sin(0.5*pi*(xt + yt +&
         zt))*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 +&
         2*aesz**2*esz**2*pi**2*ws0*cos(aesz*pi*zt)**2/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**3 +&
         aesz**2*esz*pi**2*ws0*sin(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2)/3 +&
         ros*(aesx*esx*pi*us0*ws0*sin(0.5*pi*(xt + yt + zt))**2*sin(aesx*pi*xt)/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 +&
         aesy*esy*pi*vs0*ws0*sin(aesy*pi*yt)*cos(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 -&
         aesz*esz*pi*ws0**2*cos(aesz*pi*zt)/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))**2 + 1.0*pi*us0*ws0*sin(0.5*pi*(xt&
         + yt + zt))*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) - 1.0*pi*vs0*ws0*sin(0.5*pi*(xt +&
         yt + zt))*cos(0.5*pi*(xt + yt + zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))) +&
         (apgyz*pgyz*pi*yt*cos(apgyz*pi*yt*zt) + apgz*pgz*pi*cos(apgz*pi*zt) -&
         apgzx*pgzx*pi*xt*sin(apgzx*pi*xt*zt))*(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt))
      case(8)
      !tgsrc!
        mms_source = Cpg*rog*((-atgx*pi*tgx*sin(atgx*pi*xt) -&
         atgxy*pi*tgxy*yt*sin(atgxy*pi*xt*yt) -&
         atgzx*pi*tgzx*zt*sin(atgzx*pi*xt*zt))*(-avgyz*pi*vgyz*yt*cos(avgyz*pi*yt*zt) +&
         avgz*pi*vgz*sin(avgz*pi*zt) + avgzx*pi*vgzx*xt*sin(avgzx*pi*xt*zt) +&
         awgxy*pi*wgxy*xt*cos(awgxy*pi*xt*yt) + awgy*pi*wgy*cos(awgy*pi*yt) +&
         awgyz*pi*wgyz*zt*cos(awgyz*pi*yt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (-atgxy*pi*tgxy*xt*sin(atgxy*pi*xt*yt) - atgy*pi*tgy*sin(atgy*pi*yt) +&
         atgyz*pi*tgyz*zt*cos(atgyz*pi*yt*zt))*(augyz*pi*ugyz*yt*cos(augyz*pi*yt*zt) -&
         augz*pi*ugz*sin(augz*pi*zt) - augzx*pi*ugzx*xt*sin(augzx*pi*xt*zt) +&
         awgx*pi*wgx*sin(awgx*pi*xt) - awgxy*pi*wgxy*yt*cos(awgxy*pi*xt*yt) +&
         awgzx*pi*wgzx*zt*sin(awgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) +&
         (atgyz*pi*tgyz*yt*cos(atgyz*pi*yt*zt) + atgz*pi*tgz*cos(atgz*pi*zt) -&
         atgzx*pi*tgzx*xt*sin(atgzx*pi*xt*zt))*(augxy*pi*ugxy*xt*sin(augxy*pi*xt*yt) +&
         augy*pi*ugy*sin(augy*pi*yt) - augyz*pi*ugyz*zt*cos(augyz*pi*yt*zt) +&
         avgx*pi*vgx*cos(avgx*pi*xt) - avgxy*pi*vgxy*yt*sin(avgxy*pi*xt*yt) -&
         avgzx*pi*vgzx*zt*sin(avgzx*pi*xt*zt))/(-es0 - esx*cos(aesx*pi*xt) -&
         esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0))*(-es0 - esx*cos(aesx*pi*xt)&
         - esy*cos(aesy*pi*yt) - esz*sin(aesz*pi*zt) + 1.0) -&
         kg*(-atgx**2*pi**2*tgx*cos(atgx*pi*xt) -&
         atgxy**2*pi**2*tgxy*yt**2*cos(atgxy*pi*xt*yt) -&
         atgzx**2*pi**2*tgzx*zt**2*cos(atgzx*pi*xt*zt)) -&
         kg*(-atgxy**2*pi**2*tgxy*xt**2*cos(atgxy*pi*xt*yt) -&
         atgy**2*pi**2*tgy*cos(atgy*pi*yt) -&
         atgyz**2*pi**2*tgyz*zt**2*sin(atgyz*pi*yt*zt)) -&
         kg*(-atgyz**2*pi**2*tgyz*yt**2*sin(atgyz*pi*yt*zt) -&
         atgz**2*pi**2*tgz*sin(atgz*pi*zt) -&
         atgzx**2*pi**2*tgzx*xt**2*cos(atgzx*pi*xt*zt))
      case(9)
      !tssrc!
        mms_source = Cps*ros*(us0*(-atsx*pi*tsx*sin(atsx*pi*xt) -&
         atsxy*pi*tsxy*yt*sin(atsxy*pi*xt*yt) -&
         atszx*pi*tszx*zt*sin(atszx*pi*xt*zt))*sin(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         vs0*(-atsxy*pi*tsxy*xt*sin(atsxy*pi*xt*yt) - atsy*pi*tsy*sin(atsy*pi*yt) +&
         atsyz*pi*tsyz*zt*cos(atsyz*pi*yt*zt))*cos(0.5*pi*(xt + yt + zt))**2/(es0 +&
         esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) +&
         ws0*(atsyz*pi*tsyz*yt*cos(atsyz*pi*yt*zt) + atsz*pi*tsz*cos(atsz*pi*zt) -&
         atszx*pi*tszx*xt*sin(atszx*pi*xt*zt))/(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)))*(es0 + esx*cos(aesx*pi*xt) +&
         esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt)) -&
         ks*(-atsx**2*pi**2*tsx*cos(atsx*pi*xt) -&
         atsxy**2*pi**2*tsxy*yt**2*cos(atsxy*pi*xt*yt) -&
         atszx**2*pi**2*tszx*zt**2*cos(atszx*pi*xt*zt)) -&
         ks*(-atsxy**2*pi**2*tsxy*xt**2*cos(atsxy*pi*xt*yt) -&
         atsy**2*pi**2*tsy*cos(atsy*pi*yt) -&
         atsyz**2*pi**2*tsyz*zt**2*sin(atsyz*pi*yt*zt)) -&
         ks*(-atsyz**2*pi**2*tsyz*yt**2*sin(atsyz*pi*yt*zt) -&
         atsz**2*pi**2*tsz*sin(atsz*pi*zt) -&
         atszx**2*pi**2*tszx*xt**2*cos(atszx*pi*xt*zt))
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
