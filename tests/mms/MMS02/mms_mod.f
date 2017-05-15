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
      call allocate_mms_vars

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
                esx=0.0d0, esy=0.0d0, esz=0.0d0, &
                esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
                aesx=0.5d0, aesy=0.5d0, aesz=0.5d0, &
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
        mms_function = pg0 + pgx*Cos(apgx*Pi*xt) + pgy*Cos(apgy*Pi*yt)+&
         pgxy*Cos(apgxy*Pi*xt*yt) + pgzx*Cos(apgzx*Pi*xt*zt) + &
         pgz*Sin(apgz*Pi*zt) + pgyz*Sin(apgyz*Pi*yt*zt)
      case(2)
      !ug!
        mms_function = awgy*Pi*wgy*Cos(awgy*Pi*yt) + &
         awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
         avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
         awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + &
         avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
         avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt)
      case(3)
      !vg!
        mms_function = -(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
         augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + &
         awgx*Pi*wgx*Sin(awgx*Pi*xt) - augz*Pi*ugz*Sin(augz*Pi*zt) - &
         augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
         awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt)
      case(4)
      !wg!
        mms_function = avgx*Pi*vgx*Cos(avgx*Pi*xt) - &
         augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
         augy*Pi*ugy*Sin(augy*Pi*yt) + &
         augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
         avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
         avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt)
      case(5)
      !us!
        mms_function = us0*Sin((Pi*(xt + yt + zt))/2.0d0)**2
      case(6)
      !vs!
        mms_function = vs0*Cos((Pi*(xt + yt + zt))/2.0d0)**2
      case(7)
      !ws!
        mms_function = ws0
      case(8)
      !tg!
        mms_function = Tg0 + Tgx*Cos(aTgx*Pi*xt) + Tgy*Cos(aTgy*Pi*yt)+&
        Tgxy*Cos(aTgxy*Pi*xt*yt) + Tgzx*Cos(aTgzx*Pi*xt*zt) + &
        Tgz*Sin(aTgz*Pi*zt) + Tgyz*Sin(aTgyz*Pi*yt*zt)
      case(9)
      !ts!
        mms_function = Ts0 + Tsx*Cos(aTsx*Pi*xt) + Tsy*Cos(aTsy*Pi*yt)+&
         Tsxy*Cos(aTsxy*Pi*xt*yt) + Tszx*Cos(aTszx*Pi*xt*zt) + &
         Tsz*Sin(aTsz*Pi*zt) + Tsyz*Sin(aTsyz*Pi*yt*zt)
      case(10)
      !ep_g!
        mms_function = 1.0d0 - es0 - esx*Cos(aesx*Pi*xt) - &
         esy*Cos(aesy*Pi*yt) - esxy*Cos(aesxy*Pi*xt*yt) - &
         eszx*Cos(aeszx*Pi*xt*zt) - esz*Sin(aesz*Pi*zt) - &
         esyz*Sin(aesyz*Pi*yt*zt)
      case(11)
      !rop_s!
        mms_function = ros*(es0 + esx*Cos(aesx*Pi*xt) + &
         esy*Cos(aesy*Pi*yt) + esxy*Cos(aesxy*Pi*xt*yt) + &
         eszx*Cos(aeszx*Pi*xt*zt) + esz*Sin(aesz*Pi*zt) + &
         esyz*Sin(aesyz*Pi*yt*zt))
      case(12)
      !theta_m!
        mms_function = Ths0 + Thsx*Cos(aThsx*Pi*xt) + &
         Thsy*Cos(aThsy*Pi*yt) + Thsxy*Cos(aThsxy*Pi*xt*yt) + &
         Thszx*Cos(aThszx*Pi*xt*zt) + Thsz*Sin(aThsz*Pi*zt) + &
         Thsyz*Sin(aThsyz*Pi*yt*zt)
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
                esx=0.0d0, esy=0.0d0, esz=0.0d0, &
                esxy=0.0d0, esyz=0.0d0, eszx=0.0d0, &
                aesx=0.5d0, aesy=0.5d0, aesz=0.5d0, &
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
        mms_source = -(mug*(awgxy**3*Pi**3*wgxy*xt*yt**2*&
              Cos(awgxy*Pi*xt*yt) + &
              2*awgxy**2*Pi**2*wgxy*yt*Sin(awgxy*Pi*xt*yt))) - &
         mug*(-(awgxy**3*Pi**3*wgxy*xt*yt**2*Cos(awgxy*Pi*xt*yt)) + &
            2*avgzx**2*Pi**2*vgzx*zt*Cos(avgzx*Pi*xt*zt) - &
            2*awgxy**2*Pi**2*wgxy*yt*Sin(awgxy*Pi*xt*yt) - &
            avgzx**3*Pi**3*vgzx*xt*zt**2*Sin(avgzx*Pi*xt*zt)) - &
         mug*(-2*avgzx**2*Pi**2*vgzx*zt*Cos(avgzx*Pi*xt*zt) + &
            avgzx**3*Pi**3*vgzx*xt*zt**2*Sin(avgzx*Pi*xt*zt)) + &
         (-(apgx*pgx*Pi*Sin(apgx*Pi*xt)) - apgxy*pgxy*Pi*yt*&
            Sin(apgxy*Pi*xt*yt) - &
            apgzx*pgzx*Pi*zt*Sin(apgzx*Pi*xt*zt))*&
          (1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt)) - &
         mug*(-(awgy**3*Pi**3*wgy*Cos(awgy*Pi*yt)) - &
            awgxy**3*Pi**3*wgxy*xt**3*Cos(awgxy*Pi*xt*yt) + &
            avgyz**3*Pi**3*vgyz*yt**3*Cos(avgyz*Pi*yt*zt) - &
            awgyz**3*Pi**3*wgyz*zt**3*Cos(awgyz*Pi*yt*zt) + &
            awgxy*Pi*wgxy*(-(awgxy**2*Pi**2*xt*yt**2*&
            Cos(awgxy*Pi*xt*yt)) - &
               2*awgxy*Pi*yt*Sin(awgxy*Pi*xt*yt)) - &
            avgz**3*Pi**3*vgz*Sin(avgz*Pi*zt) - &
            avgzx**3*Pi**3*vgzx*xt**3*Sin(avgzx*Pi*xt*zt) + &
            avgzx*Pi*vgzx*(2*avgzx*Pi*zt*Cos(avgzx*Pi*xt*zt) - &
               avgzx**2*Pi**2*xt*zt**2*Sin(avgzx*Pi*xt*zt)) - &
            avgyz*Pi*(-(avgyz**2*Pi**2*vgyz*yt*zt**2*&
            Cos(avgyz*Pi*yt*zt)) - &
               2*avgyz*Pi*vgyz*zt*Sin(avgyz*Pi*yt*zt)) + &
            awgyz*Pi*wgyz*(-(awgyz**2*Pi**2*yt**2*zt*&
            Cos(awgyz*Pi*yt*zt)) - &
               2*awgyz*Pi*yt*Sin(awgyz*Pi*yt*zt))) + &
         rog*(1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt))*&
          (2*(awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt) + &
               avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt) - &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) + &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt))*&
             (awgy*Pi*wgy*Cos(awgy*Pi*yt) + &
             awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + &
               avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt)) + &
            (awgy*Pi*wgy*Cos(awgy*Pi*yt) + &
            awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + &
               avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt))*&
             (-(awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) + &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) - &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            (awgy*Pi*wgy*Cos(awgy*Pi*yt) + &
            awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + &
               avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt))*&
             (-(avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt)) - &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) - &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt) + &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            (avgx*Pi*vgx*Cos(avgx*Pi*xt) - &
            augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + &
               augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt))*&
             (avgz**2*Pi**2*vgz*Cos(avgz*Pi*zt) + &
               avgzx**2*Pi**2*vgzx*xt**2*Cos(avgzx*Pi*xt*zt) + &
               awgyz*Pi*wgyz*Cos(awgyz*Pi*yt*zt) + &
               avgyz**2*Pi**2*vgyz*yt**2*Sin(avgyz*Pi*yt*zt) - &
               awgyz**2*Pi**2*wgyz*yt*zt*Sin(awgyz*Pi*yt*zt)) + &
            (-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + &
               awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - &
               augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt))*&
             (-(avgyz*Pi*vgyz*Cos(avgyz*Pi*yt*zt)) - &
               awgy**2*Pi**2*wgy*Sin(awgy*Pi*yt) - &
               awgxy**2*Pi**2*wgxy*xt**2*Sin(awgxy*Pi*xt*yt) + &
               avgyz**2*Pi**2*vgyz*yt*zt*Sin(avgyz*Pi*yt*zt) - &
               awgyz**2*Pi**2*wgyz*zt**2*Sin(awgyz*Pi*yt*zt)))
      case(3)
      !vgsrc!
        mms_source = -(mug*(-(awgxy**3*Pi**3*wgxy*xt**2*yt*&
        Cos(awgxy*Pi*xt*yt)) - &
              2*awgxy**2*Pi**2*wgxy*xt*Sin(awgxy*Pi*xt*yt))) + &
         (apgyz*pgyz*Pi*zt*Cos(apgyz*Pi*yt*zt) - &
         apgy*pgy*Pi*Sin(apgy*Pi*yt) - &
            apgxy*pgxy*Pi*xt*Sin(apgxy*Pi*xt*yt))*&
          (1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt)) - &
         mug*(awgxy**3*Pi**3*wgxy*xt**2*yt*Cos(awgxy*Pi*xt*yt) - &
            augyz**3*Pi**3*ugyz*yt*zt**2*Cos(augyz*Pi*yt*zt) + &
            2*awgxy**2*Pi**2*wgxy*xt*Sin(awgxy*Pi*xt*yt) - &
            2*augyz**2*Pi**2*ugyz*zt*Sin(augyz*Pi*yt*zt)) - &
         mug*(augyz**3*Pi**3*ugyz*yt*zt**2*Cos(augyz*Pi*yt*zt) + &
            2*augyz**2*Pi**2*ugyz*zt*Sin(augyz*Pi*yt*zt)) - &
         mug*(awgxy**3*Pi**3*wgxy*yt**3*Cos(awgxy*Pi*xt*yt) - &
            augyz**3*Pi**3*ugyz*yt**3*Cos(augyz*Pi*yt*zt) - &
            awgx**3*Pi**3*wgx*Sin(awgx*Pi*xt) - &
            awgxy*Pi*(-(awgxy**2*Pi**2*wgxy*xt**2*yt*&
            Cos(awgxy*Pi*xt*yt)) - &
               2*awgxy*Pi*wgxy*xt*Sin(awgxy*Pi*xt*yt)) + &
            augz**3*Pi**3*ugz*Sin(augz*Pi*zt) + &
            augzx**3*Pi**3*ugzx*xt**3*Sin(augzx*Pi*xt*zt) - &
            augzx*Pi*(2*augzx*Pi*ugzx*zt*Cos(augzx*Pi*xt*zt) - &
               augzx**2*Pi**2*ugzx*xt*zt**2*Sin(augzx*Pi*xt*zt)) - &
            awgzx**3*Pi**3*wgzx*zt**3*Sin(awgzx*Pi*xt*zt) + &
            awgzx*Pi*wgzx*(2*awgzx*Pi*xt*Cos(awgzx*Pi*xt*zt) - &
               awgzx**2*Pi**2*xt**2*zt*Sin(awgzx*Pi*xt*zt)) + &
            augyz*Pi*ugyz*(-(augyz**2*Pi**2*yt*zt**2*&
            Cos(augyz*Pi*yt*zt)) - &
               2*augyz*Pi*zt*Sin(augyz*Pi*yt*zt))) + &
         rog*(1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt))*&
          ((awgx**2*Pi**2*wgx*Cos(awgx*Pi*xt) - &
               augzx**2*Pi**2*ugzx*xt*zt*Cos(augzx*Pi*xt*zt) + &
               awgzx**2*Pi**2*wgzx*zt**2*Cos(awgzx*Pi*xt*zt) + &
               awgxy**2*Pi**2*wgxy*yt**2*Sin(awgxy*Pi*xt*yt) - &
               augzx*Pi*ugzx*Sin(augzx*Pi*xt*zt))*&
             (awgy*Pi*wgy*Cos(awgy*Pi*yt) + &
             awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + &
               avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt)) + &
            (awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt) + &
               avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt) - &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) + &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt))*&
             (-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + &
               awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - &
               augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt)) + &
            (avgx*Pi*vgx*Cos(avgx*Pi*xt) - &
            augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + &
               augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt))*&
             (-(augz**2*Pi**2*ugz*Cos(augz*Pi*zt)) - &
               augzx**2*Pi**2*ugzx*xt**2*Cos(augzx*Pi*xt*zt) + &
               awgzx**2*Pi**2*wgzx*xt*zt*Cos(awgzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*Sin(awgzx*Pi*xt*zt) - &
               augyz**2*Pi**2*ugyz*yt**2*Sin(augyz*Pi*yt*zt)) + &
            2*(-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + &
               awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - &
               augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt))*&
             (-(awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) + &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) - &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            (-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + &
               awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - &
               augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt))*&
             (-(avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt)) - &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) - &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt) + &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt))) 
      case(4)
      !wgsrc!
        mms_source = -(mug*(2*avgzx**2*Pi**2*vgzx*xt*&
        Cos(avgzx*Pi*xt*zt) - &
              avgzx**3*Pi**3*vgzx*xt**2*zt*Sin(avgzx*Pi*xt*zt))) + &
         (apgz*pgz*Pi*Cos(apgz*Pi*zt) + &
         apgyz*pgyz*Pi*yt*Cos(apgyz*Pi*yt*zt) - &
            apgzx*pgzx*Pi*xt*Sin(apgzx*Pi*xt*zt))*&
          (1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt)) - &
         mug*(-(augyz**3*Pi**3*ugyz*yt**2*zt*Cos(augyz*Pi*yt*zt)) - &
            2*augyz**2*Pi**2*ugyz*yt*Sin(augyz*Pi*yt*zt)) - &
         mug*(-2*avgzx**2*Pi**2*vgzx*xt*Cos(avgzx*Pi*xt*zt) + &
            augyz**3*Pi**3*ugyz*yt**2*zt*Cos(augyz*Pi*yt*zt) + &
            avgzx**3*Pi**3*vgzx*xt**2*zt*Sin(avgzx*Pi*xt*zt) + &
            2*augyz**2*Pi**2*ugyz*yt*Sin(augyz*Pi*yt*zt)) - &
         mug*(-(avgx**3*Pi**3*vgx*Cos(avgx*Pi*xt)) + &
            augyz**3*Pi**3*ugyz*zt**3*Cos(augyz*Pi*yt*zt) - &
            augy**3*Pi**3*ugy*Sin(augy*Pi*yt) - &
            augxy**3*Pi**3*ugxy*xt**3*Sin(augxy*Pi*xt*yt) + &
            augxy*Pi*ugxy*(2*augxy*Pi*yt*Cos(augxy*Pi*xt*yt) - &
               augxy**2*Pi**2*xt*yt**2*Sin(augxy*Pi*xt*yt)) + &
            avgxy**3*Pi**3*vgxy*yt**3*Sin(avgxy*Pi*xt*yt) - &
            avgxy*Pi*(2*avgxy*Pi*vgxy*xt*Cos(avgxy*Pi*xt*yt) - &
               avgxy**2*Pi**2*vgxy*xt**2*yt*Sin(avgxy*Pi*xt*yt)) + &
            avgzx**3*Pi**3*vgzx*zt**3*Sin(avgzx*Pi*xt*zt) - &
            avgzx*Pi*(2*avgzx*Pi*vgzx*xt*Cos(avgzx*Pi*xt*zt) - &
               avgzx**2*Pi**2*vgzx*xt**2*zt*Sin(avgzx*Pi*xt*zt)) - &
            augyz*Pi*(-(augyz**2*Pi**2*ugyz*yt**2*zt*&
            Cos(augyz*Pi*yt*zt)) - &
               2*augyz*Pi*ugyz*yt*Sin(augyz*Pi*yt*zt))) + &
         rog*(1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
            esxy*Cos(aesxy*Pi*xt*yt) - eszx*Cos(aeszx*Pi*xt*zt) - &
            esz*Sin(aesz*Pi*zt) - esyz*Sin(aesyz*Pi*yt*zt))*&
          ((augxy**2*Pi**2*ugxy*xt*yt*Cos(augxy*Pi*xt*yt) - &
               avgxy**2*Pi**2*vgxy*yt**2*Cos(avgxy*Pi*xt*yt) - &
               avgzx**2*Pi**2*vgzx*zt**2*Cos(avgzx*Pi*xt*zt) - &
               avgx**2*Pi**2*vgx*Sin(avgx*Pi*xt) + &
               augxy*Pi*ugxy*Sin(augxy*Pi*xt*yt))*&
             (awgy*Pi*wgy*Cos(awgy*Pi*yt) + &
             awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + &
               avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt)) + &
            (awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt) + &
               avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt) - &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) + &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt))*&
             (avgx*Pi*vgx*Cos(avgx*Pi*xt) - &
             augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + &
               augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt)) + &
            (avgx*Pi*vgx*Cos(avgx*Pi*xt) - &
            augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + &
               augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt))*&
             (-(awgxy*Pi*wgxy*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) + &
               awgxy**2*Pi**2*wgxy*xt*yt*Sin(awgxy*Pi*xt*yt) - &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            2*(avgx*Pi*vgx*Cos(avgx*Pi*xt) - &
            augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + &
               augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt))*&
             (-(avgzx**2*Pi**2*vgzx*xt*zt*Cos(avgzx*Pi*xt*zt)) - &
               augyz*Pi*ugyz*Cos(augyz*Pi*yt*zt) - &
               avgzx*Pi*vgzx*Sin(avgzx*Pi*xt*zt) + &
               augyz**2*Pi**2*ugyz*yt*zt*Sin(augyz*Pi*yt*zt)) + &
            (-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
               augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + &
               awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - &
               augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt))*&
             (augy**2*Pi**2*ugy*Cos(augy*Pi*yt) + &
               augxy**2*Pi**2*ugxy*xt**2*Cos(augxy*Pi*xt*yt) - &
               avgxy**2*Pi**2*vgxy*xt*yt*Cos(avgxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*Sin(avgxy*Pi*xt*yt) + &
               augyz**2*Pi**2*ugyz*zt**2*Sin(augyz*Pi*yt*zt)))
      case(5)
      !ussrc!
        mms_source = -(mus*Pi**2*us0*&
        Cos((Pi*(xt + yt + zt))/2.)**2)/2. + &
         (mus*Pi**2*vs0*Cos((Pi*(xt + yt + zt))/2.)**2)/2. + &
         (-(apgx*pgx*Pi*Sin(apgx*Pi*xt)) - &
         apgxy*pgxy*Pi*yt*Sin(apgxy*Pi*xt*yt) - &
            apgzx*pgzx*Pi*zt*Sin(apgzx*Pi*xt*zt))*&
          (es0 + esx*Cos(aesx*Pi*xt) + esy*Cos(aesy*Pi*yt) + &
            esxy*Cos(aesxy*Pi*xt*yt) + eszx*Cos(aeszx*Pi*xt*zt) + &
            esz*Sin(aesz*Pi*zt) + esyz*Sin(aesyz*Pi*yt*zt)) + &
         (mus*Pi**2*us0*Sin((Pi*(xt + yt + zt))/2.)**2)/2. - &
         (mus*Pi**2*vs0*Sin((Pi*(xt + yt + zt))/2.)**2)/2. - &
         3*mus*us0*((Pi**2*Cos((Pi*(xt + yt + zt))/2.)**2)/2. - &
            (Pi**2*Sin((Pi*(xt + yt + zt))/2.)**2)/2.) + &
         ros*(es0 + esx*Cos(aesx*Pi*xt) + esy*Cos(aesy*Pi*yt) + &
            esxy*Cos(aesxy*Pi*xt*yt) + eszx*Cos(aeszx*Pi*xt*zt) + &
            esz*Sin(aesz*Pi*zt) + esyz*Sin(aesyz*Pi*yt*zt))*&
          (Pi*us0*ws0*Cos((Pi*(xt + yt + zt))/2.)*&
          Sin((Pi*(xt + yt + zt))/2.) + &
            Pi*us0*vs0*Cos((Pi*(xt + yt + zt))/2.)**3*&
            Sin((Pi*(xt + yt + zt))/2.) + &
            2*Pi*us0**2*Cos((Pi*(xt + yt + zt))/2.)*&
            Sin((Pi*(xt + yt + zt))/2.)**3 - &
            Pi*us0*vs0*Cos((Pi*(xt + yt + zt))/2.)*&
            Sin((Pi*(xt + yt + zt))/2.)**3)
      case(6)
      !vssrc!
        mms_source = -(mus*Pi**2*us0*&
        Cos((Pi*(xt + yt + zt))/2.)**2)/2. + &
         (mus*Pi**2*vs0*Cos((Pi*(xt + yt + zt))/2.)**2)/2. + &
         (apgyz*pgyz*Pi*zt*Cos(apgyz*Pi*yt*zt) - &
         apgy*pgy*Pi*Sin(apgy*Pi*yt) - &
            apgxy*pgxy*Pi*xt*Sin(apgxy*Pi*xt*yt))*&
          (es0 + esx*Cos(aesx*Pi*xt) + esy*Cos(aesy*Pi*yt) + &
            esxy*Cos(aesxy*Pi*xt*yt) + eszx*Cos(aeszx*Pi*xt*zt) + &
            esz*Sin(aesz*Pi*zt) + esyz*Sin(aesyz*Pi*yt*zt)) + &
         (mus*Pi**2*us0*Sin((Pi*(xt + yt + zt))/2.)**2)/2. - &
         (mus*Pi**2*vs0*Sin((Pi*(xt + yt + zt))/2.)**2)/2. - &
         3*mus*vs0*(-(Pi**2*Cos((Pi*(xt + yt + zt))/2.)**2)/2. + &
            (Pi**2*Sin((Pi*(xt + yt + zt))/2.)**2)/2.) + &
         ros*(es0 + esx*Cos(aesx*Pi*xt) + esy*Cos(aesy*Pi*yt) + &
            esxy*Cos(aesxy*Pi*xt*yt) + eszx*Cos(aeszx*Pi*xt*zt) + &
            esz*Sin(aesz*Pi*zt) + esyz*Sin(aesyz*Pi*yt*zt))*&
          (-(Pi*vs0*ws0*Cos((Pi*(xt + yt + zt))/2.)*&
          Sin((Pi*(xt + yt + zt))/2.)) + &
            Pi*us0*vs0*Cos((Pi*(xt + yt + zt))/2.)**3*&
            Sin((Pi*(xt + yt + zt))/2.) - &
            2*Pi*vs0**2*Cos((Pi*(xt + yt + zt))/2.)**3*&
            Sin((Pi*(xt + yt + zt))/2.) - &
            Pi*us0*vs0*Cos((Pi*(xt + yt + zt))/2.)*&
            Sin((Pi*(xt + yt + zt))/2.)**3)
      case(7)
      !wssrc!
        mms_source = -(mus*Pi**2*us0*&
        Cos((Pi*(xt + yt + zt))/2.)**2)/2. + &
         (mus*Pi**2*vs0*Cos((Pi*(xt + yt + zt))/2.)**2)/2. + &
         (apgz*pgz*Pi*Cos(apgz*Pi*zt) + &
         apgyz*pgyz*Pi*yt*Cos(apgyz*Pi*yt*zt) - &
            apgzx*pgzx*Pi*xt*Sin(apgzx*Pi*xt*zt))*&
          (es0 + esx*Cos(aesx*Pi*xt) + esy*Cos(aesy*Pi*yt) + &
            esxy*Cos(aesxy*Pi*xt*yt) + eszx*Cos(aeszx*Pi*xt*zt) + &
            esz*Sin(aesz*Pi*zt) + esyz*Sin(aesyz*Pi*yt*zt)) + &
         (mus*Pi**2*us0*Sin((Pi*(xt + yt + zt))/2.)**2)/2. - &
         (mus*Pi**2*vs0*Sin((Pi*(xt + yt + zt))/2.)**2)/2. + &
         ros*(es0 + esx*Cos(aesx*Pi*xt) + esy*Cos(aesy*Pi*yt) + &
            esxy*Cos(aesxy*Pi*xt*yt) + eszx*Cos(aeszx*Pi*xt*zt) + &
            esz*Sin(aesz*Pi*zt) + esyz*Sin(aesyz*Pi*yt*zt))*&
          (Pi*us0*ws0*Cos((Pi*(xt + yt + zt))/2.)*&
          Sin((Pi*(xt + yt + zt))/2.) - &
            Pi*vs0*ws0*Cos((Pi*(xt + yt + zt))/2.)*&
            Sin((Pi*(xt + yt + zt))/2.))
      case(8)
      !tgsrc!
        mms_source = -(kg*(-(aTgx**2*Pi**2*Tgx*Cos(aTgx*Pi*xt)) - &
              aTgxy**2*Pi**2*Tgxy*yt**2*Cos(aTgxy*Pi*xt*yt) - &
              aTgzx**2*Pi**2*Tgzx*zt**2*Cos(aTgzx*Pi*xt*zt))) + &
         Cpg*rog*((-(aTgx*Pi*Tgx*Sin(aTgx*Pi*xt)) - &
         aTgxy*Pi*Tgxy*yt*Sin(aTgxy*Pi*xt*yt) - &
               aTgzx*Pi*Tgzx*zt*Sin(aTgzx*Pi*xt*zt))*&
             (awgy*Pi*wgy*Cos(awgy*Pi*yt) + &
             awgxy*Pi*wgxy*xt*Cos(awgxy*Pi*xt*yt) - &
               avgyz*Pi*vgyz*yt*Cos(avgyz*Pi*yt*zt) + &
               awgyz*Pi*wgyz*zt*Cos(awgyz*Pi*yt*zt) + &
               avgz*Pi*vgz*Sin(avgz*Pi*zt) + &
               avgzx*Pi*vgzx*xt*Sin(avgzx*Pi*xt*zt)) + &
            (aTgz*Pi*Tgz*Cos(aTgz*Pi*zt) + &
            aTgyz*Pi*Tgyz*yt*Cos(aTgyz*Pi*yt*zt) - &
               aTgzx*Pi*Tgzx*xt*Sin(aTgzx*Pi*xt*zt))*&
             (avgx*Pi*vgx*Cos(avgx*Pi*xt) - &
             augyz*Pi*ugyz*zt*Cos(augyz*Pi*yt*zt) + &
               augy*Pi*ugy*Sin(augy*Pi*yt) + &
               augxy*Pi*ugxy*xt*Sin(augxy*Pi*xt*yt) - &
               avgxy*Pi*vgxy*yt*Sin(avgxy*Pi*xt*yt) - &
               avgzx*Pi*vgzx*zt*Sin(avgzx*Pi*xt*zt)) + &
            (aTgyz*Pi*Tgyz*zt*Cos(aTgyz*Pi*yt*zt) - &
            aTgy*Pi*Tgy*Sin(aTgy*Pi*yt) - &
               aTgxy*Pi*Tgxy*xt*Sin(aTgxy*Pi*xt*yt))*&
             (-(awgxy*Pi*wgxy*yt*Cos(awgxy*Pi*xt*yt)) + &
             augyz*Pi*ugyz*yt*Cos(augyz*Pi*yt*zt) + &
               awgx*Pi*wgx*Sin(awgx*Pi*xt) - &
               augz*Pi*ugz*Sin(augz*Pi*zt) - &
               augzx*Pi*ugzx*xt*Sin(augzx*Pi*xt*zt) + &
               awgzx*Pi*wgzx*zt*Sin(awgzx*Pi*xt*zt)))*&
          (1. - es0 - esx*Cos(aesx*Pi*xt) - esy*Cos(aesy*Pi*yt) - &
          esxy*Cos(aesxy*Pi*xt*yt) - &
            eszx*Cos(aeszx*Pi*xt*zt) - esz*Sin(aesz*Pi*zt) - &
            esyz*Sin(aesyz*Pi*yt*zt)) - &
         kg*(-(aTgzx**2*Pi**2*Tgzx*xt**2*Cos(aTgzx*Pi*xt*zt)) - &
         aTgz**2*Pi**2*Tgz*Sin(aTgz*Pi*zt) - &
            aTgyz**2*Pi**2*Tgyz*yt**2*Sin(aTgyz*Pi*yt*zt)) - &
         kg*(-(aTgy**2*Pi**2*Tgy*Cos(aTgy*Pi*yt)) - &
         aTgxy**2*Pi**2*Tgxy*xt**2*Cos(aTgxy*Pi*xt*yt) - &
            aTgyz**2*Pi**2*Tgyz*zt**2*Sin(aTgyz*Pi*yt*zt))
      case(9)
      !tssrc!
        mms_source = -(ks*(-(aTsx**2*Pi**2*Tsx*Cos(aTsx*Pi*xt)) - &
              aTsxy**2*Pi**2*Tsxy*yt**2*Cos(aTsxy*Pi*xt*yt) - &
              aTszx**2*Pi**2*Tszx*zt**2*Cos(aTszx*Pi*xt*zt))) - &
         ks*(-(aTszx**2*Pi**2*Tszx*xt**2*Cos(aTszx*Pi*xt*zt)) - &
         aTsz**2*Pi**2*Tsz*Sin(aTsz*Pi*zt) - &
            aTsyz**2*Pi**2*Tsyz*yt**2*Sin(aTsyz*Pi*yt*zt)) - &
         ks*(-(aTsy**2*Pi**2*Tsy*Cos(aTsy*Pi*yt)) - &
         aTsxy**2*Pi**2*Tsxy*xt**2*Cos(aTsxy*Pi*xt*yt) - &
            aTsyz**2*Pi**2*Tsyz*zt**2*Sin(aTsyz*Pi*yt*zt)) + &
         Cps*ros*(es0 + esx*Cos(aesx*Pi*xt) + esy*Cos(aesy*Pi*yt) + &
         esxy*Cos(aesxy*Pi*xt*yt) + &
            eszx*Cos(aeszx*Pi*xt*zt) + esz*Sin(aesz*Pi*zt) + &
            esyz*Sin(aesyz*Pi*yt*zt))*&
          (vs0*Cos((Pi*(xt + yt + zt))/2.)**2*&
             (aTsyz*Pi*Tsyz*zt*Cos(aTsyz*Pi*yt*zt) - &
             aTsy*Pi*Tsy*Sin(aTsy*Pi*yt) - &
               aTsxy*Pi*Tsxy*xt*Sin(aTsxy*Pi*xt*yt)) + &
            ws0*(aTsz*Pi*Tsz*Cos(aTsz*Pi*zt) + &
            aTsyz*Pi*Tsyz*yt*Cos(aTsyz*Pi*yt*zt) - &
               aTszx*Pi*Tszx*xt*Sin(aTszx*Pi*xt*zt)) + &
            us0*(-(aTsx*Pi*Tsx*Sin(aTsx*Pi*xt)) - &
            aTsxy*Pi*Tsxy*yt*Sin(aTsxy*Pi*xt*yt) - &
               aTszx*Pi*Tszx*zt*Sin(aTszx*Pi*xt*zt))*&
               Sin((Pi*(xt + yt + zt))/2.)**2)
      case(10)
      !ropgsrc!
        mms_source = zero
      case(11)
      !ropssrc!
        mms_source = zero
      case(12)
      !thssrc!
        mms_source = -(ks*(-(aThsx**2*Pi**2*Thsx*Cos(aThsx*Pi*xt)) - &
              aThsxy**2*Pi**2*Thsxy*yt**2*Cos(aThsxy*Pi*xt*yt) - &
              aThszx**2*Pi**2*Thszx*zt**2*Cos(aThszx*Pi*xt*zt))) - &
         ks*(-(aThszx**2*Pi**2*Thszx*xt**2*Cos(aThszx*Pi*xt*zt)) - &
            aThsz**2*Pi**2*Thsz*Sin(aThsz*Pi*zt) - &
            aThsyz**2*Pi**2*Thsyz*yt**2*Sin(aThsyz*Pi*yt*zt))&
          - ks*(-(aThsy**2*Pi**2*Thsy*Cos(aThsy*Pi*yt)) - &
            aThsxy**2*Pi**2*Thsxy*xt**2*Cos(aThsxy*Pi*xt*yt) - &
            aThsyz**2*Pi**2*Thsyz*zt**2*Sin(aThsyz*Pi*yt*zt)) + &
         (3*ros*(es0 + esx*Cos(aesx*Pi*xt) + esy*Cos(aesy*Pi*yt) + &
         esxy*Cos(aesxy*Pi*xt*yt) + &
              eszx*Cos(aeszx*Pi*xt*zt) + esz*Sin(aesz*Pi*zt) + &
              esyz*Sin(aesyz*Pi*yt*zt))*&
            (vs0*Cos((Pi*(xt + yt + zt))/2.)**2*&
               (aThsyz*Pi*Thsyz*zt*Cos(aThsyz*Pi*yt*zt) - &
               aThsy*Pi*Thsy*Sin(aThsy*Pi*yt) - &
                 aThsxy*Pi*Thsxy*xt*Sin(aThsxy*Pi*xt*yt)) + &
              ws0*(aThsz*Pi*Thsz*Cos(aThsz*Pi*zt) + &
              aThsyz*Pi*Thsyz*yt*Cos(aThsyz*Pi*yt*zt) - &
                 aThszx*Pi*Thszx*xt*Sin(aThszx*Pi*xt*zt)) + &
              Pi*us0*Cos((Pi*(xt + yt + zt))/2.)*&
               (Ths0 + Thsx*Cos(aThsx*Pi*xt) + Thsy*Cos(aThsy*Pi*yt) + &
               Thsxy*Cos(aThsxy*Pi*xt*yt) + &
                 Thszx*Cos(aThszx*Pi*xt*zt) + Thsz*Sin(aThsz*Pi*zt) + &
                 Thsyz*Sin(aThsyz*Pi*yt*zt))*&
               Sin((Pi*(xt + yt + zt))/2.) - &
              Pi*vs0*Cos((Pi*(xt + yt + zt))/2.)*&
               (Ths0 + Thsx*Cos(aThsx*Pi*xt) + Thsy*Cos(aThsy*Pi*yt) + &
               Thsxy*Cos(aThsxy*Pi*xt*yt) + &
                 Thszx*Cos(aThszx*Pi*xt*zt) + Thsz*Sin(aThsz*Pi*zt) + &
                 Thsyz*Sin(aThsyz*Pi*yt*zt))*&
               Sin((Pi*(xt + yt + zt))/2.) + &
              us0*(-(aThsx*Pi*Thsx*Sin(aThsx*Pi*xt)) - &
              aThsxy*Pi*Thsxy*yt*Sin(aThsxy*Pi*xt*yt) - &
                 aThszx*Pi*Thszx*zt*Sin(aThszx*Pi*xt*zt))*&
                 Sin((Pi*(xt + yt + zt))/2.)**2))/2.
      end select


      END FUNCTION MMS_Source



      END MODULE mms
