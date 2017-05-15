!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Gw_Hw_Cw                                           C
!                                                                      C
!  Purpose: Calculate Gw, Hw, and Cw                                   C
!                                                                      C
!  Author: A. Srivastava & K. Agrawal, Princeton Univ. Date: 10-APR-98 C
!  Reviewer:                                           Date:           C
!                                                                      C
!  Modified: Sofiane Benyahia, Fluent Inc.             Date: 03-FEB-05 C
!  Purpose: Include conductivity defined by Simonin and Ahmadi         C
!           Also included Jenkins small frictional limit               C
!                                                                      C
!  Literature/Document References:                                     C
!     Johnson, P. C., and Jackson, R., Frictional-collisional          C
!        constitutive relations for granular materials, with           C
!        application to plane shearing, Journal of Fluid Mechanics,    C
!        1987, 176, pp. 67-93                                          C
!     Jenkins, J. T., and Louge, M. Y., On the flux of fluctuating     C
!        energy in a collisional grain flow at a flat frictional wall, C
!        Physics of Fluids, 1997, 9(10), pp. 2835-2840                 C
!                                                                      C
!     See calc_mu_s.f for references on frictional theory models       C
!     See calc_mu_s.f for references on Ahmadi and Simonin models      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_Gw_Hw_Cw(g0, EPS, EPG, ep_star_avg, &
            g0EP_avg, TH, Mu_g_avg, RO_g_avg, ROs_avg, &
            DP_avg, K_12_avg, Tau_12_avg, Tau_1_avg, VREL, VSLIP,&
            DEL_U, S_S, S_dd, VEL, W_VEL, M, gw, hw, cw)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE constant
      USE physprop
      USE fldvar
      USE bc
      USE run
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Radial distribution function of solids phase M with each
! other solids phase
      DOUBLE PRECISION, INTENT(IN) :: g0
! Average solids volume fraction of each solids phase
      DOUBLE PRECISION, INTENT(IN) :: EPS
! Average solids and gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPG, ep_star_avg
! Sum of eps*G_0
      DOUBLE PRECISION, INTENT(IN) :: g0EP_avg
! Average theta_m
      DOUBLE PRECISION, INTENT(INOUT) :: TH
! Average gas viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mu_g_avg
! Average gas density
      DOUBLE PRECISION, INTENT(IN) :: RO_g_avg
! Average solids density
      DOUBLE PRECISION, INTENT(IN) :: ROS_avg
! Average particle diameter of each solids phase
      DOUBLE PRECISION, INTENT(IN) :: DP_avg
! Average cross-correlation K_12 and lagrangian integral time-scale
      DOUBLE PRECISION, INTENT(IN) :: K_12_avg, Tau_12_avg, Tau_1_avg
! Magnitude of slip velocity between two phases
      DOUBLE PRECISION, INTENT(IN) :: VREL
! Slip velocity between wall and particles
      DOUBLE PRECISION, INTENT(IN) :: VSLIP
! Relevant solids velocity at wall
      DOUBLE PRECISION, INTENT(IN) :: VEL
! Relevant wall velocity
      DOUBLE PRECISION, INTENT(IN) :: W_VEL
! del.u
      DOUBLE PRECISION, INTENT(IN) :: DEL_U
! S:S
      DOUBLE PRECISION, INTENT(IN) :: S_S
! S_dd (d can be x,y or z)
      DOUBLE PRECISION, INTENT(IN) :: S_dd
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Wall momentum coefficients:
! 1st term on LHS
      DOUBLE PRECISION, INTENT(INOUT) :: Gw
! 2nd term on LHS
      DOUBLE PRECISION, INTENT(INOUT) :: Hw
! all terms appearing on RHS
      DOUBLE PRECISION, INTENT(INOUT) :: Cw
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Part of wall momentum coefficient in 1st term of LHS
      DOUBLE PRECISION :: Mu_s
! Term appearing in wall momenutm coefficient of 2nd term LHS
      DOUBLE PRECISION :: F_2
! Frictional Pressure Pf
      DOUBLE PRECISION :: Pf
! Critical pressure Pc
      DOUBLE PRECISION :: Pc
! Chi (appears in frictional boundary condition)
      DOUBLE PRECISION :: Chi, N_Pff
! Viscosity
      DOUBLE PRECISION :: Mu
! Bulk viscosity
      DOUBLE PRECISION :: Mu_b
! Viscosity corrected for interstitial fluid effects
      DOUBLE PRECISION :: Mu_star
! Reynolds number based on slip velocity
      DOUBLE PRECISION :: Re_g
! Friction Factor in drag coefficient
      DOUBLE PRECISION :: C_d
! Drag Coefficient
      DOUBLE PRECISION :: Beta, DgA
! Square root of S:S or the form suggested by Savage
      DOUBLE PRECISION :: ZETA
! Constants in Simonin model
      DOUBLE PRECISION :: Sigma_c, Tau_2_c, Tau_12_st, Nu_t
      DOUBLE PRECISION :: Tau_2, zeta_c_2, MU_2_T_Kin, Mu_2_Col
      DOUBLE PRECISION :: Tmp_Ahmadi_Const
! Other local terms
      DOUBLE PRECISION :: dpc_dphi
! Local variable for rdf
      DOUBLE PRECISION :: g0_loc
!-----------------------------------------------

! This is done here similar to bc_theta to avoid small negative values of
! Theta coming most probably from linear solver
      IF(TH .LE. ZERO)THEN
        TH = 1D-8
        IF (myPE.eq.PE_IO) THEN
           WRITE(*,*) &
              'Warning: Negative granular temp at wall set to 1e-8'
!          CALL WRITE_ERROR('THETA_HW_CW', LINE, 1)
        ENDIF
      ENDIF

      g0_loc = g0

! modify F_2 if Jenkins BC is used (sof)
      IF(JENKINS) THEN

        IF (VSLIP == ZERO) THEN
! if solids velocity field is initialized to zero, use free slip bc
           F_2 = zero

        ELSE
           IF(AHMADI) THEN
! Ahmadi model uses different solids pressure model
! the coefficient mu in Jenkins paper is defined as tan_Phi_w, that's how
! I understand it from soil mechanic papers, i.e., G.I. Tardos, powder
! Tech. 92 (1997), 61-74. See his equation (1). Define Phi_w in mfix.dat!
! here F_2 divided by VSLIP to use the same bc as Johnson&Jackson

              F_2 = tan_Phi_w*ROS_avg*EPS* &
                 ((ONE + 4.0D0*g0EP_avg) + HALF*(ONE -C_e*C_e))*TH/VSLIP

           ELSE
! Simonin or granular models use same solids pressure
              F_2 = tan_Phi_w*ROS_avg*EPS*(1d0+ 4.D0 * Eta *g0EP_avg)*TH/VSLIP
           ENDIF   ! end if(Ahmadi)

        ENDIF ! endif(vslip==0)

      ELSE    ! if(.not.jenkins)

         F_2 = (PHIP*DSQRT(3d0*TH)*Pi*ROS_avg*EPS*g0_loc)&
              /(6d0*(ONE-ep_star_avg))

      ENDIF   ! end if(Jenkins)/else


      Mu = (5d0*DSQRT(Pi*TH)*Dp_avg*ROS_avg)/96d0
      Mu_b = (256d0*Mu*EPS*g0EP_avg)/(5d0*Pi)

! This is from Wen-Yu correlation, you can put here your own single particle drag
      Re_g = EPG*RO_g_avg*Dp_avg*VREL/Mu_g_avg
      IF (Re_g.lt.1000d0) THEN
         C_d = (24.d0/(Re_g+SMALL_NUMBER))*(ONE + 0.15d0 * Re_g**0.687d0)
      ELSE
         C_d = 0.44d0
      ENDIF
      DgA = 0.75d0*C_d*Ro_g_avg*EPG*VREL/(Dp_avg*EPG**(2.65d0))
      IF(VREL == ZERO) DgA = LARGE_NUMBER
      Beta = SWITCH*EPS*DgA

! SWITCH enables us to turn on/off the modification to the
! particulate phase viscosity. If we want to simulate gas-particle
! flow then SWITCH=1 to incorporate the effect of drag on the
! particle viscosity. If we want to simulate granular flow
! without the effects of an interstitial gas, SWITCH=0.
      IF(SWITCH == ZERO .OR. Ro_g_avg == ZERO)THEN
         Mu_star = Mu
      ELSEIF(TH .LT. SMALL_NUMBER)THEN
         MU_star = ZERO
      ELSE
         Mu_star = ROS_avg*EPS* g0_loc*TH* Mu/ &
            (ROS_avg*g0EP_avg*TH + 2.0d0*SWITCH*DgA/ROS_avg* Mu)
      ENDIF

      Mu_s = ((2d0+ALPHA)/3d0)*((Mu_star/(Eta*(2d0-Eta)*&
                   g0_loc))*(ONE+1.6d0*Eta*g0EP_avg&
                   )*(ONE+1.6d0*Eta*(3d0*Eta-2d0)*&
                   g0EP_avg)+(0.6d0*Mu_b*Eta))

! particle relaxation time
      Tau_12_st = ROS_avg/(DgA+small_number)

      IF(SIMONIN) THEN !see calc_mu_s for explanation of these definitions

         Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
         Tau_2_c = DP_avg/(6.d0*EPS*g0_loc*DSQRT(16.d0*(TH+Small_number)/PI))
         zeta_c_2= 2.D0/5.D0*(ONE+ C_e)*(3.d0*C_e-ONE)
         Nu_t =  Tau_12_avg/Tau_12_st
         Tau_2 = ONE/(2.D0/Tau_12_st+Sigma_c/Tau_2_c)
         MU_2_T_Kin = (2.0D0/3.0D0*K_12_avg*Nu_t + TH * &
                     (ONE+ zeta_c_2*EPS*g0_loc))*Tau_2
         Mu_2_Col = 8.D0/5.D0*EPS*g0_loc*Eta* (MU_2_T_Kin+ &
                   Dp_avg*DSQRT(TH/PI))
         Mu_s = EPS*ROS_avg*(MU_2_T_Kin + Mu_2_Col)
      ELSEIF(AHMADI) THEN
         IF(EPS < (ONE-ep_star_avg)) THEN
            Tmp_Ahmadi_Const = &
            ONE/(ONE+ Tau_1_avg/Tau_12_st * (ONE-EPS/(ONE-ep_star_avg))**3)
         ELSE
            Tmp_Ahmadi_Const = ONE
         ENDIF
         Mu_s = Tmp_Ahmadi_Const &
             *0.1045D0*(ONE/g0_loc+3.2D0*EPS+12.1824D0*g0_loc*EPS*EPS)  &
             *Dp_avg*ROS_avg* DSQRT(TH)
      ENDIF

! Calculating frictional terms
      IF ((ONE-EPG)<= EPS_f_min) THEN
         Pf = ZERO
         Chi = ZERO
         ZETA = 1d0
      ELSE
         IF (SAVAGE.EQ.1) THEN    !form of Savage
            ZETA = ((48d0*Eta*(1d0-Eta)*ROS_avg*EPS*EPS*g0_loc*&
                    (TH**1.5d0))/&
                    (SQRT_Pi*Dp_avg*2d0*Mu_s))**0.5d0
         ELSEIF (SAVAGE.EQ.0)  THEN !S:S form
            ZETA = DSQRT(S_S)
         ELSE
            ZETA = DSQRT(S_S + (TH/(Dp_avg*Dp_avg)))
         ENDIF

         IF (EPG < ep_star_avg) THEN
            dpc_dphi = (to_SI*Fr)*((delta**5)*(2d0*(ONE-ep_star_avg-delta) - &
               2d0*eps_f_min)+((ONE-ep_star_avg-delta)-eps_f_min)&
               *(5*delta**4))/(delta**10)

            Pc = (to_SI*Fr)*(((ONE-ep_star_avg-delta) - EPS_f_min)**N_Pc)/(delta**D_Pc)
!            Pc=  1d25*(((ONE-EPG)-(ONE-ep_star_avg))**10d0)  ! this is old Pc
         ELSE
            Pc = Fr*(((ONE-EPG) - EPS_f_min)**N_Pc)/ &
               (((ONE-ep_star_avg) - (ONE-EPG) + SMALL_NUMBER)**D_Pc)
         ENDIF

         IF (DEL_U .GE. ZERO) THEN
            N_Pff = DSQRT(3d0)/(2d0*Sin_Phi) !dilatation
         ELSE
            N_Pff = N_Pf !compaction
         ENDIF

         IF ((DEL_U/(ZETA*N_Pff*DSQRT(2d0)*Sin_Phi)) .GT. 1d0) THEN
            Pf = ZERO
         ELSEIF( DEL_U == ZERO ) THEN
            Pf = Pc
         ELSE
            Pf = Pc*(1d0 - (DEL_U/(ZETA*N_Pff*DSQRT(2d0)*Sin_Phi)))**&
               (N_Pff-1d0)
         ENDIF

         Chi =  DSQRT(2d0)*Pf*Sin_Phi*(N_Pff - (N_Pff-1d0)*&
            ((Pf/Pc)**(1d0/(N_Pff-1d0))))

         IF (Chi< ZERO) THEN
            Pf = Pc*((N_Pff/(N_Pff-1d0))**(N_Pff-1d0))
            Chi = ZERO
         ENDIF

! by writing Pf & Chi in the following form, we ensure distribution
! of stresses amoung all solids phases (sof, Oct 24 2005)

         Pf = Pf * EPS/(ONE-EPG)
         Chi = Chi * EPS/(ONE-EPG)

      ENDIF

! Calculating gw, hw, cw

      Gw = ONE   ! we write this in the exact same form as non-frictional JJ BC

      IF(VSLIP == ZERO .OR. ZETA == ZERO) THEN
        Hw = ZERO
      ELSE
        Hw = F_2 + Pf*tan_Phi_w - Chi*S_dd*tan_Phi_w/ZETA
        Hw = Hw / (MU_s + Chi/(2d0*ZETA))  ! this is because Gw is set to one.
      ENDIF
      Cw = hw * W_VEL

      RETURN
      END SUBROUTINE CALC_Gw_Hw_Cw
