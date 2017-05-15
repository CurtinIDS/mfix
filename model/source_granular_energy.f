!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_GRANULAR_ENERGY                                  C
!  Purpose: Calculate the source terms in the granular energy equation C
!  Note: The center coefficient and source vector are negative.        C
!  The off-diagonal coefficients are positive.                         C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University        Date: 04-FEB-98  C
!  Reviewer: M. Syamlal                               Date:            C
!                                                                      C
!  Revision Number:1                                                   C
!  Purpose: Add Simonin and Ahmadi models                              C
!  Author: Sofiane Benyahia, Fluent Inc.               Date: 02-01-05  C
!                                                                      C
!  Literature/Document References:                                     C
!     Lun, C.K.K., S.B. Savage, D.J. Jeffrey, and N. Chepurniy,        C
!        Kinetic theories for granular flow - inelastic particles in   C
!        Couette-flow and slightly inelastic particles in a general    C
!        flow field. Journal of Fluid Mechanics, 1984. 140(MAR):       C
!        p. 223-256                                                    C
!     Gidaspow, D., Multiphase flow and fluidziation, 1994, Academic   C
!        Press Inc., California, Chapter 9.                            C
!     Koch, D. L., and Sangani, A. S., Particle pressure and marginal  C
!        stability limits for a homogeneous monodisperse gas-fluidized C
!        bed: kinetic theory and numerical simulations, Journal of     C
!        Fluid Mechanics, 1999, 400, 229-263.                          C
!                                                                      C
!     Simonin, O., 1996. Combustion and turbulence in two-phase flows, C
!        Von Karman institute for fluid dynamics, lecture series,      C
!        1996-02                                                       C
!     Balzer, G., Simonin, O., Boelle, A., and Lavieville, J., 1996,   C
!        A unifying modelling approach for the numerical prediction    C
!        of dilute and dense gas-solid two phase flow. CFB5, 5th int.  C
!        conf. on circulating fluidized beds, Beijing, China.          C
!     Cao, J. and Ahmadi, G., 1995, Gas-particle two-phase turbulent   C
!        flow in a vertical duct. Int. J. Multiphase Flow, vol. 21,    C
!        No. 6, pp. 1203-1228.                                         C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_GRANULAR_ENERGY(SOURCELHS, &
                    SOURCERHS, IJK, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE drag
      USE fldvar
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE kintheory
      USE mms
      USE parallel
      USE param
      USE param1
      USE physprop
      USE rdf
      USE run
      USE solids_pressure
      USE toleranc
      USE trace
      USE turb
      USE visc_g
      USE visc_s

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Source terms to be kept on rhs (source vector)
      DOUBLE PRECISION, INTENT(INOUT) :: sourcerhs
! Source terms to be kept on lhs (center coefficient)
      DOUBLE PRECISION, INTENT(INOUT) :: sourcelhs
! Solid phase index
      INTEGER, INTENT(IN) :: M
! Indices
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IM, JM, KM, IMJK, IJMK, IJKM, &
                 IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
! Solid phase index
      INTEGER :: MM
! Particle relaxation time
      DOUBLE PRECISION :: Tau_12_st
! Sum of eps*G_0
      DOUBLE PRECISION :: SUM_EpsGo
! Slip velocity
      DOUBLE PRECISION :: VSLIP
! coefficients in heat flux term
      DOUBLE PRECISION :: Knu_e, Knu_w, Knu_n, Knu_s, &
                          Knu_t, Knu_b
! Source terms to be kept on rhs
      DOUBLE PRECISION :: S1_rhs, S2_rhs, S3_rhs, &
                          S5a_rhs, S5b_rhs, S5c_rhs, S5_rhs, &
                          S7_rhs
! Source terms to be kept on lhs
      DOUBLE PRECISION :: S1_lhs, S3_lhs, S4_lhs, S6_lhs
!-----------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IM = Im1(I)
      JM = Jm1(J)
      KM = Km1(K)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IJKE = EAST_OF(IJK)
      IJKW = WEST_OF(IJK)
      IJKN = NORTH_OF(IJK)
      IJKS = SOUTH_OF(IJK)
      IJKT = TOP_OF(IJK)
      IJKB = BOTTOM_OF(IJK)

! initialize summation variables
      S1_rhs = ZERO
      S1_lhs = ZERO
      S2_rhs = ZERO
      S3_rhs = ZERO
      S3_lhs = ZERO
      S4_lhs = ZERO
      S6_lhs = ZERO
      S7_rhs = ZERO

! Changes needed for multitype particles, sof June 16 2005
! Sum of eps*G_0 is used instead of Eps*G_0
      SUM_EpsGo = ZERO
      DO MM = 1, MMAX
         SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,MM)*G_0(IJK,MM,MM)
      ENDDO

! Production by shear: (S:grad(vi))
! Pi_s*tr(Di)     (Lun et al. 1984)
      S1_rhs = P_S_C(IJK,M)*ZMAX(( -TRD_S_C(IJK,M) ))
      S1_lhs = P_S_C(IJK,M)*ZMAX(( TRD_S_C(IJK,M) ))

! Production by shear: (S:grad(vi))
! Mu_s*tr(Di^2)     (Lun et al. 1984)
      S2_rhs = 2.d0*MU_S_C(IJK,M)*TRD_S2(IJK,M)

! Production by shear: (S:grad(vi))
! Lambda_s*tr(Di)^2     (Lun et al. 1984)
      S3_rhs = (TRD_S_C(IJK,M)**2)*ZMAX( LAMBDA_S_C(IJK,M) )
      S3_lhs = (TRD_S_C(IJK,M)**2)*ZMAX( -LAMBDA_S_C(IJK,M) )

! Energy dissipation by collisions
      S4_lhs = (48.d0/DSQRT(PI))*ETA*(ONE-ETA)*ROP_S(IJK,M)*&
         SUM_EpsGo*DSQRT(THETA_M(IJK,M))/D_P(IJK,M)

! Need to revisit this part about granular energy MMS source terms.
      IF (USE_MMS) S4_lhs = ZERO

! Energy dissipation by viscous dampening
! Gidaspow (1994)  : addition due to role of interstitial fluid
      IF(KT_TYPE_ENUM == SIMONIN_1996 .OR. &
         KT_TYPE_ENUM == AHMADI_1995) THEN
         S6_lhs = 3.d0 * F_GS(IJK,M)
      ELSE IF(SWITCH > ZERO .AND. RO_g0 /= ZERO) THEN
         S6_lhs = SWITCH * 3.d0 * F_GS(IJK,M)
      ENDIF

! Energy production due to gas-particle slip
! Koch & Sangani (1999) : addition due to role of interstitial fluid
      IF(KT_TYPE_ENUM == SIMONIN_1996) THEN
         S7_rhs = F_GS(IJK,M)*K_12(IJK)
      ELSEIF(KT_TYPE_ENUM == AHMADI_1995) THEN
! note specific reference to F_GS of solids phase 1!
         IF(Ep_s(IJK,M) > DIL_EP_S .AND. F_GS(IJK,1) > small_number) THEN
            Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,1)
            S7_rhs = 2.d0*F_GS(IJK,M) * (ONE/(ONE+Tau_12_st/  &
               (Tau_1(ijk)+small_number)))*K_Turb_G(IJK)
         ELSE
            S7_rhs = ZERO
         ENDIF
      ELSEIF(SWITCH > ZERO .AND. RO_g0 /= ZERO) THEN
         VSLIP = DSQRT( (U_S(IJK,M)-U_G(IJK))**2 + &
            (V_S(IJK,M)-V_G(IJK))**2 + (W_S(IJK,M)-W_G(IJK))**2 )
         S7_rhs = SWITCH*81.d0*EP_S(IJK,M)*(MU_G(IJK)*VSLIP)**2 /&
            (G_0(IJK,M,M)*D_P(IJK,M)**3 * RO_S(IJK,M)*&
            DSQRT(PI)*DSQRT( THETA_M(IJK,M)+SMALL_NUMBER ) )
! Need to revisit this part about granular energy MMS source terms.
         IF (USE_MMS) S7_rhs = ZERO

      ENDIF


! The following lines are essentially commented out since Kphi_s has
! been set to zero in subroutine calc_mu_s.  To activate the feature
! activate the following lines and the lines in calc_mu_s.
! Part of Heat Flux: div (q)
! Kphi_s*grad(eps)     (Lun et al. 1984)
      IF (.FALSE.) THEN
         Knu_e = AVG_X_H(Kphi_s(IJK,M), Kphi_s(IJKE,M),I)
         Knu_w = AVG_X_H(Kphi_s(IJKW,M),Kphi_s(IJK,M),IM)
         Knu_n = AVG_Y_H(Kphi_s(IJK,M), Kphi_s(IJKN,M),J)
         Knu_s = AVG_Y_H(Kphi_s(IJKS,M),Kphi_s(IJK,M),JM)
         Knu_t = AVG_Z_H(Kphi_s(IJK,M), Kphi_s(IJKT,M),K)
         Knu_b = AVG_Z_H(Kphi_s(IJKB,M),Kphi_s(IJK,M),KM)

         S5a_rhs = Knu_e*(EP_s(IJKE,M)-EP_s(IJK,M))*&
            oDX_E(I)*AYZ(IJK) - &
            Knu_w*(EP_s(IJK,M)-EP_s(IJKW,M))*&
            oDX_E(IM)*AYZ(IMJK)
         S5b_rhs = Knu_n*(EP_s(IJKN,M)-EP_s(IJK,M))*&
            oDY_N(J)*AXZ(IJK) - &
            Knu_s*(EP_s(IJK,M)-EP_s(IJKS,M))*&
            oDY_N(JM)*AXZ(IJMK)
         S5c_rhs = Knu_t*(EP_s(IJKT,M)-EP_s(IJK,M))*&
            oX(I)*oDZ_T(K)*AXY(IJK) - &
            Knu_b*(EP_s(IJK,M)-EP_s(IJKB,M))*&
            oX(I)*oDZ_T(KM)*AXY(IJKM)
         S5_rhs = S5a_rhs + S5b_rhs + S5c_rhs
      ENDIF
      S5_rhs = ZERO


      SOURCERHS = (S1_rhs + S2_rhs + S3_rhs + S7_rhs)*VOL(IJK) + &
           S5_rhs

      SOURCELHS = ( ((S1_lhs + S3_lhs)/(THETA_M(IJK,M)+SMALL_NUMBER)) + &
          S4_lhs + S6_lhs) * VOL(IJK)

      RETURN
      END SUBROUTINE SOURCE_GRANULAR_ENERGY




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_IA_GRANULAR_ENERGY                               C
!  Purpose: Calculate the source terms of the Mth solids phase         C
!           granular energy equation                                   C
!                                                                      C
!  Author: Janine E. Galvin, Univeristy of Colorado                    C
!                                                                      C
!  Literature/Document References:                                     C
!    Iddir, Y.H., Modeling of the multiphase mixture of particles      C
!       using the kinetic theory approach, PhD Thesis, Illinois        C
!       Institute of Technology, Chicago, Illinois, 2004               C
!    Iddir, Y.H., & H. Arastoopour, Modeling of multitype particle     C
!       flow using the kinetic theory approach, AIChE J., Vol 51,      C
!       No 6, June 2005                                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_GRANULAR_ENERGY_IA(SOURCELHS, &
                    SOURCERHS, IJK, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE run
      USE drag
      USE geometry
      USE fldvar
      USE visc_g
      USE visc_s
      USE trace
      USE turb
      USE indices
      USE constant
      USE toleranc
      USE residual
      use kintheory
      USE compar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Source terms to be kept on rhs (source vector)
      DOUBLE PRECISION, INTENT(INOUT) :: sourcerhs
! Source terms to be kept on lhs (center coefficient)
      DOUBLE PRECISION, INTENT(INOUT) :: sourcelhs
! Solid phase index
      INTEGER, INTENT(IN) :: M
! Indices
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IM, JM, KM, IMJK, IJMK, IJKM,&
                 IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
! phase index
      INTEGER :: L, LM
! velocities
      DOUBLE PRECISION :: UsM_e, UsM_w, VsM_n, VsM_s, WsM_t, WsM_b,&
                          UsL_e, UsL_w, VsL_n, VsL_s, WsL_t, WsL_b,&
                          UsM_p, VsM_p, WsM_P, UsL_p, VsL_p, WsL_p
! number densities
      DOUBLE PRECISION :: NU_PM_E, NU_PM_W, NU_PM_N, NU_PM_S, NU_PM_T,&
                          NU_PM_B, NU_PM_p,&
                          NU_PL_E, NU_PL_W, NU_PL_N, NU_PL_S, NU_PL_T,&
                          NU_PL_B, NU_PL_p
! temperature of species L
      DOUBLE PRECISION :: T_PL_E, T_PL_W, T_PL_N, T_PL_S, T_PL_T,&
                          T_PL_B, T_PL_p
! particle characteristics
      DOUBLE PRECISION :: M_PM, M_PL, D_PM, D_PL
! coefficients in heat flux term
      DOUBLE PRECISION :: Knu_sL_e, Knu_sL_w, Knu_sL_n, Knu_sL_s, Knu_sL_t,&
                          Knu_sL_b, Knu_sM_e, Knu_sM_w, Knu_sM_n, Knu_sM_s,&
                          Knu_sM_t, Knu_sM_b,&
                          Kvel_s_e, Kvel_s_w, Kvel_s_n, Kvel_s_s, Kvel_s_t,&
                          Kvel_s_b,&
                          Kth_sL_e, Kth_sL_w, Kth_sL_n, Kth_sL_s, Kth_sL_t,&
                          Kth_sL_b
! Source terms to be kept on rhs
      DOUBLE PRECISION :: S10_rhs, S15_rhs, S16_rhs,&
                          S11_sum_rhs, S12_sum_rhs,&
                          S13_sum_rhs, S17_sum_rhs, S18_sum_rhs
! Source terms
      DOUBLE PRECISION :: S14a_sum, S14b_sum, S14c_sum, &
                          S9_sum, s21a_sum, s21b_sum, s21c_sum
! Source terms to be kept on lhs
      DOUBLE PRECISION :: S10_lhs, S16_lhs,&
                          S11_sum_lhs, S12_sum_lhs, S13_sum_lhs,&
                          S17_sum_lhs, S18_sum_lhs, S20_sum_lhs
!-----------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IM = Im1(I)
      JM = Jm1(J)
      KM = Km1(K)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IJKE = EAST_OF(IJK)
      IJKW = WEST_OF(IJK)
      IJKN = NORTH_OF(IJK)
      IJKS = SOUTH_OF(IJK)
      IJKT = TOP_OF(IJK)
      IJKB = BOTTOM_OF(IJK)

! initialize summation variables
      S9_sum = ZERO
      S14a_sum = ZERO
      S14b_sum = ZERO
      S14c_sum = ZERO
      S21a_sum = ZERO
      S21b_sum = ZERO
      S21c_sum = ZERO

      S11_sum_rhs = ZERO
      S12_sum_rhs = ZERO
      S13_sum_rhs = ZERO
      S17_sum_rhs = ZERO
      S18_sum_rhs = ZERO

      S11_sum_lhs = ZERO
      S12_sum_lhs = ZERO
      S13_sum_lhs = ZERO
      S17_sum_lhs = ZERO
      S18_sum_lhs = ZERO
      S20_sum_lhs = ZERO

      UsM_e = U_S(IJK,M)
      UsM_w = U_S(IMJK,M)
      VsM_n = V_S(IJK,M)
      VsM_s = V_S(IJMK,M)
      WsM_t = W_S(IJK,M)
      WsM_b = W_S(IJKM,M)
      UsM_p = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
      VsM_p = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M) )
      WsM_p = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M) )

      D_PM = D_P(IJK,M)
      M_PM = (Pi/6.d0)*D_PM**3 * RO_S(IJK,M)
      NU_PM_p = ROP_S(IJK,M)/M_PM
      NU_PM_E = ROP_S(IJKE,M)/M_PM
      NU_PM_W = ROP_S(IJKW,M)/M_PM
      NU_PM_N = ROP_S(IJKN,M)/M_PM
      NU_PM_S = ROP_S(IJKS,M)/M_PM
      NU_PM_T = ROP_S(IJKT,M)/M_PM
      NU_PM_B = ROP_S(IJKB,M)/M_PM

! Production by shear: (S:grad(vi))
! Pi_s*tr(Di)
      S10_lhs = P_S_C(IJK,M) * ZMAX(TRD_S_C(IJK,M))
      S10_rhs = P_S_C(IJK,M) * ZMAX(-TRD_S_C(IJK,M))


! Production by shear: (S:grad(vi))
! Mu_s*tr(Di^2)
      S15_rhs = 2.d0*Mu_s_c(IJK,M)*TRD_S2(IJK,M)

! Production by shear: (S:grad(vi))
! Lambda_s*tr(Di)^2
      S16_lhs = (TRD_S_C(IJK,M)**2)*ZMAX( -LAMBDA_s_c(IJK,M) )
      S16_rhs = (TRD_S_C(IJK,M)**2)*ZMAX(  LAMBDA_s_C(IJK,M) )

      DO L = 1, MMAX
          D_PL = D_P(IJK,L)
          M_PL = (Pi/6.d0)*D_PL**3 * RO_S(IJK,L)
          NU_PL_p = ROP_S(IJK,L)/M_PL
          NU_PL_E = ROP_S(IJKE,L)/M_PL
          NU_PL_W = ROP_S(IJKW,L)/M_PL
          NU_PL_N = ROP_S(IJKN,L)/M_PL
          NU_PL_S = ROP_S(IJKS,L)/M_PL
          NU_PL_T = ROP_S(IJKT,L)/M_PL
          NU_PL_B = ROP_S(IJKB,L)/M_PL

          T_PL_p = Theta_m(IJK,L)
          T_PL_E = Theta_m(IJKE,L)
          T_PL_W = Theta_m(IJKW,L)
          T_PL_N = Theta_m(IJKN,L)
          T_PL_S = Theta_m(IJKS,L)
          T_PL_T = Theta_m(IJKT,L)
          T_PL_B = Theta_m(IJKB,L)

          UsL_e = U_S(IJK,L)
          UsL_w = U_S(IMJK,L)
          VsL_n = V_S(IJK,L)
          VsL_s = V_S(IJMK,L)
          WsL_t = W_S(IJK,L)
          WsL_b = W_S(IJKM,L)
          UsL_p = AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I)
          VsL_p = AVG_Y_N(V_S(IJMK,L),V_S(IJK,L) )
          WsL_p = AVG_Z_T(W_S(IJKM,L),W_S(IJK,L) )

! Energy dissipation by collisions: Sum(Nip)
! SUM( EDT_s_ip )
          S20_sum_lhs = S20_sum_lhs + EDT_s_ip(IJK,M,L)

! Energy dissipation by collisions: SUM(Nip)
! SUM( EDvel_sL_ip* div(vp) ) !Modified by sof to include trace of V_s_L
          S11_sum_lhs = S11_sum_lhs + ZMAX(-EDvel_sL_ip(IJK,M,L)* &
               TRD_S_C(IJK,L) ) * VOL(IJK)
          S11_sum_rhs = S11_sum_rhs + ZMAX( EDvel_sL_ip(IJK,M,L)* &
               TRD_S_C(IJK,L) ) * VOL(IJK)

! Energy dissipation by collisions: Sum(Nip)
! SUM( EDvel_sM_ip* div(vi) ) !Modified by sof to include trace of V_s_M
          S12_sum_lhs = S12_sum_lhs + ZMAX(-EDvel_sM_ip(IJK,M,L)*&
               TRD_S_C(IJK,M) ) * VOL(IJK)
          S12_sum_rhs = S12_sum_rhs + ZMAX( EDvel_sM_ip(IJK,M,L)*&
               TRD_S_C(IJK,M) ) * VOL(IJK)

          IF (M .NE. L) THEN
               LM = FUNLM(L,M)

! Production by shear: (S:grad(vi))
! SUM(2*Mu_sL_ip*tr(Dk*Di) )
               S17_sum_lhs = S17_sum_lhs + 2.d0*MU_sL_ip(IJK,M,L)*&
                    ZMAX( - TRD_s2_ip(IJK,M,L) )
               S17_sum_rhs = S17_sum_rhs + 2.d0*MU_sL_ip(IJK,M,L)*&
                    ZMAX( TRD_s2_ip(IJK,M,L) )

! These two terms can be treated explicitly here by uncommenting the following
! two lines. They are currently treated by PEA algorithm in solve_granular_energy
!
!               S13_sum_lhs = S13_sum_lhs + ED_ss_ip(IJK,LM)*Theta_m(IJK,M)
!               S13_sum_rhs = S13_sum_rhs + ED_ss_ip(IJK,LM)*Theta_m(IJK,L)
!
!
! Production by shear: (S:grad(vi))
! SUM( (Xi_sL_ip-(2/3)*Mu_sL_ip)*tr(Dk)tr(Di) )
               S18_sum_lhs = S18_sum_lhs + ZMAX(-1.d0* &
                    (Xi_sL_ip(IJK,M,L)-(2.d0/3.d0)*Mu_sL_ip(IJK,M,L))*&
                    TRD_S_C(IJK,M)*TRD_S_C(IJK,L) )
               S18_sum_rhs = S18_sum_rhs + ZMAX(&
                    (Xi_sL_ip(IJK,M,L)-(2.d0/3.d0)*Mu_sL_ip(IJK,M,L))*&
                    TRD_S_C(IJK,M)*TRD_S_C(IJK,L) )

! Part of Heat Flux: div (q)
! Kth_sL_ip*[grad(Tp)]
! Note for L=M S21 terms cancel with similar term arising from grad(Ti)
               Kth_sL_e = AVG_X_S(Kth_sL_ip(IJK,M,L), Kth_sL_ip(IJKE,M,L),I)
               Kth_sL_w = AVG_X_S(Kth_sL_ip(IJKW,M,L),Kth_sL_ip(IJK,M,L), IM)
               Kth_sL_n = AVG_Y_S(Kth_sL_ip(IJK,M,L), Kth_sL_ip(IJKN,M,L),J)
               Kth_sL_s = AVG_Y_S(Kth_sL_ip(IJKS,M,L),Kth_sL_ip(IJK,M,L), JM)
               Kth_sL_t = AVG_Z_S(Kth_sL_ip(IJK,M,L), Kth_sL_ip(IJKT,M,L),K)
               Kth_sL_b = AVG_Z_S(Kth_sL_ip(IJKB,M,L),Kth_sL_ip(IJK,M,L), KM)

               S21a_sum = S21a_sum + ( (Kth_sL_e*(T_PL_E-T_PL_p) )*&
                    ODX_E(I)*AYZ(IJK) - (Kth_sL_w*(T_PL_p-T_PL_W) )*ODX_E(IM)*&
                    AYZ(IMJK) )
               S21b_sum = S21b_sum + ( (Kth_sL_n*(T_PL_N-T_PL_p) )*&
                    ODY_N(J)*AXZ(IJK) - (Kth_sL_s*(T_PL_p-T_PL_S) )*ODY_N(JM)*&
                    AXZ(IJMK) )
               S21c_sum = S21c_sum + ( (Kth_sL_t*(T_PL_T-T_PL_p) )*&
                    ODZ_T(K)*OX(I)*AXY(IJK) - (Kth_sL_b*(T_PL_p-T_PL_B) )*&
                    ODZ_T(KM)*OX(I)*AXY(IJKM) )

! Part of Heat Flux: div (q)
! Knu_s_ip*[ni*grad(np)-np*grad(ni)]
! Note S14 terms should evaluate to zero for particles from the same phase
               Knu_sL_e = AVG_X_S(Knu_sL_ip(IJK,M,L), Knu_sL_ip(IJKE,M,L),I)
               Knu_sM_e = AVG_X_S(Knu_sM_ip(IJK,M,L), Knu_sM_ip(IJKE,M,L),I)
               Knu_sL_w = AVG_X_S(Knu_sL_ip(IJKW,M,L),Knu_sL_ip(IJK,M,L), IM)
               Knu_sM_w = AVG_X_S(Knu_sM_ip(IJKW,M,L),Knu_sM_ip(IJK,M,L), IM)
               Knu_sL_n = AVG_Y_S(Knu_sL_ip(IJK,M,L), Knu_sL_ip(IJKN,M,L),J)
               Knu_sM_n = AVG_Y_S(Knu_sM_ip(IJK,M,L), Knu_sM_ip(IJKN,M,L),J)
               Knu_sL_s = AVG_Y_S(Knu_sL_ip(IJKS,M,L),Knu_sL_ip(IJK,M,L), JM)
               Knu_sM_s = AVG_Y_S(Knu_sM_ip(IJKS,M,L),Knu_sM_ip(IJK,M,L), JM)
               Knu_sL_t = AVG_Z_S(Knu_sL_ip(IJK,M,L), Knu_sL_ip(IJKT,M,L),K)
               Knu_sM_t = AVG_Z_S(Knu_sM_ip(IJK,M,L), Knu_sM_ip(IJKT,M,L),K)
               Knu_sL_b = AVG_Z_S(Knu_sL_ip(IJKB,M,L),Knu_sL_ip(IJK,M,L), KM)
               Knu_sM_b = AVG_Z_S(Knu_sM_ip(IJKB,M,L),Knu_sM_ip(IJK,M,L), KM)

               S14a_sum = S14a_sum + ( (Knu_sL_e*(NU_PL_E-NU_PL_p) - &
                    Knu_sM_e*(NU_PM_E-NU_PM_p) )*ODX_E(I)*AYZ(IJK) - (Knu_sL_w*&
                    (NU_PL_p-NU_PL_W) - Knu_sM_w*(NU_PM_p-NU_PM_W) )*ODX_E(IM)*&
                    AYZ(IMJK) )
               S14b_sum = S14b_sum + ( (Knu_sL_n*(NU_PL_N-NU_PL_p) - &
                    Knu_sM_n*(NU_PM_N-NU_PM_p) )*ODY_N(J)*AXZ(IJK) - (Knu_sL_s*&
                    (NU_PL_p-NU_PL_S) - Knu_sM_s*(NU_PM_p-NU_PM_S) )*ODY_N(JM)*&
                    AXZ(IJMK) )
               S14c_sum = S14c_sum + ( (Knu_sL_t*(NU_PL_T-NU_PL_p) - &
                    Knu_sM_t*(NU_PM_T-NU_PM_p) )*ODZ_T(K)*OX(I)*AXY(IJK) - &
                    (Knu_sL_b*(NU_PL_p-NU_PL_B) - Knu_sM_b*(NU_PM_p-NU_PM_B) )*&
                    ODZ_T(KM)*OX(I)*AXY(IJKM) )

! Part of Heat Flux: div (q)
! Kvel_s_ip*[vi-vp]
! Note S9 terms should evaluate to zero for particles from the same phase
               Kvel_s_e = AVG_X_H(Kvel_s_ip(IJK,M,L), Kvel_s_ip(IJKE,M,L),I)
               Kvel_s_w = AVG_X_H(Kvel_s_ip(IJKW,M,L),Kvel_s_ip(IJK,M,L), IM)
               Kvel_s_n = AVG_Y_H(Kvel_s_ip(IJK,M,L), Kvel_s_ip(IJKN,M,L),J)
               Kvel_s_s = AVG_Y_H(Kvel_s_ip(IJKS,M,L),Kvel_s_ip(IJK,M,L), JM)
               Kvel_s_t = AVG_Z_H(Kvel_s_ip(IJK,M,L), Kvel_s_ip(IJKT,M,L),K)
               Kvel_s_b = AVG_Z_H(Kvel_s_ip(IJKB,M,L),Kvel_s_ip(IJK,M,L), KM)

               S9_sum = S9_sum + ( Kvel_s_e*(UsM_e-UsL_e)*AYZ(IJK) - &
                    Kvel_s_w*(UsM_w-UsL_w)*AYZ(IMJK) + Kvel_s_n*(VsM_n-VsL_n)*AXZ(IJK)-&
                    Kvel_s_s*(VsM_s-VsL_s)*AXZ(IJMK) + Kvel_s_t*(WsM_t-WsL_t)*AXY(IJK)-&
                    Kvel_s_b*(WsM_b-WsL_b)*AXY(IJKM) )

          ENDIF    ! (IF M.NE.L)

      ENDDO

!  WARNING: The terms due to granular temperature gradients S21 (a,b,c) have caused
!           some converegence issues, remove them from LHS and RHS for debugging (sof).

      SOURCELHS = ( (S11_sum_lhs+S12_sum_lhs)+&
          (S10_lhs+S16_lhs+S17_sum_lhs+&
          S18_sum_lhs-S20_sum_lhs+S13_sum_lhs)*VOL(IJK) + &
          ZMAX(S21a_sum+S21b_sum+S21c_sum)+ &
          ZMAX(S14a_sum+S14b_sum+S14c_sum)+ ZMAX(S9_sum) ) / &
          Theta_m(IJK,M)

      SOURCERHS = ( S10_rhs+S15_rhs+S16_rhs+S17_sum_rhs+S18_sum_rhs+&
          S13_sum_rhs) * VOL(IJK) + S11_sum_rhs+S12_sum_rhs+ &
          ZMAX(- (S14a_sum+S14b_sum+S14c_sum) ) + ZMAX(-S9_sum) + &
          ZMAX(- (S21a_sum+S21b_sum+S21c_sum) )


      RETURN
      END SUBROUTINE SOURCE_GRANULAR_ENERGY_IA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_GRANULAR_ENERGY_GD                               C
!  Purpose: Calculate the source terms of the Mth pahse granular       C
!     energy equation                                                  C
!                                                                      C
!  Author: Janine E. Galvin                                            C
!                                                                      C
!  Literature/Document References:                                     C
!    Garzo, V., and Dufty, J., Homogeneous cooling state for a         C
!      granular mixture, Physical Review E, 1999, Vol 60 (5), 5706-    C
!      5713                                                            C
!    Garzo, Tenneti, Subramaniam, Hrenya, J. Fluid Mech., 2012, 712,   C
!      pp 129-404                                                      C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

      SUBROUTINE SOURCE_GRANULAR_ENERGY_GD(SOURCELHS, &
                    SOURCERHS, IJK, M)

!-----------------------------------------------
!  Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE drag
      USE fldvar
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE physprop
      USE rdf
      USE residual
      USE run
      USE toleranc
      USE trace
      USE turb
      USE visc_g
      USE visc_s
      use kintheory
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Source terms to be kept on rhs (source vector)
      DOUBLE PRECISION, INTENT(INOUT) :: sourcerhs
! Source terms to be kept on lhs (center coefficient)
      DOUBLE PRECISION, INTENT(INOUT) :: sourcelhs
! Solid phase index
      INTEGER, INTENT(IN) :: M
! Indices
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IM, JM, KM, IMJK, IJMK, IJKM,&
                 IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
! number densities
      DOUBLE PRECISION :: NU_PM_E, NU_PM_W, NU_PM_N, NU_PM_S, NU_PM_T,&
                          NU_PM_B, NU_PM_p
! particle characteristics
      DOUBLE PRECISION :: M_PM, D_PM
! coefficients in heat flux term
      DOUBLE PRECISION :: Knu_sM_e, Knu_sM_w, Knu_sM_n, Knu_sM_s,&
                          Knu_sM_t, Knu_sM_b
! Source terms to be kept on rhs
      DOUBLE PRECISION :: S8_rhs, S9_rhs, S10_rhs, S11_rhs, S12_rhs,&
                          S13a_rhs, S13b_rhs, S14a_rhs, S14b_rhs, S15a_rhs,&
                          S15b_rhs
! Source terms to be kept on lhs
      DOUBLE PRECISION :: S8_lhs, S10_lhs, S11_lhs, S12_lhs, S13a_lhs,&
                          S13b_lhs, S14a_lhs, S14b_lhs, S15a_lhs, S15b_lhs
! Source terms from interstitial effects
      DOUBLE PRECISION :: VSLIP, Kslip, Tslip_rhs, Tslip_lhs, Tvis_lhs
      DOUBLE PRECISION :: ep_sM
! local var. for GTSH theory
      DOUBLE PRECISION :: chi, Sgama_lhs, Spsi_rhs
!-----------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IM = Im1(I)
      JM = Jm1(J)
      KM = Km1(K)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IJKE = EAST_OF(IJK)
      IJKW = WEST_OF(IJK)
      IJKN = NORTH_OF(IJK)
      IJKS = SOUTH_OF(IJK)
      IJKT = TOP_OF(IJK)
      IJKB = BOTTOM_OF(IJK)

      S8_lhs = ZERO
      S10_lhs = ZERO
      S11_lhs = ZERO
      S12_lhs = ZERO
      S13a_lhs = ZERO
      S13b_lhs = ZERO
      S14a_lhs = ZERO
      S14b_lhs = ZERO
      S15a_lhs = ZERO
      S15b_lhs = ZERO
      S8_rhs = ZERO
      S9_rhs = ZERO
      S10_rhs = ZERO
      S11_rhs = ZERO
      S12_rhs = ZERO
      S13a_rhs = ZERO
      S13b_rhs = ZERO
      S14a_rhs = ZERO
      S14b_rhs = ZERO
      S15a_rhs = ZERO
      S15b_rhs = ZERO
      Sgama_lhs = ZERO  ! for GTSH theory
      Spsi_rhs = ZERO   ! for GTSH theory

      D_PM = D_P(IJK,M)
      M_PM = (Pi/6.d0)*D_PM**3 * RO_S(IJK,M)
      EP_SM = EP_s(IJK,M)
      CHI = G_0(IJK,M,M)

      NU_PM_p = ROP_S(IJK,M)/M_PM
      NU_PM_E = ROP_S(IJKE,M)/M_PM
      NU_PM_W = ROP_S(IJKW,M)/M_PM
      NU_PM_N = ROP_S(IJKN,M)/M_PM
      NU_PM_S = ROP_S(IJKS,M)/M_PM
      NU_PM_T = ROP_S(IJKT,M)/M_PM
      NU_PM_B = ROP_S(IJKB,M)/M_PM

! Production by shear: (S:grad(v))
! P_s*tr(D)
      S8_lhs = P_S_C(IJK,M) * ZMAX(TRD_S_C(IJK,M))
      S8_rhs = P_S_C(IJK,M) * ZMAX(-TRD_S_C(IJK,M))

! Production by shear: (S:grad(v))
! Mu_s*tr(D^2)
      S9_rhs = 2.d0*Mu_s_c(IJK,M)*TRD_S2(IJK,M)

! Production by shear: (S:grad(v))
! Lambda_s*tr(D)^2
      S10_lhs = (TRD_S_C(IJK,M)**2)*ZMAX( -LAMBDA_s_C(IJK,M) )
      S10_rhs = (TRD_S_C(IJK,M)**2)*ZMAX(  LAMBDA_s_C(IJK,M) )

! Energy dissipation by collisions: (3/2)*n*kboltz*T*zeta0
! linearized (3/2)*rop_s*T*zeta0
      S11_lhs = (3.d0/2.d0)*EDT_s_ip(IJK,M,M)
      S11_rhs = (1.d0/2.d0)*EDT_s_ip(IJK,M,M)*Theta_m(IJK,M)

! Energy dissipation by collisions: (3/2)*n*kboltz*T*zeta1
! (3/2)*rop_s*T*zeta1
      S12_lhs = ZMAX( EDvel_sM_ip(IJK,M,M) * TRD_S_C(IJK,M) )
      S12_rhs = ZMAX( -EDvel_sM_ip(IJK,M,M) * TRD_S_C(IJK,M) )*Theta_m(IJK,M)

! for GTSH theory, the dissipation terms above need to be multiplied
! by 3/2 rop_s(ijk,m)
! GTSH theory has two additional terms as source (Psi) and sink (gama)
! in eq (4.9) of GTSH JFM paper
      IF(KT_TYPE_ENUM == GTSH_2012) THEN
        S11_lhs = 1.5d0*ROP_s(IJK,M) * S11_lhs
        S11_rhs = 1.5d0*ROP_s(IJK,M) * S11_rhs
        S12_lhs = 1.5d0*ROP_s(IJK,M) * S12_lhs
        S12_rhs = 1.5d0*ROP_s(IJK,M) * S12_rhs
        Sgama_lhs = 3d0*NU_PM_p * G_gtsh(EP_SM, chi, IJK, M)
        Spsi_rhs = 1.5d0*ROP_s(IJK,M) * xsi_gtsh(ijk)
      ENDIF

! Part of Heat Flux: div (q)
! Knu_s_ip*grad(nu)
      Knu_sM_e = AVG_X_S(Kphi_s(IJK,M), Kphi_s(IJKE,M),I)
      Knu_sM_w = AVG_X_S(Kphi_s(IJKW,M),Kphi_s(IJK,M), IM)
      Knu_sM_n = AVG_Y_S(Kphi_s(IJK,M), Kphi_s(IJKN,M),J)
      Knu_sM_s = AVG_Y_S(Kphi_s(IJKS,M),Kphi_s(IJK,M), JM)
      Knu_sM_t = AVG_Z_S(Kphi_s(IJK,M), Kphi_s(IJKT,M),K)
      Knu_sM_b = AVG_Z_S(Kphi_s(IJKB,M),Kphi_s(IJK,M), KM)

      S13a_rhs = Knu_sM_e*ZMAX( (NU_PM_E-NU_PM_p) )*ODX_E(I)*AYZ(IJK)
      S13a_lhs = Knu_sM_e*ZMAX( -(NU_PM_E-NU_PM_p) )*ODX_E(I)*AYZ(IJK)

      S13b_rhs = Knu_sM_w*ZMAX( -(NU_PM_p-NU_PM_W) )*ODX_E(IM)*AYZ(IMJK)
      S13b_lhs = Knu_sM_w*ZMAX( (NU_PM_p-NU_PM_W) )*ODX_E(IM)*AYZ(IMJK)

      S14a_rhs = Knu_sM_n*ZMAX( (NU_PM_N-NU_PM_p) )*ODY_N(J)*AXZ(IJK)
      S14a_lhs = Knu_sM_n*ZMAX( -(NU_PM_N-NU_PM_p) )*ODY_N(J)*AXZ(IJK)

      S14b_rhs = Knu_sM_s*ZMAX( -(NU_PM_p-NU_PM_S) )*ODY_N(JM)*AXZ(IJMK)
      S14b_lhs = Knu_sM_s*ZMAX( (NU_PM_p-NU_PM_S) )*ODY_N(JM)*AXZ(IJMK)

      S15a_rhs = Knu_sM_t*ZMAX( (NU_PM_T-NU_PM_p) )*ODZ_T(K)*OX(I)*AXY(IJK)
      S15a_lhs = Knu_sM_t*ZMAX( -(NU_PM_T-NU_PM_p) )*ODZ_T(K)*OX(I)*AXY(IJK)

      S15b_rhs = Knu_sM_b*ZMAX( -(NU_PM_p-NU_PM_B) )*ODZ_T(KM)*OX(I)*AXY(IJKM)
      S15b_lhs = Knu_sM_b*ZMAX( (NU_PM_p-NU_PM_B) )*ODZ_T(KM)*OX(I)*AXY(IJKM)

      SOURCELHS = ( (S8_lhs+S10_lhs)*VOL(IJK) + &
          S13a_lhs+S13b_lhs+S14a_lhs+S14b_lhs+S15a_lhs+S15b_lhs)/Theta_m(IJK,M)&
          + (S11_lhs + S12_lhs + Sgama_lhs)*VOL(IJK)


      SOURCERHS = ( S8_rhs+S9_rhs+S10_rhs+S11_rhs+S12_rhs+Spsi_rhs) * VOL(IJK) + &
          S13a_rhs+S13b_rhs+S14a_rhs+S14b_rhs+S15a_rhs+S15b_rhs


! this is only done for GD_99, do not add these terms to GTSH
      IF(SWITCH > ZERO .AND. RO_g0 /= ZERO .AND. &
         (KT_TYPE_ENUM == GD_1999)) THEN

         VSLIP = (U_S(IJK,M)-U_G(IJK))**2 + (V_S(IJK,M)-V_G(IJK))**2 +&
            (W_S(IJK,M)-W_G(IJK))**2
         VSLIP = DSQRT(VSLIP)

! production by gas-particle slip: Koch & Sangani (1999)
         Kslip = SWITCH*81.d0*EP_SM*(MU_G(IJK)*VSLIP)**2.d0 / &
            (chi*D_P(IJK,M)**3.D0*RO_S(IJK,M)*DSQRT(PI))

         Tslip_rhs = 1.5d0*Kslip/( (THETA_M(IJK,M)+SMALL_NUMBER)**0.50)*VOL(IJK)
         Tslip_lhs = 0.5d0*Kslip/( (THETA_M(IJK,M)+SMALL_NUMBER)**1.50)*VOL(IJK)

! dissipation by viscous damping: Gidaspow (1994)
         Tvis_lhs = SWITCH*3d0*F_GS(IJK,M)*VOL(IJK)

         SOURCELHS = SOURCELHS + Tslip_lhs + Tvis_lhs
         SOURCERHS = SOURCERHS + Tslip_rhs
      ENDIF

      RETURN
      END SUBROUTINE SOURCE_GRANULAR_ENERGY_GD


