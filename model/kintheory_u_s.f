!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_EXPLICIT_MOM_SOURCE_S                              C
!  Purpose: Determine any additional momentum source terms that are    C
!  calculated explicitly here                                          C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_EXPLICIT_MOM_SOURCE_S()

!-----------------------------------------------
! Modules
!-----------------------------------------------

      USE param1, only: zero

! kinetic theories
      USE run, only: kt_type_enum
      USE run, only: ia_2005

! number of solids phases
      USE physprop, only: smax

! solids source term
      USE kintheory, only: ktmom_u_s, ktmom_v_s, ktmom_w_s

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Solids phase index
      INTEGER :: M
!-----------------------------------------------
! Additional interphase interaction terms that arise from kinetic theory
      DO M = 1, SMAX
         KTMOM_U_s(:,M) = ZERO
         KTMOM_V_S(:,M) = ZERO
         KTMOM_W_S(:,M) = ZERO

         IF (KT_TYPE_ENUM == IA_2005) THEN
            CALL CALC_IA_MOMSOURCE_U_S(M)
            CALL CALC_IA_MOMSOURCE_V_S(M)
            CALL CALC_IA_MOMSOURCE_W_S(M)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_EXPLICIT_MOM_SOURCE_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_IA_MOMSOURCE_U_S                                   C
!  Purpose: Determine source terms for U_S momentum equation arising   C
!           from IA kinetic theory constitutive relations for stress   C
!           and solids-solids drag (collisional momentum source)       C
!                                                                      C
!  Literature/Document References:                                     C
!    Idir, Y.H., "Modeling of the multiphase mixture of particles      C
!      using the kinetic theory approach," PhD Thesis, Illinois        C
!      Institute of Technology, Chicago, Illinois, 2004                C
!    Iddir, Y.H., & H. Arastoopour, "Modeling of multitype particle    C
!      flow using the kinetic theory approach," AIChE J., Vol 51,      C
!      No 6, June 2005                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_IA_MOMSOURCE_U_S(M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: half, zero
      USE constant, only: pi

! trace of D_s at i, j, k
      USE visc_s, only: trD_s

! number of solids phases
      USE physprop, only: mmax

! x,y,z-components of solids velocity
      USE fldvar, only: u_s, v_s, w_s
! particle diameter, bulk density, material density
      USE fldvar, only: d_p, rop_s, ro_s
! granular temperature
      USE fldvar, only: theta_m, ep_s
! dilute threshold
      USE toleranc, only: dil_ep_s

      Use kintheory

      USE geometry
      USE indices
      USE compar
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! solids phase index
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Temporary variable
      DOUBLE PRECISION :: STRESS_TERMS, DRAG_TERMS
! Indices
      INTEGER :: I, J, JM, K, KM, IJK, IJKE, IPJK, IP, IMJK, IJKN,&
                 IJKNE, IJKS, IJKSE, IPJMK, IJMK, IJKT, IJKTE,&
                 IJKB, IJKBE, IJKM, IPJKM, &
                 IJPK, IJKP, IM, IJKW
! solids phase index
      INTEGER :: L
! Viscosity values for stress calculations
      DOUBLE PRECISION :: MU_sL_pE, MU_sL_pW, MU_sL_pN, MU_sL_pS, MU_sL_pT,&
                          MU_sL_pB, Mu_sL_p
! Bulk viscosity values for calculations
      DOUBLE PRECISION :: XI_sL_pE, XI_sL_pW, LAMBDA_sL_pE, LAMBDA_sL_pW
! Variables for drag calculations
      DOUBLE PRECISION :: M_PM, M_PL, D_PM, D_PL, NU_PM_pE, NU_PM_pW, NU_PM_p, &
                          NU_PL_pE, NU_PL_pW, NU_PL_p, T_PM_pE, T_PM_pW, &
                          T_PL_pE, T_PL_pW, Fnu_s_p, FT_sM_p, FT_sL_p
! Average volume fraction
      DOUBLE PRECISION :: EPSA
! Source terms (Surface)
      DOUBLE PRECISION :: ssux, ssuy, ssuz, ssx, ssy, ssz, ssbv
! Source terms (Volumetric)
      DOUBLE PRECISION :: tauzz, DS1, DS2, DS3, DS4, DS1plusDS2
!-----------------------------------------------

! section largely based on tau_u_g:

      DO IJK = IJKSTART3, IJKEND3

! Skip walls where some values are undefined.
         IF(WALL_AT(IJK)) cycle

         D_PM = D_P(IJK,M)
         M_PM = (Pi/6d0)*D_PM**3 * RO_S(IJK,M)

         I = I_OF(IJK)
         IJKE = EAST_OF(IJK)
         EPSA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
         IF ( .NOT.SIP_AT_E(IJK) .AND. EPSA>DIL_EP_S) THEN

            IP = IP1(I)
            IM = Im1(I)
            IJKW  = WEST_OF(IJK)
            J = J_OF(IJK)
            JM = JM1(J)
            K = K_OF(IJK)
            KM = KM1(K)
            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            IPJMK = JM_OF(IPJK)
            IPJKM = IP_OF(IJKM)

            IJKN = NORTH_OF(IJK)
            IJKS = SOUTH_OF(IJK)
            IJKNE = EAST_OF(IJKN)
            IJKSE = EAST_OF(IJKS)

            IJKT = TOP_OF(IJK)
            IJKB = BOTTOM_OF(IJK)
            IJKTE = EAST_OF(IJKT)
            IJKBE = EAST_OF(IJKB)

! additional required quantities:
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)

! initialize variable
            STRESS_TERMS = ZERO
            DRAG_TERMS = ZERO

            DO L = 1,MMAX
               IF (M .ne. L) THEN

!--------------------- Sources from Stress Terms ---------------------
! Surface Forces
! standard shear stress terms (i.e. ~diffusion)
                  MU_sL_pE = MU_sL_ip(IJKE,M,L)
                  MU_sL_pW = MU_sL_ip(IJK,M,L)
                  SSUX = MU_sL_pE*(U_S(IPJK,L)-U_S(IJK,L))*AYZ_U(IJK)*ODX(IP)&
                       -MU_sL_PW*(U_S(IJK,L)-U_S(IMJK,L))*AYZ_U(IMJK)*ODX(I)

                  MU_sL_pN = AVG_X_H(AVG_Y_H(MU_sL_ip(IJK,M,L),MU_sL_ip(IJKN,M,L), J),&
                       AVG_Y_H(MU_sL_ip(IJKE,M,L),MU_sL_ip(IJKNE,M,L), J), I)
                  MU_sL_pS = AVG_X_H(AVG_Y_H(MU_sL_ip(IJKS,M,L),MU_sL_ip(IJK,M,L), JM),&
                       AVG_Y_H(MU_sL_ip(IJKSE,M,L),MU_sL_ip(IJKE,M,L), JM), I)
                  SSUY = MU_sL_pN*(U_S(IJPK,L)-U_S(IJK,L))*AXZ_U(IJK)*ODY_N(J)&
                       -MU_sL_pS*(U_S(IJK,L)-U_S(IJMK,L))*AXZ_U(IJMK)*ODY_N(JM)

                  MU_sL_pT = AVG_X_H(AVG_Z_H(MU_sL_ip(IJK,M,L),MU_sL_ip(IJKT,M,L),K),&
                       AVG_Z_H(MU_sL_ip(IJKE,M,L),MU_sL_ip(IJKTE,M,L),K),I)
                  MU_sL_pB = AVG_X_H(AVG_Z_H(MU_sL_ip(IJKB,M,L),MU_sL_ip(IJK,M,L),KM),&
                       AVG_Z_H(MU_sL_ip(IJKBE,M,L),MU_sL_ip(IJKE,M,L),KM),I)
                  SSUZ = MU_sL_pT*(U_S(IJKP,L)-U_S(IJK,L))*AXY_U(IJK)*ODZ_T(K)*OX_E(I)&
                       -MU_sL_pB*(U_S(IJK,L)-U_S(IJKM,L))*AXY_U(IJKM)*ODZ_T(KM)*OX_E(I)

! bulk viscosity term
                  XI_sL_pE = XI_sL_ip(IJKE,M,L)
                  XI_sL_pW = XI_sL_ip(IJK,M,L)
                  LAMBDA_sL_pE = -(2.d0/3.d0)*Mu_sL_pE + XI_sL_pE
                  LAMBDA_sL_pW = -(2.d0/3.d0)*Mu_sL_pW + XI_sL_pW
                  SSBV = (LAMBDA_sL_pE*TRD_S(IJKE,L)-LAMBDA_sL_pW*TRD_S(IJK,L))*AYZ(IJK)

! off diagonal shear stress terms
                  SSX = SSUX
                  SSY = MU_sL_pN*(V_S(IPJK,L)-V_S(IJK,L))*AXZ_U(IJK)*ODX_E(I)&
                       -MU_sL_pS*(V_S(IPJMK,L)-V_S(IJMK,L))*AXZ_U(IJMK)*ODX_E(I)
                  SSZ = MU_sL_pT*(W_S(IPJK,L)-W_S(IJK,L))*AXY_U(IJK)*ODX_E(I)&
                       -MU_sL_pB*(W_S(IPJKM,L)-W_S(IJKM,L))*AXY_U(IJKM)*ODX_E(I)

! special terms for cylindrical coordinates
                  IF (CYLINDRICAL) THEN

! modify Ssz: (1/x) (d/dz) (tau_zz)
!                    integral of (1/x) (d/dz) (mu*(-w/x))
!                    (normally part of tau_u_s) - explicit
                     SSZ = SSZ - (MU_sL_pT*(W_S(IPJK,L)+W_S(IJK,L))*HALF*OX_E(I)*AXY_U(IJK)&
                          -MU_sL_pB*(W_S(IPJKM,L)+W_S(IJKM,L))*HALF*OX_E(I)*AXY_U(IJKM))

! tau_zz/x terms: (tau_zz/x)
!                    integral of -(2mu/x)*(1/x)*dw/dz
!                    (normally part of tau_u_s) - explicit
                     MU_sL_p = AVG_X(MU_sL_ip(IJK,M,L),MU_sL_ip(IJKE,M,L),I)
                     tauzz = -2.d0*MU_sL_p*OX_E(I)*HALF*(&
                          ((W_S(IPJK,L)-W_S(IPJKM,L))*OX(IP)*ODZ(K))+&
                          ((W_S(IJK,L)-W_S(IJKM,L))*OX(I)*ODZ(K)) &
                          ) * VOL_U(IJK)

!                    integral of -(2mu/x)*(1/x)*u
!                    (normally part of source_u_s)
                     tauzz = tauzz + (-2.d0*MU_sL_p*OX_E(I)*OX_E(I)*&
                          U_S(IJK,L)*VOL_U(IJK))
                  ELSE
                     tauzz = ZERO
                  ENDIF
!--------------------- End Sources from Stress Term ---------------------


!--------------------- Sources from Momentum Source Term ---------------------
! Momentum source associated with the difference in the gradients in
! number density of solids phase m and all other solids phases
                  D_PL = D_P(IJK,L)
                  M_PL = (PI/6.d0)*D_PL**3 * RO_S(IJK,L)

                  NU_PM_pE = ROP_S(IJKE,M)/M_PM
                  NU_PM_pW = ROP_S(IJK,M)/M_PM
                  NU_PM_p = AVG_X(NU_PM_pW,NU_PM_pE,I)

                  NU_PL_pE = ROP_S(IJKE,L)/M_PL
                  NU_PL_pW = ROP_S(IJK,L)/M_PL
                  NU_PL_p = AVG_X(NU_PL_pW,NU_PL_pE,I)

                  Fnu_s_p = AVG_X(Fnu_s_ip(IJK,M,L),Fnu_s_ip(IJKE,M,L),I)
                  DS1 = Fnu_s_p*NU_PL_p*(NU_PM_pE-NU_PM_pW)*ODX_E(I)
                  DS2 = -Fnu_s_p*NU_PM_p*(NU_PL_pE-NU_PL_pW)*ODX_E(I)
                  DS1plusDS2 = DS1 + DS2

! Momentum source associated with the gradient in granular
! temperature of species M
                  T_PM_pE = Theta_M(IJKE,M)
                  T_PM_pW = Theta_M(IJK,M)

                  FT_sM_p = AVG_X(FT_sM_ip(IJK,M,L),FT_sM_ip(IJKE,M,L),I)
                  DS3 = FT_sM_p*(T_PM_pE-T_PM_pW)*ODX_E(I)

! Momentum source associated with the gradient in granular
! temperature of species L
                  T_PL_pE = Theta_M(IJKE,L)
                  T_PL_pW = Theta_M(IJK,L)

                  FT_sL_p = AVG_X(FT_sL_ip(IJK,M,L),FT_sL_ip(IJKE,M,L),I)
                  DS4 = FT_sL_p*(T_PL_pE-T_PL_pW)*ODX_E(I)
!--------------------- End Sources from Momentum Source Term ---------------------


! Add the terms
                  STRESS_TERMS = STRESS_TERMS + SSUX + SSUY + SSUZ + &
                       SSBV + SSX + SSY + SSZ + tauzz
                  DRAG_TERMS = DRAG_TERMS + (DS1plusDS2+DS3+DS4)*VOL_U(IJK)

               ELSE ! if m .ne. L
! for m = l all stress terms should already be handled in existing routines
! for m = l all drag terms should become zero
                  STRESS_TERMS = STRESS_TERMS + ZERO
                  DRAG_TERMS = DRAG_TERMS + ZERO

               ENDIF ! if m .ne. L
            ENDDO     ! over L

            KTMOM_U_S(IJK,M) = STRESS_TERMS + DRAG_TERMS
         ELSE
            KTMOM_U_S(IJK,M) = ZERO
         ENDIF     ! dilute
      ENDDO        ! over ijk

      RETURN
      END SUBROUTINE CALC_IA_MOMSOURCE_U_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: COLL_MOMENTUM_COEFF_IA                                  !
!  Purpose: Compute collisional momentum source terms betweem solids   !
!  phase M and solids phase L using Iddir Arastoopour (2005) kinetic   !
!  theory model that are not proportional to the relative velocity     !
!  between the two phases.  Specifically, terms proportional to the    !
!  gradient in number density and gradient in temperature              !
!                                                                      !
!  Literature/Document References:                                     !
!    Iddir, Y.H., "Modeling of the multiphase mixture of particles     !
!       using the kinetic theory approach," PhD Thesis, Illinois       !
!       Institute of Technology, Chicago, Illinois, 2004               !
!    Iddir, Y.H., & H. Arastoopour, "Modeling of Multitype particle    !
!      flow using the kinetic theory approach," AIChE J., Vol 51,      !
!      no. 6, June 2005                                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE COLL_MOMENTUM_COEFF_IA(L, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE drag
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE kintheory
      USE param1
      USE physprop
      USE rdf
      USE sendrecv

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: L
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
! Particle diameters
      DOUBLE PRECISION :: D_PM, D_PL
! Sum of particle diameters
      DOUBLE PRECISION :: DPSUM
!
      DOUBLE PRECISION :: M_PM, M_PL, MPSUM, DPSUMo2, NU_PL, NU_PM
      DOUBLE PRECISION :: Ap_lm, Dp_lm, Bp_lm
      DOUBLE PRECISION :: R0p_lm, R3p_lm, R4p_lm, R10p_lm
      DOUBLE PRECISION :: Fnus_ip, FTsM_ip, FTsL_ip, F_common_term
!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3
         IF (.NOT.WALL_AT(IJK)) THEN

            IF (M == L) THEN
               Fnu_s_ip(IJK,M,L) = ZERO
               FT_sM_ip(IJK,M,L) = ZERO
               FT_sL_ip(IJK,M,L) = ZERO

            ELSE
               D_PM = D_P(IJK,M)
               D_PL = D_P(IJK,L)
               DPSUM = D_PL + D_PM
               M_PM = (Pi/6.d0) * D_PM**3 *RO_S(IJK,M)
               M_PL = (Pi/6.d0) * D_PL**3 *RO_S(IJK,L)
               MPSUM = M_PM + M_PL
               DPSUMo2 = DPSUM/2.d0
               NU_PM = ROP_S(IJK,M)/M_PM
               NU_PL = ROP_S(IJK,L)/M_PL

               IF(Theta_m(IJK,M) > ZERO .AND. Theta_m(IJK,L) > ZERO) THEN

                  Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*Theta_m(IJK,M))/&
                        2.d0
                  Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-Theta_m(IJK,M) ))/&
                       (2.d0*MPSUM)
                  Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+M_PL*Theta_m(IJK,L) ))/&
                       (2.d0*MPSUM*MPSUM)

                  R0p_lm = ( 1.d0/( Ap_lm**1.5 * Dp_lm**2.5 ) )+ &
                            ( (15.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**2.5 * Dp_lm**3.5 ) )+&
                            ( (175.d0*(Bp_lm**4))/( 8.d0*Ap_lm**3.5 * Dp_lm**4.5 ) )

                  R3p_lm = ( 1.d0/( (Ap_lm**1.5)*(Dp_lm**3.5) ) )+&
                            ( (21.d0*Bp_lm*Bp_lm)/( 2.d0 * Ap_lm**2.5 * Dp_lm**4.5 ) )+&
                            ( (315.d0*Bp_lm**4)/( 8.d0 * Ap_lm**3.5 *Dp_lm**5.5 ) )

                  R4p_lm = ( 3.d0/( Ap_lm**2.5 * Dp_lm**3.5 ) )+&
                            ( (35.d0*Bp_lm*Bp_lm)/( 2.d0 * Ap_lm**3.5 * Dp_lm**4.5 ) )+&
                            ( (441.d0*Bp_lm**4)/( 8.d0 * Ap_lm**4.5 * Dp_lm**5.5 ) )

                  R10p_lm = ( 1.d0/( Ap_lm**2.5 * Dp_lm**2.5 ) )+&
                            ( (25.d0*Bp_lm*Bp_lm)/( 2.d0* Ap_lm**3.5 * Dp_lm**3.5 ) )+&
                            ( (1225.d0*Bp_lm**4)/( 24.d0* Ap_lm**4.5 * Dp_lm**4.5 ) )

                  F_common_term = (DPSUMo2*DPSUMo2/4.d0)*(M_PM*M_PL/MPSUM)*&
                         G_0(IJK,M,L)*(1.d0+C_E)*(M_PM*M_PL)**1.5

! Momentum source associated with the difference in the gradients in
! number density of solids phase m and all other solids phases
                  Fnus_ip = F_common_term*(PI*DPSUMo2/12.d0)*R0p_lm*&
                         (Theta_m(IJK,M)*Theta_m(IJK,L))**2.5

! Momentum source associated with the gradient in granular temperature
! of solid phase M
                  FTsM_ip = F_common_term*NU_PM*NU_PL*DPSUMo2*PI*&
                            (Theta_m(IJK,M)**1.5 * Theta_m(IJK,L)**2.5) *&
                            (  (-1.5d0/12.d0*R0p_lm)+&
                            Theta_m(IJK,L)/16.d0*(  (-M_PM*R10p_lm) - &
                            ((5.d0*M_PL*M_PL*M_PM/(192.d0*MPSUM*MPSUM))*R3p_lm)+&
                            ((5.d0*M_PM*M_PL)/(96.d0*MPSUM)*R4p_lm*Bp_lm)  )  )

! Momentum source associated with the gradient in granular temperature
! of solid phase L ! no need to recompute (sof Aug 30 2006)
                  FTsL_ip = F_common_term*NU_PM*NU_PL*DPSUMo2*PI*&
                           (Theta_m(IJK,L)**1.5 * Theta_m(IJK,M)**2.5) *&
                            (  (1.5d0/12.d0*R0p_lm)+&
                            Theta_m(IJK,M)/16.d0*(  (M_PL*R10p_lm)+&
                            (5.d0*M_PM*M_PM*M_PL/(192.d0*MPSUM*MPSUM)*R3p_lm)+&
                            (5.d0*M_PM*M_PL/(96.d0*MPSUM) *R4p_lm*Bp_lm)  )  )


                  Fnu_s_ip(IJK,M,L) = Fnus_ip

! WARNING: the following two terms have caused some convergence problems
! earlier. Set them to ZERO for debugging in case of converegence
! issues. (sof)
                  FT_sM_ip(IJK,M,L) = FTsM_ip  ! ZERO
                  FT_sL_ip(IJK,M,L) = FTsL_ip  ! ZERO
               ELSE
                  Fnu_s_ip(IJK,M,L) = ZERO
                  FT_sM_ip(IJK,M,L) = ZERO
                  FT_sL_ip(IJK,M,L) = ZERO
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE COLL_MOMENTUM_COEFF_IA

