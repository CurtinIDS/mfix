!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_IA_MOMSOURCE_V_S                                   C
!  Purpose: Determine source terms for V_S momentum equation arising   C
!           from IA kinetic theory constitutive relations for stress   C
!           and solid-solid drag                                       C
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

      SUBROUTINE CALC_IA_MOMSOURCE_V_S(M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: zero
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
      INTEGER :: I, J, K, IJK, IJKN, JP, IM,  KM, IJPK, IJMK,&
                 IJKE, IJKNE, IJKW, IJKNW, IMJPK, IMJK, IJKT,&
                 IJKTN, IJKB, IJKBN, IJKM, IJPKM,&
                 IPJK, IJKP
! Phase index
      INTEGER :: L
! Viscosity values
      DOUBLE PRECISION :: MU_sL_pE, MU_sL_pW, MU_sL_pN, MU_sL_pS, MU_sL_pT,&
                          MU_sL_pB
! Bulk viscosity values
      DOUBLE PRECISION :: XI_sL_pN, XI_sL_pS, LAMBDA_sL_pN, LAMBDA_sL_pS
! Variables for drag calculations
      DOUBLE PRECISION :: M_PM, M_PL, D_PM, D_PL, NU_PM_pN, NU_PM_pS, NU_PM_p, &
                          NU_PL_pN, NU_PL_pS, NU_PL_p, T_PM_pN, T_PM_pS, &
                          T_PL_pN, T_PL_pS, Fnu_s_p, FT_sM_p, FT_sL_p
! Average volume fraction
      DOUBLE PRECISION :: EPSA
! Source terms (Surface)
      DOUBLE PRECISION :: ssvx, ssvy, ssvz, ssx, ssy, ssz, ssbv
! Source terms (Volumetric)
      DOUBLE PRECISION :: DS1, DS2, DS3, DS4, DS1plusDS2
!-----------------------------------------------

! section largely based on tau_v_g:

      DO IJK = IJKSTART3, IJKEND3

! Skip walls where some values are undefined.
         IF(WALL_AT(IJK)) cycle

          D_PM = D_P(IJK,M)
          M_PM = (Pi/6d0) * D_PM**3 *RO_S(IJK,M)

          J = J_OF(IJK)
          IJKN = NORTH_OF(IJK)
          EPSA = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)
          IF ( .NOT.SIP_AT_N(IJK) .AND. EPSA>DIL_EP_S) THEN

               JP = JP1(J)
               I = I_OF(IJK)
               IM = IM1(I)
               K = K_OF(IJK)
               KM = KM1(K)
               IJPK = JP_OF(IJK)
               IJMK = JM_OF(IJK)
               IMJK = IM_OF(IJK)
               IJKM = KM_OF(IJK)
               IMJPK = IM_OF(IJPK)
               IJPKM = JP_OF(IJKM)

               IJKW = WEST_OF(IJK)
               IJKE = EAST_OF(IJK)
               IJKNE = EAST_OF(IJKN)
               IJKNW = NORTH_OF(IJKW)

               IJKB = BOTTOM_OF(IJK)
               IJKT = TOP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)
               IJKBN = NORTH_OF(IJKB)

! additional required quantities:
               IPJK = IP_OF(IJK)
               IJKP = KP_OF(IJK)

! initialize variable
               STRESS_TERMS = ZERO
               DRAG_TERMS = ZERO

               DO L = 1, MMAX
                    IF (M .ne. L) THEN

!--------------------- Sources from Stress Terms ---------------------
! Surface Forces
! standard shear stress terms (i.e. ~diffusion)
                         MU_sL_pE = AVG_Y_H(AVG_X_H(MU_sL_ip(IJK,M,L),MU_sL_ip(IJKE,M,L),I),&
                              AVG_X_H(MU_sL_ip(IJKN,M,L),MU_sL_ip(IJKNE,M,L),I),J)
                         MU_sL_pW = AVG_Y_H(AVG_X_H(MU_sL_ip(IJKW,M,L),MU_sL_ip(IJK,M,L),IM),&
                              AVG_X_H(MU_sL_ip(IJKNW,M,L),MU_sL_ip(IJKN,M,L),IM),J)
                         SSVX = MU_sL_pE*(V_S(IPJK,L)-V_S(IJK,L))*AYZ_V(IJK)*ODX_E(I)&
                              -MU_sL_pW*(V_S(IJK,L)-V_S(IMJK,L))*AYZ_V(IMJK)*ODX_E(IM)

                         MU_sL_pN = MU_sL_ip(IJKN,M,L)
                         MU_sL_pS = MU_sL_ip(IJK,M,L)
                         SSVY = MU_sL_pN*(V_S(IJPK,L)-V_S(IJK,L))*ODY(JP)*AXZ_V(IJK)&
                              -MU_sL_pS*(V_S(IJK,L)-V_S(IJMK,L))*ODY(J)*AXZ_V(IJMK)

                         MU_sL_pT = AVG_Y_H(AVG_Z_H(MU_sL_ip(IJK,M,L),MU_sL_ip(IJKT,M,L),K),&
                              AVG_Z_H(MU_sL_ip(IJKN,M,L),MU_sL_ip(IJKTN,M,L),K),J)
                         MU_sL_pB = AVG_Y_H(AVG_Z_H(MU_sL_ip(IJKB,M,L),MU_sL_ip(IJK,M,L),KM),&
                              AVG_Z_H(MU_sL_ip(IJKBN,M,L),MU_sL_ip(IJKN,M,L),KM),J)
                         SSVZ = MU_sL_pT*(V_S(IJKP,L)-V_S(IJK,L))*AXY_V(IJK)*ODZ_T(K)*OX(I)&
                              -MU_sL_pB*(V_S(IJK,L)-V_S(IJKM,L))*AXY_V(IJKM)*ODZ_T(KM)*OX(I)

! bulk viscosity term
                         XI_sL_pN = XI_sL_ip(IJKN,M,L)
                         XI_sL_pS = XI_sL_ip(IJK,M,L)
                         LAMBDA_sL_pN = -(2.d0/3.d0)*MU_sL_pN + XI_sL_pN
                         LAMBDA_sL_pS = -(2.d0/3.d0)*MU_sL_pS + XI_sL_pS
                         SSBV = (LAMBDA_sL_pN*TRD_S(IJKN,L)-LAMBDA_sL_pS*TRD_S(IJK,L))*AXZ(IJK)

! off diagonal shear stress terms
                         SSX = MU_sL_pE*(U_S(IJPK,L)-U_S(IJK,L))*ODY_N(J)*AYZ_V(IJK)&
                              -MU_sL_pW*(U_S(IMJPK,L)-U_S(IMJK,L))*ODY_N(J)*AYZ_V(IMJK)
                         SSY = SSVY
                         SSZ = MU_sL_pT*(W_S(IJPK,L)-W_S(IJK,L))*AXY_V(IJK)*ODY_N(J)&
                              -MU_sL_pB*(W_S(IJPKM,L)-W_S(IJKM,L))*AXY_V(IJKM)*ODY_N(J)
!--------------------- End Sources from Stress Term ---------------------


!--------------------- Sources from Momentum Source Term ---------------------
! Momentum source associated with the difference in the gradients in
! number density of solids phase m and all other solids phases
                         D_PL = D_P(IJK,L)
                         M_PL = (Pi/6d0)* D_PL**3 *RO_S(IJK,L)

                         NU_PM_pN = ROP_S(IJKN,M)/M_PM
                         NU_PM_pS = ROP_S(IJK,M)/M_PM
                         NU_PM_p = AVG_Y(NU_PM_pS,NU_PM_pN,J)

                         NU_PL_pN = ROP_S(IJKN,L)/M_PL
                         NU_PL_pS = ROP_S(IJK,L)/M_PL
                         NU_PL_p = AVG_Y(NU_PL_pS,NU_PL_pN,J)

                         Fnu_s_p = AVG_Y(Fnu_s_ip(IJK,M,L),Fnu_s_ip(IJKN,M,L),J)
                         DS1 = Fnu_s_p*NU_PL_p*(NU_PM_pN-NU_PM_pS)*ODY_N(J)
                         DS2 = -Fnu_s_p*NU_PM_p*(NU_PL_pN-NU_PL_pS)*ODY_N(J)
                         DS1plusDS2 = DS1 + DS2

! Momentum source associated with the gradient in granular
! temperature of species M
                         T_PM_pN = Theta_M(IJKN,M)
                         T_PM_pS = Theta_M(IJK,M)

                         FT_sM_p = AVG_Y(FT_sM_ip(IJK,M,L),FT_sM_ip(IJKN,M,L),J)
                         DS3 = FT_sM_p*(T_PM_pN-T_PM_pS)*ODY_N(J)

! Momentum source associated with the gradient in granular
! temperature of species L
                         T_PL_pN = Theta_M(IJKN,L)
                         T_PL_pS = Theta_M(IJK,L)

                         FT_sL_p = AVG_Y(FT_sL_ip(IJK,M,L),FT_sL_ip(IJKN,M,L),J)
                         DS4 = FT_sL_p*(T_PL_pN-T_PL_pS)*ODY_N(J)
!--------------------- End Sources from Momentum Source Term ---------------------


! Add the terms
                         STRESS_TERMS = STRESS_TERMS + SSVX + SSVY + SSVZ + &
                             SSBV + SSX + SSY + SSZ
                         DRAG_TERMS = DRAG_TERMS + (DS1plusDS2+DS3+DS4)*VOL_V(IJK)

                    ELSE ! if m .ne. L
! for m = l all stress terms should already be handled in existing routines
! for m = l all drag terms should become zero
                         STRESS_TERMS = STRESS_TERMS + ZERO
                         DRAG_TERMS = DRAG_TERMS + ZERO

                    ENDIF ! if m .ne. L
               ENDDO     ! over L

               KTMOM_V_S(IJK,M) = STRESS_TERMS + DRAG_TERMS
          ELSE
               KTMOM_V_S(IJK,M) = ZERO
          ENDIF     ! dilute
      ENDDO        ! over ijk

      RETURN
      END SUBROUTINE CALC_IA_MOMSOURCE_V_S
