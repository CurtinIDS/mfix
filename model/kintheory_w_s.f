!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_IA_MOMSOURCE_W_S                                   C
!  Purpose: Determine source terms for W_S momentum equation arising   C
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

      SUBROUTINE CALC_IA_MOMSOURCE_W_S(M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: half, zero, undefined
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
      INTEGER :: IJK, J, I, IM, IJKP, IMJK, IJKN, IJKNT, IJKS,&
                 IJKST, IJMKP, IJMK, IJKE, IJKTE, IJKW, IJKTW,&
                 IMJKP, K, IJKT, JM, KP, IJKM, IPJK,&
                 IJPK
! Phase index
      INTEGER :: L
! Viscosity values
      DOUBLE PRECISION :: MU_sL_pE, MU_sL_pW, MU_sL_pN, MU_sL_pS, MU_sL_pT,&
                          MU_sL_pB, MU_sL_p
! Bulk viscosity values
      DOUBLE PRECISION :: XI_sL_pT, XI_sL_pB, LAMBDA_sL_pT, LAMBDA_sL_pB
! Variables for drag calculations
      DOUBLE PRECISION :: M_PM, M_PL, D_PM, D_PL, NU_PM_pT, NU_PM_pB, NU_PM_p, &
                          NU_PL_pT, NU_PL_pB, NU_PL_p, T_PM_pT, T_PM_pB,&
                          T_PL_pT, T_PL_pB, Fnu_s_p, FT_sM_p, FT_sL_p
! Average velocity gradients
      DOUBLE PRECISION :: duodz
! Average volume fraction
      DOUBLE PRECISION :: EPSA
! Source terms (Surface)
      DOUBLE PRECISION :: sswx, sswy, sswz, ssx, ssy, ssz, ssbv, tauxz_x
! Source terms (Volumetric)
      DOUBLE PRECISION :: tauxz_ox, DS1, DS2, DS3, DS4, DS1plusDS2
!-----------------------------------------------

! section largely based on tau_w_g:

      DO IJK = IJKSTART3, IJKEND3

! Skip walls where some values are undefined.
         IF(WALL_AT(IJK)) cycle

          D_PM = D_P(IJK,M)
          M_PM = (Pi/6d0)*(D_PM**3.)*RO_S(IJK,M)

          K = K_OF(IJK)
          IJKT = TOP_OF(IJK)
          EPSA = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
          IF ( .NOT.SIP_AT_T(IJK) .AND. EPSA>DIL_EP_S) THEN

               J = J_OF(IJK)
               I = I_OF(IJK)
               KP = KP1(K)
               IM = IM1(I)
               JM = JM1(J)
               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IJKP = KP_OF(IJK)
               IJMKP = JM_OF(IJKP)
               IMJKP = KP_OF(IMJK)

               IJKN = NORTH_OF(IJK)
               IJKS = SOUTH_OF(IJK)
               IJKNT = TOP_OF(IJKN)
               IJKST = TOP_OF(IJKS)

               IJKE = EAST_OF(IJK)
               IJKW = WEST_OF(IJK)
               IJKTE = EAST_OF(IJKT)
               IJKTW = WEST_OF(IJKT)

! additional required quantities:
               IPJK = IP_OF(IJK)
               IJPK = JP_OF(IJK)

! initialize variable
               STRESS_TERMS = ZERO
               DRAG_TERMS = ZERO

               DO L = 1, MMAX
                    IF (M .ne. L) THEN

!--------------------- Sources from Stress Terms ---------------------
! Surface Forces
! standard shear stress terms (i.e. ~diffusion)
                         MU_sL_pE = AVG_Z_H(AVG_X_H(MU_sL_ip(IJK,M,L),MU_sL_ip(IJKE,M,L),I),&
                              AVG_X_H(MU_sL_ip(IJKT,M,L),MU_sL_ip(IJKTE,M,L),I),K)
                         MU_sL_pW = AVG_Z_H(AVG_X_H(MU_sL_ip(IJKW,M,L),MU_sL_ip(IJK,M,L),IM),&
                              AVG_X_H(MU_sL_ip(IJKTW,M,L),MU_sL_ip(IJKT,M,L),IM),K)
                         SSWX = MU_sL_pE*(W_S(IPJK,L)-W_S(IJK,L))*AYZ_W(IJK)*ODX_E(I)&
                              -MU_sL_PW*(W_S(IJK,L)-W_S(IMJK,L))*AYZ_W(IMJK)*ODX_E(IM)

                         MU_sL_pN = AVG_Z_H(AVG_Y_H(MU_sL_ip(IJK,M,L),MU_sL_ip(IJKN,M,L), J),&
                              AVG_Y_H(MU_sL_ip(IJKT,M,L),MU_sL_ip(IJKNT,M,L), J), K)
                         MU_sL_pS = AVG_Z_H(AVG_Y_H(MU_sL_ip(IJKS,M,L),MU_sL_ip(IJK,M,L), JM),&
                              AVG_Y_H(MU_sL_ip(IJKST,M,L),MU_sL_ip(IJKT,M,L), JM), K)
                         SSWY = MU_sL_pN*(W_S(IJPK,L)-W_S(IJK,L))*AXZ_W(IJK)*ODY_N(J)&
                              -MU_sL_pS*(W_S(IJK,L)-W_S(IJMK,L))*AXZ_W(IJKM)*ODY_N(JM)

                         MU_sL_pT = MU_sL_ip(IJKT,M,L)
                         MU_sL_pB = MU_sL_ip(IJK,M,L)
                         SSWZ = MU_sL_pT*(W_S(IJKP,L)-W_S(IJK,L))*AXY_W(IJK)*ODZ(KP)*OX(I)&
                              -MU_sL_pB*(W_S(IJK,L)-W_S(IJKM,L))*AXY_W(IJKM)*ODZ(K)*OX(I)

! bulk viscosity term
                         XI_sL_pT = XI_sL_ip(IJKT,M,L)
                         XI_sL_pB = XI_sL_ip(IJK,M,L)
                         LAMBDA_sL_pT = -(2d0/3d0)*MU_sL_pT + XI_sL_pT
                         LAMBDA_sL_pB = -(2d0/3d0)*MU_sL_pB + XI_sL_pB
                         SSBV = (LAMBDA_sL_pT*TRD_S(IJKT,L)-LAMBDA_sL_pB*TRD_S(IJK,L))*AXY(IJK)

! off diagonal shear stress terms
                         SSX = MU_sL_pE*(U_S(IJKP,L)-U_S(IJK,L))*OX_E(I)*AYZ_W(IJK)*ODZ_T(K)&
                              -MU_sL_pW*(U_S(IMJKP,L)-U_S(IMJK,L))*(DY(J)*HALF*(DZ(K)+&
                              DZ(KP)))*ODZ_T(K)
                         !same as oX_E(IM)*AYZ_W(IMJK), but avoids singularity

                         SSY = MU_sL_pN*(V_S(IJKP,L)-V_S(IJK,L))*AXZ_W(IJK)*ODZ_T(K)*OX(I)&
                              -MU_sL_pS*(V_S(IJMKP,L)-V_S(IJMK,L))*AXZ_W(IJMK)*ODZ_T(K)*OX(I)

                         SSZ = SSWZ

! special terms for cylindrical coordinates
                         IF (CYLINDRICAL) THEN
! tau_zz term: modify Ssz: (1/x) (d/dz) (tau_zz)
!                             integral of (1/x) (d/dz) (2mu*(u/x))
!                             (normally part of tau_w_s) - explicit
                              SSZ = SSZ +(&
                                   (MU_sL_pT*(U_S(IJKP,L)+U_S(IMJKP,L))*OX(I)*AXY_W(IJK))-&
                                   (MU_sL_pB*(U_S(IJK,L)+U_S(IMJK,L))*OX(I)*AXY_W(IJKM)))

                              MU_sL_p = AVG_Z(MU_sL_ip(IJK,M,L),MU_sL_ip(IJKT,M,L),K)
! tau_xz/x terms: (tau_xz/x)
!                             integral of (1/x)*(mu/x)*du/dz
!                             (normally part of tau_w_s) - explicit
                              IF (OX_E(IM) /= UNDEFINED) THEN
                                   DUODZ = (U_S(IMJKP,L)-U_S(IMJK,L))*OX_E(IM)*ODZ_T(K)
                              ELSE
                                   DUODZ = ZERO
                              ENDIF

                              tauxz_ox = MU_sL_p*OX(I)*HALF*(&
                                   (OX_E(I)*(U_S(IJKP,L)-U_S(IJK,L))*ODZ_T(K))+&
                                   DUODZ)

!                             integral of (1/x)*mu*dw/dx
!                             (normally part of source_w_s)
                              tauxz_ox = tauxz_ox + (MU_sL_p*OX(I)*HALF*&
                                   ( (W_S(IPJK,L)-W_S(IJK,L))*ODX_E(I) + &
                                   (W_S(IJK,L)-W_S(IMJK,L))*ODX_E(IM) ) )

!                             integral of (1/x)*(-mu/x)*w
!                             (normally part of source_w_s)
                              tauxz_ox = tauxz_ox - (OX(I)*OX(I)*MU_sL_p*W_S(IJK,L))

!                             multiply all tau_xz/x terms by volume:
                              tauxz_ox = tauxz_ox*VOL_W(IJK)

! x*tau_xz term: (1/x) (d/dz) (x*tau_xz)
!                             integral of (1/x)*d( (-x)*(mu/x)*w )/dx
!                             (normally part of source_w_s)
                              tauxz_x = -( ( MU_sL_pE*OX_E(I)*HALF*(W_S(IPJK,L)+&
                                   W_S(IJK,L))*AYZ_W(IJK) )-( MU_sL_pW*HALF*&
                                   (W_S(IJK,L)+W_S(IMJK,L))*&
                                   (DY(J)*HALF*(DZ(K)+DZ(KP)) ) ) )
                                   !same as oX_E(IM)*AYZ_W(IMJK), but avoids singularity
                         ELSE
                              tauxz_ox = ZERO
                              tauxz_x = ZERO
                         ENDIF
!--------------------- End Sources from Stress Term ---------------------


!--------------------- Sources from Momentum Source Term ---------------------
! Momentum source associated with the difference in the gradients in
! number density of solids phase m and all other solids phases
                         D_PL = D_P(IJK,L)
                         M_PL = (Pi/6d0)*(D_PL**3.)*RO_S(IJK,L)

                         NU_PM_pT = ROP_S(IJKT,M)/M_PM
                         NU_PM_pB = ROP_S(IJK,M)/M_PM
                         NU_PM_p = AVG_Z(NU_PM_pB,NU_PM_pT,K)

                         NU_PL_pT = ROP_S(IJKT,L)/M_PL
                         NU_PL_pB = ROP_S(IJK,L)/M_PL
                         NU_PL_p = AVG_Z(NU_PL_pB,NU_PL_pT,K)

                         Fnu_s_p = AVG_Z(Fnu_s_ip(IJK,M,L),Fnu_s_ip(IJKT,M,L),K)
                         DS1 = Fnu_s_p*NU_PL_p*(NU_PM_pT-NU_PM_pB)*OX(I)*ODZ_T(K)
                         DS2 = -Fnu_s_p*NU_PM_p*(NU_PL_pT-NU_PL_pB)*OX(I)*ODZ_T(K)
                         DS1plusDS2 = DS1 + DS2

! Momentum source associated with the gradient in granular
! temperature of species M
                         T_PM_pT = Theta_M(IJKT,M)
                         T_PM_pB = Theta_M(IJK,M)

                         FT_sM_p = AVG_Z(FT_sM_ip(IJK,M,L),FT_sM_ip(IJKT,M,L),K)
                         DS3 = FT_sM_p*(T_PM_pT-T_PM_pB)*OX(I)*ODZ_T(K)

! Momentum source associated with the gradient in granular
! temperature of species L
                         T_PL_pT = Theta_M(IJKT,L)
                         T_PL_pB = Theta_M(IJK,L)

                         FT_sL_p = AVG_Z(FT_sL_ip(IJK,M,L),FT_sL_ip(IJKT,M,L),K)
                         DS4 = FT_sL_p*(T_PL_pT-T_PL_pB)*OX(I)*ODZ_T(K)
!--------------------- End Sources from Momentum Source Term ---------------------


! Add the terms
                    STRESS_TERMS = STRESS_TERMS + SSWX + SSWY + SSWZ + &
                        SSBV + SSX + SSY + SSZ + tauxz_ox + tauxz_x
                    DRAG_TERMS = DRAG_TERMS + (DS1plusDS2+DS3+DS4)*VOL_W(IJK)

                    ELSE ! if m .ne. L
! for m = l all stress terms should already be handled in existing routines
! for m = l all drag terms should become zero
                         STRESS_TERMS = STRESS_TERMS + ZERO
                         DRAG_TERMS = DRAG_TERMS + ZERO
                    ENDIF ! if m .ne. L
               ENDDO     ! over L

               KTMOM_W_S(IJK,M) = STRESS_TERMS + DRAG_TERMS
          ELSE
               KTMOM_W_S(IJK,M) = ZERO
          ENDIF     ! dilute
      ENDDO     ! over ijk

      RETURN
      END SUBROUTINE CALC_IA_MOMSOURCE_W_S
