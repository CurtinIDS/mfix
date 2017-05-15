!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SOURCE_GHD_GRANULAR_ENERGY(sourcelhs,sourcerhs,IJK,IER)
!  Purpose: Calculate the source terms in the granular energy equation C
!           for GHD theory                                             C
!                                                                      C
!  Author: S. Benyahia, J. Galvin, C. Hrenya        Date: 22-JAN-09    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!     C. Hrenya handnotes and Garzo, Hrenya, Dufty papers (PRE, 2007)  C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!     Local variables: sourcelhs, sourcerhs                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_GHD_GRANULAR_ENERGY(SOURCELHS, SOURCERHS, IJK)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE run
      USE drag
      USE geometry
      USE fldvar
      USE visc_s
      USE ghdtheory
      USE trace
      USE indices
      USE constant
      USE toleranc
      USE compar        !//d
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Indices
      INTEGER          IJK, I, J, K, IM, JM, KM, IMJK, IJMK, IJKM,&
                       IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
!
!                      number densities
      DOUBLE PRECISION NiE, NiW, NiN, NiS, NiT,&
                       NiB, Nip, Ntotal
!
!                      mass of particles & non-zero gran. temp.
      DOUBLE PRECISION Mi, NonZeroTheta
!
!                      Thermal mobility-related source terms
      DOUBLE PRECISION ThermMobilityX, ThermMobilityY, ThermMobilityZ
!
!                      del.Joi and Fi.Joi terms
      DOUBLE PRECISION DelDotJoi, FiDotJoi, JoiXC, JoiYC, JoiZC
      DOUBLE PRECISION FiXC, FiYC, FiZC,ni(smax), SIGMAI(smax)
!
!                      phase index
      INTEGER          M, L
!
!                      Dufour-related source terms
      DOUBLE PRECISION DufourX, DufourY, DufourZ, DijQTerm, &
                       DijQTermE, DijQTermW, DijQTermN, DijQTermS, DijQTermT, DijQTermB,&
                       LijTermW,LijTermE,LijTermN,LijTermS,LijTermT,LijTermB

      DOUBLE PRECISION DijQTermE_H,DijQTermE_A,DijQTermW_H,DijQTermW_A,DijQTermN_H,DijQTermN_A,DijQTermS_H,&
                       DijQTermS_A,DijQTermT_H,DijQTermT_A,DijQTermB_H,DijQTermB_A,LijTermE_H,LijTermE_A, &
                       LijTermW_H,LijTermW_A,LijTermN_H,LijTermN_A,LijTermS_H,LijTermS_A,LijTermT_H, &
                       LijTermT_A,LijTermB_H,LijTermB_A
!
!                      Source terms to be kept on RHS
      DOUBLE PRECISION SOURCERHS, PressureRhs, ShearProduction, BulkViscRhs, DissDivURhs, phi_tot, &
                       SOURCE_FLUID,SINK_FLUID
!    DOUBLE PRECISION chi_ij
!
!                      Source terms to be kept on LHS
      DOUBLE PRECISION SOURCELHS, PressureLhs, CollDissipation, BulkViscLhs, DissDivULhs,VSLIP
      DOUBLE PRECISION UGC, USCM, VGC, VSCM, WGC, WSCM

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

      NonZeroTheta = MAX(THETA_M(IJK,MMAX), SMALL_NUMBER)

      Ntotal = ZERO
      phi_tot = ZERO
      DO M = 1,SMAX
          Ntotal = Ntotal + ROP_S(IJK,M)*6.d0/(Pi*D_P(IJK,M)**3 * RO_S(IJK,M))
          ni(M) = ROP_S(IJK,M)*6.d0/(Pi*D_P(IJK,M)**3 * RO_S(IJK,M))
          phi_tot = phi_tot + ROP_S(IJK,M)/RO_S(IJK,M)
          SIGMAI(M) = D_P(IJK,M)
      ENDDO

! Production by shear: (S:grad(v))
!     p_s*tr(D)
      PressureRhs = P_S_C(IJK,MMAX)*ZMAX(( -TRD_S_C(IJK,MMAX) ))
      PressureLhs = P_S_C(IJK,MMAX)*ZMAX((  TRD_S_C(IJK,MMAX) ))

!     mu_s*tr(D^2)
      ShearProduction = 2d0*MU_S_C(IJK,MMAX)*TRD_S2(IJK,MMAX)

!     lambda_s*tr(D)^2
      BulkViscRhs = ZMAX(  LAMBDA_S_C(IJK,MMAX) ) * TRD_S_C(IJK,MMAX)**2
      BulkViscLhs = ZMAX( -LAMBDA_S_C(IJK,MMAX) ) * TRD_S_C(IJK,MMAX)**2


! Energy dissipation by collisions: (3/2)*n*T*zeta
!     (3/2)*n*T*zeta0; zeroth order cooling rate term
      CollDissipation = 1.5d0*Ntotal*Zeta0(IJK)
!     (3/2)*n*T*zetaU*div(U) :
      DissDivURhs = 1.5d0*Ntotal*Theta_m(IJK,MMAX)* ZMAX( -ZetaU(IJK)*TRD_S_C(IJK,MMAX) )
      DissDivULhs = 1.5d0*Ntotal* ZMAX(  ZetaU(IJK)*TRD_S_C(IJK,MMAX) )


      DufourX = ZERO
      DufourY = ZERO
      DufourZ = ZERO
      ThermMobilityX = ZERO
      ThermMobilityY = ZERO
      ThermMobilityZ = ZERO
      DelDotJoi = ZERO
      FiDotJoi  = ZERO
      SOURCE_FLUID = ZERO
      SINK_FLUID = ZERO
      DO M = 1,SMAX

! Part of heat flux: div (q)
!     Sum_ij [ div( T^2*DijQ/nj*grad(nj)) ] -> Dufour term

          DO L = 1,SMAX
              Mi = (Pi/6.d0)*D_P(IJK,L)**3 * RO_S(IJK,L)

              Nip = ROP_S(IJK,L) /Mi
              NiE = ROP_S(IJKE,L)/Mi
              NiW = ROP_S(IJKW,L)/Mi
              NiN = ROP_S(IJKN,L)/Mi
              NiS = ROP_S(IJKS,L)/Mi

              DijQTerm  = zero
              DijQTermE = zero
              DijQTermW = zero
              DijQTermN = zero
              DijQTermS = zero

              if(ROP_S(IJK,L)/RO_S(IJK,L) > DIL_EP_S) DijQTerm = Theta_m(IJK,MMAX)**2*DijQ(IJK,M,L) / Nip
              if(ROP_S(IJKE,L)/RO_S(IJKE,L) > DIL_EP_S) DijQTermE =Theta_m(IJKE,MMAX)**2*DijQ(IJKE,M,L) / NiE
              if(ROP_S(IJKW,L)/RO_S(IJKW,L) > DIL_EP_S) DijQTermW =Theta_m(IJKW,MMAX)**2*DijQ(IJKW,M,L) / NiW
              if(ROP_S(IJKN,L)/RO_S(IJKN,L) > DIL_EP_S) DijQTermN =Theta_m(IJKN,MMAX)**2*DijQ(IJKN,M,L) / NiN
              if(ROP_S(IJKS,L)/RO_S(IJKS,L) > DIL_EP_S) DijQTermS =Theta_m(IJKS,MMAX)**2*DijQ(IJKS,M,L) / NiS

                DijQTermE_H = AVG_X_S(DijQTerm, DijQTermE, I)
                DijQTermE_A = AVG_X(DijQTerm, DijQTermE, I)
              IF(ABS(DijQTermE_H) .le. ABS(DijQTermE_A))THEN
                DijQTermE = DijQTermE_H
              ELSE
                DijQTermE = DijQTermE_A
              ENDIF

                DijQTermW_H = AVG_X_S(DijQTermW, DijQTerm, IM)
                DijQTermW_A = AVG_X(DijQTermW, DijQTerm, IM)
              IF(ABS(DijQTermW_H) .le. ABS(DijQTermW_A))THEN
                DijQTermW = DijQTermW_H
              ELSE
                DijQTermW = DijQTermW_A
              ENDIF

                DijQTermN_H = AVG_Y_S(DijQTerm, DijQTermN, J)
                DijQTermN_A = AVG_Y(DijQTerm, DijQTermN, J)
              IF(ABS(DijQTermN_H) .le. ABS(DijQTermN_A))THEN
                DijQTermN = DijQTermN_H
              ELSE
                DijQTermN = DijQTermN_A
              ENDIF

                DijQTermS_H = AVG_Y_S(DijQTermS, DijQTerm, JM)
                DijQTermS_A = AVG_Y(DijQTermS, DijQTerm, JM)
              IF(ABS(DijQTermS_H) .le. ABS(DijQTermS_A))THEN
                DijQTermS = DijQTermS_H
              ELSE
                DijQTermS = DijQTermS_A
              ENDIF

              DufourX = DufourX + ( DijQTermE*(NiE-Nip)*ODX_E(I) - &
                                    DijQTermW*(Nip-NiW)*ODX_E(IM) )*AYZ(IJK)
              DufourY = DufourY + ( DijQTermN*(NiN-Nip)*ODY_N(J) - &
                                    DijQTermS*(Nip-NiS)*ODY_N(JM) )*AXZ(IJK)

              IF(.NOT. NO_K) THEN
                 NiT = ROP_S(IJKT,L)/Mi
                 NiB = ROP_S(IJKB,L)/Mi

                 DijQTermT = zero
                 DijQTermB = zero

                 if(ROP_S(IJKT,L)/RO_S(IJKT,L) > DIL_EP_S) DijQTermT = Theta_m(IJK,MMAX)**2*DijQ(IJK,M,L) / NiT
                 if(ROP_S(IJKB,L)/RO_S(IJKB,L) > DIL_EP_S) DijQTermB = Theta_m(IJK,MMAX)**2*DijQ(IJK,M,L) / NiB

                 DijQTermT_H = AVG_Z_S(DijQTerm , DijQTermT, K)
                 DijQTermT_A = AVG_Z(DijQTerm , DijQTermT, K)

                 IF(ABS(DijQTermT_H) .le. ABS(DijQTermT_A))THEN
                   DijQTermT=DijQTermT_H
                 ELSE
                   DijQTermT=DijQTermT_A
                 ENDIF

                 DijQTermB_H = AVG_Z_S(DijQTermB, DijQTerm, KM)
                 DijQTermB_A = AVG_Z(DijQTermB, DijQTerm, KM)

                 IF(ABS(DijQTermB_H) .le. ABS(DijQTermB_A))THEN
                   DijQTermB=DijQTermB_H
                 ELSE
                   DijQTermB=DijQTermB_A
                 ENDIF

                 DufourZ = DufourZ + ( DijQTermT*(NiT-Nip)*ODZ_T(K)*OX(I) - &
                                       DijQTermB*(Nip-NiB)*ODZ_T(KM)*OX(I) )*AXY(IJK)
              ENDIF

!     Sum_ij [ div( Lij*Fj) ]; thermal mobility term
!     Where Fj = Body Force

                LijTermW_H = AVG_X_S(Lij(IJKW,M,L),Lij(IJK,M,L),IM)
                LijTermW_A = AVG_X(Lij(IJKW,M,L),Lij(IJK,M,L),IM)
              IF(ABS(LijTermW_H) .le. ABS(LijTermW_A))THEN
                LijTermW = LijTermW_H
              ELSE
                LijTermW = LijTermW_A
              ENDIF

                LijTermE_H = AVG_X_S(Lij(IJK,M,L),Lij(IJKE,M,L),I)
                LijTermE_A = AVG_X(Lij(IJK,M,L),Lij(IJKE,M,L),I)

              IF(ABS(LijTermE_H) .le. ABS(LijTermE_A))THEN
                LijTermE = LijTermE_H
              ELSE
                LijTermE = LijTermE_A
              ENDIF

                LijTermN_H = AVG_Y_S(Lij(IJK,M,L),Lij(IJKN,M,L),J)
                LijTermN_A = AVG_Y(Lij(IJK,M,L),Lij(IJKN,M,L),J)
              IF(ABS(LijTermN_H) .le. ABS(LijTermN_A))THEN
                LijTermN = LijTermN_H
              ELSE
                LijTermN = LijTermN_A
              ENDIF

                LijTermS_H = AVG_Y_S(Lij(IJKS,M,L),Lij(IJK,M,L),JM)
                LijTermS_A = AVG_Y(Lij(IJKS,M,L),Lij(IJK,M,L),JM)
              IF(ABS(LijTermS_H) .le. ABS(LijTermS_A))THEN
                LijTermS = LijTermS_H
              ELSE
                LijTermS = LijTermS_A
              ENDIF

              ThermMobilityX = ThermMobilityX + ( &
                       FiX(IJK,L) *LijTermE - FiX(IMJK,L)*LijTermW )*AYZ(IJK)

              ThermMobilityY = ThermMobilityY + ( &
                       FiY(IJK,L) *LijTermN - FiY(IJMK,L)*LijTermS )*AXZ(IJK)

              IF(.NOT. NO_K) THEN

              LijTermT_H = AVG_Z_S(Lij(IJK,M,L),Lij(IJKT,M,L),K)
              LijTermT_A = AVG_Z(Lij(IJK,M,L),Lij(IJKT,M,L),K)
              IF(ABS(LijTermT_H) .le. ABS(LijTermT_A))THEN
                LijTermT = LijTermT_H
              ELSE
                LijTermT = LijTermT_A
              ENDIF

              LijTermB_H = AVG_Z_S(Lij(IJKB,M,L),Lij(IJK,M,L),KM)
              LijTermB_A = AVG_Z(Lij(IJKB,M,L),Lij(IJK,M,L),KM)
              IF(ABS(LijTermB_H) .le. ABS(LijTermB_A))THEN
                LijTermB = LijTermB_H
              ELSE
                LijTermB = LijTermB_A
              ENDIF

              ThermMobilityZ = ThermMobilityZ + ( &
                       FiZ(IJK,L) *LijTermT  - FiZ(IJKM,L)*LijTermB )*AXY(IJK)
              ENDIF
          ENDDO ! for L = 1, smax

! Additional term arising from subtraction of 3/2*T*continuity
!     + (3/2)*T* Sum_i [ div (joi/mi) ]

          Mi = (Pi/6.d0)*D_P(IJK,M)**3 * RO_S(IJK,M)

          DelDotJoi = DelDotJoi + 1.5d0*THETA_M(IJK,MMAX)/Mi * ( &
                      (JoiX(IJK,M) - JoiX(IMJK,M))*AYZ(IJK) + &
                      (JoiY(IJK,M) - JoiY(IJMK,M))*AXZ(IJK) + &
                      (JoiZ(IJK,M) - JoiZ(IJKM,M))*AXY(IJK) )


! Species force dot species mass flux
!     Sum_i [ Fi dot joi/mi ]

!     Calculate species mass flux components at cell center
          JoiXC = AVG_X_E(JoiX(IMJK,M),JoiX(IJK,M),I)
          JoiYC = AVG_Y_N(JoiY(IJMK,M),JoiY(IJK,M))
          JoiZC = AVG_Z_T(JoiZ(IJKM,M),JoiZ(IJK,M))

! external forces evaluated at cell center
          FiXC = AVG_X_E(FiX(IMJK,M),FiX(IJK,M),I)
          FiYC = AVG_Y_N(FiY(IJMK,M),FiY(IJK,M))
          FiZC = AVG_Z_T(FiZ(IJKM,M),FiZ(IJK,M))

          FiDotJoi  = FiDotJoi  + ( JoiXC*FiXC + JoiYC*FiYC + JoiZC*FiZC ) / Mi


          IF(SWITCH > ZERO .AND. abs(RO_g0) .gt. ZERO) THEN ! do nothing for gran. flow
            UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
            VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
            WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
            USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
            VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M))
            WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

            VSLIP = (USCM-UGC)**2 + (VSCM-VGC)**2 + (WSCM-WGC)**2
            VSLIP = DSQRT(VSLIP)

! Source/Sink due to fluid do not work well with fluid-solids cases that we run
! uncomment the lines of code below to use (W. Holloway and S. Benyahia).
!
!            call chi_ij_GHD(smax,M,M,SIGMAi,phi_tot,ni,chi_ij)

!            SOURCE_FLUID = SOURCE_FLUID + (81D0*EP_S(IJK,M)*(MU_G(IJK)*VSLIP)**2D0/ &
!                   (chi_ij*(D_P(IJK,M)**3D0*RO_S(IJK,M)*THETA_M(IJK,M))**0.5D0))*VOL(IJK)

!            SINK_FLUID = SINK_FLUID + 3.d0*F_GS(IJK,M)*THETA_M(IJK,M)/Mi
          ENDIF

      ENDDO ! for M = 1, smax

      SINK_FLUID = SINK_FLUID/NonZeroTheta

      SOURCERHS = (PressureRhs + ShearProduction + BulkViscRhs + DissDivURhs)*VOL(IJK) &
                 + ZMAX(DufourX)+ZMAX(DufourY)+ZMAX(DufourZ) &
                 + ZMAX(ThermMobilityX)+ZMAX(ThermMobilityY)+ZMAX(ThermMobilityZ)      &
                 + ZMAX(DelDotJoi) + ZMAX(FiDotJoi)*VOL(IJK)

      SOURCERHS = SOURCERHS + SOURCE_FLUID

      SOURCELHS = ( (PressureLhs + BulkViscLhs)/NonZeroTheta   + &
                  (CollDissipation + DissDivULhs + SINK_FLUID) ) * VOL(IJK) + &
                  ( ZMAX(-DufourX)+ZMAX(-DufourY)+ZMAX(-DufourZ) + &
                   ZMAX(-ThermMobilityX)+ZMAX(-ThermMobilityY)+ZMAX(-ThermMobilityZ) + &
                   ZMAX(-DelDotJoi) + ZMAX(-FiDotJoi)*VOL(IJK) )/ NonZeroTheta

      RETURN
      END SUBROUTINE SOURCE_GHD_GRANULAR_ENERGY
