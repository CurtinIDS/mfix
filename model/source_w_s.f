!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_W_s                                              C
!  Purpose: Determine source terms for W_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Allow for partial-slip boundary conditions proposed by     C
!           by Johnson & Jackson (1987) if the Granular Temperature    C
!           equation is used.                                          C
!  Author: K. Agrawal, Princeton University           Date: 24-JAN-98  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_W_S(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE tau_g, only: ctau_w_g
      USE bc
      USE compar
      USE sendrecv
      use kintheory
      USE ghdtheory
      USE drag
      USE cutcell
      USE quadric
      USE mms
      USE bodyforce
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IJKT, IMJK, IJMK, IJKM, IJKP, IMJKP
      INTEGER :: IJKE, IJKW, IJKTE, IJKTW, IM, IPJK, IJPK, IJMKP
! Phase index
      INTEGER ::  M, MM, L
! Internal surface
      INTEGER :: ISV
! Pressure at top cell
      DOUBLE PRECISION :: PgT
! Average volume fraction
      DOUBLE PRECISION :: EPSA, EPStmp, epse, epsw, epsn, epss, &
                          epst, epsb, epsMix, epsMixT
      DOUBLE PRECISION :: SUM_EPS_CP
! Average density
      DOUBLE PRECISION :: ROPSA
! Average density differences
      DOUBLE PRECISION :: dro1, dro2, droa
! Average quantities
      DOUBLE PRECISION :: ugt, Cte, Ctw, cpe, cpw, MUoX
! Source terms (Surface)
      DOUBLE PRECISION :: Sdp, Sdps
! Source terms (Volumetric)
      DOUBLE PRECISION :: V0, Vmt, Vbf, Vcoa, Vcob, Vmttmp, vxza
! Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION :: Ghd_drag, avgRop
! Source terms for HYS drag relation
      DOUBLE PRECISION :: HYS_drag, avgDrag
! virtual (added) mass
      DOUBLE PRECISION :: F_vir, ROP_MA, Uge, Ugw, Vgb, Vgt, Wge, &
                          Wgw, Wgn, Wgs, Wgt, Wgb
!-----------------------------------------------

      DO M = 1, MMAX
        IF(KT_TYPE_ENUM /= GHD_2007 .OR. &
           (KT_TYPE_ENUM == GHD_2007 .AND. M==MMAX)) THEN

          IF (MOMENTUM_Z_EQ(M)) THEN


!$omp  parallel do default(shared)                                   &
!$omp  private(I, J, K, IJK, IMJK, IJMK, IJKM, IJKP, IPJK, IMJKP,    &
!$omp         IJPK, IJMKP, IJKT, IJKE, IJKW, IJKTE, IJKTW, IM,       &
!$omp         ISV, EPStmp, epsMix, epsMixT, EPSA, EPSE, EPSW, EPSN,  &
!$omp         EPSS, EPST, EPSB, SUM_EPS_CP, PGT, SDP, SDPS, ROPSA,   &
!$omp         V0, ROP_MA, MM, L, F_vir,                              &
!$omp         Uge, Ugw, Vgb, Vgt, Wge, Wgw, Wgn, Wgs, Wgt, Wgb,      &
!$omp         Ugt, VMTtmp, VMT, DRO1, DRO2, DROA, VBF, Ghd_drag,     &
!$omp         avgRop, HYS_drag, avgDrag, VCOA, VCOB, CTE, CTW,       &
!$omp         CPE, CPW, VXZA, MUOX)
            DO IJK = ijkstart3, ijkend3

! Skip walls where some values are undefined.
                IF(WALL_AT(IJK)) cycle

                I = I_OF(IJK)
                J = J_OF(IJK)
                K = K_OF(IJK)
                IMJK = IM_OF(IJK)
                IJMK = JM_OF(IJK)
                IJKM = KM_OF(IJK)
                IJKP = KP_OF(IJK)
                IPJK = IP_OF(IJK)
                IMJKP = KP_OF(IMJK)
                IJPK = JP_OF(IJK)
                IJMKP = KP_OF(IJMK)
                IJKT = TOP_OF(IJK)

                IF (KT_TYPE_ENUM == GHD_2007) THEN
! with ghd theory, m = mmax
                   EPStmp = ZERO
                   epsMix = ZERO
                   epsMixT= ZERO
                   DO L = 1, SMAX
                      EPStmp = EPStmp + AVG_Z(EP_S(IJK,L),EP_S(IJKT,L),K)
                      epsMix  = epsMix  + EP_S(IJK,L) ! epsMix, epsMixT to be used for model B
                      epsMixT = epsMixT + EP_S(IJKT,L)
                      IF(IP_AT_T(IJK)) THEN
                         W_S(IJK,L) = ZERO
                      ELSEIF(SIP_AT_T(IJK)) THEN
                         ISV = IS_ID_AT_T(IJK)
                         W_S(IJK,L) = IS_VEL_S(ISV,L)
                      ENDIF
                   ENDDO
                   EPSA = EPStmp
                ELSE
                   EPSA = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
                ENDIF

! Impermeable internal surface
                IF (IP_AT_T(IJK)) THEN
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO

! Semi-permeable internal surface
                ELSEIF (SIP_AT_T(IJK)) THEN
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  ISV = IS_ID_AT_T(IJK)
                  B_M(IJK,M) = -IS_VEL_S(ISV,M)

! Dilute flow
                ELSEIF (EPSA <= DIL_EP_S) THEN
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
                  IF (KT_TYPE_ENUM == GHD_2007) THEN
                      EPSw = ZERO
                      EPSe = ZERO
                      EPSn = ZERO
                      EPSs = ZERO
                      EPSt = ZERO
                      EPSb = ZERO
                      DO L = 1, SMAX
                        EPSw = EPSw + EP_S(WEST_OF(IJK),L)
                        EPSe = EPSe + EP_S(EAST_OF(IJK),L)
                        EPSn = EPSn + EP_S(NORTH_OF(IJK),L)
                        EPSs = EPSs + EP_S(SOUTH_OF(IJK),L)
                        EPSt = EPSt + EP_S(TOP_OF(IJK),L)
                        EPSb = EPSb + EP_S(BOTTOM_OF(IJK),L)
                      ENDDO
                  ELSE
                      EPSw = EP_S(WEST_OF(IJK),M)
                      EPSe = EP_S(EAST_OF(IJK),M)
                      EPSn = EP_S(NORTH_OF(IJK),M)
                      EPSs = EP_S(SOUTH_OF(IJK),M)
                      EPSt = EP_S(TOP_OF(IJK),M)
                      EPSb = EP_S(BOTTOM_OF(IJK),M)
                  ENDIF
! using the average boundary cell values to compute U_s (sof, Aug 23 2005)
                  IF (EPSw > DIL_EP_S .AND. .NOT.IS_AT_E(IMJK)) A_M(IJK,west,M) = ONE
                  IF (EPSe > DIL_EP_S .AND. .NOT.IS_AT_E(IJK)) A_M(IJK,east,M) = ONE
                  IF (EPSs > DIL_EP_S .AND. .NOT.IS_AT_N(IJMK)) A_M(IJK,south,M) = ONE
                  IF (EPSn > DIL_EP_S .AND. .NOT.IS_AT_N(IJK)) A_M(IJK,north,M) = ONE
                  IF (EPSb > DIL_EP_S .AND. .NOT.IS_AT_T(IJKM)) A_M(IJK,bottom,M) = ONE
                  IF (EPSt > DIL_EP_S .AND. .NOT.IS_AT_T(IJK)) A_M(IJK,top,M) = ONE
                  IF((A_M(IJK,west,M)+A_M(IJK,east,M)+A_M(IJK,south,M)+A_M(IJK,north,M)+ &
                    A_M(IJK,bottom,M)+A_M(IJK,top,M)) == ZERO) THEN
                    B_M(IJK,M) = -W_S(IJK,M)
                  ELSE
                    A_M(IJK,0,M) = -(A_M(IJK,east,M)+A_M(IJK,west,M)+A_M(IJK,north,M)+ &
                                     A_M(IJK,south,M)+A_M(IJK,top,M)+A_M(IJK,bottom,M))
                  ENDIF

! Cartesian grid implementation
               ELSEIF (BLOCKED_W_CELL_AT(IJK)) THEN
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO

! Normal case
                ELSE

! Surface forces:
! Pressure term
                  PGT = P_G(IJKT)
                  IF (CYCLIC_Z_PD) THEN
! CYCLIC_AT_T Flag is not set correctly in DMP and causes issues. This
! is avoided by using the DMP cyclic map. The flags need fixed.
!                     IF (CYCLIC_AT_T(IJK)) PGT = P_G(IJKT) - DELP_Z
                     IF (KMAP(K_OF(IJK)).EQ.KMAX1) PGT = P_G(IJKT) - DELP_Z
                  ENDIF
                  IF (MODEL_B) THEN
                     SDP = ZERO
                  ELSE
                     IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                        SDP = -P_SCALE*EPSA*(PGT - P_G(IJK))*AXY(IJK)
                     ELSE
                        SDP = -P_SCALE*EPSA*(PGT * A_WPG_T(IJK) - P_G(IJK) * A_WPG_B(IJK) )
                     ENDIF
                  ENDIF

                  IF (CLOSE_PACKED(M)) THEN
                     IF(SMAX > 1 .AND. KT_TYPE_ENUM /= GHD_2007) THEN
                        SUM_EPS_CP=0.0
                        DO MM=1,SMAX
                          IF (CLOSE_PACKED(MM))&
                            SUM_EPS_CP=SUM_EPS_CP+AVG_Z(EP_S(IJK,MM),EP_S(IJKT,MM),K)
                        ENDDO
                        SUM_EPS_CP = Max(SUM_EPS_CP, small_number)
                        SDPS = -((P_S(IJKT,M)-P_S(IJK,M))+(EPSA/SUM_EPS_CP)*&
                           (P_STAR(IJKT)-P_STAR(IJK)))*AXY(IJK)
                     ELSE
                        IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                           SDPS =-((P_S(IJKT,M)-P_S(IJK,M))+(P_STAR(IJKT)-P_STAR(IJK)))*AXY(IJK)
                        ELSE
                           SDPS =-((P_S(IJKT,M)* A_WPG_T(IJK)-P_S(IJK,M)* A_WPG_B(IJK)) &
                                 +(P_STAR(IJKT)* A_WPG_T(IJK)-P_STAR(IJK)* A_WPG_B(IJK)))
                        ENDIF
                     ENDIF
                  ELSE
                     IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                        SDPS = -(P_S(IJKT,M)-P_S(IJK,M))*AXY(IJK)
                     ELSE
                        SDPS = -(P_S(IJKT,M) * A_WPG_T(IJK) - P_S(IJK,M) * A_WPG_B(IJK))
                     ENDIF
                  ENDIF


                  IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
! Volumetric forces
                     ROPSA = AVG_Z(ROP_S(IJK,M),ROP_S(IJKT,M),K)
! Previous time step
                     V0 = AVG_Z(ROP_SO(IJK,M),ROP_SO(IJKT,M),K)*ODT
! Added mass implicit transient term {Cv eps rop_g dW/dt}
                     IF(Added_Mass .AND. M == M_AM) THEN
                       ROP_MA = AVG_Z(ROP_g(IJK)*EP_s(IJK,M),&
                                      ROP_g(IJKT)*EP_s(IJKT,M),K)
                       V0 = V0 + Cv * ROP_MA * ODT
                     ENDIF
                  ELSE
! Volumetric forces
                     ROPSA = (VOL(IJK)*ROP_S(IJK,M) + &
                        VOL(IJKT)*ROP_S(IJKT,M))/(VOL(IJK) + VOL(IJKT))
! Previous time step
                     V0 = (VOL(IJK)*ROP_SO(IJK,M) + &
                        VOL(IJKT)*ROP_SO(IJKT,M))*ODT/&
                        (VOL(IJK) + VOL(IJKT))
! Added mass implicit transient term {Cv eps rop_g dW/dt}
                     IF(Added_Mass .AND. M == M_AM) THEN
                        ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M) +&
                           VOL(IJKT)*ROP_g(IJKT)*EP_s(IJKT,M))/&
                           (VOL(IJK) + VOL(IJKT))
                        V0 = V0 + Cv * ROP_MA * ODT
                     ENDIF
                  ENDIF

! VIRTUAL MASS SECTION (explicit terms)
! adding transient term dWg/dt to virtual mass term
                  F_vir = ZERO
                  IF(Added_Mass .AND. M == M_AM .AND.&
                     (.NOT.CUT_W_TREATMENT_AT(IJK))) THEN
                     F_vir = ( (W_g(IJK) - W_gO(IJK)) )*ODT*VOL_W(IJK)

! defining gas-particles velocity at momentum cell faces (or scalar cell center)
                      Wgb = AVG_Z_T(W_g(IJKM),W_g(IJK))
                      Wgt = AVG_Z_T(W_g(IJK),W_g(IJKP))
                      Uge = AVG_Z(U_g(IJK),U_g(IJKP),K)
                      Ugw = AVG_Z(U_g(IMJK),U_g(IMJKP),K)
                      Ugt = AVG_X_E(Ugw,Uge,IP1(I))
                      Wge = AVG_X(W_g(IJK),W_g(IPJK),IP1(I))
                      Wgw = AVG_X(W_g(IMJK),W_g(IJK),I)
                      Vgb = AVG_Y_N(V_g(IJMK),V_g(IJK))
                      Vgt = AVG_Y_N(V_g(IJMKP),V_g(IJKP))
                      Wgs = AVG_Y(W_g(IJMK),W_g(IJK),JM1(J))
                      Wgn = AVG_Y(W_g(IJK),W_g(IJPK),J)

! Coriolis force
                     IF (CYLINDRICAL) F_vir = F_vir + &
                        Ugt*W_g(IJK)*OX(I)

! adding convective terms (U dW/dx + V dW/dy + W dW/dz) to virtual mass.
                     F_vir = F_vir + W_g(IJK)*OX(I) * &
                        (Wgt - Wgb)*AXY(IJK) + Ugt*(Wge - Wgw)*AYZ(IJK) + &
                        AVG_Z(Vgb,Vgt,K) * (Wgn - Wgs)*AXZ(IJK)
                     F_vir = F_vir * Cv * ROP_MA
                  ENDIF

! Interphase mass transfer
                  IF (KT_TYPE_ENUM == GHD_2007) THEN
                     VMTtmp = ZERO
                     DO L = 1,SMAX
                        VMTtmp = VMTtmp + AVG_Z(SUM_R_S(IJK,L),SUM_R_S(IJKT,L),K)
                     ENDDO
                     VMT = VMTtmp
                  ELSE
                     IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                        VMT = AVG_Z(SUM_R_S(IJK,M),SUM_R_S(IJKT,M),K)
                     ELSE
                        VMT = (VOL(IJK)*SUM_R_S(IJK,M) + VOL(IJKT)*&
                           SUM_R_S(IJKT,M))/(VOL(IJK) + VOL(IJKT))
                     ENDIF
                  ENDIF

! Body force
                  IF (MODEL_B) THEN
                     IF (KT_TYPE_ENUM == GHD_2007) THEN
                       DRO1 = ROP_S(IJK,M)  - RO_G(IJK) *epsMix
                       DRO2 = ROP_S(IJKT,M) - RO_G(IJKT)*epsMixT
                       DROA = AVG_Z(DRO1,DRO2,K)
                       VBF = DROA*BFZ_S(IJK,M)
                     ELSE
                       DRO1 = (RO_S(IJK,M)-RO_G(IJK))*EP_S(IJK,M)
                       DRO2 = (RO_S(IJK,M)-RO_G(IJKT))*EP_S(IJKT,M)
                       DROA = AVG_Z(DRO1,DRO2,K)
                       VBF = DROA*BFZ_S(IJK,M)
                     ENDIF
                  ELSE
                     VBF = ROPSA*BFZ_S(IJK,M)
                  ENDIF

! Additional force for GHD from darg force sum(beta_ig * Joi/rhop_i)
                  Ghd_drag = ZERO
                  IF (KT_TYPE_ENUM == GHD_2007) THEN
                    DO L = 1,SMAX
                      avgRop = AVG_Z(ROP_S(IJK,L),ROP_S(IJKT,L),K)
                      if(avgRop > ZERO) Ghd_drag = Ghd_drag -&
                           AVG_Z(F_GS(IJK,L),F_GS(IJKT,L),K) * JoiZ(IJK,L) / avgRop
                    ENDDO
                  ENDIF

! Additional force for HYS drag force, do not use with mixture GHD theory
                  HYS_drag = ZERO
                  IF (DRAG_TYPE_ENUM .EQ. HYS .AND. &
                      KT_TYPE_ENUM /= GHD_2007) THEN
                     DO L = 1,MMAX
                        IF (L /= M) THEN
                           avgDrag = AVG_Z(beta_ij(IJK,M,L),beta_ij(IJKT,M,L),K)
                           HYS_drag = HYS_drag - avgDrag * (W_g(ijk) - W_s(IJK,L))
                        ENDIF
                     ENDDO
                  ENDIF

! Special terms for cylindrical coordinates
                  VCOA = ZERO
                  VCOB = ZERO
                  VXZA = ZERO
                  CTE  = ZERO
                  CTW  = ZERO
                  CPE = ZERO
                  CPW = ZERO
                  IF (CYLINDRICAL) THEN
! Coriolis force
                     IMJK = IM_OF(IJK)
                     IJKP = KP_OF(IJK)
                     IMJKP = KP_OF(IMJK)
                     UGT = AVG_Z(HALF*(U_S(IJK,M)+U_S(IMJK,M)),&
                                 HALF*(U_S(IJKP,M)+U_S(IMJKP,M)),K)
                     IF (UGT > ZERO) THEN
                        VCOA = ROPSA*UGT*OX(I)
                        VCOB = ZERO
                        IF(Added_Mass .AND. M == M_AM) &
                           VCOA = VCOA + Cv*ROP_MA*UGT*OX(I)
                     ELSE
                        VCOA = ZERO
                        VCOB = -ROPSA*UGT*W_S(IJK,M)*OX(I)
                        IF(Added_Mass .AND. M == M_AM) &
                           VCOB = VCOB - Cv*ROP_MA*UGT*W_S(IJK,M)*OX(I)
                     ENDIF

! Term from tau_xz: intergral of (1/x)*(d/dx)(x*mu*(-w/x))
                     IJKE = EAST_OF(IJK)
                     IJKW = WEST_OF(IJK)
                     IJKTE = TOP_OF(IJKE)
                     IJKTW = TOP_OF(IJKW)
                     IM = IM1(I)
                     IPJK = IP_OF(IJK)

                     CTE = HALF*AVG_Z_H(AVG_X_H(EPMU_S(IJK,M),EPMU_S(IJKE,M),I),&
                                        AVG_X_H(EPMU_S(IJKT,M),EPMU_S(IJKTE,M),I),K)*&
                           OX_E(I)*AYZ_W(IJK)
                     CTW = HALF*AVG_Z_H(AVG_X_H(EPMU_S(IJKW,M),EPMU_S(IJK,M),IM),&
                                        AVG_X_H(EPMU_S(IJKTW,M),EPMU_S(IJKT,M),IM),K)*&
                           DY(J)*(HALF*(DZ(K)+DZ(KP1(K))))
                        !same as oX_E(IM)*AYZ_W(IMJK), but avoids singularity

! (mu/x)*(dw/dx) part of tau_xz/x
                     MUOX = AVG_Z(EPMU_S(IJK,M),EPMU_S(IJKT,M),K)*OX(I)
                     CPE = MUOX* HALF*MUOX*ODX_E(I)*VOL_W(IJK)
                     CPW = MUOX* HALF*MUOX*ODX_E(IM)*VOL_W(IJK)

! -(mu/x)*(w/x) part of tau_xz/x
                     VXZA = MUOX*OX(I)
                  ENDIF

! Collect the terms
                  A_M(IJK,east,M) = A_M(IJK,east,M) + CPE
                  A_M(IJK,west,M) = A_M(IJK,west,M) - CPW

                  A_M(IJK,0,M) = -(A_M(IJK,east,M)+A_M(IJK,west,M)+&
                     A_M(IJK,north,M)+A_M(IJK,south,M)+A_M(IJK,top,M)+&
                     A_M(IJK,bottom,M)+(V0+ZMAX(VMT)+VCOA+VXZA)*&
                     VOL_W(IJK)+ CTE - CTW)

                  A_M(IJK,east,M) = A_M(IJK,east,M) - CTE
                  A_M(IJK,west,M) = A_M(IJK,west,M) + CTW

                  B_M(IJK,M) = B_m(IJK, M) - (SDP + SDPS + &
                     TAU_W_S(IJK,M) + epsa*cTAU_W_G(IJK) + F_vir + &
                     ( (V0+ZMAX((-VMT)))*W_SO(IJK,M)+ &
                     VBF + VCOB + HYS_drag)*VOL_W(IJK) )
! MMS Source term.
                  IF(USE_MMS) B_M(IJK,M) = &
                     B_M(IJK,M) - MMS_W_S_SRC(IJK)*VOL_W(IJK)

                  IF (KT_TYPE_ENUM == IA_2005) THEN
                    B_M(IJK,M) = B_M(IJK,M) - KTMOM_W_S(IJK,M)
                  ELSEIF (KT_TYPE_ENUM == GHD_2007) THEN
                    B_M(IJK,M) = B_M(IJK,M) - Ghd_drag*VOL_W(IJK)
                  ENDIF

                ENDIF   ! end branching on cell type (sip/ip/dilute/block/else branches)
            ENDDO   ! end do loop over ijk
!$omp end parallel do


! modifications for cartesian grid implementation
            IF(CARTESIAN_GRID) CALL CG_SOURCE_W_S(A_M, B_M, M)
! modifications for bc
            CALL SOURCE_W_S_BC (A_M, B_M, M)
            IF(CARTESIAN_GRID) CALL CG_SOURCE_W_S_BC(A_M, B_M, M)

          ENDIF   ! end if (momentum_z_eq)
        ENDIF   ! end if for ghd theory
      ENDDO   ! end do loop over mmax

      RETURN
      END SUBROUTINE SOURCE_W_S

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_W_s_BC                                           C
!  Purpose: Determine source terms for W_s momentum eq. The terms      C
!     appear in the center coefficient and RHS vector.  The center     C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Comments: see source_w_g_bc for more detailed in-code comments      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_W_S_BC(A_M, B_M, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE bc
      USE output
      USE compar
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Solids phase index
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
!  Boundary condition
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK,&
                 IM, JM, IJKB, IJKM, IJKP
!-----------------------------------------------

! Setting the default boundary conditions
      J1 = 1
      DO K1 = kmin3,kmax3
         DO I1 = imin3,imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
! Setting the wall velocity to zero
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = -ONE
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
! Setting the wall velocity equal to the adjacent fluid velocity
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ONE
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = -ONE
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ONE
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

      I1 = 1
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = -ONE
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ONE
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

      I1 = IMAX2
      DO K1 = kmin3,kmax3
         DO J1 = jmin3,jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = -ONE
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ONE
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO
! End setting the default boundary conditions


! Setting user specified boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

! Setting wall boundary conditions
            IF (BC_TYPE_ENUM(L) == NO_SLIP_WALL) THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               IF (BC_JJ_PS(L) == 0) THEN
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           A_M(IJK,east,M) = ZERO
                           A_M(IJK,west,M) = ZERO
                           A_M(IJK,north,M) = ZERO
                           A_M(IJK,south,M) = ZERO
                           A_M(IJK,top,M) = ZERO
                           A_M(IJK,bottom,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IF (FLUID_AT(EAST_OF(IJK))) THEN
                              A_M(IJK,east,M) = -ONE
                           ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                              A_M(IJK,west,M) = -ONE
                           ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                              A_M(IJK,north,M) = -ONE
                           ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                              A_M(IJK,south,M) = -ONE
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ELSE   ! Johnson and Jackson partial slip
                  CALL JJ_BC_W_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M)
               ENDIF

            ELSEIF (BC_TYPE_ENUM(L) == FREE_SLIP_WALL) THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               IF (BC_JJ_PS(L) == 0) THEN
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           A_M(IJK,east,M) = ZERO
                           A_M(IJK,west,M) = ZERO
                           A_M(IJK,north,M) = ZERO
                           A_M(IJK,south,M) = ZERO
                           A_M(IJK,top,M) = ZERO
                           A_M(IJK,bottom,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IF (FLUID_AT(EAST_OF(IJK))) THEN
                              A_M(IJK,east,M) = ONE
                           ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                              A_M(IJK,west,M) = ONE
                           ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                              A_M(IJK,north,M) = ONE
                           ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                              A_M(IJK,south,M) = ONE
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ELSE   ! Johnson and Jackson partial slip
                  CALL JJ_BC_W_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M)
               ENDIF

            ELSEIF (BC_TYPE_ENUM(L) == PAR_SLIP_WALL) THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               IF (BC_JJ_PS(L) == 0) THEN
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IF (.NOT.WALL_AT(IJK)) CYCLE  !skip redefined cells
                           IM = IM1(I)
                           JM = JM1(J)
                           A_M(IJK,east,M) = ZERO
                           A_M(IJK,west,M) = ZERO
                           A_M(IJK,north,M) = ZERO
                           A_M(IJK,south,M) = ZERO
                           A_M(IJK,top,M) = ZERO
                           A_M(IJK,bottom,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IF (FLUID_AT(EAST_OF(IJK))) THEN
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN
                                 A_M(IJK,east,M) = -HALF
                                 A_M(IJK,0,M) = -HALF
                                 B_M(IJK,M) = -BC_WW_S(L,M)
                              ELSE
                                 IF (CYLINDRICAL) THEN
                                    A_M(IJK,0,M) = -(HALF*(BC_HW_S(L,M)-OX_E(I)&
                                       )+ODX_E(I))
                                    A_M(IJK,east,M) = -(HALF*(BC_HW_S(L,M)-OX_E(I)&
                                       )-ODX_E(I))
                                 ELSE
                                    A_M(IJK,0,M)=-(HALF*BC_HW_S(L,M)+ODX_E(I))
                                    A_M(IJK,east,M)=-(HALF*BC_HW_S(L,M)-ODX_E(I))
                                 ENDIF
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_WW_S(L,M)
                              ENDIF
                           ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN
                                 A_M(IJK,west,M) = -HALF
                                 A_M(IJK,0,M) = -HALF
                                 B_M(IJK,M) = -BC_WW_S(L,M)
                              ELSE
                                 IF (CYLINDRICAL) THEN
                                    A_M(IJK,west,M) = -(HALF*(BC_HW_S(L,M)-OX_E(IM&
                                       ))-ODX_E(IM))
                                    A_M(IJK,0,M) = -(HALF*(BC_HW_S(L,M)-OX_E(IM&
                                       ))+ODX_E(IM))
                                 ELSE
                                    A_M(IJK,west,M) = -(HALF*BC_HW_S(L,M)-ODX_E(IM&
                                       ))
                                    A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODX_E(IM&
                                       ))
                                 ENDIF
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_WW_S(L,M)
                              ENDIF
                           ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN
                                 A_M(IJK,north,M) = -HALF
                                 A_M(IJK,0,M) = -HALF
                                 B_M(IJK,M) = -BC_WW_S(L,M)
                              ELSE
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(J))
                                 A_M(IJK,north,M) = -(HALF*BC_HW_S(L,M)-ODY_N(J))
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_WW_S(L,M)
                              ENDIF
                           ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                              IF (BC_HW_S(L,M) == UNDEFINED) THEN
                                 A_M(IJK,south,M) = -HALF
                                 A_M(IJK,0,M) = -HALF
                                 B_M(IJK,M) = -BC_WW_S(L,M)
                              ELSE
                                 A_M(IJK,south,M) = -(HALF*BC_HW_S(L,M)-ODY_N(JM))
                                 A_M(IJK,0,M) = -(HALF*BC_HW_S(L,M)+ODY_N(JM))
                                 B_M(IJK,M) = -BC_HW_S(L,M)*BC_WW_S(L,M)
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ELSE   !Johnson and Jackson partial slip
                  CALL JJ_BC_W_S (I1, I2, J1, J2, K1, K2, L, M, A_M, B_M)
               ENDIF

! Setting flow boundary conditions
            ELSEIF (BC_TYPE_ENUM(L)==P_INFLOW .OR. BC_TYPE_ENUM(L)==P_OUTFLOW) THEN
               IF (BC_PLANE(L) == 'B') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           A_M(IJK,east,M) = ZERO
                           A_M(IJK,west,M) = ZERO
                           A_M(IJK,north,M) = ZERO
                           A_M(IJK,south,M) = ZERO
                           A_M(IJK,top,M) = ZERO
                           A_M(IJK,bottom,M) = ONE
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF

            ELSEIF (BC_TYPE_ENUM(L) == OUTFLOW) THEN
               IF (BC_PLANE(L) == 'B') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           A_M(IJK,east,M) = ZERO
                           A_M(IJK,west,M) = ZERO
                           A_M(IJK,north,M) = ZERO
                           A_M(IJK,south,M) = ZERO
                           A_M(IJK,top,M) = ZERO
                           A_M(IJK,bottom,M) = ONE
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IJKM = KM_OF(IJK)
                           A_M(IJKM,east,M) = ZERO
                           A_M(IJKM,west,M) = ZERO
                           A_M(IJKM,north,M) = ZERO
                           A_M(IJKM,south,M) = ZERO
                           A_M(IJKM,top,M) = ZERO
                           A_M(IJKM,bottom,M) = ONE
                           A_M(IJKM,0,M) = -ONE
                           B_M(IJKM,M) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ELSEIF (BC_PLANE(L) == 'T') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                           IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                           IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IJKP = KP_OF(IJK)
                           A_M(IJKP,east,M) = ZERO
                           A_M(IJKP,west,M) = ZERO
                           A_M(IJKP,north,M) = ZERO
                           A_M(IJKP,south,M) = ZERO
                           A_M(IJKP,top,M) = ONE
                           A_M(IJKP,bottom,M) = ZERO
                           A_M(IJKP,0,M) = -ONE
                           B_M(IJKP,M) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF

! Setting bc that are not ns, fs, psw, p_inflow, p_outflow, or outflow
            ELSE
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = -W_S(IJK,M)
                        IF (BC_PLANE(L) == 'B') THEN
                           IJKB = BOTTOM_OF(IJK)
                           A_M(IJKB,east,M) = ZERO
                           A_M(IJKB,west,M) = ZERO
                           A_M(IJKB,north,M) = ZERO
                           A_M(IJKB,south,M) = ZERO
                           A_M(IJKB,top,M) = ZERO
                           A_M(IJKB,bottom,M) = ZERO
                           A_M(IJKB,0,M) = -ONE
                           B_M(IJKB,M) = -W_S(IJKB,M)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ENDIF   ! end if (bc_type)
         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc

      RETURN
      END SUBROUTINE SOURCE_W_S_BC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: JJ_BC_W_s                                               C
!  Purpose: Implement Johnson and Jackson boundary condition           C
!                                                                      C
!  Author: K. Agrawal, A. Srivastava,                 Date: 14-APR-98  C
!          Princeton University                                        C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE JJ_BC_W_S(I1, I2, J1, J2, K1, K2, L, M, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE calc_gr_boundary
      USE compar
      USE constant
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE is
      USE output
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run
      USE rxns
      USE scales
      USE tau_s
      USE toleranc
      USE visc_s

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Boundary condition
      INTEGER, INTENT(IN) :: L
! Indices
      INTEGER, INTENT(IN) :: I1, I2, J1, J2, K1, K2
! Solids phase index
      INTEGER, INTENT(IN) ::  M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, JM, IM, IJKP
! coefficients for granular bc
      DOUBLE PRECISION :: hw, gw, cw
!-----------------------------------------------

      DO K = K1, K2
         DO J = J1, J2
            DO I = I1, I2
               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I,J,K)
               IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
               IM = IM1(I)
               JM = JM1(J)
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO

               IF (FLUID_AT(EAST_OF(IJK))) THEN
                  IJKP = KP_OF(EAST_OF(IJK))
                  IF (WALL_AT(IJKP)) CYCLE
                  IF (EP_S(EAST_OF(IJK),M) <= DIL_EP_S) THEN
                     A_M(IJK,east,M) = ONE
                  ELSE
                     IF (BC_JJ_PS(L) == 1) THEN
                        CALL CALC_GRBDRY (IJK, EAST_OF(IJK), 'E', 'W',&
                           M, L, GW, HW, CW)
                     ELSEIF (BC_JJ_PS(L) == 2) THEN
                        GW = 0D0
                        HW = 1D0
                        CW = BC_WW_S(L,M)
                     ELSE
                        GW = 1D0
                        CW = 0D0
                        HW = 0D0
                     ENDIF
                     IF (CYLINDRICAL) THEN
                        A_M(IJK,east,M) = -(HALF*(HW - OX_E(I)*GW)-ODX_E(I)*GW)
                        A_M(IJK,0,M) = -(HALF*(HW - OX_E(I)*GW)+ODX_E(I)*GW)
                     ELSE
                        A_M(IJK,east,M) = -(HALF*HW - ODX_E(I)*GW)
                        A_M(IJK,0,M) = -(HALF*HW + ODX_E(I)*GW)
                     ENDIF
                     B_M(IJK,M) = -CW
                  ENDIF

               ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                  IJKP = KP_OF(WEST_OF(IJK))
                  IF (WALL_AT(IJKP)) CYCLE
                  IF (EP_S(WEST_OF(IJK),M) <= DIL_EP_S) THEN
                     A_M(IJK,west,M) = ONE
                  ELSE
                     IF (BC_JJ_PS(L) == 1) THEN
                        CALL CALC_GRBDRY (IJK, WEST_OF(IJK), 'W', 'W',&
                           M, L, GW, HW, CW)
                     ELSEIF (BC_JJ_PS(L) == 2) THEN
                        GW = 0D0
                        HW = 1D0
                        CW = BC_WW_S(L,M)
                     ELSE
                        GW = 1D0
                        CW = 0D0
                        HW = 0D0
                     ENDIF
                     IF (CYLINDRICAL) THEN
                        A_M(IJK,west,M) = -(HALF*(HW - OX_E(IM)*GW)-ODX_E(IM)*GW)
                        A_M(IJK,0,M) = -(HALF*(HW - OX_E(IM)*GW)+ODX_E(IM)*GW)
                     ELSE
                        A_M(IJK,west,M) = -(HALF*HW - ODX_E(IM)*GW)
                        A_M(IJK,0,M) = -(HALF*HW + ODX_E(IM)*GW)
                     ENDIF
                     B_M(IJK,M) = -CW
                  ENDIF

               ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                  IJKP = KP_OF(NORTH_OF(IJK))
                  IF (WALL_AT(IJKP)) CYCLE
                  IF (EP_S(NORTH_OF(IJK),M) <= DIL_EP_S) THEN
                     A_M(IJK,north,M) = ONE
                  ELSE
                     IF (BC_JJ_PS(L) == 1) THEN
                        CALL CALC_GRBDRY (IJK, NORTH_OF(IJK), 'N', 'W',&
                           M, L, GW, HW, CW)
                     ELSEIF (BC_JJ_PS(L) == 2) THEN
                        GW = 0D0
                        HW = 1D0
                        CW = BC_WW_S(L,M)
                     ELSE
                        GW = 1D0
                        CW = 0D0
                        HW = 0D0
                     ENDIF
                     A_M(IJK,north,M) = -(HALF*HW - ODY_N(J)*GW)
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(J)*GW)
                     B_M(IJK,M) = -CW
                  ENDIF

               ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                  IJKP = KP_OF(SOUTH_OF(IJK))
                  IF (WALL_AT(IJKP)) CYCLE
                  IF (EP_S(SOUTH_OF(IJK),M) <= DIL_EP_S) THEN
                     A_M(IJK,south,M) = ONE
                  ELSE
                     IF (BC_JJ_PS(L) == 1) THEN
                        CALL CALC_GRBDRY (IJK, SOUTH_OF(IJK), 'S', 'W',&
                           M, L, GW, HW, CW)
                     ELSEIF (BC_JJ_PS(L) == 2) THEN
                        GW = 0D0
                        HW = 1D0
                        CW = BC_WW_S(L,M)
                     ELSE
                        GW = 1D0
                        CW = 0D0
                        HW = 0D0
                     ENDIF
                     A_M(IJK,south,M) = -(HALF*HW - ODY_N(JM)*GW)
                     A_M(IJK,0,M) = -(HALF*HW + ODY_N(JM)*GW)
                     B_M(IJK,M) = -CW
                  ENDIF
               ENDIF

            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE JJ_BC_W_S

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_W_S                                        C
!  Purpose: Adds point sources to the solids W-Momentum equation.      C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_W_S(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar
      use constant
      use fldvar
      use geometry
      use indices
      use param
      use param1, only: one, small_number, zero
      use physprop
      use ps
      use run
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV, M
      INTEGER :: lKT, lKB
! terms of bm expression
      DOUBLE PRECISION :: pSource
!-----------------------------------------------
      do M=1, MMAX

      PS_LP: do PSV = 1, DIMENSION_PS
         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(abs(PS_W_s(PSV,M)) < small_number) cycle PS_LP

         if(PS_W_s(PSV,M) < ZERO) then
            lKB = PS_K_B(PSV)-1
            lKT = PS_K_T(PSV)-1
         else
            lKB = PS_K_B(PSV)
            lKT = PS_K_T(PSV)
         endif

         do k = lKB, lKT
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

            ijk = funijk(i,j,k)
            if(.NOT.fluid_at(ijk)) cycle


            if(A_M(IJK,0,M) == -ONE .AND.                           &
               B_M(IJK,M) == -W_s(IJK,M)) then
               B_M(IJK,M) = -PS_W_s(PSV,M) * PS_VEL_MAG_S(PSV,M)
            else
               pSource = PS_MASSFLOW_S(PSV,M) *                     &
                  (VOL(IJK)/PS_VOLUME(PSV))

               B_M(IJK,M) = B_M(IJK,M) - pSource *                  &
                  PS_W_s(PSV,M) * PS_VEL_MAG_S(PSV,M)
            endif

         enddo
         enddo
         enddo

         enddo PS_LP

      enddo ! do M=1, MMAX

      RETURN
      END SUBROUTINE POINT_SOURCE_W_S
