!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_W_g                                              C
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
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_W_G(A_M, B_M)

! Modules
!---------------------------------------------------------------------//
      USE bodyforce, only: bfz_g
      USE bc, only: delp_z

      USE compar, only: ijkstart3, ijkend3, kmap

      USE drag, only: f_gs, beta_ij
      USE fldvar, only: p_g, ro_g, rop_g, rop_go, rop_s
      USE fldvar, only: ep_g, ep_s, epg_jfac
      USE fldvar, only: u_g, w_g, w_go, u_s, v_s, w_s, w_so

      USE fun_avg, only: avg_x, avg_z, avg_y
      USE fun_avg, only: avg_x_h, avg_z_h
      USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      USE functions, only: ip_at_t, sip_at_t, is_id_at_t
      USE functions, only: ip_of, jp_of, kp_of, im_of, jm_of, km_of
      USE functions, only: east_of, west_of, top_of, bottom_of
      USE functions, only: zmax
      USE geometry, only: kmax1, cyclic_z_pd, cylindrical
      USE geometry, only: vol, vol_w
      USE geometry, only: axy, ayz, axz, ayz_w
      USE geometry, only: ox, ox_e, dy, dz, odx_e

      USE ghdtheory, only: joiz

      USE indices, only: i_of, j_of, k_of
      USE indices, only: ip1, im1, jm1, kp1
      USE is, only: is_pc

      USE mms, only: use_mms, mms_w_g_src
      USE param
      USE param1, only: zero, one, half
      USE physprop, only: mmax, smax
      USE physprop, only: mu_g, cv
      USE run, only: momentum_z_eq
      USE run, only: model_b, added_mass, m_am
      USE run, only: kt_type_enum, drag_type_enum
      USE run, only: ghd_2007, hys
      USE run, only: odt
      USE rxns, only: sum_r_g
      USE scales, only: p_scale
      USE tau_g, only: tau_w_g
      USE toleranc, only: dil_ep_s
      USE visc_g, only: epmu_gt
      USE cutcell, only: cartesian_grid, cut_w_treatment_at
      USE cutcell, only: blocked_w_cell_at
      USE cutcell, only: a_wpg_t, a_wpg_b
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK, IJKT, IMJK, IJKP, IMJKP,&
                 IJKE, IJKW, IJKTE, IJKTW, IM, IPJK,   &
                 IJKM, IJMK, IJMKP, IJPK
! Phase index
      INTEGER :: M, L, MM
! Internal surface
      INTEGER :: ISV
! Pressure at top cell
      DOUBLE PRECISION :: PgT
! Average volume fraction
      DOUBLE PRECISION :: EPGA, EPGAJ
! Average density
      DOUBLE PRECISION :: ROPGA, ROGA
! Average viscosity
      DOUBLE PRECISION :: MUGA
! Average coefficients
      DOUBLE PRECISION :: Cte, Ctw, MUoX, cpe, cpw
! Average U_g
      DOUBLE PRECISION Ugt
! Source terms (Surface)
      DOUBLE PRECISION Sdp
! Source terms (Volumetric)
      DOUBLE PRECISION V0, Vpm, Vmt, Vbf, Vcoa, Vcob, Vxza
! Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION Ghd_drag, avgRop
! Source terms for HYS drag relation
      DOUBLE PRECISION HYS_drag, avgDrag
! virtual (added) mass
      DOUBLE PRECISION :: ROP_MA, U_se, Usw, Ust, Vsb, Vst, &
                          Wse, Wsw, Wsn, Wss, Wst, Wsb
      DOUBLE PRECISION F_vir
! jackson terms: local stress tensor quantity
      DOUBLE PRECISION :: ltau_w_g
!---------------------------------------------------------------------//

! Set reference phase to gas
      M = 0

      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN

!$omp  parallel do default(shared)                                   &
!$omp  private(I, J, K, IJK, IJKT, IJKM, IJKP, IMJK, IPJK, IJMK,     &
!$omp          IMJKP, IJPK, IJMKP, IJKTE, IJKTW, IM, IJKW, IJKE,     &
!$omp          EPGA, epgaj, PGT, SDP, ROPGA, ROGA, V0, ISV, MUGA,    &
!$omp          vpm, Vmt, Vbf, F_vir, Ghd_drag, avgRop, HYS_drag,     &
!$omp          avgdrag, MM, L, VXZA, VCOA, VCOB, CTE, CTW, UGT,      &
!$omp          cpe, cpw, MUOX, ltau_w_g)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IJKT = TOP_OF(IJK)
         IJKM = KM_OF(IJK)
         IJKP = KP_OF(IJK)
         IMJK = IM_OF(IJK)
         IPJK = IP_OF(IJK)
         IMJKP = KP_OF(IMJK)
         IJMK = JM_OF(IJK)
         IJPK = JP_OF(IJK)
         IJMKP = KP_OF(IJMK)

         EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)
! if jackson then avg ep_g otherwise 1
         EPGAJ = AVG_Z(EPG_jfac(IJK),EPG_jfac(IJKT),K)

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

! Dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN
            A_M(IJK,east,M) = ZERO
            A_M(IJK,west,M) = ZERO
            A_M(IJK,north,M) = ZERO
            A_M(IJK,south,M) = ZERO
            A_M(IJK,top,M) = ZERO
            A_M(IJK,bottom,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO
! set velocity equal to that of bottom or top cell if solids are
! present in those cells else set velocity equal to known value
            IF (EP_G(BOTTOM_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,bottom,M) = ONE
            ELSE IF (EP_G(TOP_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,top,M) = ONE
            ELSE
               B_M(IJK,M) = -W_G(IJK)
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

! Surface forces

! Pressure term
            PGT = P_G(IJKT)
            IF (CYCLIC_Z_PD) THEN
               IF (KMAP(K_OF(IJK)).EQ.KMAX1) PGT = P_G(IJKT) - DELP_Z
            ENDIF
            IF (MODEL_B) THEN
               IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*(PGT - P_G(IJK))*AXY(IJK)
               ELSE
                   SDP = -P_SCALE*(PGT * A_WPG_T(IJK) - P_G(IJK) * A_WPG_B(IJK) )
               ENDIF
            ELSE
               IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*EPGA*(PGT - P_G(IJK))*AXY(IJK)
               ELSE
                   SDP = -P_SCALE*EPGA*(PGT * A_WPG_T(IJK) - P_G(IJK) * A_WPG_B(IJK) )
               ENDIF
            ENDIF

            IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
! Volumetric forces
               ROPGA = AVG_Z(ROP_G(IJK),ROP_G(IJKT),K)
               ROGA = AVG_Z(RO_G(IJK),RO_G(IJKT),K)

! Previous time step
               V0 = AVG_Z(ROP_GO(IJK),ROP_GO(IJKT),K)*ODT

! Added mass implicit transient term {Cv eps rop_g dW/dt}
               IF(Added_Mass) THEN
                  ROP_MA = AVG_Z(ROP_g(IJK)*EP_s(IJK,M_AM),&
                     ROP_g(IJKT)*EP_s(IJKT,M_AM),K)
                  V0 = V0 + Cv * ROP_MA * ODT
               ENDIF
            ELSE
! Volumetric forces
               ROPGA = (VOL(IJK)*ROP_G(IJK) + VOL(IJKT)*ROP_G(IJKT))/&
                  (VOL(IJK) + VOL(IJKT))
               ROGA  = (VOL(IJK)*RO_G(IJK) + VOL(IJKT)*RO_G(IJKT))/&
                  (VOL(IJK) + VOL(IJKT))
! Previous time step
               V0 = (VOL(IJK)*ROP_GO(IJK) + VOL(IJKT)*ROP_GO(IJKT))*&
                  ODT/(VOL(IJK) + VOL(IJKT))
! Added mass implicit transient term {Cv eps rop_g dW/dt}
               IF(Added_Mass) THEN
                  ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M_AM) + &
                     VOL(IJKT)*ROP_g(IJKT)*EP_s(IJKT,M_AM))/&
                     (VOL(IJK) + VOL(IJKT))
                  V0 = V0 + Cv * ROP_MA * ODT
               ENDIF
            ENDIF

! VIRTUAL MASS SECTION (explicit terms)
! adding transient term dWs/dt to virtual mass term
            F_vir = ZERO
            IF(Added_Mass.AND.(.NOT.CUT_W_TREATMENT_AT(IJK))) THEN
               F_vir = ( (W_s(IJK,M_AM) - W_sO(IJK,M_AM)) )*ODT*VOL_W(IJK)

! defining gas-particles velocity at momentum cell faces (or scalar cell center)
               Wsb = AVG_Z_T(W_S(IJKM,M_AM),W_s(IJK,M_AM))
               Wst = AVG_Z_T(W_s(IJK,M_AM),W_s(IJKP,M_AM))
               U_se = AVG_Z(U_s(IJK,M_AM),U_s(IJKP,M_AM),K)
               Usw = AVG_Z(U_s(IMJK,M_AM),U_s(IMJKP,M_AM),K)
               Ust = AVG_X_E(Usw,U_se,IP1(I))
               Wse = AVG_X(W_s(IJK,M_AM),W_s(IPJK,M_AM),IP1(I))
               Wsw = AVG_X(W_s(IMJK,M_AM),W_s(IJK,M_AM),I)
               Vsb = AVG_Y_N(V_s(IJMK,M_AM),V_s(IJK,M_AM))
               Vst = AVG_Y_N(V_s(IJMKP,M_AM),V_s(IJKP,M_AM))
               Wss = AVG_Y(W_s(IJMK,M_AM),W_s(IJK,M_AM),JM1(J))
               Wsn = AVG_Y(W_s(IJK,M_AM),W_s(IJPK,M_AM),J)

! adding convective terms (U dW/dx + V dW/dy + W dW/dz) to virtual mass.
               F_vir = F_vir + W_s(IJK,M_AM)*OX(I)* &
                  (Wst - Wsb)*AXY(IJK) + Ust*(Wse - Wsw)*AYZ(IJK) + &
                  AVG_Z(Vsb,Vst,K)*(Wsn - Wss)*AXZ(IJK)

! Coriolis force
               IF (CYLINDRICAL) F_vir = F_vir + &
                  Ust*W_s(IJK,M_AM)*OX(I)
               F_vir = F_vir * Cv * ROP_MA
            ENDIF

! pressure drop through porous media
            IF (SIP_AT_T(IJK)) THEN
               ISV = IS_ID_AT_T(IJK)
               MUGA = AVG_Z(MU_G(IJK),MU_G(IJKT),K)
               VPM = MUGA/IS_PC(ISV,1)
               IF (IS_PC(ISV,2) /= ZERO) VPM = VPM + &
                  HALF*IS_PC(ISV,2)*ROPGA*ABS(W_G(IJK))
            ELSE
               VPM = ZERO
            ENDIF

! Interphase mass transfer
            IF(.NOT.CUT_W_TREATMENT_AT(IJK)) THEN
               VMT = AVG_Z(SUM_R_G(IJK),SUM_R_G(IJKT),K)
            ELSE
               VMT = (VOL(IJK)*SUM_R_G(IJK) + VOL(IJKT)*SUM_R_G(IJKT))/&
                  (VOL(IJK) + VOL(IJKT))
            ENDIF

! Body force
            IF (MODEL_B) THEN
               VBF = ROGA*BFZ_G(IJK)
            ELSE
               VBF = ROPGA*BFZ_G(IJK)
            ENDIF

! Additional force for GHD from darg force sum(beta_ig * Joi/rhop_i)
            Ghd_drag = ZERO
            IF (KT_TYPE_ENUM .EQ. GHD_2007) THEN
               DO L = 1,SMAX
                  avgRop = AVG_Z(ROP_S(IJK,L),ROP_S(IJKT,L),K)
                  if(avgRop > ZERO) Ghd_drag = Ghd_drag +&
                     AVG_Z(F_GS(IJK,L),F_GS(IJKT,L),K) * JoiZ(IJK,L) / avgRop
               ENDDO
            ENDIF

! Additional force for HYS drag force, do not use with mixture GHD theory
            avgDrag = ZERO
            HYS_drag = ZERO
            IF (DRAG_TYPE_ENUM .EQ. HYS .AND. KT_TYPE_ENUM .NE. GHD_2007) THEN
               DO MM=1,MMAX
                  DO L = 1,MMAX
                     IF (L /= MM) THEN
                        avgDrag = AVG_Z(beta_ij(IJK,MM,L),beta_ij(IJKT,MM,L),K)
                        HYS_drag = HYS_drag + avgDrag * (W_g(IJK) - W_s(IJK,L))
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

! Special terms for cylindrical coordinates
            VCOA = ZERO
            VCOB = ZERO
            VXZA = ZERO
            CTE  = ZERO
            CTW  = ZERO
            CPE  = ZERO
            CPW  = ZERO
            IF (CYLINDRICAL) THEN
! Coriolis force
               IMJK = IM_OF(IJK)
               IJKP = KP_OF(IJK)
               IMJKP = KP_OF(IMJK)
               UGT = AVG_Z(HALF*(U_G(IJK)+U_G(IMJK)),&
                  HALF*(U_G(IJKP)+U_G(IMJKP)),K)
               IF (UGT > ZERO) THEN
                  VCOA = ROPGA*UGT*OX(I)
                  VCOB = ZERO
                  IF(Added_Mass) VCOA = VCOA + Cv*ROP_MA*UGT*OX(I)
               ELSE
                  VCOA = ZERO
                  VCOB = -ROPGA*UGT*W_G(IJK)*OX(I)
                  IF(Added_Mass) VCOB = VCOB - Cv*ROP_MA*UGT*W_G(IJK)*OX(I)
               ENDIF

! if ishii, then viscosity is multiplied by void fraction otherwise by 1
! part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently
! part of (tau_xz/x + 1/x d/dx (x tau_xz)) xdxdydz =>
!         1/x d/dx (x.mu.(-w/x)) xdxdydz =>
! delta (mu/x.(-w))Ayz |E-W : at (i+1/2 - i-1/2, j, k+1/2)
               IJKE = EAST_OF(IJK)
               IJKW = WEST_OF(IJK)
               IJKTE = TOP_OF(IJKE)
               IJKTW = TOP_OF(IJKW)
               IM = IM1(I)
               IPJK = IP_OF(IJK)
               CTE = HALF*AVG_Z_H(AVG_X_H(EPMU_GT(IJK),EPMU_GT(IJKE),I),&
                                  AVG_X_H(EPMU_GT(IJKT),EPMU_GT(IJKTE),I),K)*&
                     OX_E(I)*AYZ_W(IJK)
               CTW = HALF*AVG_Z_H(AVG_X_H(EPMU_GT(IJKW),EPMU_GT(IJK),IM),&
                                  AVG_X_H(EPMU_GT(IJKTW),EPMU_GT(IJKT),IM),K)*&
                     DY(J)*(HALF*(DZ(K)+DZ(KP1(K))))
! DY(J)*HALF(DZ(k)+DZ(kp)) = oX_E(IM)*AYZ_W(IMJK), but avoids singularity

! part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently
! part of (tau_xz/x + 1/x d/dx (x tau_xz)) xdxdydz =>
!         mu/x dw/dx xdxdydz =>
! delta (mu/x.(dw/dx))Vp |p : at (i, j, k+1/2)
               MUOX = AVG_Z(EPMU_GT(IJK),EPMU_GT(IJKT),K)*OX(I)
               CPE = MUOX*HALF*VOL_W(IJK)*ODX_E(I)
               CPW = MUOX*HALF*VOL_W(IJK)*ODX_E(IM)

! part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently
! part of (tau_xz/x + 1/x d/dx (x tau_xz)) xdxdydz =>
!         1/x d/dx (x.mu.(-w/x)) xdxdydz =>
! delta (mu/x.(-w/x))Vp |p : at (i, j, k+1/2)
               VXZA = MUOX*OX(I)

               CTE = epgaj*CTE
               CTW = epgaj*CTW
               CPE = epgaj*CPE
               CPW = epgaj*CPW
               VXZA = epgaj*VXZA
            ENDIF

! if jackson, implement jackson form of governing equations (ep_g dot
! del tau_g): multiply by void fraction otherwise by 1
            ltau_w_g = epgaj*tau_w_g(ijk)

! Collect the terms
            A_M(IJK,east,M) = A_M(IJK,east,M) + CPE
            A_M(IJK,west,M) = A_M(IJK,west,M) - CPW

            A_M(IJK,0,M) = -(A_M(IJK,east,M)+A_M(IJK,west,M)+&
               A_M(IJK,north,M)+A_M(IJK,south,M)+A_M(IJK,top,M)+A_M(IJK,bottom,M)+&
               (V0+VPM+ZMAX(VMT)+VCOA + VXZA)*VOL_W(IJK) + CTE - CTW)

            A_M(IJK,east,M) = A_M(IJK,east,M) - CTE
            A_M(IJK,west,M) = A_M(IJK,west,M) + CTW

            B_M(IJK,M) = B_M(IJK,M) - ( SDP + lTAU_W_G + F_VIR + &
               ( (V0+ZMAX((-VMT)))*W_GO(IJK) + VBF + VCOB + &
               Ghd_drag+HYS_drag)*VOL_W(IJK) )

! MMS Source term.
            IF(USE_MMS) B_M(IJK,M) = &
               B_M(IJK,M) - MMS_W_G_SRC(IJK)*VOL_W(IJK)
         ENDIF   ! end branching on cell type (ip/dilute/block/else branches)
      ENDDO   ! end do loop over ijk
!$omp end parallel do

! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_W_G(A_M, B_M)
! modifications for bc
      CALL SOURCE_W_G_BC (A_M, B_M)
! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_W_G_BC(A_M, B_M)

      RETURN
      END SUBROUTINE SOURCE_W_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_W_g_BC                                           C
!  Purpose: Determine source terms for W_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_W_G_BC(A_M, B_M)

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
      USE visc_g
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE output
      USE compar
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! C_mu is constant in turbulent viscosity
      DOUBLE PRECISION, PARAMETER :: C_mu = 0.09D0
! Kappa is Von Karmen constant
      DOUBLE PRECISION, PARAMETER :: Kappa = 0.42D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK,&
                 IM, JM, IJKB, IJKM, IJKP
! Phase index
      INTEGER :: M
! Turbulent shear at walls
      DOUBLE PRECISION W_F_Slip

!-----------------------------------------------

! Set reference phase to gas
      M = 0


! Set the default boundary conditions
! The NS default setting is the where bc_type='dummy' or any default
! (i.e., bc_type=undefined) wall boundary regions are handled. Note that
! the top and bottom xy planes do not have to be explicitly addressed for
! the w-momentum equation. In this direction the velocities are defined
! at the wall (due staggered grid). They are defined as zero for a
! no penetration condition (see zero_norm_vel subroutine and code under
! the ip_at_t branch in the above source routine).
! ---------------------------------------------------------------->>>

! south xz plane
      J1 = 1
      DO K1 = kmin3,kmax3
         DO I1 = imin3,imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = -ONE
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
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

! north xz plane
      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
! Setting the wall velocity to zero (set the boundary cell value equal
! and oppostive to the adjacent fluid cell value)
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = -ONE
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
! Setting the wall velocity equal to the adjacent fluid velocity (set
! the boundary cell value equal to adjacent fluid cell value)
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

! west yz plane
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

! east yz plane
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
! ----------------------------------------------------------------<<<

! Setting user specified boundary conditions

      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

! Setting wall boundary conditions
! ---------------------------------------------------------------->>>
            IF (BC_TYPE_ENUM(L) == NO_SLIP_WALL .AND. .NOT. K_Epsilon) THEN
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

            ELSEIF (BC_TYPE_ENUM(L) == FREE_SLIP_WALL .AND. .NOT. K_Epsilon) THEN
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

            ELSEIF (BC_TYPE_ENUM(L) == PAR_SLIP_WALL .AND. .NOT. K_Epsilon) THEN
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
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,east,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_WW_G(L)
                           ELSE
                              IF (CYLINDRICAL) THEN
                                 A_M(IJK,0,M) = -( HALF*(BC_HW_G(L)-&
                                    OX_E(I)) + ODX_E(I) )
                                 A_M(IJK,east,M) = -(HALF*(BC_HW_G(L)-&
                                    OX_E(I)) - ODX_E(I))
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L)
                              ELSE
                                 A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(I))
                                 A_M(IJK,east,M) = -(HALF*BC_HW_G(L)-ODX_E(I))
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L)
                              ENDIF
                           ENDIF
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,west,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_WW_G(L)
                           ELSE
                              IF (CYLINDRICAL) THEN
                                 A_M(IJK,west,M) = -(HALF*(BC_HW_G(L)-&
                                    OX_E(IM)) - ODX_E(IM))
                                 A_M(IJK,0,M) = -(HALF*(BC_HW_G(L)-&
                                    OX_E(IM)) + ODX_E(IM))
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L)
                              ELSE
                                 A_M(IJK,west,M) = -(HALF*BC_HW_G(L)-ODX_E(IM))
                                 A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODX_E(IM))
                                 B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L)
                              ENDIF
                           ENDIF
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,north,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_WW_G(L)
                           ELSE
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(J))
                              A_M(IJK,north,M) = -(HALF*BC_HW_G(L)-ODY_N(J))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,south,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_WW_G(L)
                           ELSE
                              A_M(IJK,south,M) = -(HALF*BC_HW_G(L)-ODY_N(JM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(JM))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_WW_G(L)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

! Setting wall boundary conditions when K_EPSILON
! wall functions for V-momentum are specify in this section of the code
            ELSEIF (BC_TYPE_ENUM(L) == PAR_SLIP_WALL   .OR.  &
                    BC_TYPE_ENUM(L) == NO_SLIP_WALL    .OR.  &
                    BC_TYPE_ENUM(L) == FREE_SLIP_WALL  .AND. &
                    K_Epsilon )THEN
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
                           IF (CYLINDRICAL) THEN
                              W_F_Slip = ( ONE/&
                                 (ODX_E(I)+HALF*OX_E(I)) )*          &
                                 ( ODX_E(I) - OX_E(I) -              &
                                   RO_g(EAST_OF(IJK))*C_mu**0.25*    &
                                   SQRT(K_Turb_G((EAST_OF(IJK))))/   &
                                   MU_gT(EAST_OF(IJK))*Kappa/        &
                                   LOG(9.81D0/ODX_E(I)/(2.D0)*       &
                                   RO_g(EAST_OF(IJK))*C_mu**0.25*    &
                                   SQRT(K_Turb_G((EAST_OF(IJK))))/   &
                                   MU_g(EAST_OF(IJK))) )
                           ELSE
                              CALL Wall_Function(IJK,EAST_OF(IJK),&
                                 ODX_E(I),W_F_Slip)
                           ENDIF
                           A_M(IJK,east,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE
                           IF (BC_TYPE_ENUM(L) == PAR_SLIP_WALL) B_M(IJK,M) = -BC_WW_G(L)
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           IF (CYLINDRICAL) THEN
                              W_F_Slip =  (ONE/&
                                 (ONE*ODX_E(IM) + HALF*OX_E(IM)))*    &
                                 ( ONE*ODX_E(IM) - OX_E(IM) -         &
                                   RO_g(WEST_OF(IJK))*C_mu**0.25*     &
                                   SQRT(K_Turb_G((WEST_OF(IJK))))/    &
                                   MU_gT(WEST_OF(IJK))*Kappa/         &
                                   LOG(9.81D0/ODX_E(IM)/(2.D0)*       &
                                   RO_g(WEST_OF(IJK))*C_mu**0.25*     &
                                   SQRT(K_Turb_G((WEST_OF(IJK))))/    &
                                   MU_g(WEST_OF(IJK))) )
                           ELSE
                              CALL Wall_Function(IJK,WEST_OF(IJK),&
                                 ODX_E(IM),W_F_Slip)
                           ENDIF
                           A_M(IJK,west,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE
                           IF (BC_TYPE_ENUM(L) == PAR_SLIP_WALL) B_M(IJK,M) = -BC_WW_G(L)
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                           CALL Wall_Function(IJK,NORTH_OF(IJK),ODY_N(J),W_F_Slip)
                           A_M(IJK,north,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE
                           IF (BC_TYPE_ENUM(L) == PAR_SLIP_WALL) B_M(IJK,M) = -BC_WW_G(L)
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           CALL Wall_Function(IJK,SOUTH_OF(IJK),ODY_N(JM),W_F_Slip)
                           A_M(IJK,south,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE
                           IF (BC_TYPE_ENUM(L) == PAR_SLIP_WALL) B_M(IJK,M) = -BC_WW_G(L)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

! end setting of wall boundary conditions
! ----------------------------------------------------------------<<<

! Setting p_inflow or p_outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE_ENUM(L)==P_INFLOW .OR. BC_TYPE_ENUM(L)==P_OUTFLOW) THEN
               IF (BC_PLANE(L) == 'B') THEN
! if the fluid cell is on the bottom side of the outflow/inflow boundary
! then set the velocity in the boundary cell equal to the velocity of
! the adjacent fluid cell
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
! end setting of p_inflow or p_otuflow flow boundary conditions
! ----------------------------------------------------------------<<<

! Setting outflow flow boundary conditions
! ---------------------------------------------------------------->>>
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
! end setting of outflow flow boundary conditions
! ----------------------------------------------------------------<<<

! Setting bc that are defined but not nsw, fsw, psw, p_inflow,
! p_outflow, or outflow (at this time, this section addresses
! mass_inflow and mass_outflow type boundaries)
! ---------------------------------------------------------------->>>
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
! setting the velocity in the boundary cell equal to what is known
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = -W_G(IJK)
                        IF (BC_PLANE(L) == 'B') THEN
! if the fluid cell is on the bottom side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           IJKB = BOTTOM_OF(IJK)
                           A_M(IJKB,east,M) = ZERO
                           A_M(IJKB,west,M) = ZERO
                           A_M(IJKB,north,M) = ZERO
                           A_M(IJKB,south,M) = ZERO
                           A_M(IJKB,top,M) = ZERO
                           A_M(IJKB,bottom,M) = ZERO
                           A_M(IJKB,0,M) = -ONE
                           B_M(IJKB,M) = -W_G(IJKB)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF   ! end if/else (bc_type)
                    ! ns, fs, psw; else
                    ! p_inflow, p_outflow, or outflow; else
! end setting of 'else' flow boundary conditions
! (mass_inflow/mass_outflow)
! ----------------------------------------------------------------<<<

         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc

      RETURN
      END SUBROUTINE SOURCE_W_G_BC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_W_G                                        C
!  Purpose: Adds point sources to the gas phase W-Momentum equation.   C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_W_G(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar
      use constant
      use geometry
      use indices
      use param
      use param1, only: small_number, zero
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

! Set reference phase to gas
      M = 0

! Calculate the mass going into each IJK cell. This is done for each
! call in case the point source is time dependent.
      PS_LP: do PSV = 1, DIMENSION_PS
         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(abs(PS_W_g(PSV)) < small_number) cycle PS_LP

         if(PS_W_g(PSV) < ZERO) then
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

            pSource =  PS_MASSFLOW_G(PSV) * (VOL(IJK)/PS_VOLUME(PSV))

            B_M(IJK,M) = B_M(IJK,M) - pSource * &
               PS_W_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_W_G
