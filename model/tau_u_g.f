!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_U_g                                            C
!  Purpose: Cross terms in the gradient of stress in U_g momentum      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-DEC-96  C
!                                                                      C
!                                                                      C
!  Comments: This routine calculates the components of the gas phase   C
!  viscous stress tensor of the u-momentum equation that cannot be     C
!  cast in the form: mu.grad(u). These components are stored in the    C
!  passed variable, which is then used as a source of the u-momentum   C
!  equation.                                                           C
!                                                                      C
!  The following u component viscous stress tensor terms are           C
!  calculated here:                                                    C
!  > Combined part of 1/x d/dx (x.tau_xx) xdxdydz and                  C
!                              -tau_zz/x xdxdydz =>                    C
!    1/x d/dx (x.lambda.trcD) xdxdydz - (lambda/x.trcD) xdxdydz =>     C
!        d/dx (lambda.trcD) xdxdydz =>                                 C
!    delta (lambda.trcD)Ap|E-W                                         C
!  > part of 1/x d/dx(x.tau_xx) xdxdydz =>                             C
!            1/x d/dx (x.mu.du/dx) xdxdydz =>                          C
!    delta (mu du/dx)Ayz |E-W                                          C
!  > part of d/dy (tau_xy) xdxdydz =>                                  C
!           d/dy (mu.dv/dx) xdxdydz =>                                 C
!    delta (mu.dv/dx)Axz |N-S                                          C
!  > part of 1/x d/dz (tau_xz) xdxdydz =>                              C
!            1/x d/dz (mu.dw/dx) xdxdydz =>                            C
!    delta (mu.dw/dx)Axy |T-B                                          C
!  CYLINDRICAL TERMS:                                                  C
!  > part of 1/x d/dz (tau_xz) xdxdydz =>                              C
!            1/x d/dz (mu.(-w/x)) xdxdydz =>                           C
!    delta (mu(-w/x))Axy |T-B                                          C
!  > part of -tau_zz/x xdxdydz =>                                      C
!            -(2.mu/x).(1/x).dw/dz xdxdydz                             C
!    delta (2.mu/x.1/x.dw/dz)Vp |E-W                                   C
!                                                                      C
!  The following u-component CYLINDRICAL viscous stress tensor terms   C
!  are not calculated here. They are handled in source_u_g:            C
!  > part of the term -tau_zz/x xdxdydz =>                             C
!                     -2.mu/x. u/x xdxdydz                             C
!    delta(-2.mu.u/xp^2)Vp |p                                          C
!                                                                      C
!  To reconstitute the complete u-momentum gas phase viscous stress    C
!  tensor would require including the term calculated in source_u_g    C
!  and the 'diffusional' components (i.e., those of the form           C
!  mu.grad(u)                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TAU_U_G(lTAU_U_G, lctau_u_g)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3
      USE param1, only: zero, half
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE run
      USE toleranc, only: dil_ep_s
      USE geometry
      USE indices
      USE is
      USE compar
      USE sendrecv
      USE functions
      USE cutcell
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_U_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_U_g(DIMENSION_3)
! cTAU_U_g
      DOUBLE PRECISION, INTENT(OUT) :: lcTAU_U_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IP, JM, KM
      INTEGER :: IJK, IJKE, IJKN, IJKS, IJKT, IJKB
      INTEGER :: IJKNE, IJKSE, IJKTE, IJKBE
      INTEGER :: IPJK, IMJK, IJMK, IJKM
      INTEGER :: IPJMK, IPJKM
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average viscosity
      DOUBLE PRECISION :: MU_gte, MU_gbe, MUGA
! Average dW/Xdz
      DOUBLE PRECISION :: dWoXdz
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
! Source terms (Volumetric)
      DOUBLE PRECISION :: Vtzb
!---------------------------------------------------------------------//

      IF((.NOT.CARTESIAN_GRID).OR.(CG_SAFE_MODE(3)==1)) THEN

!$omp  parallel do default(none) &
!$omp  private(I, J, K, IP, JM, KM, IJK,                               &
!$omp          IJKE, IJKN, IJKS, IJKT, IJKB,                           &
!$omp          IJKNE, IJKSE, IJKTE, IJKBE,                             &
!$omp          IPJK, IMJK, IJMK, IJKM, IPJMK, IPJKM,                   &
!$omp          EPGA, MUGA, MU_GTE, MU_GBE, SBV, SSX, SSY, SSZ,         &
!$omp          DWOXDZ, VTZB)                                           &
!$omp  shared(ijkstart3, ijkend3, i_of, j_of, k_of, ip1, jm1, km1,     &
!$omp         do_k, cylindrical, ltau_u_g, lctau_u_g,                  &
!$omp         ep_g, epmu_gt, lambda_gt, trd_g, u_g, v_g, w_g,          &
!$omp         axy_u, axz_u, ayz, ayz_u, vol_u,                         &
!$omp         ox_e, odx_e, ox, odx, odz, eplambda_gt)
         DO IJK = IJKSTART3, IJKEND3
            I = I_OF(IJK)
            IJKE = EAST_OF(IJK)
            EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)
            IF (.NOT.IP_AT_E(IJK) .AND. EPGA>DIL_EP_S) THEN
               J = J_OF(IJK)
               K = K_OF(IJK)
               IP = IP1(I)
               JM = JM1(J)
               KM = KM1(K)

               IPJK = IP_OF(IJK)
               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IPJMK = JM_OF(IPJK)
               IPJKM = IP_OF(IJKM)

               IJKN = NORTH_OF(IJK)
               IJKNE = EAST_OF(IJKN)
               IJKS = SOUTH_OF(IJK)
               IJKSE = EAST_OF(IJKS)
               IJKT = TOP_OF(IJK)
               IJKTE = EAST_OF(IJKT)
               IJKB = BOTTOM_OF(IJK)
               IJKBE = EAST_OF(IJKB)


! Surface forces at i+1/2, j, k
! bulk viscosity term
! combines part of 1/x d/dx (x.tau_xx) xdxdydz and -tau_zz/x xdxdydz =>
! combines 1/x d/dx (x.lambda.trcD) xdxdydz - (lambda/x.trcD) xdxdydz =>
!              d/dx (lambda.trcD) xdxdydz
! delta (lambda.trcD)Ap |E-W : at (i+1 - i-1), j, k
               SBV = (EPLAMBDA_GT(IJKE)*TRD_G(IJKE)-&
                      EPLAMBDA_GT(IJK)*TRD_G(IJK))*AYZ(IJK)

! shear stress terms at i+1/2, j, k
! part of 1/x d/dx(x.tau_xx) xdxdydz =>
!         1/x d/dx (x.mu.du/dx) xdxdydz =>
! delta (mu du/dx)Ayz |E-W : at (i+1 - i-1), j, k
               SSX = EPMU_GT(IJKE)*(U_G(IPJK)-U_G(IJK))*ODX(IP)*&
                        AYZ_U(IJK) - &
                     EPMU_GT(IJK)*(U_G(IJK)-U_G(IMJK))*ODX(I)*&
                        AYZ_U(IMJK)

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dx) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at i+1/2, (j+1/2 - j-1/2), k
               SSY = AVG_X_H(AVG_Y_H(EPMU_GT(IJK),EPMU_GT(IJKN),J),&
                             AVG_Y_H(EPMU_GT(IJKE),EPMU_GT(IJKNE),J),I)*&
                        (V_G(IPJK)-V_G(IJK))*ODX_E(I)*AXZ_U(IJK) - &
                     AVG_X_H(AVG_Y_H(EPMU_GT(IJKS),EPMU_GT(IJK),JM),&
                             AVG_Y_H(EPMU_GT(IJKSE),EPMU_GT(IJKE),JM),I)*&
                        (V_G(IPJMK)-V_G(IJMK))*ODX_E(I)*AXZ_U(IJMK)

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dx) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at i+1/2, j, (k+1/2 - k-1/2)
               MU_GTE = AVG_X_H(AVG_Z_H(EPMU_GT(IJK),EPMU_GT(IJKT),K),&
                                AVG_Z_H(EPMU_GT(IJKE),EPMU_GT(IJKTE),K),I)
               MU_GBE = AVG_X_H(AVG_Z_H(EPMU_GT(IJKB),EPMU_GT(IJK),KM),&
                                AVG_Z_H(EPMU_GT(IJKBE),EPMU_GT(IJKE),KM),I)
               SSZ = MU_GTE*(W_G(IPJK)-W_G(IJK))*ODX_E(I)*&
                        AXY_U(IJK) - &
                     MU_GBE*(W_G(IPJKM)-W_G(IJKM))*ODX_E(I)*&
                        AXY_U(IJKM)


! Special terms for cylindrical coordinates
               IF (CYLINDRICAL) THEN
! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.(-w/x)) xdxdydz =>
! delta (mu(-w/x))Axy |T-B : at i+1/2, j, (k+1/2- k-1/2)
                  SSZ = SSZ - ( &
                        MU_GTE*(HALF*(W_G(IPJK)+W_G(IJK))*OX_E(I))*&
                           AXY_U(IJK)-&
                        MU_GBE*(HALF*(W_G(IPJKM)+W_G(IJKM))*OX_E(I))*&
                           AXY_U(IJKM) )

! part of -tau_zz/x xdxdydz =>
!         -(2.mu/x).(1/x).dw/dz xdxdydz
! delta (2.mu/x.1/x.dw/dz)Vp |p : at i+1/2, j, k
                  MUGA = AVG_X(EPMU_GT(IJK),EPMU_GT(IJKE),I)
                  DWOXDZ = HALF*((W_G(IJK)-W_G(IJKM))*OX(I)*ODZ(K)+&
                                 (W_G(IPJK)-W_G(IPJKM))*OX(IP)*ODZ(K))
                  VTZB = -2.d0*MUGA*OX_E(I)*DWOXDZ
               ELSE
                  VTZB = ZERO
               ENDIF

! Add the terms
               lTAU_U_G(IJK) = SBV + SSX + SSY + SSZ + VTZB*VOL_U(IJK)

! Also calculate and store the full gas phase viscous stress tensor
               CALL GET_FULL_TAU_U_G(IJK, ltau_u_g, lctau_u_g)

            ELSE
               lTAU_U_G(IJK) = ZERO
               lctau_u_g(IJK) = zero
            ENDIF   ! end if (.NOT.IP_AT_E(IJK) .AND. EPGA>DIL_EP_S)
         ENDDO   ! end do ijk
!$omp end parallel do

      ELSE
! if cartesian grid
         CALL CALC_CG_TAU_U_G(lTAU_U_G, lctau_u_g)
      ENDIF

      call send_recv(ltau_u_g,2)
      call send_recv(lctau_u_g,2)

      RETURN

    CONTAINS

      INCLUDE 'fun_avg.inc'

    END SUBROUTINE CALC_TAU_U_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_CG_Tau_U_g                                         C
!  Purpose: Cross terms in the gradient of stress in U_g momentum      C
!  based on cartesian grid cut cell.                                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_TAU_U_G(lTAu_U_G, lctau_u_g)

! Modules
!---------------------------------------------------------------------//

      USE bc
      USE compar
      USE constant
      USE cutcell
      USE fldvar
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE is
      USE param, only: dimension_3
      USE param1, only: zero, half, undefined
      USE physprop
      USE quadric
      USE run
      USE rxns
      USE toleranc, only: dil_ep_s
      USE visc_g

      IMPLICIT NONE
! Dummy arguments
!---------------------------------------------------------------------//
! TAU_U_g
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_U_g(DIMENSION_3)
! cTAU_U_g
      DOUBLE PRECISION, INTENT(OUT) :: lcTAU_U_g(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IP, JM, KM
      INTEGER :: IJK, IJKE, IJKN, IJKS, IJKT, IJKB
      INTEGER :: IJKNE, IJKSE, IJKTE, IJKBE
      INTEGER :: IPJK, IMJK, IJMK, IJKM
      INTEGER :: IPJMK, IPJKM
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average viscosity
      DOUBLE PRECISION :: MU_gte, MU_gbe
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
! Cartesian grid variables
      DOUBLE PRECISION :: DEL_H, Nx, Ny, Nz
      LOGICAL :: V_NODE_AT_NE, V_NODE_AT_NW, V_NODE_AT_SE, V_NODE_AT_SW
      LOGICAL :: W_NODE_AT_TE, W_NODE_AT_TW, W_NODE_AT_BE, W_NODE_AT_BW
      DOUBLE PRECISION :: dvdx_at_N, dvdx_at_S
      DOUBLE PRECISION :: dwdx_at_T, dwdx_at_B
      DOUBLE PRECISION :: Xi, Yi, Zi, Vi, Wi, Sx, Sy, Sz
      DOUBLE PRECISION :: MU_GT_CUT, SSY_CUT, SSZ_CUT
      DOUBLE PRECISION :: UW_g, VW_g, WW_g
      INTEGER :: BCV
      INTEGER :: BCT
!---------------------------------------------------------------------//

!$omp  parallel do default(none) &
!$omp  private(I, J, K, JM, KM, IP, IJK,                               &
!$omp          IJKE, IJKN, IJKS, IJKT, IJKB,                           &
!$omp          IJKNE, IJKSE, IJKTE, IJKBE,                             &
!$omp          IPJK, IMJK, IJMK, IJKM, IPJMK, IPJKM,                   &
!$omp          EPGA, MU_GTE, MU_GBE, SBV, SSX, SSY, SSZ,               &
!$omp          BCV, BCT, BC_TYPE_ENUM, NOC_UG, uw_g, vw_g, ww_g,       &
!$omp          del_h, nx, ny, nz, xi, yi, zi, sx, sy, sz, vi, wi,      &
!$omp          v_node_at_ne, v_node_at_sw, v_node_at_se, v_node_at_nw, &
!$omp          w_node_at_bw, w_node_at_be, w_node_at_tw, w_node_at_te, &
!$omp          dvdx_at_S, dvdx_at_N, dwdx_at_B, dwdx_at_T,             &
!$omp          mu_gt_cut, cut_tau_ug, ssy_cut, ssz_cut)                &
!$omp  shared(ijkstart3, ijkend3, i_of, j_of, k_of,  ip1, jm1, km1,    &
!$omp         do_k, ltau_u_g, lctau_u_g,                               &
!$omp         ep_g, epmu_gt, lambda_gt, trd_g, u_g, v_g, w_g,          &
!$omp         axy_u, axz_u, ayz, ayz_u, vol,                           &
!$omp         bc_type, bc_u_id, bc_uw_g, bc_vw_g, bc_ww_g, bc_hw_g,    &
!$omp         oneodx_e_u, oneodx_e_v, oneodx_e_w,                      &
!$omp         x_u, y_u, z_u, x_v, y_v, z_v, x_w, y_w, z_w,             &
!$omp         wall_v_at, wall_w_at, area_u_cut, cut_u_cell_at,         &
!$omp         blocked_v_cell_at, blocked_w_cell_at, eplambda_gt)

      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         IJKE = EAST_OF(IJK)
         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)
         IF ( .NOT.IP_AT_E(IJK) .AND. EPGA>DIL_EP_S) THEN
            J = J_OF(IJK)
            K = K_OF(IJK)
            IP = IP1(I)
            JM = JM1(J)
            KM = KM1(K)

            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            IPJMK = JM_OF(IPJK)
            IPJKM = IP_OF(IJKM)

            IJKN = NORTH_OF(IJK)
            IJKNE = EAST_OF(IJKN)
            IJKS = SOUTH_OF(IJK)
            IJKSE = EAST_OF(IJKS)
            IJKT = TOP_OF(IJK)
            IJKTE = EAST_OF(IJKT)
            IJKB = BOTTOM_OF(IJK)
            IJKBE = EAST_OF(IJKB)


! bulk viscosity term
            SBV = (EPLAMBDA_GT(IJKE)*TRD_G(IJKE)) * AYZ_U(IJK) - &
                  (EPLAMBDA_GT(IJK) *TRD_G(IJK) ) * AYZ_U(IMJK)

! shear stress terms
            IF(.NOT.CUT_U_CELL_AT(IJK)) THEN
               SSX = EPMU_GT(IJKE)*(U_G(IPJK)-U_G(IJK))*&
                        ONEoDX_E_U(IJK)*AYZ_U(IJK) - &
                     EPMU_GT(IJK)*(U_G(IJK)-U_G(IMJK))*&
                        ONEoDX_E_U(IMJK)*AYZ_U(IMJK)

               SSY = AVG_X_H(AVG_Y_H(EPMU_GT(IJK),EPMU_GT(IJKN),J),&
                             AVG_Y_H(EPMU_GT(IJKE),EPMU_GT(IJKNE),J),I)*&
                     (V_G(IPJK)-V_G(IJK))*ONEoDX_E_V(IJK)*AXZ_U(IJK) - &
                     AVG_X_H(AVG_Y_H(EPMU_GT(IJKS),EPMU_GT(IJK),JM),&
                             AVG_Y_H(EPMU_GT(IJKSE),EPMU_GT(IJKE),JM),I)*&
                     (V_G(IPJMK)-V_G(IJMK))*ONEoDX_E_V(IJMK)*AXZ_U(IJMK)

               IF(DO_K) THEN
                  MU_GTE = AVG_X_H(AVG_Z_H(EPMU_GT(IJK),EPMU_GT(IJKT),K),&
                                     AVG_Z_H(EPMU_GT(IJKE),EPMU_GT(IJKTE),K),I)
                  MU_GBE = AVG_X_H(AVG_Z_H(EPMU_GT(IJKB),EPMU_GT(IJK),KM),&
                                   AVG_Z_H(EPMU_GT(IJKBE),EPMU_GT(IJKE),KM),I)
                  SSZ = MU_GTE*(W_G(IPJK)-W_G(IJK))*&
                           ONEoDX_E_W(IJK)*AXY_U(IJK) - &
                        MU_GBE*(W_G(IPJKM)-W_G(IJKM))*&
                           ONEoDX_E_W(IJKM)*AXY_U(IJKM)
               ELSE
                  SSZ = ZERO
               ENDIF

! cut cell modifications
!---------------------------------------------------------------------//
            ELSE

               BCV = BC_U_ID(IJK)
               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               SELECT CASE (BCT)
                  CASE (CG_NSW,CG_MI)
                     CUT_TAU_UG = .TRUE.
                     NOC_UG     = .TRUE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE (CG_FSW)
                     CUT_TAU_UG = .FALSE.
                     NOC_UG     = .FALSE.
                     UW_g = ZERO
                     VW_g = ZERO
                     WW_g = ZERO
                  CASE(CG_PSW)
                     IF(BC_HW_G(BC_U_ID(IJK))==UNDEFINED) THEN   ! same as NSW
                        CUT_TAU_UG = .TRUE.
                        NOC_UG     = .TRUE.
                        UW_g = BC_UW_G(BCV)
                        VW_g = BC_VW_G(BCV)
                        WW_g = BC_WW_G(BCV)
                     ELSEIF(BC_HW_G(BC_U_ID(IJK))==ZERO) THEN   ! same as FSW
                        CUT_TAU_UG = .FALSE.
                        NOC_UG     = .FALSE.
                     ELSE                              ! partial slip
                        CUT_TAU_UG = .FALSE.
                        NOC_UG     = .FALSE.
                        UW_g = ZERO
                        VW_g = ZERO
                        WW_g = ZERO
                     ENDIF
                  CASE (NONE)
                     lTAU_U_G(IJK) = ZERO
                     CYCLE
               END SELECT

               IF(CUT_TAU_UG) THEN
                  MU_GT_CUT = (VOL(IJK)*EPMU_GT(IJK) + &
                               VOL(IPJK)*EPMU_GT(IJKE))/&
                              (VOL(IJK) + VOL(IPJK))
               ELSE
                  MU_GT_CUT = ZERO
               ENDIF

! SSX:
               CALL GET_DEL_H(IJK,'U_MOMENTUM', X_U(IJK), &
                  Y_U(IJK), Z_U(IJK), Del_H, Nx, Ny, Nz)

               SSX = EPMU_GT(IJKE)*(U_G(IPJK)-U_G(IJK))*&
                        ONEoDX_E_U(IJK)*AYZ_U(IJK) - &
                     EPMU_GT(IJK)*(U_G(IJK)-U_G(IMJK))*&
                        ONEoDX_E_U(IMJK)*AYZ_U(IMJK) - &
                     MU_GT_CUT * (U_g(IJK) - UW_g) / DEL_H * &
                        (Nx**2) * Area_U_CUT(IJK)

! SSY:
               V_NODE_AT_NE = ((.NOT.BLOCKED_V_CELL_AT(IPJK)).AND.&
                               (.NOT.WALL_V_AT(IPJK)))
               V_NODE_AT_NW = ((.NOT.BLOCKED_V_CELL_AT(IJK)).AND.&
                               (.NOT.WALL_V_AT(IJK)))
               V_NODE_AT_SE = ((.NOT.BLOCKED_V_CELL_AT(IPJMK)).AND.&
                               (.NOT.WALL_V_AT(IPJMK)))
               V_NODE_AT_SW = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.&
                               (.NOT.WALL_V_AT(IJMK)))

               IF(V_NODE_AT_NE.AND.V_NODE_AT_NW) THEN
                  Vi = HALF * (V_G(IPJK) + V_G(IJK))
                  Xi = HALF * (X_V(IPJK) + X_V(IJK))
                  Yi = HALF * (Y_V(IPJK) + Y_V(IJK))
                  Zi = HALF * (Z_V(IPJK) + Z_V(IJK))
                  Sx = X_V(IPJK) - X_V(IJK)
                  Sy = Y_V(IPJK) - Y_V(IJK)
                  Sz = Z_V(IPJK) - Z_V(IJK)
                  CALL GET_DEL_H(IJK,'U_MOMENTUM',Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dvdx_at_N = (V_G(IPJK)-V_G(IJK)) * ONEoDX_E_V(IJK)
                  IF(NOC_UG) dvdx_at_N = dvdx_at_N - ( (Vi-VW_g)*&
                     ONEoDX_E_V(IJK)/DEL_H * (Sy*Ny+Sz*Nz) )
               ELSE
                  dvdx_at_N =  ZERO
               ENDIF

               IF(V_NODE_AT_SE.AND.V_NODE_AT_SW) THEN
                  Vi = HALF * (V_G(IPJMK) + V_G(IJMK))
                  Xi = HALF * (X_V(IPJMK) + X_V(IJMK))
                  Yi = HALF * (Y_V(IPJMK) + Y_V(IJMK))
                  Zi = HALF * (Z_V(IPJMK) + Z_V(IJMK))
                  Sx = X_V(IPJMK) - X_V(IJMK)
                  Sy = Y_V(IPJMK) - Y_V(IJMK)
                  Sz = Z_V(IPJMK) - Z_V(IJMK)
                  CALL GET_DEL_H(IJK,'U_MOMENTUM', Xi, Yi, Zi, &
                     Del_H, Nx, Ny, Nz)

                  dvdx_at_S = (V_G(IPJMK)-V_G(IJMK)) * ONEoDX_E_V(IJMK)
                  IF(NOC_UG) dvdx_at_S = dvdx_at_S - ( (Vi-VW_g)*&
                     ONEoDX_E_V(IJMK)/DEL_H * (Sy*Ny+Sz*Nz) )
               ELSE
                  dvdx_at_S =  ZERO
               ENDIF

               IF(V_NODE_AT_NW) THEN
                  CALL GET_DEL_H(IJK,'U_MOMENTUM', X_V(IJK), &
                     Y_V(IJK), Z_V(IJK), Del_H, Nx, Ny, Nz)
                  SSY_CUT = -MU_GT_CUT * (V_G(IJK)-VW_g)/DEL_H * &
                     (Nx*Ny) * Area_U_CUT(IJK)
               ELSE
                  SSY_CUT =  ZERO
               ENDIF

               SSY = AVG_X_H(AVG_Y_H(EPMU_GT(IJK),EPMU_GT(IJKN),J),&
                             AVG_Y_H(EPMU_GT(IJKE),EPMU_GT(IJKNE),J),I)*&
                        dvdx_at_N*AXZ_U(IJK) - &
                     AVG_X_H(AVG_Y_H(EPMU_GT(IJKS),EPMU_GT(IJK),JM),&
                             AVG_Y_H(EPMU_GT(IJKSE),EPMU_GT(IJKE),JM),I)*&
                        dvdx_at_S*AXZ_U(IJMK) + SSY_CUT

! SSZ:
               IF(DO_K) THEN
                  W_NODE_AT_TE = ((.NOT.BLOCKED_W_CELL_AT(IPJK)).AND.&
                                  (.NOT.WALL_W_AT(IPJK)))
                  W_NODE_AT_TW = ((.NOT.BLOCKED_W_CELL_AT(IJK)).AND.&
                                  (.NOT.WALL_W_AT(IJK)))
                  W_NODE_AT_BE = ((.NOT.BLOCKED_W_CELL_AT(IPJKM)).AND.&
                                  (.NOT.WALL_W_AT(IPJKM)))
                  W_NODE_AT_BW = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.&
                                  (.NOT.WALL_W_AT(IJKM)))

                  IF(W_NODE_AT_TE.AND.W_NODE_AT_TW) THEN
                     Wi = HALF * (W_G(IPJK) + W_G(IJK))
                     Xi = HALF * (X_W(IPJK) + X_W(IJK))
                     Yi = HALF * (Y_W(IPJK) + Y_W(IJK))
                     Zi = HALF * (Z_W(IPJK) + Z_W(IJK))
                     Sx = X_W(IPJK) - X_W(IJK)
                     Sy = Y_W(IPJK) - Y_W(IJK)
                     Sz = Z_W(IPJK) - Z_W(IJK)
                     CALL GET_DEL_H(IJK, 'U_MOMENTUM', Xi, Yi, Zi, &
                        Del_H, Nx, Ny, Nz)

                     dwdx_at_T = (W_G(IPJK)-W_G(IJK)) * ONEoDX_E_W(IJK)
                     IF(NOC_UG) dwdx_at_T = dwdx_at_T - ( (Wi-WW_g) * &
                        ONEoDX_E_W(IJK)/DEL_H * (Sy*Ny+Sz*Nz) )
                  ELSE
                     dwdx_at_T =  ZERO
                  ENDIF

                  IF(W_NODE_AT_BE.AND.W_NODE_AT_BW) THEN
                     Wi = HALF * (W_G(IPJKM) + W_G(IJKM))
                     Xi = HALF * (X_W(IPJKM) + X_W(IJKM))
                     Yi = HALF * (Y_W(IPJKM) + Y_W(IJKM))
                     Zi = HALF * (Z_W(IPJKM) + Z_W(IJKM))
                     Sx = X_W(IPJKM) - X_W(IJKM)
                     Sy = Y_W(IPJKM) - Y_W(IJKM)
                     Sz = Z_W(IPJKM) - Z_W(IJKM)

                     CALL GET_DEL_H(IJK,'U_MOMENTUM', Xi, Yi, Zi,&
                        Del_H, Nx, Ny, Nz)

                     dwdx_at_B = (W_G(IPJKM)-W_G(IJKM)) * ONEoDX_E_W(IJKM)
                     IF(NOC_UG) dwdx_at_B = dwdx_at_B - ( (Wi-WW_g) *&
                        ONEoDX_E_W(IJKM)/DEL_H * (Sy*Ny+Sz*Nz))
                  ELSE
                     dwdx_at_B =  ZERO
                  ENDIF

                  IF(W_NODE_AT_TW) THEN
                     CALL GET_DEL_H(IJK,'U_MOMENTUM', X_W(IJK), &
                        Y_W(IJK), Z_W(IJK), Del_H, Nx, Ny, Nz)
                     SSZ_CUT = -MU_GT_CUT * (W_G(IJK)-WW_g)/DEL_H * &
                        (Nx*Nz) * Area_U_CUT(IJK)
                  ELSE
                     SSZ_CUT =  ZERO
                  ENDIF

                  MU_GTE = AVG_X_H(AVG_Z_H(EPMU_GT(IJK),EPMU_GT(IJKT),K),&
                                   AVG_Z_H(EPMU_GT(IJKE),EPMU_GT(IJKTE),K),I)
                  MU_GBE = AVG_X_H(AVG_Z_H(EPMU_GT(IJKB),EPMU_GT(IJK),KM),&
                                   AVG_Z_H(EPMU_GT(IJKBE),EPMU_GT(IJKE),KM),I)
                  SSZ = MU_GTE*dwdx_at_T*AXY_U(IJK) - &
                        MU_GBE*dwdx_at_B*AXY_U(IJKM)  + SSZ_CUT
               ELSE
                 SSZ = ZERO
               ENDIF ! DO_K

            ENDIF  ! end if/else cut_u_cell_at

! Add the terms
            lTAU_U_G(IJK) = SBV + SSX + SSY + SSZ

! Also calculate and store the full gas phase viscous stress tensor
            CALL GET_FULL_TAU_U_G(IJK, ltau_u_g, lctau_u_g)
         ELSE
            lTAU_U_G(IJK) = ZERO
            lctau_u_g(IJK) = zero
         ENDIF
      ENDDO
!$omp end parallel do

      RETURN
      END SUBROUTINE CALC_CG_TAU_U_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GET_FULL_Tau_U_g                                        C
!  Purpose: Calculate the divergence of the complete gas phase stress  C
!  tensor including all terms in the u-direction.                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_FULL_TAU_U_G(IJK, ltau_u_g, lctau_u_g)

! Modules
!---------------------------------------------------------------------//
      USE fldvar, only: u_g

      USE functions, only: east_of, flow_at_e
      USE functions, only: im_of, jm_of, km_of
      USE functions, only: ip_of, jp_of, kp_of

      USE fun_avg, only: avg_x

      USE geometry, only: ox_e, do_k, cylindrical, vol_u
      USE indices, only: i_of

      USE param
      USE param1, only: zero
      USE visc_g, only: mu_gt, df_gu
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! ijk index
      INTEGER, INTENT(IN) :: IJK
! TAU_U_g
      DOUBLE PRECISION, INTENT(IN) :: lTAU_U_g(DIMENSION_3)
! cTAU_U_g
      DOUBLE PRECISION, INTENT(INOUT) :: lcTAU_U_g(DIMENSION_3)

! Local Variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, IJKE
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
! average viscosity
      DOUBLE PRECISION :: MUGA
! source terms
      DOUBLE PRECISION :: SSX, SSY, SSZ
! cylindrical terms
      DOUBLE PRECISION :: vtza
!---------------------------------------------------------------------//

      IPJK = IP_OF(IJK)
      IJPK = JP_OF(IJK)
      IJKP = KP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

! additional terms (found in source_u_g)
      VTZA = ZERO
      IF (CYLINDRICAL) THEN
         I = I_OF(IJK)
         IJKE = EAST_OF(IJK)
! part of -tau_zz/x xdxdydz =>
!         -(2mu/x)*(u/x) xdxdydz =>
! delta(-2.mu.u/x^2)V |p
         MUGA = AVG_X(MU_GT(IJK),MU_GT(IJKE),I)
         VTZA = 2.d0*MUGA*OX_E(I)*OX_E(I)*VOL_U(IJK)*U_G(IJK)
      ENDIF

! convection terms: see conv_dif_u_g
      SSX = ZERO
      SSY = ZERO
      SSZ = ZERO
      IF (FLOW_AT_E(IJK)) THEN
! part of 1/x d/dx (x.tau_xx) xdxdydz =>
!         1/x d/dx (x.mu.du/dx) xdxdydz =>
! delta (mu.du/dx.Ayz) |E-W : at (i+1 - i-1), j, k
         SSX = DF_GU(IJK,east)*(U_G(IPJK) - U_G(IJK)) - &
               DF_GU(IJK,west)*(U_G(IJK) - U_G(IJKM))

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.du/dy) xdxdydz =>
! delta (mu.du/dy.Axz) |N-S : at (i+1/2, j+1/2 - j-1/2, k)
         SSY = DF_GU(IJK,north)*(U_G(IJPK)-U_G(IJK)) - &
               DF_GU(IJK,south)*(U_G(IJK)-U_G(IJMK))

         IF (DO_K) THEN
! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu/x.du/dz) xdxdydz =>
! delta (mu/x.du/dz.Axy) |T-B : at (i+1/2, j, k+1/2 - k-1/2)
            SSZ = DF_GU(IJK,top)*(U_G(IJKP)-U_G(IJK)) - &
                  DF_GU(IJK,bottom)*(U_G(IJK)-U_G(IJKM))
         ENDIF
      ENDIF   ! end if flow_at_e

! Add the terms
      lctau_u_g(IJK) = (lTAU_U_G(IJK) + SSX + SSY + SSZ + VTZA)

      RETURN
      END SUBROUTINE GET_FULL_TAU_U_G
