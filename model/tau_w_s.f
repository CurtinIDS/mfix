!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_W_s                                            C
!  Purpose: Cross terms in the gradient of stress in W_s momentum      c
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!                                                                      C
!                                                                      C
!  Comments: This routine calculates the components of the solids      C
!  phase viscous stress tensor of the w-momentum equation that cannot  C
!  be cast in the form: mu.grad(w). These components are stored in the C
!  passed variable, which is then used as a source of the w-momentum   C
!  equation.                                                           C
!  For greater details see calc_tau_w_g.                               C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TAU_W_S(lTAU_W_S)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3, dimension_m
      USE param1, only: zero
      USE physprop, only: smax, mmax
      USE fldvar, only: ep_s
      USE toleranc, only: dil_ep_s
      USE run, only: kt_type_enum, ghd_2007
      USE cutcell, only: cartesian_grid, cg_safe_mode

      USE functions
      USE fun_avg
      USE compar
      USE geometry
      USE indices
      USE sendrecv
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_W_s
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_W_s(DIMENSION_3, DIMENSION_M)

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, K, IJKT
! Phase index
      INTEGER :: M, L
! Average volume fraction
      DOUBLE PRECISION :: EPSA, EPStmp
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
! Cylindrical source terms (Volumetric)
      DOUBLE PRECISION :: Vxz
!---------------------------------------------------------------------//

      DO M = 1, MMAX
         IF(KT_TYPE_ENUM == GHD_2007 .AND. M /= MMAX) CYCLE

!!$omp  parallel do &
!!$omp  private(IJK, K, IJKT,                                         &
!!$omp&         EPSA, EPStmp, SSX, SSY, SSZ, SBV, VXZ)                &
!!$omp& schedule(static)

         DO IJK = IJKSTART3, IJKEND3

! Skip walls where some values are undefined.
            IF(WALL_AT(IJK)) cycle

            K = K_OF(IJK)
            IJKT = TOP_OF(IJK)

            IF (KT_TYPE_ENUM == GHD_2007) THEN
               EPStmp = ZERO
               DO L = 1, SMAX
                  EPStmp = EPStmp + AVG_Z(EP_S(IJK,L),EP_S(IJKT,L),K)
               ENDDO
               EPSA = EPStmp
            ELSE
               EPSA = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
            ENDIF

            IF ( .NOT.SIP_AT_T(IJK) .AND. EPSA>DIL_EP_S) THEN

               IF((.NOT.CARTESIAN_GRID).OR.(CG_SAFE_MODE(5)==1)) THEN
                   CALL CALC_REG_TAU_W_S(IJK, M, SSX, SSY, SSZ, SBV, VXZ)
               ELSE
                   VXZ = ZERO
                   CALL CALC_CG_TAU_W_S(IJK, M, SSX, SSY, SSZ, SBV)
               ENDIF

! Add the terms
               lTAU_W_S(IJK,M) = SBV + SSX + SSY + SSZ + VXZ*VOL_W(IJK)
            ELSE
               lTAU_W_S(IJK,M) = ZERO
            ENDIF   ! end if/else .not.sip_at_t .and. epsa>dil_ep_s

         ENDDO   ! end do ijk
!!$omp end parallel do

      ENDDO   ! end do m=1,mmax

      call send_recv(ltau_w_s,2)
      RETURN
      END SUBROUTINE CALC_TAU_W_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: calc_reg_tau_w_s                                        C
!  Purpose: Cross terms in the gradient of stress in W_s momentum      C
!  based on NON cartesian grid cut cell.                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_REG_TAU_W_S(IJK, M, SSX, SSY, SSZ, SBV, VXZ)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, half, one, undefined
      USE fldvar, only: u_s, v_s, w_s
      USE visc_s, only: epmu_s, eplambda_s
      USE visc_s, only: trd_s

      USE functions
      USE fun_avg
      USE compar
      USE geometry
      USE indices
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! current ijk index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M
! Source terms (surface)
      DOUBLE PRECISION, INTENT(OUT) :: SSX
      DOUBLE PRECISION, INTENT(OUT) :: SSY
      DOUBLE PRECISION, INTENT(OUT) :: SSZ
      DOUBLE PRECISION, INTENT(OUT) :: SBV
! Cylindrical source terms
      DOUBLE PRECISION, INTENT(OUT) :: VXZ

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JM, KP
      INTEGER :: IJKP, IMJK, IJMK, IJKM
      INTEGER :: IJKE, IJKW, IJKN, IJKS, IJKT
      INTEGER :: IJKNT, IJKST, IJKTE, IJKTW
      INTEGER :: IJMKP, IMJKP
! Average velocity gradients
      DOUBLE PRECISION :: duodz
!---------------------------------------------------------------------//

      K = K_OF(IJK)
      J = J_OF(IJK)
      I = I_OF(IJK)
      IM = IM1(I)
      JM = JM1(J)
      KP = KP1(K)

      IJKP = KP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IJMKP = JM_OF(IJKP)
      IMJKP = KP_OF(IMJK)

      IJKT = TOP_OF(IJK)
      IJKN = NORTH_OF(IJK)
      IJKNT = TOP_OF(IJKN)
      IJKS = SOUTH_OF(IJK)
      IJKST = TOP_OF(IJKS)
      IJKE = EAST_OF(IJK)
      IJKTE = EAST_OF(IJKT)
      IJKW = WEST_OF(IJK)
      IJKTW = WEST_OF(IJKT)

! Surface forces

! bulk viscosity term
! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (lambda.trcD) xdxdydz=>
! delta (lambda.trcD)Ap |T-B : at (i, j, k+1 - k-1)
      SBV = (EPLAMBDA_S(IJKT,M)*TRD_S(IJKT,M)-&
             EPLAMBDA_S(IJK,M)*TRD_S(IJK,M))*AXY(IJK)

! shear stress terms
! part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently
! part of (tau_xz/x + 1/x d/dx (x tau_xz) ) xdxdydz =>
!         1/x d/dx(mu.du/dz) xdxdydz =>
! delta (mu/x du/dz)Ayz |E-W : at (i+1/2-i-1/2, j, k+1/2)
      SSX = AVG_Z_H(AVG_X_H(EPMU_S(IJK,M),EPMU_S(IJKE,M),I),&
                    AVG_X_H(EPMU_S(IJKT,M),EPMU_S(IJKTE,M),I),K)*&
               (U_S(IJKP,M)-U_S(IJK,M))*OX_E(I)*ODZ_T(K)*AYZ_W(IJK) - &
            AVG_Z_H(AVG_X_H(EPMU_S(IJKW,M),EPMU_S(IJK,M),IM),&
                    AVG_X_H(EPMU_S(IJKTW,M),EPMU_S(IJKT,M),IM),K)*&
               (U_S(IMJKP,M)-U_S(IMJK,M))*ODZ_T(K)*DY(J)*(HALF*(DZ(K)+DZ(KP)))
! Same as oX_E(IM)*AYZ_W(IMJK), but avoids singularity

! part of d/dy (tau_zy) xdxdydz =>
!         d/dy (mu/x dv/dz) xdxdydz =>
! delta (mu/x dv/dz)Axz |N-S : at (i, j+1/2 - j-1/2, k+1/2)
      SSY = AVG_Z_H(AVG_Y_H(EPMU_S(IJK,M),EPMU_S(IJKN,M),J),&
                    AVG_Y_H(EPMU_S(IJKT,M),EPMU_S(IJKNT,M),J),K)*&
               (V_S(IJKP,M)-V_S(IJK,M))*OX(I)*ODZ_T(K)*AXZ_W(IJK) - &
            AVG_Z_H(AVG_Y_H(EPMU_S(IJKS,M),EPMU_S(IJK,M),JM),&
                    AVG_Y_H(EPMU_S(IJKST,M),EPMU_S(IJKT,M),JM),K)*&
               (V_S(IJMKP,M)-V_S(IJMK,M))*OX(I)*ODZ_T(K)*AXZ_W(IJMK)

! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (mu/x dw/dz) xdxdydz =>
! delta (mu/x dw/dz)Axy |T-B : at (i, j, k+1 - k-1)
      SSZ = EPMU_S(IJKT,M)*(W_S(IJKP,M)-W_S(IJK,M))*&
               OX(I)*ODZ(KP)*AXY_W(IJK) - &
            EPMU_S(IJK,M)*(W_S(IJK,M)-W_S(IJKM,M))*&
               OX(I)*ODZ(K)* AXY_W(IJKM)

! Special terms for cylindrical coordinates
      VXZ = ZERO
      IF (CYLINDRICAL) THEN

! part of 1/x d/dz (tau_zz) xdxdydz =>
!         1/x d/dz (2.mu/x u) xdxdydz =>
! delta (2.mu/x u)Axy |T-B : at (i, j, k+1 - k-1)
         SSZ = SSZ + &
            EPMU_S(IJKT,M)*(U_S(IJKP,M)+U_S(IMJKP,M))*OX(I)*AXY_W(IJK) - &
            EPMU_S(IJK,M)*(U_S(IJK,M)+U_S(IMJK,M))*OX(I)*AXY_W(IJKM)

! part of 1/x^2 d/dx (x^2 tau_xz) xdxdydz => or equivalently
! part of (tau_xz/x + 1/x d/dx (x tau_xz)) xdxdydz =>
!         1/x (mu/x du/dz) xdxdydz =>
! delta (1/x mu/x du/dz)Vp : at (i, j, k+1/2)
         IF (OX_E(IM) /= UNDEFINED) THEN
            DUODZ = (U_S(IMJKP,M)-U_S(IMJK,M))*OX_E(IM)*ODZ_T(K)
         ELSE
            DUODZ = ZERO
         ENDIF
         VXZ = AVG_Z(EPMU_S(IJK,M),EPMU_S(IJKT,M),K)*OX(I)*HALF*&
               ((U_S(IJKP,M)-U_S(IJK,M))*OX_E(I)*ODZ_T(K)+DUODZ)
      ENDIF

      RETURN
      END SUBROUTINE CALC_REG_TAU_W_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: calc_cg_tau_w_s                                         C
!  Purpose: Cross terms in the gradient of stress in W_s momentum      C
!  based on cartesian grid cut cell.                                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_TAU_W_S(IJK, M, SSX, SSY, SSZ, SBV)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, half, one, undefined
      USE fldvar, only: u_s, v_s, w_s
      USE visc_s, only: epmu_s, eplambda_s
      USE visc_s, only: trd_s

      USE functions
      USE fun_avg
      USE compar
      USE geometry
      USE indices
      USE sendrecv

! for cutcell:
! wall velocity for slip conditions
      USE bc, only: bc_hw_s, bc_uw_s, bc_vw_s, bc_ww_s
      USE bc
      USE quadric
      USE cutcell
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! current ijk index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M
! Source terms (surface)
      DOUBLE PRECISION, INTENT(OUT) :: SSX
      DOUBLE PRECISION, INTENT(OUT) :: SSY
      DOUBLE PRECISION, INTENT(OUT) :: SSZ
      DOUBLE PRECISION, INTENT(OUT) :: SBV

! Local Variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JM, KP
      INTEGER :: IJKP, IMJK, IJMK, IJKM
      INTEGER :: IJKE, IJKW, IJKN, IJKS, IJKT
      INTEGER :: IJKNT, IJKST, IJKTE, IJKTW
      INTEGER :: IJMKP, IMJKP

! for cartesian grid implementation:
      DOUBLE PRECISION :: DEL_H, Nx, Ny, Nz
      LOGICAL :: U_NODE_AT_ET, U_NODE_AT_EB, U_NODE_AT_WT, U_NODE_AT_WB
      LOGICAL :: V_NODE_AT_NT, V_NODE_AT_NB, V_NODE_AT_ST, V_NODE_AT_SB
      DOUBLE PRECISION :: dudz_at_E, dudz_at_W
      DOUBLE PRECISION :: dvdz_at_N, dvdz_at_S
      DOUBLE PRECISION :: MU_S_CUT, SSX_CUT, SSY_CUT
      DOUBLE PRECISION :: Xi, Yi, Zi, Ui, Vi, Sx, Sy, Sz
      DOUBLE PRECISION :: UW_s, VW_s, WW_s
      INTEGER :: BCV
      INTEGER :: BCT
!---------------------------------------------------------------------//

      K = K_OF(IJK)
      J = J_OF(IJK)
      I = I_OF(IJK)
      IM = IM1(I)
      JM = JM1(J)
      KP = KP1(K)

      IJKP = KP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IJMKP = JM_OF(IJKP)
      IMJKP = KP_OF(IMJK)

      IJKT = TOP_OF(IJK)
      IJKN = NORTH_OF(IJK)
      IJKNT = TOP_OF(IJKN)
      IJKS = SOUTH_OF(IJK)
      IJKST = TOP_OF(IJKS)
      IJKE = EAST_OF(IJK)
      IJKTE = EAST_OF(IJKT)
      IJKW = WEST_OF(IJK)
      IJKTW = WEST_OF(IJKT)

! Surface forces
! bulk viscosity term
      SBV =  (EPLAMBDA_S(IJKT,M)*TRD_S(IJKT,M)) * AXY_W(IJK) - &
             (EPLAMBDA_S(IJK,M) *TRD_S(IJK,M) ) * AXY_W(IJKM)

! shear stress terms
      IF(.NOT.CUT_W_CELL_AT(IJK)) THEN
         SSX = AVG_Z_H(AVG_X_H(EPMU_S(IJK,M),EPMU_S(IJKE,M),I),&
                       AVG_X_H(EPMU_S(IJKT,M),EPMU_S(IJKTE,M),I),K)*&
                  (U_S(IJKP,M)-U_S(IJK,M))*ONEoDZ_T_U(IJK)*AYZ_W(IJK) - &
               AVG_Z_H(AVG_X_H(EPMU_S(IJKW,M),EPMU_S(IJK,M),IM),&
                       AVG_X_H(EPMU_S(IJKTW,M),EPMU_S(IJKT,M),IM),K)*&
                  (U_S(IMJKP,M)-U_S(IMJK,M))*ONEoDZ_T_U(IMJK)*AYZ_W(IMJK)
         SSY = AVG_Z_H(AVG_Y_H(EPMU_S(IJK,M),EPMU_S(IJKN,M),J),&
                       AVG_Y_H(EPMU_S(IJKT,M),EPMU_S(IJKNT,M),J),K)*&
                  (V_S(IJKP,M)-V_S(IJK,M))*ONEoDZ_T_V(IJK)*AXZ_W(IJK) - &
               AVG_Z_H(AVG_Y_H(EPMU_S(IJKS,M),EPMU_S(IJK,M),JM),&
                       AVG_Y_H(EPMU_S(IJKST,M),EPMU_S(IJKT,M),JM),K)*&
                  (V_S(IJMKP,M)-V_S(IJMK,M))*ONEoDZ_T_V(IJMK)*AXZ_W(IJMK)
         SSZ = EPMU_S(IJKT,M)*(W_S(IJKP,M)-W_S(IJK,M))*&
                  ONEoDZ_T_W(IJK)*AXY_W(IJK) - &
               EPMU_S(IJK,M)*(W_S(IJK,M)-W_S(IJKM,M))*&
                  ONEoDZ_T_W(IJKM)*AXY_W(IJKM)

! cut cell modifications
!---------------------------------------------------------------------//
      ELSE

         BCV = BC_W_ID(IJK)
         IF(BCV > 0 ) THEN
            BCT = BC_TYPE_ENUM(BCV)
         ELSE
            BCT = NONE
         ENDIF

         SELECT CASE (BCT)
            CASE (CG_NSW,CG_MI)
               CUT_TAU_WS = .TRUE.
               NOC_WS     = .TRUE.
               UW_s = ZERO
               VW_s = ZERO
               WW_s = ZERO
            CASE (CG_FSW)
               CUT_TAU_WS = .FALSE.
               NOC_WS     = .FALSE.
               UW_s = ZERO
               VW_s = ZERO
               WW_s = ZERO
            CASE(CG_PSW)
               IF(BC_HW_S(BC_W_ID(IJK),M)==UNDEFINED) THEN   ! same as NSW
                  CUT_TAU_WS = .TRUE.
                  NOC_WS     = .TRUE.
                  UW_s = BC_UW_S(BCV,M)
                  VW_s = BC_VW_S(BCV,M)
                  WW_s = BC_WW_S(BCV,M)
               ELSEIF(BC_HW_S(BC_W_ID(IJK),M)==ZERO) THEN   ! same as FSW
                  CUT_TAU_WS = .FALSE.
                  NOC_WS     = .FALSE.
                  UW_s = ZERO
                  VW_s = ZERO
                  WW_s = ZERO
               ELSE                              ! partial slip
                  CUT_TAU_WS = .FALSE.
                  NOC_WS     = .FALSE.
               ENDIF
            CASE (NONE)
! equivalent of setting tau_w_s to zero at this ijk cell
               SSX = ZERO
               SSY = ZERO
               SSZ = ZERO
               SBV = ZERO
               RETURN
         END SELECT

         IF(CUT_TAU_WS) THEN
            MU_S_CUT = (VOL(IJK)*EPMU_S(IJK,M) + &
                        VOL(IJKP)*EPMU_S(IJKT,M))/(VOL(IJK) + VOL(IJKP))
         ELSE
            MU_S_CUT = ZERO
         ENDIF

! SSX:
         U_NODE_AT_ET = ((.NOT.BLOCKED_U_CELL_AT(IJKP)).AND.&
                         (.NOT.WALL_U_AT(IJKP)))
         U_NODE_AT_EB = ((.NOT.BLOCKED_U_CELL_AT(IJK)).AND.&
                         (.NOT.WALL_U_AT(IJK)))
         U_NODE_AT_WT = ((.NOT.BLOCKED_U_CELL_AT(IMJKP)).AND.&
                         (.NOT.WALL_U_AT(IMJKP)))
         U_NODE_AT_WB = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.&
                         (.NOT.WALL_U_AT(IMJK)))

         IF(U_NODE_AT_ET.AND.U_NODE_AT_EB) THEN
            Ui = HALF * (U_S(IJKP,M) + U_S(IJK,M))
            Xi = HALF * (X_U(IJKP) + X_U(IJK))
            Yi = HALF * (Y_U(IJKP) + Y_U(IJK))
            Zi = HALF * (Z_U(IJKP) + Z_U(IJK))
            Sx = X_U(IJKP) - X_U(IJK)
            Sy = Y_U(IJKP) - Y_U(IJK)
            Sz = Z_U(IJKP) - Z_U(IJK)
            CALL GET_DEL_H(IJK, 'W_MOMENTUM', Xi, Yi, Zi, Del_H, &
               Nx, Ny, Nz)
            dudz_at_E = (U_S(IJKP,M) - U_S(IJK,M)) * ONEoDZ_T_U(IJK)
            IF(NOC_WS) dudz_at_E = dudz_at_E - ((Ui-UW_s) * &
                          ONEoDZ_T_U(IJK)/DEL_H*(Sx*Nx+Sy*Ny))
         ELSE
            dudz_at_E =  ZERO
         ENDIF

         IF(U_NODE_AT_WT.AND.U_NODE_AT_WB) THEN
            Ui = HALF * (U_S(IMJKP,M) + U_S(IMJK,M))
            Xi = HALF * (X_U(IMJKP) + X_U(IMJK))
            Yi = HALF * (Y_U(IMJKP) + Y_U(IMJK))
            Zi = HALF * (Z_U(IMJKP) + Z_U(IMJK))
            Sx = X_U(IMJKP) - X_U(IMJK)
            Sy = Y_U(IMJKP) - Y_U(IMJK)
            Sz = Z_U(IMJKP) - Z_U(IMJK)
            CALL GET_DEL_H(IJK, 'W_MOMENTUM', Xi, Yi, Zi, Del_H, &
               Nx, Ny, Nz)
            dudz_at_W = (U_S(IMJKP,M) - U_S(IMJK,M)) * ONEoDZ_T_U(IMJK)
            IF(NOC_WS) dudz_at_W = dudz_at_W - ((Ui-UW_s) * &
                          ONEoDZ_T_U(IMJK)/DEL_H*(Sx*Nx+Sy*Ny))
         ELSE
            dudz_at_W =  ZERO
         ENDIF

         IF(U_NODE_AT_EB) THEN
            CALL GET_DEL_H(IJK, 'W_MOMENTUM', X_U(IJK), Y_U(IJK), &
               Z_U(IJK), Del_H, Nx, Ny, Nz)
            SSX_CUT = -MU_S_CUT * (U_S(IJK,M) - UW_s) / DEL_H * &
                      (Nz*Nx) * Area_W_CUT(IJK)
         ELSE
            SSX_CUT =  ZERO
         ENDIF

         SSX = AVG_Z_H(AVG_X_H(EPMU_S(IJK,M),EPMU_S(IJKE,M),I),&
                       AVG_X_H(EPMU_S(IJKT,M),EPMU_S(IJKTE,M),I),K)*&
                  dudz_at_E*AYZ_W(IJK) - &
               AVG_Z_H(AVG_X_H(EPMU_S(IJKW,M),EPMU_S(IJK,M),IM),&
                       AVG_X_H(EPMU_S(IJKTW,M),EPMU_S(IJKT,M),IM),K)*&
                  dudz_at_W*AYZ_W(IMJK) + SSX_CUT

! SSY:
         V_NODE_AT_NT = ((.NOT.BLOCKED_V_CELL_AT(IJKP)).AND.&
                         (.NOT.WALL_V_AT(IJKP)))
         V_NODE_AT_NB = ((.NOT.BLOCKED_V_CELL_AT(IJK)).AND.&
                         (.NOT.WALL_V_AT(IJK)))
         V_NODE_AT_ST = ((.NOT.BLOCKED_V_CELL_AT(IJMKP)).AND.&
                         (.NOT.WALL_V_AT(IJMKP)))
         V_NODE_AT_SB = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.&
                         (.NOT.WALL_V_AT(IJMK)))

         IF(V_NODE_AT_NT.AND.V_NODE_AT_NB) THEN
            Vi = HALF * (V_S(IJKP,M) + V_S(IJK,M))
            Xi = HALF * (X_V(IJKP) + X_V(IJK))
            Yi = HALF * (Y_V(IJKP) + Y_V(IJK))
            Zi = HALF * (Z_V(IJKP) + Z_V(IJK))
            Sx = X_V(IJKP) - X_V(IJK)
            Sy = Y_V(IJKP) - Y_V(IJK)
            Sz = Z_V(IJKP) - Z_V(IJK)
            CALL GET_DEL_H(IJK, 'W_MOMENTUM', Xi, Yi, Zi, Del_H, &
               Nx, Ny, Nz)
            dvdz_at_N = (V_S(IJKP,M) - V_S(IJK,M)) * ONEoDZ_T_V(IJK)
            IF(NOC_WS) dvdz_at_N = dvdz_at_N - ((Vi-VW_s) * &
                          ONEoDZ_T_V(IJK)/DEL_H*(Sx*Nx+Sy*Ny))
         ELSE
            dvdz_at_N =  ZERO
         ENDIF

         IF(V_NODE_AT_ST.AND.V_NODE_AT_SB) THEN
            Vi = HALF * (V_S(IJMKP,M) + V_S(IJMK,M))
            Xi = HALF * (X_V(IJMKP) + X_V(IJMK))
            Yi = HALF * (Y_V(IJMKP) + Y_V(IJMK))
            Zi = HALF * (Z_V(IJMKP) + Z_V(IJMK))
            Sx = X_V(IJMKP) - X_V(IJMK)
            Sy = Y_V(IJMKP) - Y_V(IJMK)
            Sz = Z_V(IJMKP) - Z_V(IJMK)
            CALL GET_DEL_H(IJK, 'W_MOMENTUM', Xi, Yi, Zi, Del_H, &
               Nx, Ny, Nz)
            dvdz_at_S = (V_S(IJMKP,M) - V_S(IJMK,M)) * ONEoDZ_T_V(IJMK)
            IF(NOC_WS) dvdz_at_S = dvdz_at_S - ((Vi-VW_s) * &
                          ONEoDZ_T_V(IJMK)/DEL_H*(Sx*Nx+Sy*Ny))
         ELSE
            dvdz_at_S =  ZERO
         ENDIF

         IF(V_NODE_AT_NB) THEN
            CALL GET_DEL_H(IJK, 'W_MOMENTUM', X_V(IJK), Y_V(IJK),&
               Z_V(IJK), Del_H, Nx, Ny, Nz)
            SSY_CUT = -MU_S_CUT * (V_S(IJK,M) - VW_s) / DEL_H * &
                      (Nz*Ny) * Area_W_CUT(IJK)
         ELSE
            SSY_CUT =  ZERO
         ENDIF

         SSY = AVG_Z_H(AVG_Y_H(EPMU_S(IJK,M),EPMU_S(IJKN,M),J),&
                       AVG_Y_H(EPMU_S(IJKT,M),EPMU_S(IJKNT,M),J),K)*&
                  dvdz_at_N*AXZ_W(IJK) - &
               AVG_Z_H(AVG_Y_H(EPMU_S(IJKS,M),EPMU_S(IJK,M),JM),&
                       AVG_Y_H(EPMU_S(IJKST,M),EPMU_S(IJKT,M),JM),K)*&
                  dvdz_at_S*AXZ_W(IJMK) + SSY_CUT

! SSZ:
         CALL GET_DEL_H(IJK, 'W_MOMENTUM', X_W(IJK), Y_W(IJK), &
            Z_W(IJK), Del_H, Nx, Ny, Nz)
         SSZ = EPMU_S(IJKT,M)*(W_S(IJKP,M)-W_S(IJK,M))*&
                  ONEoDZ_T_W(IJK)*AXY_W(IJK) - &
               EPMU_S(IJK,M)*(W_S(IJK,M)-W_S(IJKM,M))*&
                  OX(I)*ONEoDZ_T_W(IJKM)*AXY_W(IJKM) &
               -MU_S_CUT * (W_S(IJK,M) - WW_s) / DEL_H * &
                  (Nz**2) * Area_W_CUT(IJK)

      ENDIF  ! end if/else cut_cell

      RETURN
      END SUBROUTINE CALC_CG_TAU_W_S
