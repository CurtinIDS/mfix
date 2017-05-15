!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_Tau_U_s                                            C
!  Purpose: Cross terms in the gradient of stress in U_s momentum      C
!                                                                      C
!  Comments: This routine calculates the components of the solids      C
!  phase viscous stress tensor of the u-momentum equation that cannot  C
!  be cast in the form: mu.grad(u). These components are stored in the C
!  passed variable, which is then used as a source of the u-momentum   C
!  equation.                                                           C
!  For greater details see calc_tau_u_g.                               C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TAU_U_S(lTAU_U_S)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3, dimension_m
      USE param1, only: zero
      USE physprop, only: smax, mmax
      USE fldvar, only: v_s
      USE fldvar, only: ep_s
      USE toleranc, only: dil_ep_s
      USE run, only: kt_type_enum, ghd_2007
      USE run, only: shear
      USE vshear, only: VSH
      USE cutcell, only: cartesian_grid, cg_safe_mode

      USE compar
      USE geometry
      USE indices
      USE functions
      USE sendrecv
      USE fun_avg
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! TAU_U_s
      DOUBLE PRECISION, INTENT(OUT) :: lTAU_U_s(DIMENSION_3, DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, I, IJKE
! Phase index
      INTEGER :: M, L
! Average volume fraction
      DOUBLE PRECISION :: EPSA, EPStmp
! Source terms (Surface)
      DOUBLE PRECISION :: Sbv, Ssx, Ssy, Ssz
! Cylindrical source terms
      DOUBLE PRECISION :: Vtzb
!---------------------------------------------------------------------//


      DO M = 1, MMAX
         IF(KT_TYPE_ENUM == GHD_2007 .AND. M /= MMAX) CYCLE

         IF (SHEAR) THEN
!!$omp  parallel do private(IJK)
            DO IJK = IJKSTART3, IJKEND3
               IF (FLUID_AT(IJK)) THEN
                  V_S(IJK,m)=V_S(IJK,m)+VSH(IJK)
               ENDIF
            ENDDO
         ENDIF

!!$omp  parallel do &
!!$omp  private(I, IJK, IJKE,                                         &
!!$omp&         EPSA, EPStmp, SSX, SSY, SSZ, SBV, VTZB)               &
!!$omp&  schedule(static)
         DO IJK = IJKSTART3, IJKEND3

! Skip walls where some values are undefined.
            IF(WALL_AT(IJK)) cycle
            I = I_OF(IJK)
            IJKE = EAST_OF(IJK)

            IF (KT_TYPE_ENUM == GHD_2007) THEN
               EPStmp = ZERO
               DO L = 1, SMAX
                 EPStmp = EPStmp + AVG_X(EP_S(IJK,L),EP_S(IJKE,L),I)
               ENDDO
               EPSA = EPStmp
            ELSE
               EPSA = AVG_X(EP_S(IJK,M),EP_S(IJKE,M),I)
            ENDIF

            IF ( .NOT.SIP_AT_E(IJK) .AND. EPSA>DIL_EP_S) THEN

               IF ((.NOT.CARTESIAN_GRID) .OR. CG_SAFE_MODE(3) ==1) THEN
                  CALL CALC_REG_TAU_U_S(IJK, M, SSX, SSY, SSZ, SBV, VTZB)
               ELSE
                  VTZB = ZERO   ! no cylindrical terms in CG
                  CALL CALC_CG_TAU_U_S(IJK, M, SSX, SSY, SSZ, SBV)
               ENDIF

! Add the terms
               lTAU_U_S(IJK,M) = SBV + SSX + SSY + SSZ + VTZB*VOL_U(IJK)
            ELSE
               lTAU_U_S(IJK,M) = ZERO
            ENDIF   ! end if/else .not.sip_at_e .and. epsa>dil_ep_s
         ENDDO   ! end do ijk
!!$omp end parallel do

         IF (SHEAR) THEN
!!$omp  parallel do private(IJK)
            DO IJK = IJKSTART3, IJKEND3
               IF (FLUID_AT(IJK)) THEN
                  V_S(IJK,m)=V_S(IJK,m)-VSH(IJK)
               ENDIF
            ENDDO
         ENDIF

      ENDDO   ! end do m=1,mmax

      call send_recv(ltau_u_s,2)
      RETURN
      END SUBROUTINE CALC_TAU_U_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: calc_reg_tau_u_s                                        C
!  Purpose: Cross terms in the gradient of stress in U_s momentum      C
!  based on NON cartesian grid cut cell.                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_REG_TAU_U_S(IJK, M, SSX, SSY, SSZ, SBV, VTZB)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, half
      USE fldvar, only: u_s, v_s, w_s
      USE visc_s, only: epmu_s, eplambda_s
      USE visc_s, only: trd_s

      USE compar
      USE geometry
      USE indices
      USE functions
      USE fun_avg

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
      DOUBLE PRECISION, INTENT(OUT) :: VTZB

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, JM, KM, IP
      INTEGER :: IJKE, IJKN, IJKS, IJKT, IJKB
      INTEGER :: IJKNE, IJKSE, IJKTE, IJKBE
      INTEGER :: IPJK, IMJK, IJKM, IJMK
      INTEGER :: IPJMK, IPJKM
! Average viscosity
      DOUBLE PRECISION :: MU_ste, MU_sbe, MUSA
! Average dW/Xdz
      DOUBLE PRECISION :: dWoXdz
!---------------------------------------------------------------------//

      I = I_OF(IJK)
      IP = IP1(I)
      J = J_OF(IJK)
      JM = JM1(J)
      K = K_OF(IJK)
      KM = KM1(K)

      IPJK = IP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IPJKM = IP_OF(IJKM)
      IPJMK = JM_OF(IPJK)

      IJKE = EAST_OF(IJK)
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
      SBV = (EPLAMBDA_S(IJKE,M)*TRD_S(IJKE,M)-&
             EPLAMBDA_S(IJK,M)*TRD_S(IJK,M))*AYZ(IJK)

! shear stress terms at i+1/2, j, k
! part of 1/x d/dx(x.tau_xx) xdxdydz =>
!         1/x d/dx (x.mu.du/dx) xdxdydz =>
! delta (mu du/dx)Ayz |E-W : at (i+1 - i-1), j, k
      SSX = EPMU_S(IJKE,M)*(U_S(IPJK,M)-U_S(IJK,M))*&
               ODX(IP)*AYZ_U(IJK) - &
            EPMU_S(IJK,M)*(U_S(IJK,M)-U_S(IMJK,M))*&
               ODX(I)*AYZ_U(IMJK)

! part of d/dy (tau_xy) xdxdydz =>
!         d/dy (mu.dv/dx) xdxdydz =>
! delta (mu.dv/dx)Axz |N-S : at i+1/2, (j+1/2 - j-1/2), k
      SSY = AVG_X_H(AVG_Y_H(EPMU_S(IJK,M),EPMU_S(IJKN,M),J),&
                    AVG_Y_H(EPMU_S(IJKE,M),EPMU_S(IJKNE,M),J),I)*&
               (V_S(IPJK,M)-V_S(IJK,M))*ODX_E(I)*AXZ_U(IJK) - &
            AVG_X_H(AVG_Y_H(EPMU_S(IJKS,M),EPMU_S(IJK,M),JM),&
                    AVG_Y_H(EPMU_S(IJKSE,M),EPMU_S(IJKE,M),JM),I)*&
               (V_S(IPJMK,M)-V_S(IJMK,M))*ODX_E(I)*AXZ_U(IJMK)

! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.dw/dx) xdxdydz =>
! delta (mu.dw/dx)Axy |T-B : at i+1/2, j, (k+1/2 - k-1/2)
      MU_STE = AVG_X_H(AVG_Z_H(EPMU_S(IJK,M),EPMU_S(IJKT,M),K),&
                       AVG_Z_H(EPMU_S(IJKE,M),EPMU_S(IJKTE,M),K),I)
      MU_SBE = AVG_X_H(AVG_Z_H(EPMU_S(IJKB,M),EPMU_S(IJK,M),KM),&
                       AVG_Z_H(EPMU_S(IJKBE,M),EPMU_S(IJKE,M),KM),I)
      SSZ = MU_STE*(W_S(IPJK,M)-W_S(IJK,M))*ODX_E(I)*&
                        AXY_U(IJK) - &
            MU_SBE*(W_S(IPJKM,M)-W_S(IJKM,M))*ODX_E(I)*&
                        AXY_U(IJKM)


! Special terms for cylindrical coordinates
      IF (CYLINDRICAL) THEN
! part of 1/x d/dz (tau_xz) xdxdydz =>
!         1/x d/dz (mu.(-w/x)) xdxdydz =>
! delta (mu(-w/x))Axy |T-B : at i+1/2, j, (k+1/2- k-1/2)
         SSZ = SSZ - (&
               MU_STE*(HALF*(W_S(IPJK,M)+W_S(IJK,M))*&
                  OX_E(I))*AXY_U(IJK)-&
               MU_SBE*(HALF*(W_S(IPJKM,M)+W_S(IJKM,M))*&
                  OX_E(I))*AXY_U(IJKM))

! part of -tau_zz/x xdxdydz =>
!         -(2.mu/x).(1/x).dw/dz xdxdydz
! delta (2.mu/x.1/x.dw/dz)Vp |p : at i+1/2, j, k
         MUSA = AVG_X(EPMU_S(IJK,M),EPMU_S(IJKE,M),I)
         DWOXDZ = HALF*&
                  ((W_S(IJK,M)-W_S(IJKM,M))*OX(I)*ODZ(K)+&
                   (W_S(IPJK,M)-W_S(IPJKM,M))*OX(IP)*ODZ(K))
         VTZB = -2.d0*MUSA*OX_E(I)*DWOXDZ
      ELSE
         VTZB = ZERO
      ENDIF

      RETURN
      END SUBROUTINE CALC_REG_TAU_U_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: calc_cg_tau_u_s                                         C
!  Purpose: Cross terms in the gradient of stress in U_s momentum      C
!  based on cartesian grid cut cell.                                   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_TAU_U_S(IJK, M, SSX, SSY, SSZ, SBV)

! Modules
!---------------------------------------------------------------------//

      USE compar
      USE fldvar, only: u_s, v_s, w_s
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE param1, only: zero, half, undefined
      USE visc_s, only: epmu_s, eplambda_s
      USE visc_s, only: trd_s

! for cutcell:
! wall velocity for slip conditions
      USE bc, only: bc_hw_s, bc_uw_s, bc_vw_s, bc_ww_s, cg_nsw, cg_fsw, cg_psw, cg_mi, none, bc_type_enum
      USE quadric
      USE cutcell

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

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, JM, KM, IP
      INTEGER :: IJKE, IJKN, IJKS, IJKT, IJKB
      INTEGER :: IJKNE, IJKSE, IJKTE, IJKBE
      INTEGER :: IPJK, IMJK, IJKM, IJMK
      INTEGER :: IPJMK, IPJKM
! Average viscosity
      DOUBLE PRECISION :: MU_ste, MU_sbe

! for cartesian grid implementation:
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      LOGICAL :: V_NODE_AT_NE, V_NODE_AT_NW, V_NODE_AT_SE, V_NODE_AT_SW
      LOGICAL :: W_NODE_AT_TE, W_NODE_AT_TW, W_NODE_AT_BE, W_NODE_AT_BW
      DOUBLE PRECISION :: dvdx_at_N, dvdx_at_S
      DOUBLE PRECISION :: dwdx_at_T, dwdx_at_B
      DOUBLE PRECISION :: Xi, Yi, Zi, Vi, Wi, Sx, Sy, Sz
      DOUBLE PRECISION :: MU_S_CUT, SSY_CUT, SSZ_CUT
      DOUBLE PRECISION :: UW_s, VW_s, WW_s
      INTEGER :: BCV
      INTEGER :: BCT
!---------------------------------------------------------------------//

      I = I_OF(IJK)
      IP = IP1(I)
      J = J_OF(IJK)
      JM = JM1(J)
      K = K_OF(IJK)
      KM = KM1(K)

      IPJK = IP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IPJKM = IP_OF(IJKM)
      IPJMK = JM_OF(IPJK)

      IJKE = EAST_OF(IJK)
      IJKN = NORTH_OF(IJK)
      IJKNE = EAST_OF(IJKN)
      IJKS = SOUTH_OF(IJK)
      IJKSE = EAST_OF(IJKS)
      IJKT = TOP_OF(IJK)
      IJKTE = EAST_OF(IJKT)
      IJKB = BOTTOM_OF(IJK)
      IJKBE = EAST_OF(IJKB)


! bulk viscosity term
      SBV = (EPLAMBDA_S(IJKE,M)*TRD_S(IJKE,M)) * AYZ_U(IJK) - &
            (EPLAMBDA_S(IJK,M) *TRD_S(IJK,M) ) * AYZ_U(IMJK)

! shear stress terms
      IF(.NOT.CUT_U_CELL_AT(IJK))   THEN
         SSX = EPMU_S(IJKE,M)*(U_S(IPJK,M)-U_S(IJK,M))*&
                  ONEoDX_E_U(IJK)*AYZ_U(IJK) - &
               EPMU_S(IJK,M)*(U_S(IJK,M)-U_S(IMJK,M))*&
                  ONEoDX_E_U(IMJK)*AYZ_U(IMJK)
         SSY = AVG_X_H(AVG_Y_H(EPMU_S(IJK,M),EPMU_S(IJKN,M),J),&
                       AVG_Y_H(EPMU_S(IJKE,M),EPMU_S(IJKNE,M),J),I)*&
                  (V_S(IPJK,M)-V_S(IJK,M))*&
                  ONEoDX_E_V(IJK)*AXZ_U(IJK) - &
               AVG_X_H(AVG_Y_H(EPMU_S(IJKS,M),EPMU_S(IJK,M),JM),&
                       AVG_Y_H(EPMU_S(IJKSE,M),EPMU_S(IJKE,M),JM),I)*&
                  (V_S(IPJMK,M)-V_S(IJMK,M))*&
                  ONEoDX_E_V(IJMK)*AXZ_U(IJMK)
         IF(DO_K) THEN
            MU_STE = AVG_X_H(AVG_Z_H(EPMU_S(IJK,M),EPMU_S(IJKT,M),K),&
                             AVG_Z_H(EPMU_S(IJKE,M),EPMU_S(IJKTE,M),K),I)
            MU_SBE = AVG_X_H(AVG_Z_H(EPMU_S(IJKB,M),EPMU_S(IJK,M),KM),&
                             AVG_Z_H(EPMU_S(IJKBE,M),EPMU_S(IJKE,M),KM),I)
            SSZ = MU_STE*(W_S(IPJK,M)-W_S(IJK,M))*&
                     ONEoDX_E_W(IJK)*AXY_U(IJK) - &
                  MU_SBE*(W_S(IPJKM,M)-W_S(IJKM,M))*&
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
               CUT_TAU_US = .TRUE.
               NOC_US     = .TRUE.
               UW_s = ZERO
               VW_s = ZERO
               WW_s = ZERO
            CASE (CG_FSW)
               CUT_TAU_US = .FALSE.
               NOC_US     = .FALSE.
               UW_s = ZERO
               VW_s = ZERO
               WW_s = ZERO
            CASE(CG_PSW)
               IF(BC_HW_S(BC_U_ID(IJK),M)==UNDEFINED) THEN   ! same as NSW
                  CUT_TAU_US = .TRUE.
                  NOC_US     = .TRUE.
                  UW_s = BC_UW_S(BCV,M)
                  VW_s = BC_VW_S(BCV,M)
                  WW_s = BC_WW_S(BCV,M)
               ELSEIF(BC_HW_S(BC_U_ID(IJK),M)==ZERO) THEN   ! same as FSW
                  CUT_TAU_US = .FALSE.
                  NOC_US     = .FALSE.
                  UW_s = ZERO
                  VW_s = ZERO
                  WW_s = ZERO
               ELSE                              ! partial slip
                  CUT_TAU_US = .FALSE.
                  NOC_US     = .FALSE.
               ENDIF
            CASE (NONE)
! equivalent of setting tau_u_s to zero at this ijk cell
               SSX = ZERO
               SSY = ZERO
               SSZ = ZERO
               SBV = ZERO
               RETURN
         END SELECT

         IF(CUT_TAU_US) THEN
            MU_S_CUT = (VOL(IJK)*EPMU_S(IJK,M) + &
                        VOL(IPJK)*EPMU_S(IJKE,M))/&
                       (VOL(IJK) + VOL(IPJK))
         ELSE
            MU_S_CUT = ZERO
         ENDIF


! SSX:
         CALL GET_DEL_H(IJK, 'U_MOMENTUM', X_U(IJK), &
            Y_U(IJK), Z_U(IJK), Del_H, Nx, Ny, Nz)
         SSX = EPMU_S(IJKE,M)*(U_S(IPJK,M)-U_S(IJK,M))*&
                  ONEoDX_E_U(IJK)*AYZ_U(IJK) - &
               EPMU_S(IJK,M)*(U_S(IJK,M)-U_S(IMJK,M))*&
                  ONEoDX_E_U(IMJK)*AYZ_U(IMJK) - &
               MU_S_CUT * (U_S(IJK,M) - UW_s) / DEL_H * &
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
            Vi = HALF * (V_S(IPJK,M) + V_S(IJK,M))
            Xi = HALF * (X_V(IPJK) + X_V(IJK))
            Yi = HALF * (Y_V(IPJK) + Y_V(IJK))
            Zi = HALF * (Z_V(IPJK) + Z_V(IJK))
            Sx = X_V(IPJK) - X_V(IJK)
            Sy = Y_V(IPJK) - Y_V(IJK)
            Sz = Z_V(IPJK) - Z_V(IJK)
            CALL GET_DEL_H(IJK, 'U_MOMENTUM', Xi, Yi, Zi, &
               Del_H, Nx, Ny, Nz)
            dvdx_at_N = (V_S(IPJK,M) - V_S(IJK,M)) * &
                         ONEoDX_E_V(IJK)
            IF(NOC_US) dvdx_at_N = dvdx_at_N - ((Vi-VW_s) *&
                          ONEoDX_E_V(IJK) /DEL_H*(Sy*Ny+Sz*Nz))
         ELSE
            dvdx_at_N =  ZERO
         ENDIF

         IF(V_NODE_AT_SE.AND.V_NODE_AT_SW) THEN
            Vi = HALF * (V_S(IPJMK,M) + V_S(IJMK,M))
            Xi = HALF * (X_V(IPJMK) + X_V(IJMK))
            Yi = HALF * (Y_V(IPJMK) + Y_V(IJMK))
            Zi = HALF * (Z_V(IPJMK) + Z_V(IJMK))
            Sx = X_V(IPJMK) - X_V(IJMK)
            Sy = Y_V(IPJMK) - Y_V(IJMK)
            Sz = Z_V(IPJMK) - Z_V(IJMK)
            CALL GET_DEL_H(IJK, 'U_MOMENTUM', Xi, Yi, Zi, &
               Del_H, Nx, Ny, Nz)
            dvdx_at_S = (V_S(IPJMK,M)-V_S(IJMK,M))*ONEoDX_E_V(IJMK)
            IF(NOC_US) dvdx_at_S = dvdx_at_S - ((Vi-VW_s) * &
                          ONEoDX_E_V(IJMK)/DEL_H*(Sy*Ny+Sz*Nz))
         ELSE
            dvdx_at_S =  ZERO
         ENDIF

         IF(V_NODE_AT_NW) THEN
            CALL GET_DEL_H(IJK, 'U_MOMENTUM', X_V(IJK), &
               Y_V(IJK), Z_V(IJK), Del_H, Nx, Ny, Nz)
            SSY_CUT = -MU_S_CUT * (V_S(IJK,M) - VW_s)/DEL_H*&
                       (Nx*Ny) * Area_U_CUT(IJK)
         ELSE
            SSY_CUT =  ZERO
         ENDIF

         SSY = AVG_X_H(AVG_Y_H(EPMU_S(IJK,M),EPMU_S(IJKN,M),J),&
                       AVG_Y_H(EPMU_S(IJKE,M),EPMU_S(IJKNE,M),J),I)*&
                  dvdx_at_N*AXZ_U(IJK) - &
               AVG_X_H(AVG_Y_H(EPMU_S(IJKS,M),EPMU_S(IJK,M),JM),&
                       AVG_Y_H(EPMU_S(IJKSE,M),EPMU_S(IJKE,M),JM),I)*&
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
               Wi = HALF * (W_S(IPJK,M) + W_S(IJK,M))
               Xi = HALF * (X_W(IPJK) + X_W(IJK))
               Yi = HALF * (Y_W(IPJK) + Y_W(IJK))
               Zi = HALF * (Z_W(IPJK) + Z_W(IJK))
               Sx = X_W(IPJK) - X_W(IJK)
               Sy = Y_W(IPJK) - Y_W(IJK)
               Sz = Z_W(IPJK) - Z_W(IJK)
               CALL GET_DEL_H(IJK, 'U_MOMENTUM', Xi, Yi, Zi, &
                  Del_H, Nx, Ny, Nz)
               dwdx_at_T = (W_S(IPJK,M) - W_S(IJK,M)) * &
                           ONEoDX_E_W(IJK)
               IF(NOC_US) dwdx_at_T = dwdx_at_T - ((Wi-WW_s)*&
                             ONEoDX_E_W(IJK)/DEL_H*(Sy*Ny+Sz*Nz))
            ELSE
               dwdx_at_T =  ZERO
            ENDIF

            IF(W_NODE_AT_BE.AND.W_NODE_AT_BW) THEN
               Wi = HALF * (W_S(IPJKM,M) + W_S(IJKM,M))
               Xi = HALF * (X_W(IPJKM) + X_W(IJKM))
               Yi = HALF * (Y_W(IPJKM) + Y_W(IJKM))
               Zi = HALF * (Z_W(IPJKM) + Z_W(IJKM))
               Sx = X_W(IPJKM) - X_W(IJKM)
               Sy = Y_W(IPJKM) - Y_W(IJKM)
               Sz = Z_W(IPJKM) - Z_W(IJKM)
               CALL GET_DEL_H(IJK, 'U_MOMENTUM', Xi, Yi, &
                  Zi, Del_H, Nx, Ny, Nz)
               dwdx_at_B = (W_S(IPJKM,M) - W_S(IJKM,M)) * &
                           ONEoDX_E_W(IJKM)
               IF(NOC_US) dwdx_at_B = dwdx_at_B - ((Wi-WW_s)*&
                             ONEoDX_E_W(IJKM)/DEL_H*(Sy*Ny+Sz*Nz))
            ELSE
               dwdx_at_B =  ZERO
            ENDIF

            IF(W_NODE_AT_TW) THEN
               CALL GET_DEL_H(IJK, 'U_MOMENTUM', X_W(IJK), &
                  Y_W(IJK), Z_W(IJK), Del_H, Nx, Ny, Nz)
               SSZ_CUT = -MU_S_CUT * (W_S(IJK,M) - WW_s) / &
                         DEL_H * (Nx*Nz) * Area_U_CUT(IJK)
            ELSE
               SSZ_CUT =  ZERO
            ENDIF

            MU_STE = AVG_X_H(AVG_Z_H(EPMU_S(IJK,M),EPMU_S(IJKT,M),K),&
                             AVG_Z_H(EPMU_S(IJKE,M),EPMU_S(IJKTE,M),K),I)
            MU_SBE = AVG_X_H(AVG_Z_H(EPMU_S(IJKB,M),EPMU_S(IJK,M),KM),&
                             AVG_Z_H(EPMU_S(IJKBE,M),EPMU_S(IJKE,M),KM),I)
            SSZ = MU_STE*dwdx_at_T*AXY_U(IJK) -  &
                  MU_SBE*dwdx_at_B*AXY_U(IJKM) + SSZ_CUT
         ELSE
            SSZ = ZERO
         ENDIF  ! end if do_k

      ENDIF  ! end if/else cut_cell

      RETURN
      END SUBROUTINE CALC_CG_TAU_U_S
