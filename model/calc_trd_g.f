!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_trD_g                                              C
!  Purpose: Calculate the trace of the gas phase rate of strain        C
!  tensor at i, j, k                                                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TRD_G(lTRD_G)

! Modules
!-----------------------------------------------
      USE param, only: DIMENSION_3
      USE param1, only: ZERO
      USE parallel
      USE geometry
      USE fldvar
      USE indices
      USE compar
      USE sendrecv
      USE functions
      USE cutcell
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! Strain rate tensor components for mth solids phase
      DOUBLE PRECISION, INTENT(OUT) :: ltrD_g(DIMENSION_3)

! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IMJK, IJMK, IJKM, IM

!-----------------------------------------------


!!$omp  parallel do private( IJK, I, J, K, IM, IMJK, IJMK, IJKM ) &
!!$omp& schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3
         IF (.NOT.WALL_AT(IJK)) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IM = IM1(I)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

            IF (.NOT. CUT_CELL_AT(IJK)) THEN
! at i, j, k:
! du/dx + dv/dy + 1/x dw/dz + u/x =
! 1/x d(xu)/dx + dv/dy + 1/x dw/dz  =
               lTRD_G(IJK) = &
                  (X_E(I)*U_G(IJK)-X_E(IM)*U_G(IMJK))*OX(I)*ODX(I) +&
                  (V_G(IJK)-V_G(IJMK))*ODY(J) + &
                  (W_G(IJK)-W_G(IJKM))*(OX(I)*ODZ(K))
            ELSE
               CALL CALC_CG_TRD_G(IJK,lTRD_G)
            ENDIF

         ELSE
            lTRD_G(IJK) = ZERO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_TRD_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_CG_trD_g                                           C
!  Purpose: Calculate the trace of the gas phase rate of strain        C
!  tensor at i, j, k with cartesian grid modifications                 C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_TRD_G(IJK,lTRD_G)

! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE cutcell
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE parallel
      USE param, only: DIMENSION_3
      USE param1, only: ZERO, HALF, UNDEFINED
      USE quadric
      USE sendrecv
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! index
      INTEGER, INTENT(IN) :: IJK
! Strain rate tensor components for mth solids phase
      DOUBLE PRECISION, INTENT(INOUT) :: ltrD_g(DIMENSION_3)

! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IMJK, IJMK, IJKM, IM

      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      DOUBLE PRECISION :: dudx,dvdy,dwdz
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: UW_g,VW_g,WW_g
      LOGICAL :: U_NODE_AT_E, U_NODE_AT_W
      LOGICAL :: V_NODE_AT_N, V_NODE_AT_S
      LOGICAL :: W_NODE_AT_T, W_NODE_AT_B
      INTEGER :: BCV
      INTEGER :: BCT
!-----------------------------------------------


      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IM = IM1(I)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

      IF(FLOW_AT(IJK)) THEN
         lTRD_G(IJK) = ZERO
         RETURN
      ENDIF

      BCV = BC_ID(IJK)
      IF(BCV > 0 ) THEN
         BCT = BC_TYPE_ENUM(BCV)
      ELSE
         BCT = NONE
      ENDIF
      SELECT CASE (BCT)
         CASE (CG_NSW)
            NOC_TRDG = .TRUE.
            UW_g = ZERO
            VW_g = ZERO
            WW_g = ZERO
         CASE (CG_FSW)
            NOC_TRDG = .FALSE.
            UW_g = ZERO
            VW_g = ZERO
            WW_g = ZERO
         CASE(CG_PSW)
            IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
               NOC_TRDG = .TRUE.
               UW_g = BC_UW_G(BCV)
               VW_g = BC_VW_G(BCV)
               WW_g = BC_WW_G(BCV)
            ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
               NOC_TRDG = .FALSE.
               UW_g = ZERO
               VW_g = ZERO
               WW_g = ZERO
            ELSE                              ! partial slip
               NOC_TRDG = .FALSE.
            ENDIF
         CASE (CG_MI)
            lTRD_G(IJK) = ZERO
            RETURN
         CASE (CG_PO)
            lTRD_G(IJK) = ZERO
            RETURN
         CASE (NONE)
            lTRD_G(IJK) = ZERO
            RETURN
      END SELECT


! du/dx
!=======================================================================
      U_NODE_AT_E = ((.NOT.BLOCKED_U_CELL_AT(IJK)) .AND.&
                     (.NOT.WALL_U_AT(IJK)))
      U_NODE_AT_W = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.&
                     (.NOT.WALL_U_AT(IMJK)))
      IF(U_NODE_AT_E.AND.U_NODE_AT_W) THEN
         Ui = HALF * (U_G(IJK) + U_G(IMJK))
         Xi = HALF * (X_U(IJK) + X_U(IMJK))
         Yi = HALF * (Y_U(IJK) + Y_U(IMJK))
         Zi = HALF * (Z_U(IJK) + Z_U(IMJK))
         Sx = X_U(IJK) - X_U(IMJK)
         Sy = Y_U(IJK) - Y_U(IMJK)
         Sz = Z_U(IJK) - Z_U(IMJK)
         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
         IF(abs(Sx) > ZERO) THEN
            dudx =  (U_G(IJK) - U_G(IMJK))/Sx
            IF(NOC_TRDG) dudx = dudx - ((Ui-UW_g)/(Sx*DEL_H)*&
               (Sy*Ny+Sz*Nz))
         ELSE
            dudx = ZERO
         ENDIF
      ELSE
         dudx = ZERO
      ENDIF

! dv/dy
!=======================================================================
      V_NODE_AT_N = ((.NOT.BLOCKED_V_CELL_AT(IJK)) .AND.&
                     (.NOT.WALL_V_AT(IJK)))
      V_NODE_AT_S = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.&
                     (.NOT.WALL_V_AT(IJMK)))
      IF(V_NODE_AT_N.AND.V_NODE_AT_S) THEN
         Vi = HALF * (V_G(IJK) + V_G(IJMK))
         Xi = HALF * (X_V(IJK) + X_V(IJMK))
         Yi = HALF * (Y_V(IJK) + Y_V(IJMK))
         Zi = HALF * (Z_V(IJK) + Z_V(IJMK))
         Sx = X_V(IJK) - X_V(IJMK)
         Sy = Y_V(IJK) - Y_V(IJMK)
         Sz = Z_V(IJK) - Z_V(IJMK)
         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
         IF(abs(Sy) > ZERO) THEN
            dvdy =  (V_G(IJK) - V_G(IJMK))/Sy
            IF(NOC_TRDG) dvdy = dvdy - ((Vi-VW_g)/(Sy*DEL_H)*&
               (Sx*Nx+Sz*Nz))
         ELSE
            dvdy =  ZERO
         ENDIF
      ELSEIF (V_NODE_AT_N.AND.(.NOT.V_NODE_AT_S).AND.NOC_TRDG) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_V(IJK),Y_V(IJK),&
                        Z_V(IJK),DEL_H,Nx,Ny,Nz)
         dvdy = (V_g(IJK) - VW_g) / DEL_H * Ny
      ELSEIF ((.NOT.V_NODE_AT_N).AND.V_NODE_AT_S.AND.NOC_TRDG) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_V(IJMK),Y_V(IJMK),&
                        Z_V(IJMK),DEL_H,Nx,Ny,Nz)
         dvdy = (V_g(IJMK) - VW_g) / DEL_H * Ny
      ELSE
         dvdy = ZERO
      ENDIF

! dw/dz
!=======================================================================
      IF(NO_K) THEN
         dwdz = ZERO
      ELSE
         W_NODE_AT_T = ((.NOT.BLOCKED_W_CELL_AT(IJK)) .AND.&
            (.NOT.WALL_W_AT(IJK)))
         W_NODE_AT_B = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.&
            (.NOT.WALL_W_AT(IJKM)))
         IF(W_NODE_AT_T.AND.W_NODE_AT_B) THEN
            Wi = HALF * (W_G(IJK) + W_G(IJKM))
            Xi = HALF * (X_W(IJK) + X_W(IJKM))
            Yi = HALF * (Y_W(IJK) + Y_W(IJKM))
            Zi = HALF * (Z_W(IJK) + Z_W(IJKM))
            Sx = X_W(IJK) - X_W(IJKM)
            Sy = Y_W(IJK) - Y_W(IJKM)
            Sz = Z_W(IJK) - Z_W(IJKM)
            CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
            IF(abs(Sz) > ZERO) THEN
               dwdz =  (W_G(IJK) - W_G(IJKM))/Sz
               IF(NOC_TRDG) dwdz = dwdz - ((Wi-WW_g)/(Sz*DEL_H)*&
                  (Sx*Nx+Sy*Ny))
            ELSE
               dwdz = ZERO
            ENDIF
         ELSEIF (W_NODE_AT_T.AND.(.NOT.W_NODE_AT_B).AND.NOC_TRDG) THEN
            CALL GET_DEL_H(IJK,'SCALAR',X_W(IJK),Y_W(IJK),&
                           Z_W(IJK),DEL_H,Nx,Ny,Nz)
            dwdz = (W_g(IJK) - WW_g) / DEL_H * Nz
         ELSEIF ((.NOT.W_NODE_AT_T).AND.W_NODE_AT_B.AND.NOC_TRDG) THEN
            CALL GET_DEL_H(IJK,'SCALAR',X_W(IJKM),Y_W(IJKM),&
                           Z_W(IJKM),DEL_H,Nx,Ny,Nz)
            dwdz = (W_g(IJKM) - WW_g) / DEL_H * Nz
         ELSE
            dwdz = ZERO
         ENDIF
      ENDIF   ! end if/else no_k

      lTRD_G(IJK) = dudx + dvdy + dwdz

      RETURN
      END SUBROUTINE CALC_CG_TRD_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GET_FACE_VEL_GAS                                        C
!  Purpose: Evaluate the velocity components at each of the faces of   C
!  a scalar cell                                                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_FACE_VEL_GAS(IJK, uge, ugw, ugn, ugs, ugt, ugb, ugcc, &
                                       vge, vgw, vgn, vgs, vgt, vgb, &
                                       wge, wgw, wgn, wgs, wgt, wgb, wgcc)

! Modules
!-----------------------------------------------------------------------
      use fldvar, only: u_g, v_g, w_g
      use functions, only: im_of, jm_of, km_of
      use functions, only: ip_of, jp_of, kp_of
      use fun_avg, only: avg_x, avg_x_e
      use fun_avg, only: avg_y, avg_y_n
      use fun_avg, only: avg_z, avg_z_t
      use geometry, only: cylindrical
      use indices, only: i_of, j_of, k_of
      use indices, only: im1, jm1, km1
      use param1, only: zero
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------------------------------
! index
      INTEGER, INTENT(IN) :: IJK
! U_g at the east (i+1/2, j, k) and west face (i-1/2, j, k)
      DOUBLE PRECISION, INTENT(OUT) :: UgE, UgW
! U_g at the north (i, j+1/2, k) and south face (i, j-1/2, k) 
      DOUBLE PRECISION, INTENT(OUT) :: UgN, UgS
! U_g at the top (i, j, k+1/2) and bottom face (i, j, k-1/2)
      DOUBLE PRECISION, INTENT(OUT) :: UgT, UgB
! U_g at the center of a scalar cell (i, j, k)
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION, INTENT(OUT) :: UgcC

! V_g at the east (i+1/2, j, k) and west face (i-1/2, j, k)
      DOUBLE PRECISION, INTENT(OUT) :: VgE, VgW
! V_g at the north (i, j+1/2, k) and south face (i, j-1/2, k)
      DOUBLE PRECISION, INTENT(OUT) :: VgN, VgS
! V_g at the top (i, j, k+1/2) and bottom face (i, j, k-1/2)
      DOUBLE PRECISION, INTENT(OUT) :: VgT, VgB

! W_g at the east (i+1/2, j, k) and west face (i-1/2, j, k)
      DOUBLE PRECISION, INTENT(OUT) :: WgE, WgW
! W_g at the north (i, j+1/2, k) and south face (i, j-1/2, k)
      DOUBLE PRECISION, INTENT(OUT) :: WgN, WgS
! W_g at the top (i, j, k+1/2) and bottom face (i, j, k-1/2)
      DOUBLE PRECISION, INTENT(OUT) :: WgT, WgB
! W_g at the center of a scalar cell (i, j, k).
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION, INTENT(OUT) :: WgcC

! Local variables
!-----------------------------------------------------------------------
! Cell indices
      INTEGER :: I, J, K, IM, JM, KM
      INTEGER :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK
      INTEGER :: IJMKP, IJMKM, IJPKM
!-----------------------------------------------------------------------

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IM = IM1(I)
      JM = JM1(J)
      KM = KM1(K)
      IMJK = IM_OF(IJK)
      IPJK = IP_OF(IJK)
      IJMK = JM_OF(IJK)
      IJPK = JP_OF(IJK)
      IJKM = KM_OF(IJK)
      IJKP = KP_OF(IJK)
      IMJPK = IM_OF(IJPK)
      IMJMK = IM_OF(IJMK)
      IMJKP = IM_OF(IJKP)
      IMJKM = IM_OF(IJKM)
      IPJKM = IP_OF(IJKM)
      IPJMK = IP_OF(IJMK)
      IJMKP = JM_OF(IJKP)
      IJMKM = JM_OF(IJKM)
      IJPKM = JP_OF(IJKM)

! Find fluid velocity values at faces of the cell
      UgN = AVG_Y(AVG_X_E(U_G(IMJK),U_G(IJK),I), &
                  AVG_X_E(U_G(IMJPK),U_G(IJPK),I),J)     !i, j+1/2, k
      UgS = AVG_Y(AVG_X_E(U_G(IMJMK),U_G(IJMK),I),&
                  AVG_X_E(U_G(IMJK),U_G(IJK),I),JM)      !i, j-1/2, k
      UgT = AVG_Z(AVG_X_E(U_G(IMJK),U_G(IJK),I),&
                  AVG_X_E(U_G(IMJKP),U_G(IJKP),I),K)     !i, j, k+1/2
      UgB = AVG_Z(AVG_X_E(U_G(IMJKM),U_G(IJKM),I),&
                  AVG_X_E(U_G(IMJK),U_G(IJK),I),KM)      !i, j, k-1/2
      UgE = U_G(IJK)
      UgW = U_G(IMJK)

      VgE = AVG_X(AVG_Y_N(V_G(IJMK),V_G(IJK)),&
                  AVG_Y_N(V_G(IPJMK),V_G(IPJK)),I)       !i+1/2, j, k
      VgW = AVG_X(AVG_Y_N(V_G(IMJMK),V_G(IMJK)),&
                  AVG_Y_N(V_G(IJMK),V_G(IJK)),IM)        !i-1/2, j, k
      VgT = AVG_Z(AVG_Y_N(V_G(IJMK),V_G(IJK)),&
                  AVG_Y_N(V_G(IJMKP),V_G(IJKP)),K)       !i, j, k+1/2
      VgB = AVG_Z(AVG_Y_N(V_G(IJMKM),V_G(IJKM)),&
                  AVG_Y_N(V_G(IJMK),V_G(IJK)),KM)        !i, j, k-1/2
      VgN = V_G(IJK)
      VgS = V_G(IJMK)

      WgN = AVG_Y(AVG_Z_T(W_G(IJKM),W_G(IJK)),&
                  AVG_Z_T(W_G(IJPKM),W_G(IJPK)),J)       !i, j+1/2, k
      WgS = AVG_Y(AVG_Z_T(W_G(IJMKM),W_G(IJMK)),&
                  AVG_Z_T(W_G(IJKM),W_G(IJK)),JM)        !i, j-1/2, k
      WgE = AVG_X(AVG_Z_T(W_G(IJKM),W_G(IJK)),&
                  AVG_Z_T(W_G(IPJKM),W_G(IPJK)),I)       !i+1/2, j, k
      WgW = AVG_X(AVG_Z_T(W_G(IMJKM),W_G(IMJK)),&
                  AVG_Z_T(W_G(IJKM),W_G(IJK)),IM)        !i-1/2, j, k
      WgT = W_G(IJK)
      WgB = W_G(IJKM)

      IF (CYLINDRICAL) THEN
         UgcC = AVG_X_E(U_G(IMJK),U_G(IJK),I)             !i, j, k
         WgcC = AVG_Z_T(W_G(IJKM),W_G(IJK))               !i, j, k
      ELSE
         UgcC = ZERO
         WgcC = ZERO
      ENDIF

      RETURN
      END SUBROUTINE GET_FACE_VEL_GAS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DERIV_VEL_GAS                                      C
!  Purpose: Calculate the gradient of the gas phase velocity and       C
!  the rate of strain tensor at i, j, k                                C
!                                                                      C
!         |  du/dx    du/dy   du/dz  |                                 C
! DELV =  |  dv/dx    dv/dy   dv/dz  |  =  dVi/dxj                     C
!         |  dw/dx    dw/dy   dw/dz  |                                 C
!                                                                      C
! 1/2(DELV+DELV^T) = 1/2(dVi/dxj+dVj/dxi)                              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DERIV_VEL_GAS(IJK, lVelGradG, lRateStrainG)

! Modules
!-----------------------------------------------------------------------
      use cutcell, only: cut_cell_at
      use geometry, only: odx, ody, odz
      use geometry, only: ox
      use indices, only: i_of, j_of, k_of
      use param1, only: half
! Dummy arguments
!-----------------------------------------------------------------------
! index
      INTEGER, INTENT(IN) :: IJK
! gradient of velocity
      DOUBLE PRECISION, INTENT(OUT) :: lVelGradG(3,3)    ! delV
! rate of strain tensor
      DOUBLE PRECISION, INTENT(OUT) :: lRateStrainG(3,3) ! D_g

! Local variables
!-----------------------------------------------------------------------
! Cell indices
      INTEGER :: I, J, K
! face values of velocity
      DOUBLE PRECISION :: uge, ugw, ugn, ugs, ugt, ugb, ugcc
      DOUBLE PRECISION :: vge, vgw, vgn, vgs, vgt, vgb
      DOUBLE PRECISION :: wge, wgw, wgn, wgs, wgt, wgb, wgcc
!-----------------------------------------------------------------------

      IF(.NOT.CUT_CELL_AT(IJK)) THEN
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

! Set the velocity components at the scalar faces
         CALL GET_FACE_VEL_GAS(IJK, uge, ugw, ugn, ugs, ugt, ugb, ugcc, &
                          vge, vgw, vgn, vgs, vgt, vgb, &
                          wge, wgw, wgn, wgs, wgt, wgb, wgcc)

! Gradient of gas velocity  at cell center (i, j, k)
         lVelGradG(1,1) = (UgE-UgW)*ODX(I)
         lVelGradG(1,2) = (UgN-UgS)*ODY(J)
         lVelGradG(1,3) = (UgT-UgB)*(OX(I)*ODZ(K))-WgcC*OX(I)

         lVelGradG(2,1) = (VgE-VgW)*ODX(I)
         lVelGradG(2,2) = (VgN-VgS)*ODY(J)
         lVelGradG(2,3) = (VgT-VgB)*(OX(I)*ODZ(K))

         lVelGradG(3,1) = (WgE-WgW)*ODX(I)
         lVelGradG(3,2) = (WgN-WgS)*ODY(J)
         lVelGradG(3,3) = (WgT-WgB)*(OX(I)*ODZ(K)) + UgcC*OX(I)


! Strain rate tensor, D_G, at cell center (i, j, k)
! or calculated as 0.5*(DelG+DelG^T)
         lRateStrainG(1,1) = (UgE-UgW)*ODX(I)
         lRateStrainG(1,2) = HALF*((UgN-UgS)*ODY(J)+&
                                   (VgE-VgW)*ODX(I))
         lRateStrainG(1,3) = HALF*((WgE-WgW)*ODX(I)+&
                                   (UgT-UgB)*(OX(I)*ODZ(K))-&
                             WgcC*OX(I))

         lRateStrainG(2,1) = lRateStrainG(1,2)
         lRateStrainG(2,2) = (VgN-VgS)*ODY(J)
         lRateStrainG(2,3) = HALF*((VgT-VgB)*(OX(I)*ODZ(K))+&
                                   (WgN-WgS)*ODY(J))

         lRateStrainG(3,1) = lRateStrainG(1,3)
         lRateStrainG(3,2) = lRateStrainG(2,3)
         lRateStrainG(3,3) = (WgT-WgB)*(OX(I)*ODZ(K)) +&
                             UgcC*OX(I)

      ELSE   ! cut cell
         CALL CALC_CG_DERIV_VEL_GAS(IJK,lVelGradG,lRateStrainG)
      ENDIF

      RETURN
      END SUBROUTINE CALC_DERIV_VEL_GAS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_CG_DERIV_VEL_GAS                                  C
!  Purpose: Calculate velocity derivatives in scalar cut-cell          C
!           gas phase                                                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 25-JAN-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_DERIV_VEL_GAS(IJK, lVelGradG, lRateStrainG)

! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE fldvar
      USE indices
      USE compar
      USE sendrecv
      USE bc
      USE cutcell
      USE quadric
      USE functions
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! index
      INTEGER, INTENT(IN) :: IJK
! velocity gradient gas phase
!         |  du/dx    du/dy   du/dz  |
! DELV =  |  dv/dx    dv/dy   dv/dz  |  =  dUi/dxj
!         |  dw/dx    dw/dy   dw/dz  |
      DOUBLE PRECISION, INTENT(OUT) :: lVelGradG(3,3)
! rate of strain tensor gas phase
      DOUBLE PRECISION, INTENT(OUT) :: lRateStrainG(3,3)

! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IMJK, IJMK, IJKM
! loop counters
      INTEGER :: P, Q
!
      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      DOUBLE PRECISION :: dudx,dudy,dudz
      DOUBLE PRECISION :: dvdx,dvdy,dvdz
      DOUBLE PRECISION :: dwdx,dwdy,dwdz
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: UW_g,VW_g,WW_g

      LOGICAL :: U_NODE_AT_E, U_NODE_AT_W
      LOGICAL :: V_NODE_AT_N, V_NODE_AT_S
      LOGICAL :: W_NODE_AT_T, W_NODE_AT_B
      INTEGER :: BCV
      INTEGER :: BCT
!-----------------------------------------------

! initialize
      lVelGradG = ZERO
      lRateStrainG = ZERO

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

      IF(FLOW_AT(IJK)) THEN
         lVelGradG = ZERO
         RETURN
      ENDIF

      BCV = BC_ID(IJK)
      IF(BCV > 0 ) THEN
         BCT = BC_TYPE_ENUM(BCV)
      ELSE
         BCT = NONE
      ENDIF
      SELECT CASE (BCT)
         CASE (CG_NSW)
            NOC_TRDG = .TRUE.
            UW_g = ZERO
            VW_g = ZERO
            WW_g = ZERO
         CASE (CG_FSW)
            NOC_TRDG = .FALSE.
            UW_g = ZERO
            VW_g = ZERO
            WW_g = ZERO
         CASE(CG_PSW)
            IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
               NOC_TRDG = .TRUE.
               UW_g = BC_UW_G(BCV)
               VW_g = BC_VW_G(BCV)
               WW_g = BC_WW_G(BCV)
            ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
               NOC_TRDG = .FALSE.
               UW_g = ZERO
               VW_g = ZERO
               WW_g = ZERO
            ELSE                              ! partial slip
               NOC_TRDG = .FALSE.
            ENDIF
         CASE (CG_MI)
            lVelGradG = ZERO
            RETURN
         CASE (CG_PO)
            lVelGradG = ZERO
            RETURN
         CASE (NONE)
            lVelGradG = ZERO
            RETURN
      END SELECT


! du/dx, du/dy, du/dz
!=======================================================================
      U_NODE_AT_E = ((.NOT.BLOCKED_U_CELL_AT(IJK)) .AND.&
                     (.NOT.WALL_U_AT(IJK)))
      U_NODE_AT_W = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.&
                     (.NOT.WALL_U_AT(IMJK)))

      IF(U_NODE_AT_E.AND.U_NODE_AT_W) THEN
         Ui = HALF * (U_G(IJK) + U_G(IMJK))
         Xi = HALF * (X_U(IJK) + X_U(IMJK))
         Yi = HALF * (Y_U(IJK) + Y_U(IMJK))
         Zi = HALF * (Z_U(IJK) + Z_U(IMJK))
         Sx = X_U(IJK) - X_U(IMJK)
         Sy = Y_U(IJK) - Y_U(IMJK)
         Sz = Z_U(IJK) - Z_U(IMJK)
         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
         IF(abs(Sx) > ZERO) THEN
            dudx =  (U_G(IJK) - U_G(IMJK))/Sx
            dudy =  ZERO
            dudz =  ZERO
            IF(NOC_TRDG) THEN
               dudx = dudx - ((Ui-UW_g)/(Sx*DEL_H)*(Sy*Ny+Sz*Nz))
               dudy = (Ui-UW_g) / DEL_H * Ny
               dudz = (Ui-UW_g) / DEL_H * Nz
            ENDIF
         ELSE
            dudx = ZERO
            dudy = ZERO
            dudz = ZERO
         ENDIF
      ELSEIF (U_NODE_AT_E.AND.(.NOT.U_NODE_AT_W).AND.NOC_TRDG) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_U(IJK),Y_U(IJK),Z_U(IJK),&
                        DEL_H,Nx,Ny,Nz)
         dudx = (U_g(IJK) - UW_g) / DEL_H * Nx
         dudy = (U_g(IJK) - UW_g) / DEL_H * Ny
         dudz = (U_g(IJK) - UW_g) / DEL_H * Nz
      ELSEIF ((.NOT.U_NODE_AT_E).AND.U_NODE_AT_W.AND.NOC_TRDG) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_U(IMJK),Y_U(IMJK),Z_U(IMJK),&
                        DEL_H,Nx,Ny,Nz)
         dudx = (U_g(IMJK) - UW_g) / DEL_H * Nx
         dudy = (U_g(IMJK) - UW_g) / DEL_H * Ny
         dudz = (U_g(IMJK) - UW_g) / DEL_H * Nz
      ELSE
         dudx = ZERO
         dudy = ZERO
         dudz = ZERO
      ENDIF
      lVelGradG(1,1) = dudx
      lVelGradG(1,2) = dudy
      lVelGradG(1,3) = dudz

! dv/dx, dv/dy, dv/dz
!=======================================================================
      V_NODE_AT_N = ((.NOT.BLOCKED_V_CELL_AT(IJK)) .AND.&
         (.NOT.WALL_V_AT(IJK)))
      V_NODE_AT_S = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.&
         (.NOT.WALL_V_AT(IJMK)))

      IF(V_NODE_AT_N.AND.V_NODE_AT_S) THEN
         Vi = HALF * (V_G(IJK) + V_G(IJMK))
         Xi = HALF * (X_V(IJK) + X_V(IJMK))
         Yi = HALF * (Y_V(IJK) + Y_V(IJMK))
         Zi = HALF * (Z_V(IJK) + Z_V(IJMK))
         Sx = X_V(IJK) - X_V(IJMK)
         Sy = Y_V(IJK) - Y_V(IJMK)
         Sz = Z_V(IJK) - Z_V(IJMK)
         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
         IF(abs(Sy) > ZERO) THEN
            dvdx =  ZERO
            dvdy =  (V_G(IJK) - V_G(IJMK))/Sy
            dvdz =  ZERO
            IF(NOC_TRDG) THEN
               dvdx = (Vi-VW_g) / DEL_H * Nx
               dvdy = dvdy - ((Vi-VW_g)/(Sy*DEL_H)*(Sx*Nx+Sz*Nz))
               dvdz = (Vi-VW_g) / DEL_H * Nz
            ENDIF
         ELSE
            dvdx =  ZERO
            dvdy =  ZERO
            dvdz =  ZERO
         ENDIF
      ELSEIF (V_NODE_AT_N.AND.(.NOT.V_NODE_AT_S).AND.NOC_TRDG) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_V(IJK),Y_V(IJK),Z_V(IJK),&
                        DEL_H,Nx,Ny,Nz)
         dvdx = (V_g(IJK) - VW_g) / DEL_H * Nx
         dvdy = (V_g(IJK) - VW_g) / DEL_H * Ny
         dvdz = (V_g(IJK) - VW_g) / DEL_H * Nz
      ELSEIF ((.NOT.V_NODE_AT_N).AND.V_NODE_AT_S.AND.NOC_TRDG) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_V(IJMK),Y_V(IJMK),Z_V(IJMK),&
                        DEL_H,Nx,Ny,Nz)
         dvdx = (V_g(IJMK) - VW_g) / DEL_H * Nx
         dvdy = (V_g(IJMK) - VW_g) / DEL_H * Ny
         dvdz = (V_g(IJMK) - VW_g) / DEL_H * Nz
      ELSE
         dvdx =  ZERO
         dvdy =  ZERO
         dvdz =  ZERO
      ENDIF
      lVelGradG(2,1) = dvdx
      lVelGradG(2,2) = dvdy
      lVelGradG(2,3) = dvdz

! dw/dx, dw/dy, dw/dz
!=======================================================================
      IF(NO_K) THEN
         dwdx = ZERO
         dwdy = ZERO
         dwdz = ZERO
      ELSE
         W_NODE_AT_T = ((.NOT.BLOCKED_W_CELL_AT(IJK)) .AND.&
            (.NOT.WALL_W_AT(IJK)))
         W_NODE_AT_B = ((.NOT.BLOCKED_W_CELL_AT(IJKM)).AND.&
            (.NOT.WALL_W_AT(IJKM)))
         IF(W_NODE_AT_T.AND.W_NODE_AT_B) THEN
            Wi = HALF * (W_G(IJK) + W_G(IJKM))
            Xi = HALF * (X_W(IJK) + X_W(IJKM))
            Yi = HALF * (Y_W(IJK) + Y_W(IJKM))
            Zi = HALF * (Z_W(IJK) + Z_W(IJKM))
            Sx = X_W(IJK) - X_W(IJKM)
            Sy = Y_W(IJK) - Y_W(IJKM)
            Sz = Z_W(IJK) - Z_W(IJKM)
            CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
            IF(abs(Sz) > ZERO) THEN
               dwdx =  ZERO
               dwdy =  ZERO
               dwdz =  (W_G(IJK) - W_G(IJKM))/Sz
               IF(NOC_TRDG) THEN
                  dwdx = (Wi-WW_g) / DEL_H * Nx
                  dwdy = (Wi-WW_g) / DEL_H * Ny
                  dwdz = dwdz - ((Wi-WW_g)/(Sz*DEL_H)*(Sx*Nx+Sy*Ny))
               ENDIF
            ELSE
               dwdx = ZERO
               dwdy = ZERO
               dwdz = ZERO
            ENDIF
         ELSEIF (W_NODE_AT_T.AND.(.NOT.W_NODE_AT_B).AND.NOC_TRDG) THEN
            CALL GET_DEL_H(IJK,'SCALAR',X_W(IJK),Y_W(IJK),Z_W(IJK),&
                           DEL_H,Nx,Ny,Nz)
            dwdx = (W_g(IJK) - WW_g) / DEL_H * Nx
            dwdy = (W_g(IJK) - WW_g) / DEL_H * Ny
            dwdz = (W_g(IJK) - WW_g) / DEL_H * Nz
         ELSEIF ((.NOT.W_NODE_AT_T).AND.W_NODE_AT_B.AND.NOC_TRDG) THEN
            CALL GET_DEL_H(IJK,'SCALAR',X_W(IJKM),Y_W(IJKM),Z_W(IJKM),&
                           DEL_H,Nx,Ny,Nz)
            dwdx = (W_g(IJKM) - WW_g) / DEL_H * Nx
            dwdy = (W_g(IJKM) - WW_g) / DEL_H * Ny
            dwdz = (W_g(IJKM) - WW_g) / DEL_H * Nz
         ELSE
            dwdx = ZERO
            dwdy = ZERO
            dwdz = ZERO
         ENDIF
      ENDIF
      lVelGradG(3,1) = dwdx
      lVelGradG(3,2) = dwdy
      lVelGradG(3,3) = dwdz

      DO P = 1,3
         DO Q = 1,3
            lRateStrainG(P,Q) = HALF * (lVelGradG(P,Q)+lVelGradG(Q,P))
         ENDDO
      ENDDO


      RETURN
      END SUBROUTINE CALC_CG_DERIV_VEL_GAS
