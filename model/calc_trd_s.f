!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_trD_s                                              C
!  Purpose: Calculate the trace of the solids phase rate of strain     C
!  tensor at i, j, k                                                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-DEC-96  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_TRD_S(lTRD_S)

! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE fldvar
      USE indices
      USE physprop
      USE compar
      USE sendrecv
      USE functions
      USE cutcell
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! Strain rate tensor components for mth solids phase
      DOUBLE PRECISION, INTENT(OUT) :: ltrD_s(DIMENSION_3, DIMENSION_M)

! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K
      INTEGER :: IJK, IMJK, IJMK, IJKM
      INTEGER :: IM, M

!-----------------------------------------------

      DO M = 1, MMAX
!!$omp    parallel do private(ijk,i,j,k,im,imjk,ijmk,ijkm)
         DO IJK = ijkstart3, ijkend3
            IF (.NOT.WALL_AT(IJK)) THEN
               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)
               IM = IM1(I)
               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)

               IF(.NOT.CUT_CELL_AT(IJK)) THEN
! at i, j, k:
! du/dx + dv/dy + 1/x dw/dz + u/x =
! 1/x d(xu)/dx + dv/dy + 1/x dw/dz =
                  lTRD_S(IJK,M) = (X_E(I)*U_S(IJK,M)-&
                     X_E(IM)*U_S(IMJK,M))*OX(I)*ODX(I) + &
                     (V_S(IJK,M)-V_S(IJMK,M))*ODY(J) + &
                     (W_S(IJK,M)-W_S(IJKM,M))*(OX(I)*ODZ(K))
               ELSE
                  CALL CALC_CG_TRD_S(IJK,M,ltrD_s)
               ENDIF

            ELSE
               ltrD_s(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

      END SUBROUTINE CALC_TRD_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_CG_trD_s                                           C
!  Purpose: Calculate the trace of the solids phase rate of strain     C
!  tensor at i, j, k with cartesian grid modifications                 C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_TRD_S(IJK,M,ltrD_s)

! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE geometry
      USE fldvar
      USE indices
      USE physprop
      USE compar
      USE sendrecv
      USE functions
      USE cutcell

      USE bc
      USE cutcell
      USE quadric
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! index
      INTEGER, INTENT(IN) :: IJK
! phase index
      INTEGER, INTENT(IN) :: M
! Strain rate tensor components for mth solids phase
      DOUBLE PRECISION, INTENT(INOUT) :: ltrD_s(DIMENSION_3, DIMENSION_M)

! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K
      INTEGER :: IMJK, IJMK, IJKM
      INTEGER :: IM

      DOUBLE PRECISION :: DEL_H,Nx,Ny,Nz
      DOUBLE PRECISION :: dudx,dvdy,dwdz
      DOUBLE PRECISION :: Xi,Yi,Zi,Ui,Vi,Wi,Sx,Sy,Sz
      DOUBLE PRECISION :: UW_s,VW_s,WW_s

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
          lTRD_S(IJK,M) = ZERO
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
            NOC_TRDS = .TRUE.
            UW_s = ZERO
            VW_s = ZERO
            WW_s = ZERO
         CASE (CG_FSW)
            NOC_TRDS = .FALSE.
            UW_s = ZERO
            VW_s = ZERO
            WW_s = ZERO
         CASE(CG_PSW)
            IF(BC_JJ_PS(BCV)==1) THEN   ! Johnson-Jackson partial slip bc
               NOC_TRDS = .FALSE.
               UW_s = BC_UW_S(BCV,M)
               VW_s = BC_VW_S(BCV,M)
               WW_s = BC_WW_S(BCV,M)
            ELSEIF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
               NOC_TRDS = .TRUE.
               UW_s = BC_UW_S(BCV,M)
               VW_s = BC_VW_S(BCV,M)
               WW_s = BC_WW_S(BCV,M)
            ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN    ! same as FSW
               NOC_TRDS = .FALSE.
            ELSE                                 ! partial slip
               NOC_TRDS = .FALSE.
               UW_s = ZERO
               VW_s = ZERO
               WW_s = ZERO
            ENDIF
         CASE (CG_MI)
            lTRD_S(IJK,M) = ZERO
            RETURN
         CASE (CG_PO)
            lTRD_S(IJK,M) = ZERO
            RETURN
         CASE (NONE)
            lTRD_S(IJK,M) = ZERO
            RETURN
      END SELECT


! du/dx
!=======================================================================
      U_NODE_AT_E = ((.NOT.BLOCKED_U_CELL_AT(IJK)) .AND.&
                     (.NOT.WALL_U_AT(IJK)))
      U_NODE_AT_W = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.&
                     (.NOT.WALL_U_AT(IMJK)))
      IF(U_NODE_AT_E.AND.U_NODE_AT_W) THEN
         Ui = HALF * (U_S(IJK,M) + U_S(IMJK,M))
         Xi = HALF * (X_U(IJK) + X_U(IMJK))
         Yi = HALF * (Y_U(IJK) + Y_U(IMJK))
         Zi = HALF * (Z_U(IJK) + Z_U(IMJK))
         Sx = X_U(IJK) - X_U(IMJK)
         Sy = Y_U(IJK) - Y_U(IMJK)
         Sz = Z_U(IJK) - Z_U(IMJK)
         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
         IF(Sx /= ZERO) THEN
            dudx =  (U_S(IJK,M) - U_S(IMJK,M))/Sx
            IF(NOC_TRDS) dudx = dudx - ((Ui-UW_s)/(Sx*DEL_H)*&
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
         Vi = HALF * (V_S(IJK,M) + V_S(IJMK,M))
         Xi = HALF * (X_V(IJK) + X_V(IJMK))
         Yi = HALF * (Y_V(IJK) + Y_V(IJMK))
         Zi = HALF * (Z_V(IJK) + Z_V(IJMK))
         Sx = X_V(IJK) - X_V(IJMK)
         Sy = Y_V(IJK) - Y_V(IJMK)
         Sz = Z_V(IJK) - Z_V(IJMK)
         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
         IF(Sy /= ZERO) THEN
            dvdy =  (V_S(IJK,M) - V_S(IJMK,M))/Sy
            IF(NOC_TRDS) dvdy = dvdy - ((Vi-VW_s)/(Sy*DEL_H)*&
               (Sx*Nx+Sz*Nz))
         ELSE
            dvdy =  ZERO
         ENDIF
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
            Wi = HALF * (W_S(IJK,M) + W_S(IJKM,M))
            Xi = HALF * (X_W(IJK) + X_W(IJKM))
            Yi = HALF * (Y_W(IJK) + Y_W(IJKM))
            Zi = HALF * (Z_W(IJK) + Z_W(IJKM))
            Sx = X_W(IJK) - X_W(IJKM)
            Sy = Y_W(IJK) - Y_W(IJKM)
            Sz = Z_W(IJK) - Z_W(IJKM)
            CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
            IF(Sz /= ZERO) THEN
               dwdz =  (W_S(IJK,M) - W_S(IJKM,M))/Sz
               IF(NOC_TRDS) dwdz = dwdz - ((Wi-WW_s)/(Sz*DEL_H)*&
                  (Sx*Nx+Sy*Ny))
            ELSE
               dwdz = ZERO
            ENDIF
         ELSE
            dwdz = ZERO
         ENDIF
      ENDIF  ! NO_K

      lTRD_S(IJK,M) = dudx + dvdy + dwdz

      RETURN
      END SUBROUTINE CALC_CG_TRD_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GET_FACE_VEL_SOLIDS                                     C
!  Purpose: Evaluate the velocity components at each of the faces of   C
!  a scalar cell                                                       C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_FACE_VEL_SOLIDS(IJK, M, us_e, usw, usn, uss, ust, usb, uscc, &
                                       vse, vsw, vsn, vss, vst, vsb, &
                                       wse, wsw, wsn, wss, wst, wsb, wscc)

! Modules
!-----------------------------------------------------------------------
      use fldvar, only: u_s, v_s, w_s
      use functions, only: im_of, jm_of, km_of
      use functions, only: ip_of, jp_of, kp_of
      use functions, only: is_at_e, is_at_n, is_at_t
      use functions, only: wall_at
      use fun_avg, only: avg_x, avg_x_e
      use fun_avg, only: avg_y, avg_y_n
      use fun_avg, only: avg_z, avg_z_t
      use geometry, only: cylindrical
      use geometry, only: xlength, odx_e
      use indices, only: i_of, j_of, k_of
      use indices, only: im1, jm1, km1
      USE is, only: any_is_defined
      use param1, only: zero, one
      use run, only: shear, V_SH
      Use vshear, only: VSH
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------------------------------
! index
      INTEGER, INTENT(IN) :: IJK
! solids index
      INTEGER, INTENT(IN) :: M
! U_s at the east (i+1/2, j, k) and west face (i-1/2, j, k)
      DOUBLE PRECISION, INTENT(OUT) :: Us_E, UsW
! U_s at the north (i, j+1/2, k) and south face (i, j-1/2, k)
      DOUBLE PRECISION, INTENT(OUT) :: UsN, UsS
! U_s at the top (i, j, k+1/2) and bottom face (i, j, k-1/2)
      DOUBLE PRECISION, INTENT(OUT) :: UsT, UsB
! U_s at the center of a scalar cell (i, j, k)
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION, INTENT(OUT) :: UscC

! V_s at the east (i+1/2, j, k) and west face (i-1/2, j, k)
      DOUBLE PRECISION, INTENT(OUT) :: VsE, VsW
! V_s at the north (i, j+1/2, k) and south face (i, j-1/2, k)
      DOUBLE PRECISION, INTENT(OUT) :: VsN, VsS
! V_s at the top (i, j, k+1/2) and bottom face (i, j, k-1/2)
      DOUBLE PRECISION, INTENT(OUT) :: VsT, VsB

! W_s at the east (i+1/2, j, k) and west face (i-1/2, j, k)
      DOUBLE PRECISION, INTENT(OUT) :: WsE, WsW
! W_s at the north (i, j+1/2, k) and south face (i, j-1/2, k)
      DOUBLE PRECISION, INTENT(OUT) :: WsN, WsS
! W_s at the top (i, j, k+1/2) and bottom face (i, j, k-1/2)
      DOUBLE PRECISION, INTENT(OUT) :: WsT, WsB
! W_s at the center of a scalar cell (i, j, k).
! Calculated for Cylindrical coordinates only.
      DOUBLE PRECISION, INTENT(OUT) :: WscC

! Local variables
!-----------------------------------------------------------------------
! Cell indices
      INTEGER :: I, J, K, IM, JM, KM
      INTEGER :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: IMJPK, IMJMK, IMJKP, IMJKM, IPJKM, IPJMK
      INTEGER :: IJMKP, IJMKM, IJPKM
! shear value calculation
      DOUBLE PRECISION :: SRT
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

      UsN = AVG_Y(AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I),&
                  AVG_X_E(U_s(IMJPK, M), U_s(IJPK, M), I), J)   !i, j+1/2, k
      UsS = AVG_Y(AVG_X_E(U_s(IMJMK, M), U_s(IJMK, M), I),&
                  AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I), JM)    !i, j-1/2, k
      UsT = AVG_Z(AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I),&
                  AVG_X_E(U_s(IMJKP, M), U_s(IJKP, M), I), K)   !i, j, k+1/2
      UsB = AVG_Z(AVG_X_E(U_s(IMJKM, M), U_s(IJKM, M), I),&
                  AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I), KM)    !i, j, k-1/2
      Us_E = U_S(IJK,M)
      UsW = U_S(IMJK,M)

      IF (SHEAR)  THEN
         SRT=(2d0*V_sh/XLENGTH)
         VsE = AVG_X(AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
                     AVG_Y_N((V_s(IPJMK, M) - VSH(IPJMK) + &
                              VSH(IJMK) + SRT*ONE/oDX_E(I)), &
                             (V_s(IPJK, M) - VSH(IPJK) + &
                              VSH(IJK) + SRT*ONE/oDX_E(I))), I) !i+1/2, j, k
         VsW = AVG_X(AVG_Y_N((V_s(IMJMK, M) - VSH(IMJMK) + &
                              VSH(IJMK) - SRT*ONE/oDX_E(IM)),&
                             (V_s(IMJK, M) - VSH(IMJK) + &
                              VSH(IJK) - SRT*ONE/oDX_E(IM)) ),&
                     AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), IM)    !i-1/2, j, k
      ELSE
         VsE = AVG_X(AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
                     AVG_Y_N(V_s(IPJMK, M), V_s(IPJK, M)), I )  !i+1/2, j, k
         VsW = AVG_X(AVG_Y_N(V_s(IMJMK, M), V_s(IMJK, M)),&
                     AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), IM )   !i-1/2, j, k
      ENDIF
      VsT = AVG_Z(AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)),&
                  AVG_Y_N(V_s(IJMKP, M), V_s(IJKP, M)), K )     !i, j, k+1/2
      VsB = AVG_Z(AVG_Y_N(V_s(IJMKM, M), V_s(IJKM, M)),&
                  AVG_Y_N(V_s(IJMK, M), V_s(IJK, M)), KM )      !i, j, k-1/2
      VsN = V_S(IJK,M)
      VsS = V_S(IJMK,M)

      WsN = AVG_Y(AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)),&
                  AVG_Z_T(W_s(IJPKM, M), W_s(IJPK, M)), J )     !i, j+1/2, k
      WsS = AVG_Y(AVG_Z_T(W_s(IJMKM, M), W_s(IJMK, M)),&
                  AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)), JM )      !i, j-1/2, k
      WsE = AVG_X(AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)),&
                  AVG_Z_T(W_s(IPJKM, M), W_s(IPJK, M)), I)      !i+1/2, j, k
      WsW = AVG_X(AVG_Z_T(W_s(IMJKM, M), W_s(IMJK, M)),&
                  AVG_Z_T(W_s(IJKM, M), W_s(IJK, M)), IM )      !i-1/2, j, k
      WsT = W_S(IJK,M)
      WsB = W_S(IJKM,M)

      IF(CYLINDRICAL) THEN
         UscC = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)          !i, j, k
         WscC = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))             !i, j, k
      ELSE
         UscC = ZERO
         WscC = ZERO
      ENDIF

! Check for IS surfaces and modify solids velocity-comp accordingly
      IF(ANY_IS_DEFINED) THEN
         IF(IS_AT_N(IJK)  .AND. .NOT.WALL_AT(IJPK)) &
            UsN = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
         IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) &
            UsS = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
         IF(IS_AT_T(IJK)  .AND. .NOT.WALL_AT(IJKP)) &
            UsT = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)
         IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) &
            UsB = AVG_X_E(U_s(IMJK, M), U_s(IJK, M), I)

         IF(IS_AT_E(IJK)  .AND. .NOT.WALL_AT(IPJK)) &
            VsE = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
         IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) &
            VsW = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
         IF(IS_AT_T(IJK)  .AND. .NOT.WALL_AT(IJKP)) &
            VsT = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))
         IF(IS_AT_T(IJKM) .AND. .NOT.WALL_AT(IJKM)) &
            VsB = AVG_Y_N(V_s(IJMK, M), V_s(IJK, M))

         IF(IS_AT_N(IJK)  .AND. .NOT.WALL_AT(IJPK)) &
            WsN = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
         IF(IS_AT_N(IJMK) .AND. .NOT.WALL_AT(IJMK)) &
            WsS = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
         IF(IS_AT_E(IJK)  .AND. .NOT.WALL_AT(IPJK)) &
            WsE = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
         IF(IS_AT_E(IMJK) .AND. .NOT.WALL_AT(IMJK)) &
            WsW = AVG_Z_T(W_s(IJKM, M), W_s(IJK, M))
      ENDIF

      RETURN
      END SUBROUTINE GET_FACE_VEL_SOLIDS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DERIV_VEL_SOLIDS                                   C
!  Purpose: Calculate the gradient of the solids phase velocity and    C
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
      SUBROUTINE CALC_DERIV_VEL_SOLIDS(IJK, M, lVelGradS, lRateStrainS)

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
! solids phase index
      INTEGER, INTENT(IN) :: M
! gradient of velocity
      DOUBLE PRECISION, INTENT(OUT) :: lVelGradS(3,3)    ! delV
! rate of strain tensor
      DOUBLE PRECISION, INTENT(OUT) :: lRateStrainS(3,3) ! D_s

! Local variables
!-----------------------------------------------------------------------
! Cell indices
      INTEGER :: I, J, K
! face values of velocity
      DOUBLE PRECISION :: us_e, usw, usn, uss, ust, usb, uscc
      DOUBLE PRECISION :: vse, vsw, vsn, vss, vst, vsb
      DOUBLE PRECISION :: wse, wsw, wsn, wss, wst, wsb, wscc
!-----------------------------------------------------------------------

      IF(.NOT.CUT_CELL_AT(IJK)) THEN
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

! Set the velocity components at the scalar faces
         CALL GET_FACE_VEL_SOLIDS(IJK, M, &
            us_e, usw, usn, uss, ust, usb, uscc, &
             vse, vsw, vsn, vss, vst, vsb, &
             wse, wsw, wsn, wss, wst, wsb, wscc)

! Gradient of gas velocity  at cell center (i, j, k)
         lVelGradS(1,1) = (Us_E-UsW)*ODX(I)
         lVelGradS(1,2) = (UsN-UsS)*ODY(J)
         lVelGradS(1,3) = (UsT-UsB)*(OX(I)*ODZ(K))-WscC*OX(I)

         lVelGradS(2,1) = (VsE-VsW)*ODX(I)
         lVelGradS(2,2) = (VsN-VsS)*ODY(J)
         lVelGradS(2,3) = (VsT-VsB)*(OX(I)*ODZ(K))

         lVelGradS(3,1) = (WsE-WsW)*ODX(I)
         lVelGradS(3,2) = (WsN-WsS)*ODY(J)
         lVelGradS(3,3) = (WsT-WsB)*(OX(I)*ODZ(K)) + UscC*OX(I)


! Strain rate tensor, D_S, at cell center (i, j, k)
! or calculated as 0.5*(DelS+DelS^T)
         lRateStrainS(1,1) = (Us_E-UsW)*ODX(I)
         lRateStrainS(1,2) = HALF*((UsN-UsS)*ODY(J)+&
                                   (VsE-VsW)*ODX(I))
         lRateStrainS(1,3) = HALF*((WsE-WsW)*ODX(I)+&
                                   (UsT-UsB)*(OX(I)*ODZ(K))-&
                             WscC*OX(I))

         lRateStrainS(2,1) = lRateStrainS(1,2)
         lRateStrainS(2,2) = (VsN-VsS)*ODY(J)
         lRateStrainS(2,3) = HALF*((VsT-VsB)*(OX(I)*ODZ(K))+&
                                   (WsN-WsS)*ODY(J))

         lRateStrainS(3,1) = lRateStrainS(1,3)
         lRateStrainS(3,2) = lRateStrainS(2,3)
         lRateStrainS(3,3) = (WsT-WsB)*(OX(I)*ODZ(K)) +&
                             UscC*OX(I)

      ELSE   ! cut cell
         CALL CALC_CG_DERIV_VEL_SOLIDS(IJK,M,lVelGradS,lRateStrainS)
      ENDIF

      RETURN
      END SUBROUTINE CALC_DERIV_VEL_SOLIDS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_CG_DERIV_VEL_SOLIDS                               C
!  Purpose: Calculate velocity derivatives in scalar cut-cell          C
!           Solids phase                                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 25-JAN-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_CG_DERIV_VEL_SOLIDS(IJK, M, lVelGradS, lRateStrainS)

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
! solids phase index
      INTEGER, INTENT(IN) :: M
! velocity gradient of solids phase M
!         |  du/dx    du/dy   du/dz  |
! DELV =  |  dv/dx    dv/dy   dv/dz  |  =  dUi/dxj
!         |  dw/dx    dw/dy   dw/dz  |
      DOUBLE PRECISION, INTENT(OUT) :: lVelGradS(3,3)
! rate of strain tensor solids phase
      DOUBLE PRECISION, INTENT(OUT) :: lRateStrainS(3,3)

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
      DOUBLE PRECISION :: UW_s,VW_s,WW_s

      LOGICAL :: U_NODE_AT_E, U_NODE_AT_W
      LOGICAL :: V_NODE_AT_N, V_NODE_AT_S
      LOGICAL :: W_NODE_AT_T, W_NODE_AT_B
      INTEGER :: BCV
      INTEGER :: BCT
!-----------------------------------------------

!!$omp  parallel do private(IJK, I, J, K, IMJK, IJMK, IJKM) &
!!!!$omp& schedule(dynamic,chunk_size)

! initialize
      lVelGradS = ZERO
      lRateStrainS = ZERO

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

      IF(FLOW_AT(IJK)) THEN
         lVelGradS = ZERO
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
            NOC_TRDS = .TRUE.
            UW_s = ZERO
            VW_s = ZERO
            WW_s = ZERO
         CASE (CG_FSW)
            NOC_TRDS = .FALSE.
            UW_s = ZERO
            VW_s = ZERO
            WW_s = ZERO
            RETURN
         CASE(CG_PSW)
            IF(BC_JJ_PS(BCV)==1) THEN   ! Johnson-Jackson partial slip bc
               NOC_TRDS = .FALSE.
               UW_s = BC_UW_S(BCV,M)
               VW_s = BC_VW_S(BCV,M)
               WW_s = BC_WW_S(BCV,M)
            ELSEIF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
               NOC_TRDS = .TRUE.
               UW_s = BC_UW_S(BCV,M)
               VW_s = BC_VW_S(BCV,M)
               WW_s = BC_WW_S(BCV,M)
            ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
               NOC_TRDS = .FALSE.
               UW_s = ZERO
               VW_s = ZERO
               WW_s = ZERO
            ELSE                              ! partial slip
               NOC_TRDS = .FALSE.
            ENDIF
         CASE (CG_MI)
            lVelGradS = ZERO
            RETURN
         CASE (CG_PO)
            lVelGradS = ZERO
            RETURN
         CASE (NONE)
            lVelGradS = ZERO
            RETURN
      END SELECT



! du/dx, du/dy, du/dz
!=======================================================================
      U_NODE_AT_E = ((.NOT.BLOCKED_U_CELL_AT(IJK)) .AND.&
                     (.NOT.WALL_U_AT(IJK)))
      U_NODE_AT_W = ((.NOT.BLOCKED_U_CELL_AT(IMJK)).AND.&
                     (.NOT.WALL_U_AT(IMJK)))

      IF(U_NODE_AT_E.AND.U_NODE_AT_W) THEN
         Ui = HALF * (U_S(IJK,M) + U_S(IMJK,M))
         Xi = HALF * (X_U(IJK) + X_U(IMJK))
         Yi = HALF * (Y_U(IJK) + Y_U(IMJK))
         Zi = HALF * (Z_U(IJK) + Z_U(IMJK))
         Sx = X_U(IJK) - X_U(IMJK)
         Sy = Y_U(IJK) - Y_U(IMJK)
         Sz = Z_U(IJK) - Z_U(IMJK)
         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
         IF(abs(Sx) > ZERO) THEN
            dudx =  (U_S(IJK,M) - U_S(IMJK,M))/Sx
            dudy =  ZERO
            dudz =  ZERO
            IF(NOC_TRDS) THEN
               dudx = dudx - ((Ui-UW_s)/(Sx*DEL_H)*(Sy*Ny+Sz*Nz))
               dudy = (Ui-UW_s) / DEL_H * Ny
               dudz = (Ui-UW_s) / DEL_H * Nz
            ENDIF
         ELSE
            dudx = ZERO
            dudy = ZERO
            dudz = ZERO
         ENDIF
      ELSEIF (U_NODE_AT_E.AND.(.NOT.U_NODE_AT_W).AND.NOC_TRDS) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_U(IJK),Y_U(IJK),Z_U(IJK),&
                        DEL_H,Nx,Ny,Nz)
         dudx = (U_s(IJK,M) - UW_s) / DEL_H * Nx
         dudy = (U_s(IJK,M) - UW_s) / DEL_H * Ny
         dudz = (U_s(IJK,M) - UW_s) / DEL_H * Nz
      ELSEIF ((.NOT.U_NODE_AT_E).AND.U_NODE_AT_W.AND.NOC_TRDS) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_U(IMJK),Y_U(IMJK),Z_U(IMJK),&
                        DEL_H,Nx,Ny,Nz)
         dudx = (U_s(IMJK,M) - UW_s) / DEL_H * Nx
         dudy = (U_s(IMJK,M) - UW_s) / DEL_H * Ny
         dudz = (U_s(IMJK,M) - UW_s) / DEL_H * Nz
      ELSE
         dudx = ZERO
         dudy = ZERO
         dudz = ZERO
      ENDIF
      lVelGradS(1,1) = dudx
      lVelGradS(1,2) = dudy
      lVelGradS(1,3) = dudz

! dv/dx, dv/dy, dv/dz
!=======================================================================
      V_NODE_AT_N = ((.NOT.BLOCKED_V_CELL_AT(IJK)) .AND.&
                     (.NOT.WALL_V_AT(IJK)))
      V_NODE_AT_S = ((.NOT.BLOCKED_V_CELL_AT(IJMK)).AND.&
                     (.NOT.WALL_V_AT(IJMK)))

      IF(V_NODE_AT_N.AND.V_NODE_AT_S) THEN
         Vi = HALF * (V_S(IJK,M) + V_S(IJMK,M))
         Xi = HALF * (X_V(IJK) + X_V(IJMK))
         Yi = HALF * (Y_V(IJK) + Y_V(IJMK))
         Zi = HALF * (Z_V(IJK) + Z_V(IJMK))
         Sx = X_V(IJK) - X_V(IJMK)
         Sy = Y_V(IJK) - Y_V(IJMK)
         Sz = Z_V(IJK) - Z_V(IJMK)
         CALL GET_DEL_H(IJK,'SCALAR',Xi,Yi,Zi,Del_H,Nx,Ny,Nz)
         IF(abs(Sy) > ZERO) THEN
            dvdx =  ZERO
            dvdy =  (V_S(IJK,M) - V_S(IJMK,M))/Sy
            dvdz =  ZERO
            IF(NOC_TRDS) THEN
               dvdx = (Vi-VW_s) / DEL_H * Nx
               dvdy = dvdy - ((Vi-VW_s)/(Sy*DEL_H)*(Sx*Nx+Sz*Nz))
               dvdz = (Vi-VW_s) / DEL_H * Nz
            ENDIF
         ELSE
            dvdx =  ZERO
            dvdy =  ZERO
            dvdz =  ZERO
         ENDIF
      ELSEIF (V_NODE_AT_N.AND.(.NOT.V_NODE_AT_S).AND.NOC_TRDS) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_V(IJK),Y_V(IJK),Z_V(IJK),&
                        DEL_H,Nx,Ny,Nz)
         dvdx = (V_s(IJK,M) - VW_s) / DEL_H * Nx
         dvdy = (V_s(IJK,M) - VW_s) / DEL_H * Ny
         dvdz = (V_s(IJK,M) - VW_s) / DEL_H * Nz
      ELSEIF ((.NOT.V_NODE_AT_N).AND.V_NODE_AT_S.AND.NOC_TRDS) THEN
         CALL GET_DEL_H(IJK,'SCALAR',X_V(IJMK),Y_V(IJMK),Z_V(IJMK),&
                        DEL_H,Nx,Ny,Nz)
         dvdx = (V_s(IJMK,M) - VW_s) / DEL_H * Nx
         dvdy = (V_s(IJMK,M) - VW_s) / DEL_H * Ny
         dvdz = (V_s(IJMK,M) - VW_s) / DEL_H * Nz
      ELSE
         dvdx =  ZERO
         dvdy =  ZERO
         dvdz =  ZERO
      ENDIF
      lVelGradS(2,1) = dvdx
      lVelGradS(2,2) = dvdy
      lVelGradS(2,3) = dvdz

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
            Wi = HALF * (W_S(IJK,M) + W_S(IJKM,M))
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
               dwdz =  (W_S(IJK,M) - W_S(IJKM,M))/Sz
               IF(NOC_TRDS) THEN
                  dwdx = (Wi-WW_s) / DEL_H * Nx
                  dwdy = (Wi-WW_s) / DEL_H * Ny
                  dwdz = dwdz - ((Wi-WW_s)/(Sz*DEL_H)*(Sx*Nx+Sy*Ny))
               ENDIF
            ELSE
               dwdx = ZERO
               dwdy = ZERO
               dwdz = ZERO
            ENDIF
         ELSEIF (W_NODE_AT_T.AND.(.NOT.W_NODE_AT_B).AND.NOC_TRDS) THEN
            CALL GET_DEL_H(IJK,'SCALAR',X_W(IJK),Y_W(IJK),Z_W(IJK),&
                           DEL_H,Nx,Ny,Nz)
            dwdx = (W_s(IJK,M) - WW_s) / DEL_H * Nx
            dwdy = (W_s(IJK,M) - WW_s) / DEL_H * Ny
            dwdz = (W_s(IJK,M) - WW_s) / DEL_H * Nz
         ELSEIF ((.NOT.W_NODE_AT_T).AND.W_NODE_AT_B.AND.NOC_TRDS) THEN
            CALL GET_DEL_H(IJK,'SCALAR',X_W(IJKM),Y_W(IJKM),Z_W(IJKM),&
                           DEL_H,Nx,Ny,Nz)
            dwdx = (W_s(IJKM,M) - WW_s) / DEL_H * Nx
            dwdy = (W_s(IJKM,M) - WW_s) / DEL_H * Ny
            dwdz = (W_s(IJKM,M) - WW_s) / DEL_H * Nz
         ELSE
            dwdx = ZERO
            dwdy = ZERO
            dwdz = ZERO
         ENDIF
      ENDIF
      lVelGradS(3,1) = dwdx
      lVelGradS(3,2) = dwdy
      lVelGradS(3,3) = dwdz

      DO P = 1,3
         DO Q = 1,3
            lRateStrainS(P,Q) = HALF * (lVelGradS(P,Q)+lVelGradS(Q,P))
         ENDDO
      ENDDO


      RETURN
      END SUBROUTINE CALC_CG_DERIV_VEL_SOLIDS
