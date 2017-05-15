!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_3D_ALPHA_U_CUT_CELL                                C
!  Purpose: Calculate the correction term alpha_U                      C
!           for U-momentum cut cells,                                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_3D_ALPHA_U_CUT_CELL

      USE bc
      USE compar, ONLY: iend1, jend1, kend1, IJKSTART3, IJKEND3, mype, pe_io
      USE cutcell
      USE exit, only: mfix_exit
      USE functions, ONLY: FUNIJK
      USE geometry, ONLY: DO_K, NO_K, ayz, ayz_u, flag_e
      USE indices, ONLY: I_OF, J_OF, K_OF
      USE param1, ONLY: half, one, undefined, zero

      IMPLICIT NONE
      DOUBLE PRECISION:: Xe,Ye,Ze,Xn,Yn,Zn,Xt,Yt,Zt
      INTEGER :: I,J,K,IM,JM,KM,IP,JP,KP,IJK,IMJK,IJMK,IPJK,IJPK,IJKM,IJKP
      DOUBLE PRECISION :: Xi
      DOUBLE PRECISION :: Sx,Sy,Sz
      DOUBLE PRECISION :: Nx,Ny,Nz
      DOUBLE PRECISION :: DELH_ec , DELH_e
      DOUBLE PRECISION :: DELH_nc , DELH_n
      DOUBLE PRECISION :: DELH_tc , DELH_t
      LOGICAL :: V_NODE_AT_NE,V_NODE_AT_NW
      LOGICAL :: W_NODE_AT_TE,W_NODE_AT_TW
      INTEGER :: BCV
      INTEGER ::BCT

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'COMPUTING INTERPOLATION FACTORS IN U-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)

      DELH_U = UNDEFINED

      Theta_Ue  = HALF
      Theta_Ue_bar = HALF
      Theta_U_ne  = HALF
      Theta_U_nw  = HALF
      Theta_U_te  = HALF
      Theta_U_tw  = HALF
      ALPHA_Ue_c  = ONE
      NOC_U_E  = ZERO
      Theta_Un  = HALF
      Theta_Un_bar = HALF
      ALPHA_Un_c  = ONE
      NOC_U_N  = ZERO
      Theta_Ut  = HALF
      Theta_Ut_bar = HALF
      ALPHA_Ut_c  = ONE
      NOC_U_T  = ZERO

      A_UPG_E = AYZ_U
      A_UPG_W = AYZ_U


      DO IJK = IJKSTART3, IJKEND3
         IF(INTERIOR_CELL_AT(IJK)) THEN

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IM = I - 1
            JM = J - 1
            KM = K - 1

            IP = I + 1
            JP = J + 1
            KP = K + 1

            IMJK = FUNIJK(IM,J,K)
            IPJK = FUNIJK(IP,J,K)
            IJMK = FUNIJK(I,JM,K)
            IJPK = FUNIJK(I,JP,K)
            IJKM = FUNIJK(I,J,KM)
            IJKP = FUNIJK(I,J,KP)

            CALL GET_CELL_NODE_COORDINATES(IJK,'U_MOMENTUM')

!======================================================================
!  Get Interpolation factors at East face
!======================================================================

            Theta_Ue(IJK) = DELX_Ue(IJK) / (DELX_Ue(IJK) + DELX_Uw(IPJK))
            Theta_Ue_bar(IJK) = ONE - Theta_Ue(IJK)

!======================================================================
!  Get Interpolation factors at North face
!======================================================================

            V_NODE_AT_NW = ((.NOT.BLOCKED_V_CELL_AT(IJK)).AND.(.NOT.WALL_V_AT(IJK)))
            V_NODE_AT_NE = ((.NOT.BLOCKED_V_CELL_AT(IPJK)).AND.(.NOT.WALL_V_AT(IPJK)))

            IF(V_NODE_AT_NW.AND.V_NODE_AT_NE) THEN
                Theta_U_ne(IJK) = (X_U_nc(IJK) - X_V(IJK) ) / (X_V(IPJK) - X_V(IJK))
                Theta_U_nw(IJK) = ONE - Theta_U_ne(IJK)
            ELSE IF (V_NODE_AT_NE.AND.(.NOT.V_NODE_AT_NW)) THEN
               IF(NO_K) THEN
                  Xi = Xn_U_int(IJK)
               ELSE
                  Xi = HALF * (Xn_U_int(IJK) + Xn_U_int(IJKM))
               ENDIF
               Theta_U_ne(IJK) = (X_U_nc(IJK) - Xi) / (X_V(IPJK) - Xi)
               Theta_U_nw(IJK) = ONE - Theta_U_ne(IJK)
            ELSE IF ((.NOT.V_NODE_AT_NE).AND.V_NODE_AT_NW) THEN
               IF(NO_K) THEN
                  Xi = Xn_U_int(IJK)
               ELSE
                  Xi = HALF * (Xn_U_int(IJK) + Xn_U_int(IJKM))
               ENDIF
               Theta_U_ne(IJK) = (X_U_nc(IJK) - X_V(IJK) ) / (Xi - X_V(IJK))
               Theta_U_nw(IJK) = ONE - Theta_U_ne(IJK)
            ELSE
               Theta_U_ne(IJK) = ZERO
               Theta_U_nw(IJK) = ZERO
            ENDIF


            IF ((Theta_U_ne(IJK)>=ONE).OR.(Theta_U_ne(IJK)<=ZERO)) THEN
               Theta_U_ne(IJK) = HALF
               Theta_U_nw(IJK) = HALF
            ENDIF

!======================================================================
!  Get Interpolation factors at Top face
!======================================================================

            IF(DO_K) THEN

               W_NODE_AT_TW = ((.NOT.BLOCKED_W_CELL_AT(IJK)).AND.(.NOT.WALL_W_AT(IJK)))
               W_NODE_AT_TE = ((.NOT.BLOCKED_W_CELL_AT(IPJK)).AND.(.NOT.WALL_W_AT(IPJK)))

               IF(W_NODE_AT_TW.AND.W_NODE_AT_TE) THEN
                  Theta_U_te(IJK) = (X_U_tc(IJK) - X_W(IJK) ) / (X_W(IPJK) - X_W(IJK))
                  Theta_U_tw(IJK) = ONE - Theta_U_te(IJK)

               ELSE IF (W_NODE_AT_TE.AND.(.NOT.W_NODE_AT_TW)) THEN
                  Xi = HALF * (Xn_U_int(IJK) + Xn_U_int(IJMK))
                  Theta_U_te(IJK) = (X_U_tc(IJK) - Xi) / (X_W(IPJK) - Xi)
                  Theta_U_tw(IJK) = ONE - Theta_U_te(IJK)
               ELSE IF ((.NOT.W_NODE_AT_TE).AND.W_NODE_AT_TW) THEN
                  Xi = HALF * (Xn_U_int(IJK) + Xn_U_int(IJMK))
                  Theta_U_te(IJK) = (X_U_tc(IJK) - X_W(IJK) ) / (Xi - X_W(IJK))
                  Theta_U_tw(IJK) = ONE - Theta_U_te(IJK)
               ELSE
                  Theta_U_te(IJK) = ZERO
                  Theta_U_tw(IJK) = ZERO
               ENDIF

               IF ((Theta_U_te(IJK)>=ONE).OR.(Theta_U_te(IJK)<=ZERO)) THEN
                  Theta_U_te(IJK) = HALF
                  Theta_U_tw(IJK) = HALF
               ENDIF
            ENDIF

            BCV = BC_U_ID(IJK)

            IF(BCV>0) THEN
               BCT = BC_TYPE_ENUM(BCV)
            ELSE
               BCT = BLANK
            ENDIF

            IF(BCT ==CG_NSW.OR.BCT ==CG_PSW) THEN

               CALL GET_DEL_H(IJK,'U_MOMENTUM',X_U(IJK),Y_U(IJK),Z_U(IJK),DELH_U(IJK),Nx,Ny,Nz)

               IF(DELH_U(IJK)<ZERO.AND.(.NOT.WALL_U_AT(IJK))) THEN
                  WRITE(*,*) 'NEGATIVE DELH-U AT XYZ=', X_U(IJK),Y_U(IJK),Z_U(IJK)
                  WRITE(*,*) 'DELH_U=', DELH_U(IJK)
                  WRITE(*,*) 'AYZ=', AYZ(IJK)
                  WRITE(*,*) 'MFIX WILL EXIT NOW.'
                  CALL MFIX_EXIT(MYPE)
               ENDIF

!======================================================================
!  Get Alpha Correction factors at East face
!======================================================================

! Location of the interpolated velocity along east face
! Xe has not moved based on definition of Theta_Ue

               Xe = X_NODE(8)
               Ye = Theta_Ue_bar(IJK) * Y_U(IJK) + Theta_Ue(IJK) * Y_U(IPJK)
               Ze = Theta_Ue_bar(IJK) * Z_U(IJK) + Theta_Ue(IJK) * Z_U(IPJK)

               CALL GET_DEL_H(IJK,'U_MOMENTUM',X_U_ec(IJK),Y_U_ec(IJK),Z_U_ec(IJK),DELH_ec,Nx,Ny,Nz)
               CALL GET_DEL_H(IJK,'U_MOMENTUM',Xe,Ye,Ze,DELH_e,Nx,Ny,Nz)

               IF((DELH_ec == UNDEFINED).OR.(DELH_e == UNDEFINED)) THEN
                  ALPHA_Ue_c(IJK) = ONE
               ELSE
                  ALPHA_Ue_c(IJK) = DMIN1(ALPHA_MAX , DELH_ec / DELH_e )
               ENDIF


               IF(BLOCKED_U_CELL_AT(IPJK)) ALPHA_Ue_c(IJK) = ZERO
               IF(WALL_U_AT(IJK).AND.WALL_U_AT(IPJK)) ALPHA_Ue_c(IJK) = ZERO
               IF(I == IEND1) ALPHA_Ue_c(IJK) = ONE

               IF(ALPHA_Ue_c(IJK)<ZERO) THEN
                 WRITE(*,*) 'NEGATIVE ALPHA_Ue_c at IJK=',IJK
                 WRITE(*,*) 'MFIX WILL EXIT NOW.'
                 CALL MFIX_EXIT(MYPE)
               ENDIF

!======================================================================
!  Get Non-Ortogonality correction factor at East face
!======================================================================

               Sx = X_U(IPJK) - X_U(IJK)
               Sy = Y_U(IPJK) - Y_U(IJK)
               Sz = Z_U(IPJK) - Z_U(IJK)

               NOC_U_E(IJK) = (Sy * Ny + Sz * Nz)/(Sx * DELH_e)

               IF(BLOCKED_U_CELL_AT(IPJK)) NOC_U_E(IJK) = ZERO
               IF(WALL_U_AT(IJK).AND.WALL_U_AT(IPJK)) NOC_U_E(IJK) = ZERO
               IF(I == IEND1) NOC_U_E(IJK) = ZERO

!======================================================================
!  Get Alpha Correction factors at North face
!======================================================================

               Theta_Un(IJK) = DELY_Un(IJK) / (DELY_Un(IJK) + DELY_Us(IJPK))
               Theta_Un_bar(IJK) = ONE - Theta_Un(IJK)

               IF ((Theta_Un(IJK)>=ONE).OR.(Theta_Un(IJK)<=ZERO)) THEN
                  Theta_Un(IJK) = HALF
                  Theta_Un_bar(IJK) = HALF
               ENDIF


! Location of the interpolated velocity along North face
! Yn has not moved based on definition of Theta_Un


               Xn = Theta_Un_bar(IJK) * X_U(IJK) + Theta_Un(IJK) * X_U(IJPK)
               Yn = Y_NODE(8)
               Zn = Theta_Un_bar(IJK) * Z_U(IJK) + Theta_Un(IJK) * Z_U(IJPK)

               CALL GET_DEL_H(IJK,'U_MOMENTUM',X_U_nc(IJK),Y_U_nc(IJK),Z_U_nc(IJK),DELH_nc,Nx,Ny,Nz)
               CALL GET_DEL_H(IJK,'U_MOMENTUM',Xn,Yn,Zn,DELH_n,Nx,Ny,Nz)

               IF((DELH_nc == UNDEFINED).OR.(DELH_n == UNDEFINED)) THEN
                  ALPHA_Un_c(IJK) = ONE
               ELSE
                  ALPHA_Un_c(IJK) = DMIN1(ALPHA_MAX , DELH_nc / DELH_n )
               ENDIF


               IF(BLOCKED_U_CELL_AT(IJPK)) ALPHA_Un_c(IJK) = ZERO
               IF(WALL_U_AT(IJK).AND.WALL_U_AT(IJPK)) ALPHA_Un_c(IJK) = ZERO
               IF(J == JEND1) ALPHA_Un_c(IJK) = ONE

               IF(ALPHA_Un_c(IJK)<ZERO) ALPHA_Un_c(IJK) = ONE

!======================================================================
!  Get Non-Ortogonality correction factor at North face
!======================================================================

               Sx = X_U(IJPK) - X_U(IJK)
               Sy = Y_U(IJPK) - Y_U(IJK)
               Sz = Z_U(IJPK) - Z_U(IJK)

               NOC_U_N(IJK) = (Sx * Nx + Sz * Nz)/(Sy * DELH_n)

               IF(BLOCKED_U_CELL_AT(IJPK)) NOC_U_N(IJK) = ZERO
               IF(WALL_U_AT(IJK).AND.WALL_U_AT(IJPK)) NOC_U_N(IJK) = ZERO
               IF(J == JEND1) NOC_U_N(IJK) = ZERO


               IF(DO_K) THEN
!======================================================================
!  Get Alpha Correction factors at Top face
!======================================================================

                  Theta_Ut(IJK) = DELZ_Ut(IJK) / (DELZ_Ut(IJK) + DELZ_Ub(IJKP))
                  Theta_Ut_bar(IJK) = ONE - Theta_Ut(IJK)

                  IF ((Theta_Ut(IJK)>=ONE).OR.(Theta_Ut(IJK)<=ZERO)) THEN
                     Theta_Ut(IJK) = HALF
                     Theta_Ut_bar(IJK) = HALF
                  ENDIF


! Location of the interpolated velocity along Top face
! Zt has not moved based on definition of Theta_Ut

                  Xt = Theta_Ut_bar(IJK) * X_U(IJK) + Theta_Ut(IJK) * X_U(IJKP)
                  Yt = Theta_Ut_bar(IJK) * Y_U(IJK) + Theta_Ut(IJK) * Y_U(IJKP)
                  Zt = Z_NODE(8)

                  CALL GET_DEL_H(IJK,'U_MOMENTUM',X_U_tc(IJK),Y_U_tc(IJK),Z_U_tc(IJK),DELH_tc,Nx,Ny,Nz)
                  CALL GET_DEL_H(IJK,'U_MOMENTUM',Xt,Yt,Zt,DELH_t,Nx,Ny,Nz)

                  IF((DELH_tc == UNDEFINED).OR.(DELH_t == UNDEFINED)) THEN
                     ALPHA_Ut_c(IJK) = ONE
                  ELSE
                     ALPHA_Ut_c(IJK) = DMIN1(ALPHA_MAX , DELH_tc / DELH_t )
                  ENDIF

                  IF(BLOCKED_U_CELL_AT(IJKP)) ALPHA_Ut_c(IJK) = ZERO
                  IF(WALL_U_AT(IJK).AND.WALL_U_AT(IJKP)) ALPHA_Ut_c(IJK) = ZERO
                  IF(K == KEND1) ALPHA_Ut_c(IJK) = ONE

                  IF(ALPHA_Ut_c(IJK)<ZERO) ALPHA_Ut_c(IJK) = ONE

!======================================================================
!  Get Non-Ortogonality correction factor at Top face
!======================================================================

                  Sx = X_U(IJKP) - X_U(IJK)
                  Sy = Y_U(IJKP) - Y_U(IJK)
                  Sz = Z_U(IJKP) - Z_U(IJK)

                  NOC_U_T(IJK) = (Sx * Nx + Sy * Ny)/(Sz * DELH_t)


                  IF(BLOCKED_U_CELL_AT(IJKP)) NOC_U_T(IJK) = ZERO
                  IF(WALL_U_AT(IJK).AND.WALL_U_AT(IJKP)) NOC_U_T(IJK) = ZERO
                  IF(K == KEND1) NOC_U_T(IJK) = ZERO

               ENDIF

            ENDIF



!======================================================================
!  Get Surfaces used to compute pressure gradient
!======================================================================

            IF (  CUT_U_CELL_AT(IJK)  )  THEN

               SELECT CASE (PG_OPTION)

                  CASE(0)                               ! do nothing
                                                        ! will be treated below

                  CASE(1)

                     A_UPG_E(IJK) = dmax1(AYZ_U(IJK),AYZ_U(IMJK))
                     A_UPG_W(IJK) = A_UPG_E(IJK)

                  CASE(2)

                     A_UPG_E(IJK) = AYZ_U(IJK)
                     A_UPG_W(IJK) = AYZ_U(IMJK)

                  CASE DEFAULT

                     WRITE(*,*)'INVALID PRESSURE GRADIENT OPTION:',PG_OPTION
                     WRITE(*,*)'PG_OPTION SHOULD BE SET EQUAL TO 0, 1 OR 2.'
                     WRITE(*,*)'PLEASE VERIFY MFIX.DAT FILE AND TRY AGAIN.'
                     WRITE(*,*) 'MFIX WILL EXIT NOW.'
                     CALL MFIX_EXIT(myPE)

               END SELECT


               IF (BLOCKED_CELL_AT(IJK).OR.BLOCKED_CELL_AT(IPJK)) THEN
                  IF(.NOT.WALL_U_AT(IJK)) THEN
                     WALL_U_AT(IJK) = .TRUE.
                     FLAG_E(IJK) = 0
                     IF(PRINT_WARNINGS) THEN
                        WRITE(*,*) 'WARNING: ONLY ONE PRESSURE NODE DETECTED IN U-MOMENTUM CELL IJK =',IJK
                        WRITE(*,*) 'RESETTING U-MOMENTUM CELL AS U_WALL CELL.'
                     ENDIF
!                     write(*,*) 'MFiX will exit now.'
!                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF

            ELSE   ! Regular cell

               A_UPG_E(IJK) = AYZ_U(IJK)
               A_UPG_W(IJK) = AYZ_U(IJK)

            ENDIF

         ENDIF
      END DO

      IF(PG_OPTION==0) THEN
         A_UPG_E = AYZ
         A_UPG_W = AYZ
      ENDIF

      RETURN

      END SUBROUTINE GET_3D_ALPHA_U_CUT_CELL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_3D_ALPHA_V_CUT_CELL                                C
!  Purpose: Calculate the correction term alpha_V                      C
!           for V-momentum cut cells,                                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_3D_ALPHA_V_CUT_CELL

      USE bc
      USE compar, ONLY: iend1, jend1, kend1, IJKSTART3, IJKEND3, mype, pe_io
      USE cutcell
      USE functions, ONLY: FUNIJK
      USE exit, only: mfix_exit
      USE geometry, ONLY: DO_K, NO_K, axz, axz_v, flag_n
      USE indices, ONLY: I_OF, J_OF, K_OF
      USE param1, ONLY: half, one, undefined, zero

      IMPLICIT NONE
      DOUBLE PRECISION:: Xe,Ye,Ze,Xn,Yn,Zn,Xt,Yt,Zt
      INTEGER :: I,J,K,IM,JM,KM,IP,JP,KP,IJK,IMJK,IJMK,IPJK,IJPK,IJKM,IJKP
      DOUBLE PRECISION :: Yi
      DOUBLE PRECISION :: Sx,Sy,Sz
      DOUBLE PRECISION :: Nx,Ny,Nz
      DOUBLE PRECISION :: DELH_ec , DELH_e
      DOUBLE PRECISION :: DELH_nc , DELH_n
      DOUBLE PRECISION :: DELH_tc , DELH_t
      LOGICAL :: U_NODE_AT_NE, U_NODE_AT_SE
      LOGICAL :: W_NODE_AT_NT, W_NODE_AT_ST
      INTEGER :: BCV
      INTEGER ::BCT

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'COMPUTING INTERPOLATION FACTORS IN V-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)

      DELH_V = UNDEFINED

      Theta_V_ne = HALF
      Theta_V_se  = HALF
      Theta_Vn  = HALF
      Theta_Vn_bar = HALF
      Theta_V_nt  = HALF
      Theta_V_st = HALF
      Theta_Ve  = HALF
      Theta_Ve_bar = HALF
      ALPHA_Ve_c  = ONE
      NOC_V_E  = ZERO
      ALPHA_Vn_c  = ONE
      NOC_V_N  = ZERO
      Theta_Vt  = HALF
      Theta_Vt_bar = HALF
      ALPHA_Vt_c  = ONE
      NOC_V_T  = ZERO

      A_VPG_N = AXZ_V
      A_VPG_S = AXZ_V


      DO IJK = IJKSTART3, IJKEND3
         IF(INTERIOR_CELL_AT(IJK)) THEN

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IM = I - 1
            JM = J - 1
            KM = K - 1

            IP = I + 1
            JP = J + 1
            KP = K + 1

            IMJK = FUNIJK(IM,J,K)
            IPJK = FUNIJK(IP,J,K)
            IJMK = FUNIJK(I,JM,K)
            IJPK = FUNIJK(I,JP,K)
            IJKM = FUNIJK(I,J,KM)
            IJKP = FUNIJK(I,J,KP)

            CALL GET_CELL_NODE_COORDINATES(IJK,'V_MOMENTUM')

!======================================================================
!  Get Interpolation factors at East face
!======================================================================

            U_NODE_AT_NE = ((.NOT.BLOCKED_U_CELL_AT(IJPK)).AND.(.NOT.WALL_U_AT(IJPK)))
            U_NODE_AT_SE = ((.NOT.BLOCKED_U_CELL_AT(IJK)).AND.(.NOT.WALL_U_AT(IJK)))

            IF(U_NODE_AT_SE.AND.U_NODE_AT_NE) THEN
               Theta_V_ne(IJK) = (Y_V_ec(IJK) - Y_U(IJK)    ) / (Y_U(IJPK) - Y_U(IJK))
               Theta_V_se(IJK) = ONE - Theta_V_ne(IJK)
            ELSE IF (U_NODE_AT_SE.AND.(.NOT.U_NODE_AT_NE)) THEN
               IF(NO_K) THEN
                  Yi = Ye_V_int(IJK)
               ELSE
                  Yi = HALF * (Ye_V_int(IJK) + Ye_V_int(IJKM))
               ENDIF
               Theta_V_ne(IJK) = (Y_V_ec(IJK) - Y_U(IJK)    ) / (Yi - Y_U(IJK))
               Theta_V_se(IJK) = ONE - Theta_V_ne(IJK)
            ELSE IF ((.NOT.U_NODE_AT_SE).AND.U_NODE_AT_NE) THEN
               IF(NO_K) THEN
                  Yi = Ye_V_int(IJK)
               ELSE
                  Yi = HALF * (Ye_V_int(IJK) + Ye_V_int(IJKM))
               ENDIF
               Theta_V_ne(IJK) = (Y_V_ec(IJK) - Yi    ) / (Y_U(IJPK) - Yi)
               Theta_V_se(IJK) = ONE - Theta_V_ne(IJK)
            ELSE
               Theta_V_ne(IJK) = ZERO
               Theta_V_se(IJK) = ZERO
            ENDIF

            IF ((Theta_V_ne(IJK)>=ONE).OR.(Theta_V_ne(IJK)<=ZERO)) THEN
               Theta_V_ne(IJK) = HALF
               Theta_V_se(IJK) = HALF
            ENDIF

!======================================================================
!  Get Interpolation factors at North face
!======================================================================

            Theta_Vn(IJK) = DELY_Vn(IJK) / (DELY_Vn(IJK) + DELY_Vs(IJPK))
            Theta_Vn_bar(IJK) = ONE - Theta_Vn(IJK)

            IF ((Theta_Vn(IJK)>=ONE).OR.(Theta_Vn(IJK)<=ZERO)) THEN
               Theta_Vn(IJK) = HALF
               Theta_Vn_bar(IJK) = HALF
            ENDIF

            IF(DO_K) THEN
!======================================================================
!  Get Interpolation factors at Top face
!======================================================================

               W_NODE_AT_NT = ((.NOT.BLOCKED_W_CELL_AT(IJPK)).AND.(.NOT.WALL_W_AT(IJPK)))
               W_NODE_AT_ST = ((.NOT.BLOCKED_W_CELL_AT(IJK)).AND.(.NOT.WALL_W_AT(IJK)))

               IF(W_NODE_AT_ST.AND.W_NODE_AT_NT) THEN
                  Theta_V_nt(IJK) = (Y_V_tc(IJK) - Y_W(IJK)    ) / (Y_W(IJPK) - Y_W(IJK))
                  Theta_V_st(IJK) = ONE - Theta_V_nt(IJK)
               ELSE IF (W_NODE_AT_ST.AND.(.NOT.W_NODE_AT_NT)) THEN
                  Yi = HALF * (Ye_V_int(IJK) + Ye_V_int(IMJK))
                  Theta_V_nt(IJK) = (Y_V_tc(IJK) - Y_W(IJK)    ) / (Yi - Y_W(IJK))
                  Theta_V_st(IJK) = ONE - Theta_V_nt(IJK)
               ELSE IF ((.NOT.W_NODE_AT_ST).AND.W_NODE_AT_NT) THEN
                  Yi = HALF * (Ye_V_int(IJK) + Ye_V_int(IMJK))
                  Theta_V_nt(IJK) = (Y_V_tc(IJK) - Yi ) / (Y_W(IJPK) - Yi)
                  Theta_V_st(IJK) = ONE - Theta_V_nt(IJK)
               ELSE
                  Theta_V_nt(IJK) = ZERO
                  Theta_V_st(IJK) = ZERO
               ENDIF

               IF ((Theta_V_nt(IJK)>=ONE).OR.(Theta_V_nt(IJK)<=ZERO)) THEN
                  Theta_V_nt(IJK) = HALF
                  Theta_V_st(IJK) = HALF
               ENDIF

            ENDIF

            BCV = BC_V_ID(IJK)
            IF(BCV>0) THEN
               BCT = BC_TYPE_ENUM(BCV)
            ELSE
               BCT = BLANK
            ENDIF

            IF(BCT ==CG_NSW.OR.BCT ==CG_PSW) THEN

               CALL GET_DEL_H(IJK,'V_MOMENTUM',X_V(IJK),Y_V(IJK),Z_V(IJK),DELH_V(IJK),Nx,Ny,Nz)

               IF(DELH_V(IJK)<ZERO.AND.(.NOT.WALL_V_AT(IJK))) THEN
                 WRITE(*,*) 'NEGATIVE DELH-V AT XYZ=', X_V(IJK),Y_V(IJK),Z_V(IJK)
                 WRITE(*,*) 'DELH_V=', DELH_V(IJK)
                 WRITE(*,*) 'AYZ=', AXZ(IJK)
                 WRITE(*,*) 'MFIX WILL EXIT NOW.'
                 CALL MFIX_EXIT(MYPE)
              ENDIF

!======================================================================
!  Get Alpha Correction factors at East face
!======================================================================

               Theta_Ve(IJK) = DELX_Ve(IJK) / (DELX_Ve(IJK) + DELX_Vw(IPJK))
               Theta_Ve_bar(IJK) = ONE - Theta_Ve(IJK)

            IF ((Theta_Ve(IJK)>=ONE).OR.(Theta_Ve(IJK)<=ZERO)) THEN
               Theta_Ve(IJK) = HALF
               Theta_Ve_bar(IJK) = HALF
            ENDIF


! Location of the interpolated velocity along East face
! Xe has not moved based on definition of Theta_Ve
               Xe = X_NODE(8)
               Ye = Theta_Ve_bar(IJK) * Y_V(IJK) + Theta_Ve(IJK) * Y_V(IPJK)
               Ze = Theta_Ve_bar(IJK) * Z_V(IJK) + Theta_Ve(IJK) * Z_V(IPJK)


               CALL GET_DEL_H(IJK,'V_MOMENTUM',X_V_ec(IJK),Y_V_ec(IJK),Z_V_ec(IJK),DELH_ec,Nx,Ny,Nz)

               CALL GET_DEL_H(IJK,'V_MOMENTUM',Xe,Ye,Ze,DELH_e,Nx,Ny,Nz)

               IF((DELH_ec == UNDEFINED).OR.(DELH_e == UNDEFINED)) THEN
                  ALPHA_Ve_c(IJK) = ONE
               ELSE
                  ALPHA_Ve_c(IJK) = DMIN1(ALPHA_MAX , DELH_ec / DELH_e )
               ENDIF


               IF(BLOCKED_V_CELL_AT(IPJK)) ALPHA_Ve_c(IJK) = ZERO
               IF(WALL_V_AT(IJK).AND.WALL_V_AT(IPJK)) ALPHA_Ve_c(IJK) = ZERO
               IF(I == IEND1) ALPHA_Ve_c(IJK) = ONE

               IF(ALPHA_Ve_c(IJK)<ZERO) ALPHA_Ve_c(IJK) = ONE

!======================================================================
!  Get Non-Ortogonality correction factor at East face
!======================================================================

               Sx = X_V(IPJK) - X_V(IJK)
               Sy = Y_V(IPJK) - Y_V(IJK)
               Sz = Z_V(IPJK) - Z_V(IJK)

               NOC_V_E(IJK) = (Sy * Ny + Sz * Nz)/(Sx * DELH_e)

               IF(BLOCKED_V_CELL_AT(IPJK)) NOC_V_E(IJK) = ZERO
               IF(WALL_V_AT(IJK).AND.WALL_V_AT(IPJK)) NOC_V_E(IJK) = ZERO
               IF(I == IEND1) NOC_V_E(IJK) = ZERO

!======================================================================
!  Get Alpha Correction factors at North face
!======================================================================

! Location of the interpolated velocity along north face
! Yn has not moved based on definition of Theta_V

               Xn = Theta_Vn_bar(IJK) * X_V(IJK) + Theta_Vn(IJK) * X_V(IJPK)
               Yn = Y_NODE(8)
               Zn = Theta_Vn_bar(IJK) * Z_V(IJK) + Theta_Vn(IJK) * Z_V(IJPK)

               CALL GET_DEL_H(IJK,'V_MOMENTUM',X_V_nc(IJK),Y_V_nc(IJK),Z_V_nc(IJK),DELH_nc,Nx,Ny,Nz)

               CALL GET_DEL_H(IJK,'V_MOMENTUM',Xn,Yn,Zn,DELH_n,Nx,Ny,Nz)

               IF((DELH_nc == UNDEFINED).OR.(DELH_n == UNDEFINED)) THEN
                  ALPHA_Vn_c(IJK) = ONE
               ELSE
                  ALPHA_Vn_c(IJK) = DMIN1(ALPHA_MAX , DELH_nc / DELH_n )
               ENDIF

               IF(BLOCKED_V_CELL_AT(IJPK)) ALPHA_Vn_c(IJK) = ZERO
               IF(WALL_V_AT(IJK).AND.WALL_V_AT(IJPK)) ALPHA_Vn_c(IJK) = ZERO
               IF(J == JEND1) ALPHA_Vn_c(IJK) = ONE

               IF(ALPHA_Vn_c(IJK)<ZERO) THEN
                 WRITE(*,*) 'NEGATIVE ALPHA_Vn_c at IJK=',IJK
                 WRITE(*,*) 'MFIX WILL EXIT NOW.'
                 CALL MFIX_EXIT(MYPE)
               ENDIF

!======================================================================
!  Get Non-Ortogonality correction factor at North face
!======================================================================

               Sx = X_V(IJPK) - X_V(IJK)
               Sy = Y_V(IJPK) - Y_V(IJK)
               Sz = Z_V(IJPK) - Z_V(IJK)

               NOC_V_N(IJK) = (Sx * Nx + Sz * Nz)/(Sy * DELH_n)

               IF(BLOCKED_V_CELL_AT(IJPK)) NOC_V_N(IJK) = ZERO
               IF(WALL_V_AT(IJK).AND.WALL_V_AT(IJPK)) NOC_V_N(IJK) = ZERO
               IF(J == JEND1) NOC_V_N(IJK) = ZERO

               IF(DO_K) THEN

!======================================================================
!  Get Alpha Correction factors at Top face
!======================================================================

                  Theta_Vt(IJK) = DELZ_Vt(IJK) / (DELZ_Vt(IJK) + DELZ_Vb(IJKP))
                  Theta_Vt_bar(IJK) = ONE - Theta_Vt(IJK)

                  IF ((Theta_Vt(IJK)>=ONE).OR.(Theta_Vt(IJK)<=ZERO)) THEN
                     Theta_Vt(IJK) = HALF
                     Theta_Vt_bar(IJK) = HALF
                  ENDIF

! Location of the interpolated velocity along Top face
! Zt has not moved based on definition of Theta_Ut

                  Xt = Theta_Vt_bar(IJK) * X_V(IJK) + Theta_Vt(IJK) * X_V(IJKP)
                  Yt = Theta_Vt_bar(IJK) * Y_V(IJK) + Theta_Vt(IJK) * Y_V(IJKP)
                  Zt = Z_NODE(8)

                  CALL GET_DEL_H(IJK,'V_MOMENTUM',X_V_tc(IJK),Y_V_tc(IJK),Z_V_tc(IJK),DELH_tc,Nx,Ny,Nz)

                  CALL GET_DEL_H(IJK,'V_MOMENTUM',Xt,Yt,Zt,DELH_t,Nx,Ny,Nz)

                  IF((DELH_tc == UNDEFINED).OR.(DELH_t == UNDEFINED)) THEN
                     ALPHA_Vt_c(IJK) = ONE
                  ELSE
                     ALPHA_Vt_c(IJK) = DMIN1(ALPHA_MAX , DELH_tc / DELH_t )
                  ENDIF

                  IF(BLOCKED_V_CELL_AT(IJKP)) ALPHA_Vt_c(IJK) = ZERO
                  IF(WALL_V_AT(IJK).AND.WALL_V_AT(IJKP)) ALPHA_Vt_c(IJK) = ZERO
                  IF(K == KEND1) ALPHA_Vt_c(IJK) = ONE

                  IF(ALPHA_Vt_c(IJK)<ZERO) ALPHA_Vt_c(IJK) = ONE

!======================================================================
!  Get Non-Ortogonality correction factor at Top face
!======================================================================

                  Sx = X_V(IJKP) - X_V(IJK)
                  Sy = Y_V(IJKP) - Y_V(IJK)
                  Sz = Z_V(IJKP) - Z_V(IJK)

                  NOC_V_T(IJK) = (Sx * Nx + Sy * Ny)/(Sz * DELH_t)


                  IF(BLOCKED_V_CELL_AT(IJKP)) NOC_V_T(IJK) = ZERO
                  IF(WALL_V_AT(IJK).AND.WALL_V_AT(IJKP)) NOC_V_T(IJK) = ZERO
                  IF(K == KEND1) NOC_V_T(IJK) = ZERO

               ENDIF

            ENDIF

!======================================================================
!  Get Surfaces used to compute pressure gradient
!======================================================================

            IF (  CUT_V_CELL_AT(IJK)  )  THEN

               SELECT CASE (PG_OPTION)

                  CASE(0)                               ! do nothing
                                                     ! will be treated below

                  CASE(1)

                     A_VPG_N(IJK) = dmax1(AXZ_V(IJK),AXZ_V(IJMK))
                     A_VPG_S(IJK) = A_VPG_N(IJK)

                  CASE(2)

                     A_VPG_N(IJK) = AXZ_V(IJK)
                     A_VPG_S(IJK) = AXZ_V(IJMK)


                  CASE DEFAULT

                     WRITE(*,*)'INVALID PRESSURE GRADIENT OPTION:',PG_OPTION
                     WRITE(*,*)'PG_OPTION SHOULD BE SET EQUAL TO 0, 1 OR 2.'
                     WRITE(*,*)'PLEASE VERIFY MFIX.DAT FILE AND TRY AGAIN.'
                     WRITE(*,*) 'MFIX WILL EXIT NOW.'
                     CALL MFIX_EXIT(myPE)

               END SELECT


               IF (BLOCKED_CELL_AT(IJK).OR.BLOCKED_CELL_AT(IJPK)) THEN
                  IF(.NOT.WALL_V_AT(IJK)) THEN
                     WALL_V_AT(IJK) = .TRUE.
                     FLAG_N(IJK) = 0
                     IF(PRINT_WARNINGS) THEN
                        WRITE(*,*) 'WARNING: ONLY ONE PRESSURE NODE DETECTED IN V-MOMENTUM CELL IJK =',IJK
                        WRITE(*,*) 'RESETTING U-MOMENTUM CELL AS V_WALL CELL.'
                     ENDIF
!                     write(*,*) 'MFiX will exit now.'
!                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF

            ELSE

               A_VPG_N(IJK) = AXZ_V(IJK)
               A_VPG_S(IJK) = AXZ_V(IJK)

            ENDIF

         ENDIF
      END DO

      IF(PG_OPTION==0) THEN
         A_VPG_N = AXZ
         A_VPG_S = AXZ
      ENDIF

      RETURN

      END SUBROUTINE GET_3D_ALPHA_V_CUT_CELL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_3D_ALPHA_W_CUT_CELL                                C
!  Purpose: Calculate the correction term alpha_W                      C
!           for W-momentum cut cells,                                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_3D_ALPHA_W_CUT_CELL

      USE bc
      USE compar, ONLY: iend1, jend1, kend1, IJKSTART3, IJKEND3, mype, pe_io
      USE cutcell
      USE exit, only: mfix_exit
      USE functions, ONLY: FUNIJK
      USE geometry, ONLY: axy, axy_w, flag_t
      USE indices, ONLY: I_OF, J_OF, K_OF
      USE param1, ONLY: half, one, undefined, zero

      IMPLICIT NONE
      DOUBLE PRECISION:: Xe,Ye,Ze,Xn,Yn,Zn,Xt,Yt,Zt
      INTEGER :: I,J,K,IM,JM,KM,IP,JP,KP,IJK,IMJK,IJMK,IPJK,IJPK,IJKM,IJKP
      DOUBLE PRECISION :: Zi
      DOUBLE PRECISION :: Sx,Sy,Sz
      DOUBLE PRECISION :: Nx,Ny,Nz
      DOUBLE PRECISION :: DELH_ec , DELH_e
      DOUBLE PRECISION :: DELH_nc , DELH_n
      DOUBLE PRECISION :: DELH_tc , DELH_t
      LOGICAL :: U_NODE_AT_TE, U_NODE_AT_BE
      LOGICAL :: V_NODE_AT_TN, V_NODE_AT_BN
      INTEGER :: BCV
      INTEGER ::BCT

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'COMPUTING INTERPOLATION FACTORS IN W-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)

      DELH_W = UNDEFINED

      Theta_W_te = HALF
      Theta_W_be = HALF
      Theta_W_tn = HALF
      Theta_W_bn = HALF
      Theta_Wt = HALF
      Theta_Wt_bar = HALF
      Theta_We = HALF
      Theta_We_bar = HALF
      ALPHA_We_c = ONE
      NOC_W_E = ZERO
      Theta_Wn = HALF
      Theta_Wn_bar = HALF
      ALPHA_Wn_c = ONE
      NOC_W_N = ZERO
      ALPHA_Wt_c = ONE
      NOC_W_T = ZERO
      A_WPG_T = AXY_W
      A_WPG_B = AXY_W

      DO IJK = IJKSTART3, IJKEND3
         IF(INTERIOR_CELL_AT(IJK)) THEN

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IM = I - 1
            JM = J - 1
            KM = K - 1

            IP = I + 1
            JP = J + 1
            KP = K + 1

            IMJK = FUNIJK(IM,J,K)
            IPJK = FUNIJK(IP,J,K)
            IJMK = FUNIJK(I,JM,K)
            IJPK = FUNIJK(I,JP,K)
            IJKM = FUNIJK(I,J,KM)
            IJKP = FUNIJK(I,J,KP)

            CALL GET_CELL_NODE_COORDINATES(IJK,'W_MOMENTUM')

!======================================================================
!  Get Interpolation factors at East face
!======================================================================

            U_NODE_AT_TE = ((.NOT.BLOCKED_U_CELL_AT(IJKP)).AND.(.NOT.WALL_U_AT(IJKP)))
            U_NODE_AT_BE = ((.NOT.BLOCKED_U_CELL_AT(IJK)).AND.(.NOT.WALL_U_AT(IJK)))


            IF(U_NODE_AT_TE.AND.U_NODE_AT_BE) THEN
               Theta_W_te(IJK) = (Z_W_ec(IJK) - Z_U(IJK)    ) / (Z_U(IJKP) - Z_U(IJK))
               Theta_W_be(IJK) = ONE - Theta_W_te(IJK)
            ELSE IF (U_NODE_AT_BE.AND.(.NOT.U_NODE_AT_TE)) THEN
               Zi = HALF * (Zt_W_int(IJK) + Zt_W_int(IJMK))
               Theta_W_te(IJK) = (Z_W_ec(IJK) - Z_U(IJK)    ) / (Zi - Z_U(IJK))
               Theta_W_be(IJK) = ONE - Theta_W_te(IJK)
            ELSE IF ((.NOT.U_NODE_AT_BE).AND.U_NODE_AT_TE) THEN
               Zi = HALF * (Zt_W_int(IJK) + Zt_W_int(IJMK))
               Theta_W_te(IJK) = (Z_W_ec(IJK) - Zi) / (Z_U(IJKP) - Zi)
               Theta_W_be(IJK) = ONE - Theta_W_te(IJK)
            ELSE
               Theta_W_te(IJK) = ZERO
               Theta_W_be(IJK) = ZERO
            ENDIF

            IF ((Theta_W_te(IJK)>=ONE).OR.(Theta_W_te(IJK)<=ZERO)) THEN
               Theta_W_te(IJK) = HALF
               Theta_W_be(IJK) = HALF
            ENDIF

!======================================================================
!  Get Interpolation factors at North face
!======================================================================

            V_NODE_AT_TN = ((.NOT.BLOCKED_V_CELL_AT(IJKP)).AND.(.NOT.WALL_V_AT(IJKP)))
            V_NODE_AT_BN = ((.NOT.BLOCKED_V_CELL_AT(IJK)).AND.(.NOT.WALL_V_AT(IJK)))


            IF(V_NODE_AT_TN.AND.V_NODE_AT_BN) THEN
               Theta_W_tn(IJK) = (Z_W_nc(IJK) - Z_V(IJK)    ) / (Z_V(IJKP) - Z_V(IJK))
               Theta_W_bn(IJK) = ONE - Theta_W_tn(IJK)
            ELSE IF (V_NODE_AT_BN.AND.(.NOT.V_NODE_AT_TN)) THEN
               Zi = HALF * (Zt_W_int(IJK) + Zt_W_int(IMJK))
               Theta_W_tn(IJK) = (Z_W_nc(IJK) - Z_V(IJK)    ) / (Zi - Z_V(IJK))
               Theta_W_bn(IJK) = ONE - Theta_W_bn(IJK)
            ELSE IF ((.NOT.V_NODE_AT_BN).AND.V_NODE_AT_TN) THEN
               Zi = HALF * (Zt_W_int(IJK) + Zt_W_int(IMJK))
               Theta_W_tn(IJK) = (Z_W_nc(IJK) - Zi) / (Z_V(IJKP) - Zi)
               Theta_W_bn(IJK) = ONE - Theta_W_bn(IJK)
            ELSE
               Theta_W_tn(IJK) = ZERO
               Theta_W_bn(IJK) = ZERO
            ENDIF

            IF ((Theta_W_tn(IJK)>=ONE).OR.(Theta_W_tn(IJK)<=ZERO)) THEN
               Theta_W_tn(IJK) = HALF
               Theta_W_bn(IJK) = HALF
            ENDIF

!======================================================================
!  Get Interpolation factors at Top face
!======================================================================

            Theta_Wt(IJK) = DELZ_Wt(IJK) / (DELZ_Wt(IJK) + DELZ_Wb(IJKP))
            Theta_Wt_bar(IJK) = ONE - Theta_Wt(IJK)

            IF ((Theta_Wt(IJK)>=ONE).OR.(Theta_Wt(IJK)<=ZERO)) THEN
               Theta_Wt(IJK) = HALF
               Theta_Wt_bar(IJK) = HALF
            ENDIF

            BCV = BC_W_ID(IJK)
            IF(BCV>0) THEN
               BCT = BC_TYPE_ENUM(BCV)
            ELSE
               BCT = BLANK
            ENDIF

            IF(BCT ==CG_NSW.OR.BCT ==CG_PSW) THEN

               CALL GET_DEL_H(IJK,'W_MOMENTUM',X_W(IJK),Y_W(IJK),Z_W(IJK),DELH_W(IJK),Nx,Ny,Nz)

               IF(DELH_W(IJK)<ZERO.AND.(.NOT.WALL_W_AT(IJK))) THEN
                  WRITE(*,*) 'NEGATIVE DELH-W AT XYZ=', X_W(IJK),Y_W(IJK),Z_W(IJK)
                  WRITE(*,*) 'DELH_W=', DELH_W(IJK)
                  WRITE(*,*) 'AXY=', AXY(IJK)
                  WRITE(*,*) 'MFIX WILL EXIT NOW.'
                  CALL MFIX_EXIT(MYPE)
               ENDIF

!======================================================================
!  Get Alpha Correction factors at East face
!======================================================================

               Theta_We(IJK) = DELX_We(IJK) / (DELX_We(IJK) + DELX_Ww(IPJK))
               Theta_We_bar(IJK) = ONE - Theta_We(IJK)

               IF ((Theta_We(IJK)>=ONE).OR.(Theta_We(IJK)<=ZERO)) THEN
                  Theta_We(IJK) = HALF
                  Theta_We_bar(IJK) = HALF
               ENDIF

! Location of the interpolated velocity along East face
! Xe has not moved based on definition of Theta_Ve
               Xe = X_NODE(8)
               Ye = Theta_We_bar(IJK) * Y_W(IJK) +Theta_We(IJK) * Y_W(IPJK)
               Ze = Theta_We_bar(IJK) * Z_W(IJK) +Theta_We(IJK) * Z_W(IPJK)

               CALL GET_DEL_H(IJK,'W_MOMENTUM',X_W_ec(IJK),Y_W_ec(IJK),Z_W_ec(IJK),DELH_ec,Nx,Ny,Nz)

               CALL GET_DEL_H(IJK,'W_MOMENTUM',Xe,Ye,Ze,DELH_e,Nx,Ny,Nz)


               IF((DELH_ec == UNDEFINED).OR.(DELH_e == UNDEFINED)) THEN
                  ALPHA_We_c(IJK) = ONE
               ELSE
                  ALPHA_We_c(IJK) = DMIN1(ALPHA_MAX , DELH_ec / DELH_e )
               ENDIF

               IF(BLOCKED_W_CELL_AT(IPJK)) ALPHA_We_c(IJK) = ZERO
               IF(WALL_W_AT(IJK).AND.WALL_W_AT(IPJK)) ALPHA_We_c(IJK) = ZERO
               IF(I == IEND1) ALPHA_We_c(IJK) = ONE

               IF(ALPHA_We_c(IJK)<ZERO) ALPHA_We_c(IJK) = ONE

!======================================================================
!  Get Non-Ortogonality correction factor at East face
!======================================================================

               Sx = X_W(IPJK) - X_W(IJK)
               Sy = Y_W(IPJK) - Y_W(IJK)
               Sz = Z_W(IPJK) - Z_W(IJK)

               NOC_W_E(IJK) = (Sy * Ny + Sz * Nz)/(Sx * DELH_e)


               IF(BLOCKED_W_CELL_AT(IPJK)) NOC_W_E(IJK) = ZERO
               IF(WALL_W_AT(IJK).AND.WALL_W_AT(IPJK)) NOC_W_E(IJK) = ZERO
               IF(I == IEND1) NOC_W_E(IJK) = ZERO

!======================================================================
!  Get Alpha Correction factors at North face
!======================================================================

               Theta_Wn(IJK) = DELY_Wn(IJK) / (DELY_Wn(IJK) + DELY_Ws(IJPK))
               Theta_Wn_bar(IJK) = ONE - Theta_Wn(IJK)

               IF ((Theta_Wn(IJK)>=ONE).OR.(Theta_Wn(IJK)<=ZERO)) THEN
                  Theta_Wn(IJK) = HALF
                  Theta_Wn_bar(IJK) = HALF
               ENDIF


! Location of the interpolated velocity along north face
! Yn has not moved based on definition of Theta_V

               Xn = Theta_Wn_bar(IJK) * X_W(IJK) + Theta_Wn(IJK) * X_W(IJPK)
               Yn = Y_NODE(8)
               Zn = Theta_Wn_bar(IJK) * Z_W(IJK) + Theta_Wn(IJK) * Z_W(IJPK)

               CALL GET_DEL_H(IJK,'W_MOMENTUM',X_W_nc(IJK),Y_W_nc(IJK),Z_W_nc(IJK),DELH_nc,Nx,Ny,Nz)

               CALL GET_DEL_H(IJK,'W_MOMENTUM',Xn,Yn,Zn,DELH_n,Nx,Ny,Nz)

               IF((DELH_nc == UNDEFINED).OR.(DELH_n == UNDEFINED)) THEN
                  ALPHA_Wn_c(IJK) = ONE
               ELSE
                  ALPHA_Wn_c(IJK) = DMIN1(ALPHA_MAX , DELH_nc / DELH_n )
               ENDIF

               IF(BLOCKED_W_CELL_AT(IJPK)) ALPHA_Wn_c(IJK) = ZERO
               IF(WALL_W_AT(IJK).AND.WALL_W_AT(IJPK)) ALPHA_Wn_c(IJK) = ZERO
               IF(J == JEND1) ALPHA_Wn_c(IJK) = ONE

               IF(ALPHA_Wn_c(IJK)<ZERO) ALPHA_Wn_c(IJK) = ONE

!======================================================================
!  Get Non-Ortogonality correction factor at North face
!======================================================================

               Sx = X_W(IJPK) - X_W(IJK)
               Sy = Y_W(IJPK) - Y_W(IJK)
               Sz = Z_W(IJPK) - Z_W(IJK)

               NOC_W_N(IJK) = (Sx * Nx + Sz * Nz)/(Sy * DELH_n)

               IF(BLOCKED_W_CELL_AT(IJPK)) NOC_W_N(IJK) = ZERO
               IF(WALL_W_AT(IJK).AND.WALL_W_AT(IJPK)) NOC_W_N(IJK) = ZERO
               IF(J == JEND1) NOC_W_N(IJK) = ZERO

!======================================================================
!  Get Alpha Correction factors at Top face
!======================================================================

! Location of the interpolated velocity along Top face
! Zt has not moved based on definition of Theta_Ut

               Xt = Theta_Wt_bar(IJK) * X_W(IJK) + Theta_Wt(IJK) * X_W(IJKP)
               Yt = Theta_Wt_bar(IJK) * Y_W(IJK) + Theta_Wt(IJK) * Y_W(IJKP)
               Zt = Z_NODE(8)

               CALL GET_DEL_H(IJK,'W_MOMENTUM',X_W_tc(IJK),Y_W_tc(IJK),Z_W_tc(IJK),DELH_tc,Nx,Ny,Nz)

               CALL GET_DEL_H(IJK,'W_MOMENTUM',Xt,Yt,Zt,DELH_t,Nx,Ny,Nz)

               IF((DELH_tc == UNDEFINED).OR.(DELH_t == UNDEFINED)) THEN
                  ALPHA_Wt_c(IJK) = ONE
               ELSE
                  ALPHA_Wt_c(IJK) = DMIN1(ALPHA_MAX , DELH_tc / DELH_t )
               ENDIF

               IF(BLOCKED_W_CELL_AT(IJKP)) ALPHA_Wt_c(IJK) = ZERO
               IF(WALL_W_AT(IJK).AND.WALL_W_AT(IJKP)) ALPHA_Wt_c(IJK) = ZERO
               IF(K == KEND1) ALPHA_Wt_c(IJK) = ONE

               IF(ALPHA_Wt_c(IJK)<ZERO) THEN
                  WRITE(*,*) 'NEGATIVE ALPHA_Wt_c at IJK=',IJK
                  WRITE(*,*) 'MFIX WILL EXIT NOW.'
                  CALL MFIX_EXIT(MYPE)
               ENDIF

!======================================================================
!  Get Non-Ortogonality correction factor at Top face
!======================================================================

               Sx = X_W(IJKP) - X_W(IJK)
               Sy = Y_W(IJKP) - Y_W(IJK)
               Sz = Z_W(IJKP) - Z_W(IJK)

               NOC_W_T(IJK) = (Sx * Nx + Sy * Ny)/(Sz * DELH_t)

               IF(BLOCKED_W_CELL_AT(IJKP)) NOC_W_T(IJK) = ZERO
               IF(WALL_W_AT(IJK).AND.WALL_W_AT(IJKP)) NOC_W_T(IJK) = ZERO
               IF(K == KEND1) NOC_W_T(IJK) = ZERO

            ENDIF

!======================================================================
!  Get Surfaces used to compute pressure gradient
!======================================================================

            IF (  CUT_W_CELL_AT(IJK)  )  THEN

               SELECT CASE (PG_OPTION)

                  CASE(0)                               ! do nothing
                                                     ! will be treated below

                  CASE(1)

                     A_WPG_T(IJK) = dmax1(AXY_W(IJK),AXY_W(IJKM))
                     A_WPG_B(IJK) = A_WPG_T(IJK)

                  CASE(2)

                     A_WPG_T(IJK) = AXY_W(IJK)
                     A_WPG_B(IJK) = AXY_W(IJKM)

                  CASE DEFAULT

                     WRITE(*,*)'INVALID PRESSURE GRADIENT OPTION:',PG_OPTION
                     WRITE(*,*)'PG_OPTION SHOULD BE SET EQUAL TO 0, 1 OR 2.'
                     WRITE(*,*)'PLEASE VERIFY MFIX.DAT FILE AND TRY AGAIN.'
                     WRITE(*,*) 'MFIX WILL EXIT NOW.'
                     CALL MFIX_EXIT(myPE)

               END SELECT

               IF (BLOCKED_CELL_AT(IJK).OR.BLOCKED_CELL_AT(IJKP)) THEN
                  IF(.NOT.WALL_W_AT(IJK)) THEN
                     WALL_W_AT(IJK) = .TRUE.
                     FLAG_T(IJK) = 0
                     IF(PRINT_WARNINGS) THEN
                        WRITE(*,*) 'WARNING: ONLY ONE PRESSURE NODE DETECTED IN W-MOMENTUM CELL IJK =',IJK
                        WRITE(*,*) 'RESETTING U-MOMENTUM CELL AS W_WALL CELL.'
                     ENDIF
!                     write(*,*) 'MFiX will exit now.'
!                     CALL MFIX_EXIT(myPE)
                  ENDIF
               ENDIF

            ELSE

               A_WPG_T(IJK) = AXY_W(IJK)
               A_WPG_B(IJK) = AXY_W(IJK)

            ENDIF

         ENDIF
      END DO

      IF(PG_OPTION==0) THEN
         A_WPG_T = AXY
         A_WPG_B = AXY
      ENDIF

      RETURN

      END SUBROUTINE GET_3D_ALPHA_W_CUT_CELL
