!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_W_g(A_m, B_m, IER)                           C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for W_g momentum eq.                                                C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CG_SOURCE_W_G(A_M, B_M)

! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE fldvar
      USE fun_avg
      use geometry, only: vol, vol_w
      use geometry, only: axy, axz, ayz
      use geometry, only: axy_w, axz_w, ayz_w
      USE indices
      USE param
      USE param1
      USE physprop
      USE run
      USE toleranc
      USE visc_g
      USE cutcell
      USE quadric
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION :: B_m(DIMENSION_3, 0:DIMENSION_M)

! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IJKT, IMJK, IJKP, IMJKP
      INTEGER :: IJKE, IJKW, IJKTE, IM, IPJK
! Phase index
      INTEGER :: M
! Average volume fraction
      DOUBLE PRECISION :: EPGA
      INTEGER :: JM,IP,JP,IJMK,IJPK
      INTEGER :: IJKC,IJKN,IJKNE,IJKS,IJKSE,IPJMK,IJKM,KM,KP,IJMKP
      INTEGER :: IJKTN,IJKWT,IJKST
      DOUBLE PRECISION :: We,Ww,Wn,Ws,Wt,Wb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_GT_E,MU_GT_W,MU_GT_N,MU_GT_S
      DOUBLE PRECISION :: MU_GT_T,MU_GT_B,MU_GT_CUT
      DOUBLE PRECISION :: WW_g
      INTEGER :: BCV
      INTEGER :: BCT
! virtual (added) mass
      DOUBLE PRECISION :: F_vir, ROP_MA
      DOUBLE PRECISION :: U_se, Usw, Wse, Wsw, Wsn, Wss
      DOUBLE PRECISION :: Wst, Wsb, Usc,Vsc,Vsn,Vss
! Wall function
      DOUBLE PRECISION :: W_F_Slip
!-----------------------------------------------
      IF(CG_SAFE_MODE(5)==1) RETURN

      M = 0
      IF (.NOT.MOMENTUM_Z_EQ(0)) RETURN
!!!$omp  parallel do private( I, J, K, IJK, IJKT, ISV, Sdp, V0, Vpm, Vmt, Vbf, &
!!!$omp&  PGT, ROGA, IMJK, IJKP, IMJKP, IJKW, IJKTE, IJKTW, IM, IPJK,  &
!!!$omp&  CTE, CTW, SXZB, EPMUOX, VXZA, VXZB, UGT, VCOA, VCOB, IJKE,&
!!!$omp&  MUGA, ROPGA, EPGA, LINE)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IJKT = TOP_OF(IJK)
         EPGA = AVG_Z(EP_G(IJK),EP_G(IJKT),K)
         IF (IP_AT_T(IJK)) THEN
! do nothing

! dilute flow
         ELSE IF (EPGA <= DIL_EP_S) THEN
! do nothing

         ELSEIF(INTERIOR_CELL_AT(IJK)) THEN
            BCV = BC_W_ID(IJK)

            IF(BCV > 0 ) THEN
               BCT = BC_TYPE_ENUM(BCV)
            ELSE
               BCT = NONE
            ENDIF

            SELECT CASE (BCT)
               CASE (CG_NSW)
                  NOC_WG = .TRUE.
                  WW_g = ZERO
                  MU_GT_CUT =  (VOL(IJK)*MU_GT(IJK) + VOL(IJKT)*MU_GT(IJKT))/(VOL(IJK) + VOL(IJKT))

                  IF(.NOT.K_EPSILON) THEN
                     A_M(IJK,0,M) = A_M(IJK,0,M) - MU_GT_CUT * Area_W_CUT(IJK)/DELH_W(IJK)
                  ELSE
                     CALL Wall_Function(IJK,IJK,ONE/DELH_W(IJK),W_F_Slip)
                     A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_GT_CUT * Area_W_CUT(IJK)*(ONE-W_F_Slip)/DELH_W(IJK)
                  ENDIF

               CASE (CG_FSW)
                  NOC_WG = .FALSE.
                  WW_g = ZERO
               CASE(CG_PSW)
                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     NOC_WG = .TRUE.
                     WW_g = BC_WW_G(BCV)
                     MU_GT_CUT = (VOL(IJK)*MU_GT(IJK) + VOL(IJKT)*MU_GT(IJKT))/(VOL(IJK) + VOL(IJKT))
                     A_M(IJK,0,M) = A_M(IJK,0,M) - MU_GT_CUT * Area_W_CUT(IJK)/DELH_W(IJK)
                     B_M(IJK,M) = B_M(IJK,M) - MU_GT_CUT * WW_g * Area_W_CUT(IJK)/DELH_W(IJK)
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     NOC_WG = .FALSE.
                     WW_g = ZERO
                  ELSE                              ! partial slip
                     NOC_WG = .FALSE.
                     WW_g = BC_WW_G(BCV)
                     MU_GT_CUT = (VOL(IJK)*MU_GT(IJK) + VOL(IJKT)*MU_GT(IJKT))/(VOL(IJK) + VOL(IJKT))
                     A_M(IJK,0,M) = A_M(IJK,0,M) - MU_GT_CUT * Area_W_CUT(IJK)*(BC_HW_G(BCV))
                     B_M(IJK,M) = B_M(IJK,M) - MU_GT_CUT * WW_g * Area_W_CUT(IJK)*(BC_HW_G(BCV))

                  ENDIF
               CASE (NONE,CG_MI)
                  NOC_WG = .FALSE.
            END SELECT

            IF(NOC_WG) THEN

               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IPJK = IP_OF(IJK)
               IJPK = JP_OF(IJK)
               IJKP = KP_OF(IJK)

               We = Theta_We_bar(IJK)  * W_g(IJK)  + Theta_We(IJK)  * W_g(IPJK)
               Ww = Theta_We_bar(IMJK) * W_g(IMJK) + Theta_We(IMJK) * W_g(IJK)

               Wn = Theta_Wn_bar(IJK)  * W_g(IJK)  + Theta_Wn(IJK)  * W_g(IJPK)
               Ws = Theta_Wn_bar(IJMK) * W_g(IJMK) + Theta_Wn(IJMK) * W_g(IJK)

               Wt = Theta_Wt_bar(IJK)  * W_g(IJK)  + Theta_Wt(IJK)  * W_g(IJKP)
               Wb = Theta_Wt_bar(IJKM) * W_g(IJKM) + Theta_Wt(IJKM) * W_g(IJK)

               IJKE = EAST_OF(IJK)

               ijkt = top_of(ijk)

               IF (WALL_AT(IJK)) THEN
                  IJKC = IJKT
               ELSE
                  IJKC = IJK
               ENDIF

               IP = IP1(I)
               IM = IM1(I)
               IJKN = NORTH_OF(IJK)
               IJKNE = EAST_OF(IJKN)

               JM = JM1(J)
               IPJMK = IP_OF(IJMK)
               IJKS = SOUTH_OF(IJK)
               IJKSE = EAST_OF(IJKS)

               KP = KP1(K)
               IJKT = TOP_OF(IJK)
               IJKE = EAST_OF(IJK)
               IJKP = KP_OF(IJK)
               IJKTN = NORTH_OF(IJKT)
               IJKTE = EAST_OF(IJKT)
               IJKW = WEST_OF(IJK)
               IJKWT = TOP_OF(IJKW)
               IJKS = SOUTH_OF(IJK)
               IJKST = TOP_OF(IJKS)

               MU_GT_E = AVG_Z_H(AVG_X_H(MU_GT(IJKC),MU_GT(IJKE),I),&
                         AVG_X_H(MU_GT(IJKT),MU_GT(IJKTE),I),K)

               MU_GT_W = AVG_Z_H(AVG_X_H(MU_GT(IJKW),MU_GT(IJKC),IM),&
                         AVG_X_H(MU_GT(IJKWT),MU_GT(IJKT),IM),K)

               MU_GT_N = AVG_Z_H(AVG_Y_H(MU_GT(IJKC),MU_GT(IJKN),J),&
                         AVG_Y_H(MU_GT(IJKT),MU_GT(IJKTN),J),K)

               MU_GT_S = AVG_Z_H(AVG_Y_H(MU_GT(IJKS),MU_GT(IJKC),JM),&
                         AVG_Y_H(MU_GT(IJKST),MU_GT(IJKT),JM),K)

               MU_GT_T = MU_GT(IJKT)
               MU_GT_B = MU_GT(IJKC)

               B_NOC =     MU_GT_E * Ayz_W(IJK)  * (We-WW_g) * NOC_W_E(IJK)  &
                       -   MU_GT_W * Ayz_W(IMJK) * (Ww-WW_g) * NOC_W_E(IMJK) &
                       +   MU_GT_N * Axz_W(IJK)  * (Wn-WW_g) * NOC_W_N(IJK)  &
                       -   MU_GT_S * Axz_W(IJMK) * (Ws-WW_g) * NOC_W_N(IJMK) &
                       +   MU_GT_T * Axy_W(IJK)  * (Wt-WW_g) * NOC_W_T(IJK)  *2.0d0&
                       -   MU_GT_B * Axy_W(IJKM) * (Wb-WW_g) * NOC_W_T(IJKM) *2.0D0

               B_M(IJK,M) = B_M(IJK,M)   +  B_NOC
            ENDIF

            IF(CUT_W_TREATMENT_AT(IJK)) THEN

! BEGIN VIRTUAL MASS SECTION (explicit terms)
! adding transient term  dWs/dt to virtual mass term
               F_vir = ZERO
               IF(Added_Mass) THEN

                  F_vir = ( (W_s(IJK,M_AM) - W_sO(IJK,M_AM)) )*ODT*VOL_W(IJK)
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
                  IJMK = FUNIJK(I,JM,K)
                  IPJK = FUNIJK(IP,J,K)
                  IJPK = FUNIJK(I,JP,K)
                  IJKP = FUNIJK(I,J,KP)
                  IJKM = FUNIJK(I,J,KM)
                  IMJKP = KP_OF(IMJK)
                  IJMKP = KP_OF(IJMK)
                  IJKE = EAST_OF(IJK)

! defining gas-particles velocity at momentum cell faces (or scalar cell center)
                  Wse = Theta_We_bar(IJK) * W_s(IJK,M_AM) + Theta_We(IJK) * W_s(IPJK,M_AM)
                  Wsw = Theta_We_bar(IMJK) * W_s(IMJK,M_AM) + Theta_We(IMJK) * W_s(IJK,M_AM)
                  U_se = Theta_W_te(IJK) * U_s(IJK,M_AM) + Theta_W_be(IJK) * U_s(IJKP,M_AM)
                  Usw = Theta_W_te(IMJK) * U_s(IMJK,M_AM) + Theta_W_be(IMJK) * U_s(IMJKP,M_AM)
                  Usc = (DELX_we(IJK) * Usw + DELX_ww(IJK) * U_se) / (DELX_we(IJK) + DELX_ww(IJK))
                  Wsn = Theta_Wn_bar(IJK) * W_s(IJK,M_AM) + Theta_Wn(IJK) * W_s(IJPK,M_AM)
                  Wss = Theta_Wn_bar(IJMK) * W_s(IJMK,M_AM) + Theta_Wn(IJMK) * W_s(IJK,M_AM)
                  Vsn =  Theta_W_tn(IJK)  * V_s(IJK,M_AM)  + Theta_W_bn(IJK)  * V_s(IJKP,M_AM)
                  Vss =  Theta_W_tn(IJMK) * V_s(IJMK,M_AM) + Theta_W_bn(IJMK) * V_s(IJMKP,M_AM)
                  Vsc = (DELY_wn(IJK) * Vss + DELY_ws(IJK) * Vsn) / (DELY_wn(IJK) + DELY_ws(IJK))
                  Wst = Theta_Wt_bar(IJK)  * W_s(IJK,M_AM)  + Theta_Wt(IJK)  * W_s(IJKP,M_AM)
                  Wsb = Theta_Wt_bar(IMJK) * W_s(IMJK,M_AM) + Theta_Wt(IMJK) * W_s(IMJKP,M_AM)

! adding convective terms (U dW/dx + V dW/dy + W dW/dz) to virtual mass
                  F_vir = F_vir +  Usc * (Wse - Wsw)*AYZ(IJK)    + &
                                   Vsc * (Wsn - Wss)*AXZ(IJK)  + &
                                   W_s(IJK,M_AM)*(Wst - Wsb)*AXY(IJK)


                  ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M_AM) + &
                            VOL(IJKT)*ROP_g(IJKT)*EP_s(IJKT,M_AM))/(VOL(IJK) + VOL(IJKT))

                  F_vir = F_vir * Cv * ROP_MA
                  B_M(IJK,M) = B_M(IJK,M) - F_vir ! explicit part of virtual mass force

               ENDIF
! END VIRTUAL MASS SECTION

            ENDIF
         ENDIF
      ENDDO

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CG_SOURCE_W_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CG_SOURCE_W_g_BC(A_m, B_m, IER)                        C
!  Purpose: Determine contribution of cut-cell to source terms         C
!  for W_g momentum eq.                                                C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CG_SOURCE_W_G_BC(A_M, B_M)

! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE fldvar
      USE fun_avg
      USE indices
      USE param
      USE param1
      USE physprop
      USE run
      USE toleranc
      USE visc_g
      USE cutcell
      USE quadric
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION :: B_m(DIMENSION_3, 0:DIMENSION_M)

! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKB
! Solids phase
      INTEGER :: M
      INTEGER :: BCV
      INTEGER :: BCT

!-----------------------------------------------

      IF(CG_SAFE_MODE(5)==1) RETURN
      M = 0

      DO IJK = ijkstart3, ijkend3

         BCV = BC_W_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE_ENUM(BCV)
         ELSE
            BCT = NONE
         ENDIF

         SELECT CASE (BCT)

            CASE (CG_NSW)

               IF(WALL_W_AT(IJK)) THEN

                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE

                  B_M(IJK,M) = ZERO

               ENDIF

            CASE (CG_FSW)

               IF(WALL_W_AT(IJK)) THEN

                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE

!                  B_M(IJK,M) = - W_g(W_MASTER_OF(IJK))  ! Velocity of master node

                  B_M(IJK,M) = ZERO

                  IF(DABS(NORMAL_W(IJK,3))/=ONE) THEN

                     IF (W_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                        A_M(IJK,east,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                        A_M(IJK,west,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                        A_M(IJK,north,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                        A_M(IJK,south,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                        A_M(IJK,top,M) = ONE
                     ELSEIF (W_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                        A_M(IJK,bottom,M) = ONE
                     ENDIF

                  ENDIF

               ENDIF

            CASE (CG_PSW)

               IF(WALL_W_AT(IJK)) THEN

                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE


                  IF(BC_HW_G(BCV)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = -BC_WW_G(BCV)
                  ELSEIF(BC_HW_G(BCV)==ZERO) THEN   ! same as FSW
                     B_M(IJK,M) = ZERO

                     IF(DABS(NORMAL_W(IJK,3))/=ONE) THEN

                        IF (W_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                           A_M(IJK,east,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                           A_M(IJK,west,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                           A_M(IJK,north,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                           A_M(IJK,south,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                           A_M(IJK,top,M) = ONE
                        ELSEIF (W_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                           A_M(IJK,bottom,M) = ONE
                        ENDIF

                     ENDIF

                  ELSE                              ! partial slip

                  ENDIF

               ENDIF


            CASE (CG_MI)

               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE

               IF(BC_W_g(BCV)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_W_g(BCV)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_g(BCV)*NORMAL_W(IJK,3)
               ENDIF


               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJKB,east,M) = ZERO
                  A_M(IJKB,west,M) = ZERO
                  A_M(IJKB,north,M) = ZERO
                  A_M(IJKB,south,M) = ZERO
                  A_M(IJKB,top,M) = ZERO
                  A_M(IJKB,bottom,M) = ZERO
                  A_M(IJKB,0,M) = -ONE

                  IF(BC_W_g(BCV)/=UNDEFINED) THEN
                     B_M(IJKB,M) = - BC_W_g(BCV)
                  ELSE
                     B_M(IJKB,M) = - BC_VELMAG_g(BCV)*NORMAL_W(IJK,3)
                  ENDIF


               ENDIF

            CASE (CG_PO)

               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO

               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJK,bottom,M) = ONE
                  A_M(IJK,0,M) = -ONE

               ENDIF

         END SELECT

         BCV = BC_ID(IJK)

         IF(BCV > 0 ) THEN
            BCT = BC_TYPE_ENUM(BCV)
         ELSE
            BCT = NONE
         ENDIF

         SELECT CASE (BCT)

            CASE (CG_MI)

               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE

               IF(BC_W_g(BCV)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_W_g(BCV)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_g(BCV)*NORMAL_S(IJK,3)
               ENDIF


               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJKB,east,M) = ZERO
                  A_M(IJKB,west,M) = ZERO
                  A_M(IJKB,north,M) = ZERO
                  A_M(IJKB,south,M) = ZERO
                  A_M(IJKB,top,M) = ZERO
                  A_M(IJKB,bottom,M) = ZERO
                  A_M(IJKB,0,M) = -ONE

                  IF(BC_W_g(BCV)/=UNDEFINED) THEN
                     B_M(IJKB,M) = - BC_W_g(BCV)
                  ELSE
                     B_M(IJKB,M) = - BC_VELMAG_g(BCV)*NORMAL_S(IJK,3)
                  ENDIF


               ENDIF

            CASE (CG_PO)

               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO

               IJKB = BOTTOM_OF(IJK)
               IF(FLUID_AT(IJKB)) THEN

                  A_M(IJK,bottom,M) = ONE
                  A_M(IJK,0,M) = -ONE

               ENDIF

         END SELECT


      ENDDO

      RETURN


    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CG_SOURCE_W_G_BC
