!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CG_SOURCE_W_s                                           C
!  Purpose: Determine contribution of cut-cell to source terms         C
!           for W_s momentum eq.                                       C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CG_SOURCE_W_S(A_M, B_M, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param, only: dimension_3, dimension_m
      USE param1, only: zero, undefined

! number of solids phases
      USE physprop, only: smax, mmax
! virtual (added) mass coefficient Cv
      USE physprop, only: cv

      USE bc
      USE calc_gr_boundary
      USE compar
      USE cutcell
      USE fldvar
      USE fun_avg
      USE geometry
      USE indices
      USE quadric
      USE run
      USE toleranc, only: dil_ep_s
      USE visc_s, only: mu_s

      IMPLICIT NONE

!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Solids phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IJKT, IMJK, IJMK, IJKM, IJKP, IMJKP
      INTEGER :: IJKE, IJKW, IJKTE, IM, IPJK
! Solids phase index
      INTEGER :: L
! Average volume fraction
      DOUBLE PRECISION :: EPSA, EPStmp
! virtual (added) mass
      DOUBLE PRECISION :: F_vir, ROP_MA
      DOUBLE PRECISION :: Uge, Ugw, Wge, Wgw, Wgn, &
                          Wgs, Wgt, Wgb, Ugc, Vgc, Vgn, Vgs
!
      DOUBLE PRECISION :: F_2
! for cartesian grid:
      INTEGER :: JM,IP,JP,IJPK,IJKC,IJKN,IJKNE,IJKS,IJKSE,IPJMK,KM,KP,IJMKP
      INTEGER :: IJKTN,IJKWT,IJKST
      DOUBLE PRECISION :: We,Ww,Wn,Ws,Wt,Wb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_S_E,MU_S_W,MU_S_N,MU_S_S,MU_S_T,MU_S_B,MU_S_CUT
      DOUBLE PRECISION :: WW_s
      INTEGER :: BCV
      INTEGER :: BCT
!-----------------------------------------------

      IF(CG_SAFE_MODE(5)==1) RETURN


      IF(KT_TYPE_ENUM /= GHD_2007 .OR. &
         (KT_TYPE_ENUM == GHD_2007 .AND. M==MMAX)) THEN

         IF (MOMENTUM_Z_EQ(M)) THEN

!!$omp  parallel do &
!!$omp& private(IJK, I, J, K, IJKT, EPSA, EPStmp, &
!!$omp&         IMJK,IJKP,IMJKP, &
!!$omp&         IJKE,IJKW,IJKTE,IJKTW,IM,IPJK, &
!!$omp&         ROP_MA,F_vir)
            DO IJK = ijkstart3, ijkend3

! Skip walls where some values are undefined.
               IF(WALL_AT(IJK)) cycle

               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)
               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IJKT = TOP_OF(IJK)

               IF (KT_TYPE_ENUM == GHD_2007) THEN   ! with ghd theory, m = mmax
                  EPStmp = ZERO
                  DO L = 1, SMAX
                    EPStmp = EPStmp + AVG_Z(EP_S(IJK,L),EP_S(IJKT,L),K)
                  ENDDO
                  EPSA = EPStmp
                ELSE
                  EPSA = AVG_Z(EP_S(IJK,M),EP_S(IJKT,M),K)
                ENDIF

! Impermeable internal surface
                IF (IP_AT_T(IJK)) THEN   ! do nothing

! Semi-permeable internal surface
                ELSEIF (SIP_AT_T(IJK)) THEN   ! do nothing

! dilute flow
                ELSEIF (EPSA <= DIL_EP_S) THEN   ! do nothing


                ELSEIF(INTERIOR_CELL_AT(IJK)) THEN
                  BCV = BC_W_ID(IJK)
                  IF(BCV > 0 ) THEN
                     BCT = BC_TYPE_ENUM(BCV)
                  ELSE
                     BCT = NONE
                  ENDIF

                  SELECT CASE (BCT)
                     CASE (CG_NSW)
                        NOC_WS = .TRUE.
                        WW_s = ZERO
                        MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKT)*MU_S(IJKT,M))/(VOL(IJK) + VOL(IJKT))
                        A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_W_CUT(IJK)/DELH_W(IJK)
                     CASE (CG_FSW)
                        NOC_WS = .FALSE.
                        WW_s = ZERO
                     CASE(CG_PSW)
                        IF(BC_JJ_PS(BCV)==1) THEN   ! Johnson-Jackson partial slip bc
                           NOC_WS = .FALSE.
                           WW_s = BC_WW_S(BCV,M)
                           CALL CG_CALC_GRBDRY(IJK, 'W_MOMENTUM', M, BCV, F_2)
                           A_M(IJK,0,M) = A_M(IJK,0,M) - Area_W_CUT(IJK)*F_2
                           B_M(IJK,M) = B_M(IJK,M) - Area_W_CUT(IJK)*F_2*WW_s
                        ELSEIF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                           NOC_WS = .TRUE.
                           WW_s = BC_WW_S(BCV,M)
                           MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKT)*MU_S(IJKT,M))/(VOL(IJK) + VOL(IJKT))
                           A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_W_CUT(IJK)/DELH_W(IJK)
                           B_M(IJK,M) = B_M(IJK,M) - MU_S_CUT * WW_s * Area_W_CUT(IJK)/DELH_W(IJK)
                        ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                           NOC_WS = .FALSE.
                           WW_s = ZERO
                        ELSE                              ! partial slip
                           NOC_WS = .FALSE.
                           WW_s = BC_WW_S(BCV,M)
                           MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKT)*MU_S(IJKT,M))/(VOL(IJK) + VOL(IJKT))
                           A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_W_CUT(IJK)*(BC_HW_S(BCV,M))
                           B_M(IJK,M) = B_M(IJK,M) - MU_S_CUT * WW_s * Area_W_CUT(IJK)*(BC_HW_S(BCV,M))
                        ENDIF
                     CASE (NONE,CG_MI)
                        NOC_WS = .FALSE.

                  END SELECT


                  IF(NOC_WS) THEN
                     IMJK = IM_OF(IJK)
                     IJMK = JM_OF(IJK)
                     IJKM = KM_OF(IJK)
                     IPJK = IP_OF(IJK)
                     IJPK = JP_OF(IJK)
                     IJKP = KP_OF(IJK)

                     We = Theta_We_bar(IJK)  * W_S(IJK,M)  + Theta_We(IJK)  * W_S(IPJK,M)
                     Ww = Theta_We_bar(IMJK) * W_S(IMJK,M) + Theta_We(IMJK) * W_S(IJK,M)
                     Wn = Theta_Wn_bar(IJK)  * W_S(IJK,M)  + Theta_Wn(IJK)  * W_S(IJPK,M)
                     Ws = Theta_Wn_bar(IJMK) * W_S(IJMK,M) + Theta_Wn(IJMK) * W_S(IJK,M)
                     Wt = Theta_Wt_bar(IJK)  * W_S(IJK,M)  + Theta_Wt(IJK)  * W_S(IJKP,M)
                     Wb = Theta_Wt_bar(IJKM) * W_S(IJKM,M) + Theta_Wt(IJKM) * W_S(IJK,M)

                     IF (WALL_AT(IJK)) THEN
                        IJKC = IJKT
                     ELSE
                        IJKC = IJK
                     ENDIF

                     IJKE = EAST_OF(IJK)
                     ijkt = top_of(ijk)
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

                     MU_S_E = AVG_Z_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                              AVG_X_H(MU_S(IJKT,M),MU_S(IJKTE,M),I),K)
                     MU_S_W = AVG_Z_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                              AVG_X_H(MU_S(IJKWT,M),MU_S(IJKT,M),IM),K)
                     MU_S_N = AVG_Z_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                              AVG_Y_H(MU_S(IJKT,M),MU_S(IJKTN,M),J),K)
                     MU_S_S = AVG_Z_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                              AVG_Y_H(MU_S(IJKST,M),MU_S(IJKT,M),JM),K)
                     MU_S_T = MU_S(IJKT,M)
                     MU_S_B = MU_S(IJKC,M)

                     B_NOC =     MU_S_E * Ayz_W(IJK)  * (We-WW_s) * NOC_W_E(IJK)  &
                             -   MU_S_W * Ayz_W(IMJK) * (Ww-WW_s) * NOC_W_E(IMJK) &
                             +   MU_S_N * Axz_W(IJK)  * (Wn-WW_s) * NOC_W_N(IJK)  &
                             -   MU_S_S * Axz_W(IJMK) * (Ws-WW_s) * NOC_W_N(IJMK) &
                             +   MU_S_T * Axy_W(IJK)  * (Wt-WW_s) * NOC_W_T(IJK)  *2.0d0&
                             -   MU_S_B * Axy_W(IJKM) * (Wb-WW_s) * NOC_W_T(IJKM) *2.0D0

                     B_M(IJK,M) = B_M(IJK,M)   +  B_NOC
                  ENDIF


                  IF(CUT_W_TREATMENT_AT(IJK)) THEN
! VIRTUAL MASS SECTION (explicit terms)
! adding transient term  dWg/dt to virtual mass term
                     F_vir = ZERO
                     IF(Added_Mass.AND. M == M_AM ) THEN
                        F_vir = ( (W_g(IJK) - W_gO(IJK)) )*ODT*VOL_W(IJK)
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
                        Wge = Theta_We_bar(IJK) * W_g(IJK) + Theta_We(IJK) * W_g(IPJK)
                        Wgw = Theta_We_bar(IMJK) * W_g(IMJK) + Theta_We(IMJK) * W_g(IJK)
                        Uge = Theta_W_te(IJK) * U_g(IJK) + Theta_W_be(IJK) * U_g(IJKP)
                        Ugw = Theta_W_te(IMJK) * U_g(IMJK) + Theta_W_be(IMJK) * U_g(IMJKP)
                        Ugc = (DELX_we(IJK) * Ugw + DELX_ww(IJK) * Uge) / (DELX_we(IJK) + DELX_ww(IJK))
                        Wgn = Theta_Wn_bar(IJK) * W_g(IJK) + Theta_Wn(IJK) * W_g(IJPK)
                        Wgs = Theta_Wn_bar(IJMK) * W_g(IJMK) + Theta_Wn(IJMK) * W_g(IJK)
                        Vgn =  Theta_W_tn(IJK)  * V_g(IJK)  + Theta_W_bn(IJK)  * V_g(IJKP)
                        Vgs =  Theta_W_tn(IJMK) * V_g(IJMK) + Theta_W_bn(IJMK) * V_g(IJMKP)
                        Vgc = (DELY_wn(IJK) * Vgs + DELY_ws(IJK) * Vgn) / (DELY_wn(IJK) + DELY_ws(IJK))
                        Wgt = Theta_Wt_bar(IJK)  * W_g(IJK)  + Theta_Wt(IJK)  * W_g(IJKP)
                        Wgb = Theta_Wt_bar(IMJK) * W_g(IMJK) + Theta_Wt(IMJK) * W_g(IMJKP)

! adding convective terms (U dW/dx + V dW/dy + W dW/dz) to virtual mass
                        F_vir = F_vir +  Ugc * (Wge - Wgw)*AYZ(IJK)    + &
                                         Vgc * (Wgn - Wgs)*AXZ(IJK)  + &
                                         W_g(IJK)*(Wgt - Wgb)*AXY(IJK)
                        ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M) + &
                                  VOL(IJKT)*ROP_g(IJKT)*EP_s(IJKT,M))/&
                                 (VOL(IJK) + VOL(IJKT))
                        F_vir = F_vir * Cv * ROP_MA
                        B_M(IJK,M) = B_M(IJK,M) - F_vir ! explicit part of virtual mass force
                     ENDIF   ! end if added_mass .and. m==m_am)
                  ENDIF   ! end if cut_w_treatment_at


               ENDIF ! end if sip or ip or dilute flow branch
            ENDDO   ! end do ijk
         ENDIF   ! end if momentum_z_eq
      ENDIF   ! end if for GHD Theory

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CG_SOURCE_W_S

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine name: CG_SOURCE_W_s_BC                                   C
!  Purpose: Determine contribution of cut-cell to source terms         C
!           for W_s momentum eq.                                       C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CG_SOURCE_W_S_BC(A_M, B_M, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1, only: zero, one, undefined

      USE fldvar
      USE bc

      USE compar
      USE geometry
      USE indices

      USE cutcell
      USE quadric
      USE fun_avg

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Solids phase
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKB
! for cartesian grid:
      INTEGER :: BCV
      INTEGER :: BCT
!-----------------------------------------------

      IF(CG_SAFE_MODE(5)==1) RETURN

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
!                  B_M(IJK,M) = - W_s(W_MASTER_OF(IJK),M)  ! Velocity of master node
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
                  IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = -BC_WW_S(BCV,M)
                  ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
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
                  ELSE                              ! partial slip  (WARNING:currently same as FSW)
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
               ENDIF


            CASE (CG_MI)
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               IF(BC_W_s(BCV,M)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_W_s(BCV,M)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_s(BCV,M)*NORMAL_W(IJK,3)
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
                  IF(BC_W_s(BCV,M)/=UNDEFINED) THEN
                     B_M(IJKB,M) = - BC_W_s(BCV,M)
                  ELSE
                     B_M(IJKB,M) = - BC_VELMAG_s(BCV,M)*NORMAL_W(IJK,3)
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


! this is straight repeat of section of above code...?
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
               IF(BC_W_s(BCV,M)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_W_s(BCV,M)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_s(BCV,M)*NORMAL_S(IJK,3)
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
                  IF(BC_W_s(BCV,M)/=UNDEFINED) THEN
                     B_M(IJKB,M) = - BC_W_s(BCV,M)
                  ELSE
                     B_M(IJKB,M) = - BC_VELMAG_s(BCV,M)*NORMAL_S(IJK,3)
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

      END SUBROUTINE CG_SOURCE_W_S_BC
