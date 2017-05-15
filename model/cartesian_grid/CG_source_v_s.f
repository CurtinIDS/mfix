!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CG_SOURCE_V_s                                           C
!  Purpose: Determine contribution of cut-cell to source terms         C
!           for V_s momentum eq.                                       C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CG_SOURCE_V_S(A_M, B_M, M)

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
      INTEGER :: I, J, K, IJK,IMJK, IJMK, IJKM, IJKN
! Solids phase index
      INTEGER :: L
! average volume fraction
      DOUBLE PRECISION EPSA, EPStmp
! virtual (added) mass
      DOUBLE PRECISION :: F_vir, ROP_MA
      DOUBLE PRECISION :: Vgn, Vgs, Uge, Ugw, Ugc, Vge, &
                          Vgw, Wgt, Wgb, Wgc,Vgt, Vgb
!
      DOUBLE PRECISION :: F_2
! for cartesian grid
      INTEGER :: IM,JM,IP,JP,KM,KP
      INTEGER :: IPJK,IJPK,IJKP,IJKC,IJKE,IJKNE,IJKW,IJKWN,IMJPK,IJPKM
      INTEGER :: IJKT,IJKTN,IJKB,IJKBN
      DOUBLE PRECISION :: Vn,Vs,Ve,Vw, Vt,Vb
      DOUBLE PRECISION :: B_NOC
      DOUBLE PRECISION :: MU_S_E,MU_S_W,MU_S_N,MU_S_S,MU_S_T,MU_S_B,MU_S_CUT
      DOUBLE PRECISION :: VW_s
      INTEGER :: BCV
      INTEGER :: BCT
!-----------------------------------------------

      IF(CG_SAFE_MODE(4)==1) RETURN

      IF(KT_TYPE_ENUM /= GHD_2007 .OR. &
         (KT_TYPE_ENUM == GHD_2007 .AND. M==MMAX)) THEN

         IF (MOMENTUM_Y_EQ(M)) THEN

!!$omp  parallel do private(I, J, K, IJK, IJKN, &
!!$omp&                     EPSA,EPStmp,ROP_MA,F_vir) &
!!$omp&  schedule(static)
            DO IJK = ijkstart3, ijkend3

! Skip walls where some values are undefined.
               IF(WALL_AT(IJK)) cycle

               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)
               IMJK = IM_OF(IJK)
               IJMK = JM_OF(IJK)
               IJKM = KM_OF(IJK)
               IJKN = NORTH_OF(IJK)

               IF (KT_TYPE_ENUM == GHD_2007) THEN   ! with ghd theory, m = mmax
                  EPStmp = ZERO
                  DO L = 1, SMAX
                     EPStmp = EPStmp + AVG_Y(EP_S(IJK,L),EP_S(IJKN,L),J)
                  ENDDO
                  EPSA = EPStmp
               ELSE
                  EPSA = AVG_Y(EP_S(IJK,M),EP_S(IJKN,M),J)
               ENDIF


! Impermeable internal surface
               IF (IP_AT_N(IJK)) THEN   ! do nothing

! Semi-permeable internal surface
               ELSEIF (SIP_AT_N(IJK)) THEN   ! do nothing

! dilute flow
               ELSEIF (EPSA <= DIL_EP_S) THEN   ! do nothing


               ELSEIF(INTERIOR_CELL_AT(IJK)) THEN
                  BCV = BC_V_ID(IJK)
                  IF(BCV > 0 ) THEN
                     BCT = BC_TYPE_ENUM(BCV)
                  ELSE
                     BCT = NONE
                  ENDIF

                  SELECT CASE (BCT)
                     CASE (CG_NSW)
                        NOC_VS = .TRUE.
                        VW_s = ZERO
                        MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKN)*MU_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
                        A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_S_CUT * Area_V_CUT(IJK)/DELH_V(IJK)
                     CASE (CG_FSW)
                        NOC_VS = .FALSE.
                        VW_s = ZERO
                     CASE(CG_PSW)
                        IF(BC_JJ_PS(BCV)==1) THEN   ! Johnson-Jackson partial slip bc
                           NOC_VS = .FALSE.
                           VW_s = BC_VW_S(BCV,M)
                           CALL CG_CALC_GRBDRY(IJK, 'V_MOMENTUM', M, BCV, F_2)
                           A_M(IJK,0,M) = A_M(IJK,0,M) - Area_V_CUT(IJK)*F_2
                           B_M(IJK,M) = B_M(IJK,M) - Area_V_CUT(IJK)*F_2*VW_s
                        ELSEIF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                           NOC_VS = .TRUE.
                           VW_s = BC_VW_S(BCV,M)
                           MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKN)*MU_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
                           A_M(IJK,0,M) = A_M(IJK,0,M)  - MU_S_CUT * Area_V_CUT(IJK)/DELH_V(IJK)
                           B_M(IJK,M) = B_M(IJK,M) - MU_S_CUT * VW_s * Area_V_CUT(IJK)/DELH_V(IJK)
                        ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                           NOC_VS = .FALSE.
                           VW_s = ZERO
                        ELSE                              ! partial slip
                           NOC_VS = .FALSE.
                           VW_s = BC_VW_S(BCV,M)
                           MU_S_CUT = (VOL(IJK)*MU_S(IJK,M) + VOL(IJKN)*MU_S(IJKN,M))/(VOL(IJK) + VOL(IJKN))
                           A_M(IJK,0,M) = A_M(IJK,0,M) - MU_S_CUT * Area_V_CUT(IJK)*(BC_HW_S(BCV,M))
                           B_M(IJK,M) = B_M(IJK,M) - MU_S_CUT * VW_s * Area_V_CUT(IJK)*(BC_HW_S(BCV,M))
                        ENDIF
                     CASE (NONE, CG_MI)
                        NOC_VS = .FALSE.

                  END SELECT


                  IF(NOC_VS) THEN
                     IMJK = IM_OF(IJK)
                     IJMK = JM_OF(IJK)
                     IJKM = KM_OF(IJK)
                     IPJK = IP_OF(IJK)
                     IJPK = JP_OF(IJK)
                     IJKP = KP_OF(IJK)

                     Vn = Theta_Vn_bar(IJK)  * V_S(IJK,M)  + Theta_Vn(IJK)  * V_S(IJPK,M)
                     Vs = Theta_Vn_bar(IJMK) * V_S(IJMK,M) + Theta_Vn(IJMK) * V_S(IJK,M)
                     Ve = Theta_Ve_bar(IJK)  * V_S(IJK,M)  + Theta_Ve(IJK)  * V_S(IPJK,M)
                     Vw = Theta_Ve_bar(IMJK) * V_S(IMJK,M) + Theta_Ve(IMJK) * V_S(IJK,M)

                     IJKN = NORTH_OF(IJK)
                     IF (WALL_AT(IJK)) THEN
                        IJKC = IJKN
                     ELSE
                        IJKC = IJK
                     ENDIF

                     JP = JP1(J)
                     IJKE = EAST_OF(IJK)
                     IJKNE = EAST_OF(IJKN)
                     IM = IM1(I)
                     IJKW = WEST_OF(IJK)
                     IJKWN = NORTH_OF(IJKW)
                     IMJPK = JP_OF(IMJK)
                     KM = KM1(K)
                     IJKT = TOP_OF(IJK)
                     IJKTN = NORTH_OF(IJKT)
                     IJKB = BOTTOM_OF(IJK)
                     IJKBN = NORTH_OF(IJKB)

                     MU_S_E = AVG_Y_H(AVG_X_H(MU_S(IJKC,M),MU_S(IJKE,M),I),&
                              AVG_X_H(MU_S(IJKN,M),MU_S(IJKNE,M),I),J)
                     MU_S_W = AVG_Y_H(AVG_X_H(MU_S(IJKW,M),MU_S(IJKC,M),IM),&
                              AVG_X_H(MU_S(IJKWN,M),MU_S(IJKN,M),IM),J)
                     MU_S_N = MU_S(IJKN,M)
                     MU_S_S = MU_S(IJKC,M)

                     B_NOC =   MU_S_N * Axz_V(IJK)  * (Vn-VW_s) * NOC_V_N(IJK)   *2.0d0&
                             - MU_S_S * Axz_V(IJMK) * (Vs-VW_s) * NOC_V_N(IJMK)  *2.0d0&
                             + MU_S_E * Ayz_V(IJK)  * (Ve-VW_s) * NOC_V_E(IJK)   &
                             - MU_S_W * Ayz_V(IMJK) * (Vw-VW_s) * NOC_V_E(IMJK)

                     IF(DO_K) THEN

                        Vt = Theta_Vt_bar(IJK)  * V_S(IJK,M)  + Theta_Vt(IJK)  * V_S(IJKP,M)
                        Vb = Theta_Vt_bar(IJKM) * V_S(IJKM,M) + Theta_Vt(IJKM) * V_S(IJK,M)

                        MU_S_T = AVG_Y_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                                 AVG_Z_H(MU_S(IJKN,M),MU_S(IJKTN,M),K),J)
                        MU_S_B = AVG_Y_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                                 AVG_Z_H(MU_S(IJKBN,M),MU_S(IJKN,M),KM),J)

                        B_NOC = B_NOC + MU_S_T * Axy_V(IJK)  * (Vt-VW_s) * NOC_V_T(IJK)   &
                                      - MU_S_B * Axy_V(IJKM) * (Vb-VW_s) * NOC_V_T(IJKM)

                     ENDIF
                     B_M(IJK,M) = B_M(IJK,M) + B_NOC
                  ENDIF   ! end if noc_vs


                  IF(CUT_V_TREATMENT_AT(IJK)) THEN
! VIRTUAL MASS SECTION (explicit terms)
! adding transient term dVg/dt to virtual mass term
                     F_vir = ZERO
                     IF(Added_Mass.AND. M==M_AM ) THEN
                        F_vir = ( (V_g(IJK) - V_gO(IJK)) )*ODT*VOL_V(IJK)
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
                        IMJPK = IM_OF(IJPK)
                        IJKN = NORTH_OF(IJK)

! defining gas-particles velocity at momentum cell faces (or scalar cell center)
                        Vge = Theta_Ve_bar(IJK)  * V_g(IJK)  + Theta_Ve(IJK)  * V_g(IPJK)
                        Vgw = Theta_Ve_bar(IMJK) * V_g(IMJK) + Theta_Ve(IMJK) * V_g(IJK)
                        Uge =  Theta_V_ne(IJK)  * U_g(IJK)  + Theta_V_se(IJK)  * U_g(IJPK)
                        Ugw =  Theta_V_ne(IMJK) * U_g(IMJK) + Theta_V_se(IMJK) * U_g(IMJPK)
                        Ugc = (DELX_ve(IJK) * Ugw + DELX_vw(IJK) * Uge) / (DELX_ve(IJK) + DELX_vw(IJK))
                        Vgn = Theta_Vn_bar(IJK)  * V_g(IJK)  + Theta_Vn(IJK)  * V_g(IJPK)
                        Vgs = Theta_Vn_bar(IJMK) * V_g(IJMK) + Theta_Vn(IJMK) * V_g(IJK)
                        IF(DO_K) THEN
                           IJPKM = KM_OF(IJPK)
                           Vgt = Theta_Vt_bar(IJK)  * V_g(IJK)  + Theta_Vt(IJK)  * V_g(IJKP)
                           Vgb = Theta_Vt_bar(IJKM) * V_g(IJKM) + Theta_Vt(IJKM) * V_g(IJK)
                           Wgt = Theta_V_nt(IJK)  * W_g(IJK)  + Theta_V_st(IJK)  * W_g(IJPK)
                           Wgb = Theta_V_nt(IJKM) * W_g(IJKM) + Theta_V_st(IJKM) * W_g(IJPKM)
                           Wgc = (DELZ_vt(IJK) * Wgb + DELZ_vb(IJK) * Wgt) / (DELZ_vt(IJK) + DELZ_vb(IJK))
                           F_vir = F_vir +  Wgc* (Vgt - Vgb)*AXY(IJK)
                        ENDIF

! adding convective terms (U dV/dx + V dV/dy) to virtual mass; W dV/dz added above.
                        F_vir = F_vir + V_g(IJK)*(Vgn - Vgs)*AXZ(IJK) + &
                                Ugc*(Vge - Vgw)*AYZ(IJK)
                        ROP_MA =  (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M) + &
                                   VOL(IJKN)*ROP_g(IJKN)*EP_s(IJKN,M))/&
                                  (VOL(IJK) + VOL(IJKN))
                        F_vir = F_vir * Cv * ROP_MA
                        B_M(IJK,M) = B_M(IJK,M) - F_vir ! adding explicit-part of virtual mass force.
                     ENDIF   ! end if added_mass .and. m==m_am)
                  ENDIF   ! end if cut_v_treatment_at


               ENDIF ! end if sip or ip or dilute flow branch
            ENDDO   ! end do ijk
         ENDIF   ! end if momentum_y_eq
      ENDIF   ! end if for GHD Theory

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CG_SOURCE_V_S

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CG_SOURCE_V_s_BC                                        C
!  Purpose: Determine contribution of cut-cell to source terms         C
!           for V_s momentum eq.                                       C
!                                                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 01-MAY-09  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CG_SOURCE_V_S_BC(A_M, B_M, M)

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
      INTEGER :: IJK, IJKS
! for cartesian grid:
      INTEGER :: BCV
      INTEGER :: BCT
!-----------------------------------------------

      IF(CG_SAFE_MODE(4)==1) RETURN

      DO IJK = ijkstart3, ijkend3

         BCV = BC_V_ID(IJK)
         IF(BCV > 0 ) THEN
            BCT = BC_TYPE_ENUM(BCV)
         ELSE
            BCT = NONE
         ENDIF

         SELECT CASE (BCT)
            CASE (CG_NSW)
               IF(WALL_V_AT(IJK)) THEN
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
               IF(WALL_V_AT(IJK)) THEN
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
!                  B_M(IJK,M) = - V_s(V_MASTER_OF(IJK),M)    ! Velocity of master node
                  B_M(IJK,M) = ZERO
                  IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN
                     IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                        A_M(IJK,east,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                        A_M(IJK,west,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                        A_M(IJK,north,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                        A_M(IJK,south,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                        A_M(IJK,top,M) = ONE
                     ELSEIF (V_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                        A_M(IJK,bottom,M) = ONE
                     ENDIF
                  ENDIF
               ENDIF


            CASE (CG_PSW)
               IF(WALL_V_AT(IJK)) THEN
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  IF(BC_HW_S(BCV,M)==UNDEFINED) THEN   ! same as NSW
                     B_M(IJK,M) = -BC_VW_S(BCV,M)
                  ELSEIF(BC_HW_S(BCV,M)==ZERO) THEN   ! same as FSW
                     B_M(IJK,M) = ZERO
                     IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN
                        IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                           A_M(IJK,east,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                           A_M(IJK,west,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                           A_M(IJK,north,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                           A_M(IJK,south,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                           A_M(IJK,top,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
                           A_M(IJK,bottom,M) = ONE
                        ENDIF
                     ENDIF
                  ELSE                              ! partial slip  (WARNING:currently same as FSW)
                     B_M(IJK,M) = ZERO
                     IF(DABS(NORMAL_V(IJK,2))/=ONE) THEN
                        IF (V_MASTER_OF(IJK) == EAST_OF(IJK)) THEN
                           A_M(IJK,east,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == WEST_OF(IJK)) THEN
                           A_M(IJK,west,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == NORTH_OF(IJK)) THEN
                           A_M(IJK,north,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == SOUTH_OF(IJK)) THEN
                           A_M(IJK,south,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == TOP_OF(IJK)) THEN
                           A_M(IJK,top,M) = ONE
                        ELSEIF (V_MASTER_OF(IJK) == BOTTOM_OF(IJK)) THEN
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
               IF(BC_V_s(BCV,M)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_V_s(BCV,M)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_s(BCV,M)*NORMAL_V(IJK,2)
               ENDIF
               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN
                  A_M(IJKS,east,M) = ZERO
                  A_M(IJKS,west,M) = ZERO
                  A_M(IJKS,north,M) = ZERO
                  A_M(IJKS,south,M) = ZERO
                  A_M(IJKS,top,M) = ZERO
                  A_M(IJKS,bottom,M) = ZERO
                  A_M(IJKS,0,M) = -ONE
                  IF(BC_V_s(BCV,M)/=UNDEFINED) THEN
                     B_M(IJKS,M) = - BC_V_s(BCV,M)
                  ELSE
                     B_M(IJKS,M) = - BC_VELMAG_s(BCV,M)*NORMAL_V(IJK,2)
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

               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN

                  A_M(IJK,south,M) = ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF

         END SELECT   ! end select bc_type


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
               IF(BC_V_s(BCV,M)/=UNDEFINED) THEN
                  B_M(IJK,M) = - BC_V_s(BCV,M)
               ELSE
                  B_M(IJK,M) = - BC_VELMAG_s(BCV,M)*NORMAL_S(IJK,2)
               ENDIF
               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN
                  A_M(IJKS,east,M) = ZERO
                  A_M(IJKS,west,M) = ZERO
                  A_M(IJKS,north,M) = ZERO
                  A_M(IJKS,south,M) = ZERO
                  A_M(IJKS,top,M) = ZERO
                  A_M(IJKS,bottom,M) = ZERO
                  A_M(IJKS,0,M) = -ONE
                  IF(BC_V_s(BCV,M)/=UNDEFINED) THEN
                     B_M(IJKS,M) = - BC_V_s(BCV,M)
                  ELSE
                     B_M(IJKS,M) = - BC_VELMAG_s(BCV,M)*NORMAL_S(IJK,2)
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
               IJKS = SOUTH_OF(IJK)
               IF(FLUID_AT(IJKS)) THEN
                  A_M(IJK,south,M) = ONE
                  A_M(IJK,0,M) = -ONE
               ENDIF

         END SELECT

      ENDDO

      RETURN

    CONTAINS

      INCLUDE 'functions.inc'

      END SUBROUTINE CG_SOURCE_V_S_BC
