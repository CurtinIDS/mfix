!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SEND_RECEIVE_CUT_CELL_VARIABLES                        C
!  Purpose: Send/receive all relevant cut cell related variables       C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SEND_RECEIVE_CUT_CELL_VARIABLES

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell
      USE functions

      USE mpi_utility      !//d pnicol : for gather


      IMPLICIT NONE

      call SEND_RECEIVE_1D_LOGICAL(WALL_U_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(WALL_V_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(WALL_W_AT,2)

      call send_recv(area_cut, 2)
      call send_recv(Area_U_CUT,2)
      call send_recv(Area_V_CUT,2)
      call send_recv(Area_W_CUT,2)

      call send_recv(DELX_Ue,2)
      call send_recv(DELX_Uw,2)
      call send_recv(DELY_Un,2)
      call send_recv(DELY_Us,2)
      call send_recv(DELZ_Ut,2)
      call send_recv(DELZ_Ub,2)

      call send_recv(DELX_Ve,2)
      call send_recv(DELX_Vw,2)
      call send_recv(DELY_Vn,2)
      call send_recv(DELY_Vs,2)
      call send_recv(DELZ_Vt,2)
      call send_recv(DELZ_Vb,2)

      call send_recv(DELX_We,2)
      call send_recv(DELX_Ww,2)
      call send_recv(DELY_Wn,2)
      call send_recv(DELY_Ws,2)
      call send_recv(DELZ_Wt,2)
      call send_recv(DELZ_Wb,2)

      call send_recv(X_U,2)
      call send_recv(Y_U,2)
      call send_recv(Z_U,2)

      call send_recv(X_U_ec,2)
      call send_recv(Y_U_ec,2)
      call send_recv(Z_U_ec,2)

      call send_recv(X_U_nc,2)
      call send_recv(Y_U_nc,2)
      call send_recv(Z_U_nc,2)

      call send_recv(X_U_tc,2)
      call send_recv(Y_U_tc,2)
      call send_recv(Z_U_tc,2)

      call send_recv(X_V,2)
      call send_recv(Y_V,2)
      call send_recv(Z_V,2)

      call send_recv(X_V_ec,2)
      call send_recv(Y_V_ec,2)
      call send_recv(Z_V_ec,2)

      call send_recv(X_V_nc,2)
      call send_recv(Y_V_nc,2)
      call send_recv(Z_V_nc,2)

      call send_recv(X_V_tc,2)
      call send_recv(Y_V_tc,2)
      call send_recv(Z_V_tc,2)

      call send_recv(X_W,2)
      call send_recv(Y_W,2)
      call send_recv(Z_W,2)

      call send_recv(X_W_ec,2)
      call send_recv(Y_W_ec,2)
      call send_recv(Z_W_ec,2)

      call send_recv(X_W_nc,2)
      call send_recv(Y_W_nc,2)
      call send_recv(Z_W_nc,2)

      call send_recv(X_W_tc,2)
      call send_recv(Y_W_tc,2)
      call send_recv(Z_W_tc,2)

      call send_recv(DELH_Scalar,2)
      call send_recv(DELH_U,2)
      call send_recv(DELH_V,2)
      call send_recv(DELH_W,2)

      call send_recv(NORMAL_S,2)
      call send_recv(NORMAL_U,2)
      call send_recv(NORMAL_V,2)
      call send_recv(NORMAL_W,2)

      call send_recv(REFP_S,2)
      call send_recv(REFP_U,2)
      call send_recv(REFP_V,2)
      call send_recv(REFP_W,2)

      call send_recv(Theta_Ue,2)
      call send_recv(Theta_Ue_bar,2)
      call send_recv(Theta_U_ne,2)
      call send_recv(Theta_U_nw,2)
      call send_recv(Theta_U_te,2)
      call send_recv(Theta_U_tw,2)
      call send_recv(ALPHA_Ue_c,2)
      call send_recv(NOC_U_E,2)
      call send_recv(Theta_Un,2)
      call send_recv(Theta_Un_bar,2)
      call send_recv(ALPHA_Un_c,2)
      call send_recv(NOC_U_N,2)
      call send_recv(Theta_Ut,2)
      call send_recv(Theta_Ut_bar,2)
      call send_recv(ALPHA_Ut_c,2)
      call send_recv(NOC_U_T,2)
      call send_recv(A_UPG_E,2)
      call send_recv(A_UPG_W,2)

      call send_recv(Theta_V_ne,2)
      call send_recv(Theta_V_se,2)
      call send_recv(Theta_Vn,2)
      call send_recv(Theta_Vn_bar,2)
      call send_recv(Theta_V_nt,2)
      call send_recv(Theta_V_st,2)
      call send_recv(Theta_Ve,2)
      call send_recv(Theta_Ve_bar,2)
      call send_recv(ALPHA_Ve_c,2)
      call send_recv(NOC_V_E,2)
      call send_recv(ALPHA_Vn_c,2)
      call send_recv(NOC_V_N,2)
      call send_recv(Theta_Vt,2)
      call send_recv(Theta_Vt_bar,2)
      call send_recv(ALPHA_Vt_c,2)
      call send_recv(NOC_V_T,2)
      call send_recv(A_VPG_N,2)
      call send_recv(A_VPG_S,2)

      call send_recv(Theta_W_te,2)
      call send_recv(Theta_W_be,2)
      call send_recv(Theta_W_tn,2)
      call send_recv(Theta_W_bn,2)
      call send_recv(Theta_Wt,2)
      call send_recv(Theta_Wt_bar,2)
      call send_recv(Theta_We,2)
      call send_recv(Theta_We_bar,2)
      call send_recv(ALPHA_We_c,2)
      call send_recv(NOC_W_E,2)
      call send_recv(Theta_Wn,2)
      call send_recv(Theta_Wn_bar,2)
      call send_recv(ALPHA_Wn_c,2)
      call send_recv(NOC_W_N,2)
      call send_recv(ALPHA_Wt_c,2)
      call send_recv(NOC_W_T,2)
      call send_recv(A_WPG_T,2)
      call send_recv(A_WPG_B,2)

      call send_recv(ONEoDX_E_U,2)
      call send_recv(ONEoDY_N_U,2)
      call send_recv(ONEoDZ_T_U,2)

      call send_recv(ONEoDX_E_V,2)
      call send_recv(ONEoDY_N_V,2)
      call send_recv(ONEoDZ_T_V,2)

      call send_recv(ONEoDX_E_W,2)
      call send_recv(ONEoDY_N_W,2)
      call send_recv(ONEoDZ_T_W,2)

      call SEND_RECEIVE_1D_LOGICAL(CUT_TREATMENT_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_U_TREATMENT_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_V_TREATMENT_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_W_TREATMENT_AT,2)

      call SEND_RECEIVE_1D_LOGICAL(CUT_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_U_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_V_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(CUT_W_CELL_AT,2)

      call SEND_RECEIVE_1D_LOGICAL(SMALL_CELL_AT,2)
      call send_recv(SMALL_CELL_FLAG,2)


      call SEND_RECEIVE_1D_LOGICAL(BLOCKED_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(BLOCKED_U_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(BLOCKED_V_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(BLOCKED_W_CELL_AT,2)

      call SEND_RECEIVE_1D_LOGICAL(STANDARD_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(STANDARD_U_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(STANDARD_V_CELL_AT,2)
      call SEND_RECEIVE_1D_LOGICAL(STANDARD_W_CELL_AT,2)

      call send_recv(U_MASTER_OF,2)
      call send_recv(V_MASTER_OF,2)
      call send_recv(W_MASTER_OF,2)

      call send_recv(BC_ID,2)
      call send_recv(BC_U_ID,2)
      call send_recv(BC_V_ID,2)
      call send_recv(BC_W_ID,2)

      call send_recv(FLAG,2)
      call send_recv(FLAG_E,2)
      call send_recv(FLAG_N,2)
      call send_recv(FLAG_T,2)

      call send_recv(AYZ,2)
      call send_recv(AXZ,2)
      call send_recv(AXY,2)
      call send_recv(VOL,2)

      call send_recv(AYZ_U,2)
      call send_recv(AXZ_U,2)
      call send_recv(AXY_U,2)
      call send_recv(VOL_U,2)

      call send_recv(AYZ_V,2)
      call send_recv(AXZ_V,2)
      call send_recv(AXY_V,2)
      call send_recv(VOL_V,2)

      call send_recv(AYZ_W,2)
      call send_recv(AXZ_W,2)
      call send_recv(AXY_W,2)
      call send_recv(VOL_W,2)

      RETURN


      END SUBROUTINE SEND_RECEIVE_CUT_CELL_VARIABLES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SEND_RECEIVE_2D_LOGICAL                                C
!  Purpose: Emulates send/receive for a 2D logical array               C
!  using temporary integer                                             C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SEND_RECEIVE_1D_LOGICAL(L1D,NLAYERS)

      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell

      USE mpi_utility      !//d pnicol : for gather

      IMPLICIT NONE
      INTEGER :: IJK,NLAYERS
      INTEGER, DIMENSION(DIMENSION_3) :: I1D
      LOGICAL, DIMENSION(DIMENSION_3) :: L1D

      IF((NLAYERS/=1).AND.(NLAYERS/=2)) THEN
         WRITE(*,*)' NLAYERS=',NLAYERS
         WRITE(*,*)' SEND_RECEIVE_1D_LOGICAL ERROR: NLAYER MUST BE EQUAL TO 1 or 2'
         CALL MFIX_EXIT(MYPE)
      ENDIF

      DO IJK = IJKSTART3, IJKEND3
         IF(L1D(IJK)) THEN
            I1D(IJK) = 1
         ELSE
            I1D(IJK) = 0
         ENDIF
      ENDDO

      call send_recv(I1D,NLAYERS)

      DO IJK = IJKSTART3, IJKEND3
         IF(I1D(IJK)==1) THEN
            L1D(IJK) = .TRUE.
         ELSE
            L1D(IJK) = .FALSE.
         ENDIF
      ENDDO


      RETURN


      END SUBROUTINE SEND_RECEIVE_1D_LOGICAL











