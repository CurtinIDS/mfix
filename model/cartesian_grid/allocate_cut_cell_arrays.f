
      SUBROUTINE ALLOCATE_CUT_CELL_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: ALLOCATE_ARRAYS
!  Purpose: allocate arrays
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      Use indices

      USE cutcell
      USE stl
      USE discretelement

      IMPLICIT NONE

      DIMENSION_MAX_CUT_CELL = INT(FAC_DIM_MAX_CUT_CELL*DIMENSION_3G)

      Allocate(  INTERIOR_CELL_AT  (DIMENSION_3) )

      Allocate( XG_E(0:DIMENSION_I) )
      Allocate( YG_N(0:DIMENSION_J) )
      Allocate( ZG_T(0:DIMENSION_K) )

      Allocate(  X_U (DIMENSION_3) )
      Allocate(  Y_U (DIMENSION_3) )
      Allocate(  Z_U (DIMENSION_3) )

      Allocate(  X_V (DIMENSION_3) )
      Allocate(  Y_V (DIMENSION_3) )
      Allocate(  Z_V (DIMENSION_3) )

      Allocate(  X_W (DIMENSION_3) )
      Allocate(  Y_W (DIMENSION_3) )
      Allocate(  Z_W (DIMENSION_3) )

      Allocate(  INTERSECT_X  (DIMENSION_3) )
      Allocate(  INTERSECT_Y  (DIMENSION_3) )
      Allocate(  INTERSECT_Z  (DIMENSION_3) )

      Allocate(  X_int (DIMENSION_3) )
      Allocate(  Y_int (DIMENSION_3) )
      Allocate(  Z_int (DIMENSION_3) )

      Allocate(  X_NEW_POINT  (DIMENSION_MAX_CUT_CELL) )
      Allocate(  Y_NEW_POINT  (DIMENSION_MAX_CUT_CELL) )
      Allocate(  Z_NEW_POINT  (DIMENSION_MAX_CUT_CELL) )

      Allocate(  X_NEW_U_POINT  (DIMENSION_MAX_CUT_CELL) )
      Allocate(  Y_NEW_U_POINT  (DIMENSION_MAX_CUT_CELL) )
      Allocate(  Z_NEW_U_POINT  (DIMENSION_MAX_CUT_CELL) )

      Allocate(  X_NEW_V_POINT  (DIMENSION_MAX_CUT_CELL) )
      Allocate(  Y_NEW_V_POINT  (DIMENSION_MAX_CUT_CELL) )
      Allocate(  Z_NEW_V_POINT  (DIMENSION_MAX_CUT_CELL) )

      Allocate(  X_NEW_W_POINT  (DIMENSION_MAX_CUT_CELL) )
      Allocate(  Y_NEW_W_POINT  (DIMENSION_MAX_CUT_CELL) )
      Allocate(  Z_NEW_W_POINT  (DIMENSION_MAX_CUT_CELL) )

      Allocate(  NUMBER_OF_NODES  (DIMENSION_3) )
      Allocate(  NUMBER_OF_U_NODES  (DIMENSION_3) )
      Allocate(  NUMBER_OF_V_NODES  (DIMENSION_3) )
      Allocate(  NUMBER_OF_W_NODES  (DIMENSION_3) )

     NUMBER_OF_NODES   = 0
     NUMBER_OF_U_NODES = 0
     NUMBER_OF_V_NODES = 0
     NUMBER_OF_W_NODES = 0

      Allocate(  CONNECTIVITY  (DIMENSION_3,15) )
      Allocate(  CONNECTIVITY_U  (DIMENSION_3,15) )
      Allocate(  CONNECTIVITY_V  (DIMENSION_3,15) )
      Allocate(  CONNECTIVITY_W  (DIMENSION_3,15) )

      Allocate(  PARTITION  (DIMENSION_3) )

      Allocate(  WALL_U_AT (DIMENSION_3) )
      Allocate(  WALL_V_AT (DIMENSION_3) )
      Allocate(  WALL_W_AT (DIMENSION_3) )

      WALL_U_AT = .FALSE.
      WALL_V_AT = .FALSE.
      WALL_W_AT = .FALSE.

      Allocate( Area_CUT  (DIMENSION_3) )
      Allocate( Area_U_CUT  (DIMENSION_3) )
      Allocate( Area_V_CUT  (DIMENSION_3) )
      Allocate( Area_W_CUT  (DIMENSION_3) )


      Allocate( DELX_Ue  (DIMENSION_3) )
      Allocate( DELX_Uw  (DIMENSION_3) )
      Allocate( DELY_Un  (DIMENSION_3) )
      Allocate( DELY_Us  (DIMENSION_3) )
      Allocate( DELZ_Ut  (DIMENSION_3) )
      Allocate( DELZ_Ub  (DIMENSION_3) )

      Allocate( DELX_Ve  (DIMENSION_3) )
      Allocate( DELX_Vw  (DIMENSION_3) )
      Allocate( DELY_Vn  (DIMENSION_3) )
      Allocate( DELY_Vs  (DIMENSION_3) )
      Allocate( DELZ_Vt  (DIMENSION_3) )
      Allocate( DELZ_Vb  (DIMENSION_3) )

      Allocate( DELX_We  (DIMENSION_3) )
      Allocate( DELX_Ww  (DIMENSION_3) )
      Allocate( DELY_Wn  (DIMENSION_3) )
      Allocate( DELY_Ws  (DIMENSION_3) )
      Allocate( DELZ_Wt  (DIMENSION_3) )
      Allocate( DELZ_Wb  (DIMENSION_3) )

      Allocate( X_U_ec  (DIMENSION_3) )
      Allocate( Y_U_ec  (DIMENSION_3) )
      Allocate( Z_U_ec  (DIMENSION_3) )
      Allocate( X_U_nc  (DIMENSION_3) )
      Allocate( Y_U_nc  (DIMENSION_3) )
      Allocate( Z_U_nc  (DIMENSION_3) )
      Allocate( X_U_tc  (DIMENSION_3) )
      Allocate( Y_U_tc  (DIMENSION_3) )
      Allocate( Z_U_tc  (DIMENSION_3) )

      Allocate( X_V_ec  (DIMENSION_3) )
      Allocate( Y_V_ec  (DIMENSION_3) )
      Allocate( Z_V_ec  (DIMENSION_3) )
      Allocate( X_V_nc  (DIMENSION_3) )
      Allocate( Y_V_nc  (DIMENSION_3) )
      Allocate( Z_V_nc  (DIMENSION_3) )
      Allocate( X_V_tc  (DIMENSION_3) )
      Allocate( Y_V_tc  (DIMENSION_3) )
      Allocate( Z_V_tc  (DIMENSION_3) )

      Allocate( X_W_ec  (DIMENSION_3) )
      Allocate( Y_W_ec  (DIMENSION_3) )
      Allocate( Z_W_ec  (DIMENSION_3) )
      Allocate( X_W_nc  (DIMENSION_3) )
      Allocate( Y_W_nc  (DIMENSION_3) )
      Allocate( Z_W_nc  (DIMENSION_3) )
      Allocate( X_W_tc  (DIMENSION_3) )
      Allocate( Y_W_tc  (DIMENSION_3) )
      Allocate( Z_W_tc  (DIMENSION_3) )

      Allocate( DELH_Scalar  (DIMENSION_3) )

      Allocate( DELH_U  (DIMENSION_3) )
      Allocate( Theta_Ue  (DIMENSION_3) )
      Allocate( Theta_Ue_bar (DIMENSION_3) )
      Allocate( Theta_U_ne  (DIMENSION_3) )
      Allocate( Theta_U_nw  (DIMENSION_3) )
      Allocate( Theta_U_te  (DIMENSION_3) )
      Allocate( Theta_U_tw  (DIMENSION_3) )
      Allocate( ALPHA_Ue_c  (DIMENSION_3) )
      Allocate( NOC_U_E  (DIMENSION_3) )
      Allocate( Theta_Un  (DIMENSION_3) )
      Allocate( Theta_Un_bar (DIMENSION_3) )
      Allocate( ALPHA_Un_c  (DIMENSION_3) )
      Allocate( NOC_U_N  (DIMENSION_3) )
      Allocate( Theta_Ut  (DIMENSION_3) )
      Allocate( Theta_Ut_bar (DIMENSION_3) )
      Allocate( ALPHA_Ut_c  (DIMENSION_3) )
      Allocate( NOC_U_T  (DIMENSION_3) )
      Allocate( A_UPG_E (DIMENSION_3) )
      Allocate( A_UPG_W (DIMENSION_3) )

      Allocate( DELH_V  (DIMENSION_3) )
      Allocate( Theta_V_ne  (DIMENSION_3) )
      Allocate( Theta_V_se  (DIMENSION_3) )
      Allocate( Theta_Vn  (DIMENSION_3) )
      Allocate( Theta_Vn_bar (DIMENSION_3) )
      Allocate( Theta_V_nt  (DIMENSION_3) )
      Allocate( Theta_V_st (DIMENSION_3) )
      Allocate( Theta_Ve  (DIMENSION_3) )
      Allocate( Theta_Ve_bar (DIMENSION_3) )
      Allocate( ALPHA_Ve_c  (DIMENSION_3) )
      Allocate( NOC_V_E  (DIMENSION_3) )
      Allocate( ALPHA_Vn_c  (DIMENSION_3) )
      Allocate( NOC_V_N  (DIMENSION_3) )
      Allocate( Theta_Vt  (DIMENSION_3) )
      Allocate( Theta_Vt_bar (DIMENSION_3) )
      Allocate( ALPHA_Vt_c  (DIMENSION_3) )
      Allocate( NOC_V_T  (DIMENSION_3) )
      Allocate( A_VPG_N (DIMENSION_3) )
      Allocate( A_VPG_S (DIMENSION_3) )

      Allocate( DELH_W (DIMENSION_3) )
      Allocate( Theta_W_te (DIMENSION_3) )
      Allocate( Theta_W_be (DIMENSION_3) )
      Allocate( Theta_W_tn (DIMENSION_3) )
      Allocate( Theta_W_bn (DIMENSION_3) )
      Allocate( Theta_Wt (DIMENSION_3) )
      Allocate( Theta_Wt_bar (DIMENSION_3) )
      Allocate( Theta_We (DIMENSION_3) )
      Allocate( Theta_We_bar (DIMENSION_3) )
      Allocate( ALPHA_We_c (DIMENSION_3) )
      Allocate( NOC_W_E (DIMENSION_3) )
      Allocate( Theta_Wn (DIMENSION_3) )
      Allocate( Theta_Wn_bar (DIMENSION_3) )
      Allocate( ALPHA_Wn_c (DIMENSION_3) )
      Allocate( NOC_W_N (DIMENSION_3) )
      Allocate( ALPHA_Wt_c (DIMENSION_3) )
      Allocate( NOC_W_T (DIMENSION_3) )
      Allocate( A_WPG_T (DIMENSION_3) )
      Allocate( A_WPG_B (DIMENSION_3) )


      Allocate( NORMAL_S (DIMENSION_3,3) )
      Allocate( NORMAL_U (DIMENSION_3,3) )
      Allocate( NORMAL_V (DIMENSION_3,3) )
      Allocate( NORMAL_W (DIMENSION_3,3) )

      Allocate( REFP_S (DIMENSION_3,3) )
      Allocate( REFP_U (DIMENSION_3,3) )
      Allocate( REFP_V (DIMENSION_3,3) )
      Allocate( REFP_W (DIMENSION_3,3) )

      Allocate(  ONEoDX_E_U (DIMENSION_3) )
      Allocate(  ONEoDY_N_U (DIMENSION_3) )
      Allocate(  ONEoDZ_T_U (DIMENSION_3) )

      Allocate(  ONEoDX_E_V (DIMENSION_3) )
      Allocate(  ONEoDY_N_V (DIMENSION_3) )
      Allocate(  ONEoDZ_T_V (DIMENSION_3) )

      Allocate(  ONEoDX_E_W (DIMENSION_3) )
      Allocate(  ONEoDY_N_W (DIMENSION_3) )
      Allocate(  ONEoDZ_T_W (DIMENSION_3) )

      Allocate(  Xn_int (DIMENSION_3) )
      Allocate(  Xn_U_int (DIMENSION_3) )
      Allocate(  Xn_V_int (DIMENSION_3) )
      Allocate(  Xn_W_int (DIMENSION_3) )

      Allocate(  Ye_int (DIMENSION_3) )
      Allocate(  Ye_U_int (DIMENSION_3) )
      Allocate(  Ye_V_int (DIMENSION_3) )
      Allocate(  Ye_W_int (DIMENSION_3) )

      Allocate(  Zt_int (DIMENSION_3) )
      Allocate(  Zt_U_int (DIMENSION_3) )
      Allocate(  Zt_V_int (DIMENSION_3) )
      Allocate(  Zt_W_int (DIMENSION_3) )

      Allocate(  SNAP (DIMENSION_3) )

      SNAP = .FALSE.

      Allocate(  CUT_TREATMENT_AT (DIMENSION_3) )
      Allocate(  CUT_U_TREATMENT_AT (DIMENSION_3) )
      Allocate(  CUT_V_TREATMENT_AT (DIMENSION_3) )
      Allocate(  CUT_W_TREATMENT_AT (DIMENSION_3) )


      CUT_TREATMENT_AT = .FALSE.
      CUT_U_TREATMENT_AT = .FALSE.
      CUT_V_TREATMENT_AT = .FALSE.
      CUT_W_TREATMENT_AT = .FALSE.

      Allocate(  CUT_CELL_AT (DIMENSION_3) )
      Allocate(  CUT_U_CELL_AT (DIMENSION_3) )
      Allocate(  CUT_V_CELL_AT (DIMENSION_3) )
      Allocate(  CUT_W_CELL_AT (DIMENSION_3) )

      CUT_CELL_AT = .FALSE.
      CUT_U_CELL_AT = .FALSE.
      CUT_V_CELL_AT = .FALSE.
      CUT_W_CELL_AT = .FALSE.

      Allocate( SMALL_CELL_AT  (DIMENSION_3) )
      SMALL_CELL_AT = .FALSE.

      Allocate( SMALL_CELL_FLAG  (DIMENSION_3) )
      SMALL_CELL_FLAG = 0

      Allocate(  BLOCKED_CELL_AT (DIMENSION_3) )
      Allocate(  BLOCKED_U_CELL_AT (DIMENSION_3) )
      Allocate(  BLOCKED_V_CELL_AT (DIMENSION_3) )
      Allocate(  BLOCKED_W_CELL_AT (DIMENSION_3) )

      BLOCKED_CELL_AT   = .FALSE.
      BLOCKED_U_CELL_AT = .FALSE.
      BLOCKED_V_CELL_AT = .FALSE.
      BLOCKED_W_CELL_AT = .FALSE.

      Allocate(  STANDARD_CELL_AT (DIMENSION_3) )
      Allocate(  STANDARD_U_CELL_AT (DIMENSION_3) )
      Allocate(  STANDARD_V_CELL_AT (DIMENSION_3) )
      Allocate(  STANDARD_W_CELL_AT (DIMENSION_3) )


      Allocate(  VORTICITY (DIMENSION_3) )
      Allocate(  LAMBDA2 (DIMENSION_3) )

      Allocate(  TRD_G_OUT (DIMENSION_3) )
      Allocate(  PP_G_OUT (DIMENSION_3) )
      Allocate(  EPP_OUT (DIMENSION_3) )

      Allocate(  dudx_OUT (DIMENSION_3) )
      Allocate(  dvdy_OUT (DIMENSION_3) )
      Allocate(  delv_OUT (DIMENSION_3) )

      Allocate(  U_MASTER_OF (DIMENSION_3) )
      Allocate(  V_MASTER_OF (DIMENSION_3) )
      Allocate(  W_MASTER_OF (DIMENSION_3) )

      Allocate(  BC_ID (DIMENSION_3) )
      Allocate(  BC_U_ID (DIMENSION_3) )
      Allocate(  BC_V_ID (DIMENSION_3) )
      Allocate(  BC_W_ID (DIMENSION_3) )

      BC_ID   = 0
      BC_U_ID = 0
      BC_V_ID = 0
      BC_W_ID = 0

      Allocate(  DEBUG_CG (DIMENSION_3,15) )

      Allocate(  U_g_CC (DIMENSION_3) )
      Allocate(  V_g_CC (DIMENSION_3) )
      Allocate(  W_g_CC (DIMENSION_3) )

      Allocate(  U_s_CC (DIMENSION_3, DIMENSION_M) )
      Allocate(  V_s_CC (DIMENSION_3, DIMENSION_M) )
      Allocate(  W_s_CC (DIMENSION_3, DIMENSION_M) )

      ALLOCATE(N_FACET_AT(DIMENSION_3))
      N_FACET_AT = 0

      ALLOCATE(LIST_FACET_AT(DIMENSION_3,DIM_FACETS_PER_CELL))

      ALLOCATE(POTENTIAL_CUT_CELL_AT(DIMENSION_3))


      Allocate(  F_AT (DIMENSION_3) )

      Allocate(  DWALL (DIMENSION_3) )


      RETURN
      END SUBROUTINE ALLOCATE_CUT_CELL_ARRAYS


