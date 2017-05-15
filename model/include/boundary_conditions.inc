!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: BOUNDARY_CONDITIONS_LOCKED                                !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / BOUNDARY_CONDITIONS_LOCKED /&

      BC_TYPE, &
      BC_X_w, BC_X_e, BC_Y_s, BC_Y_n, BC_Z_b, BC_Z_t, &
      BC_I_w, BC_I_e, BC_J_s, BC_J_n, BC_K_b, BC_K_t, &

      BC_PO_APPLY_TO_DES, BC_PIC_MI_CONST_NPC, BC_PIC_MI_CONST_STATWT, &
      U_g0, U_s0, V_g0, V_s0, W_g0, W_s0, BC_ROP_g, BC_RO_g


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: BOUNDARY_CONDITIONS_UNLOCKED                              !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / BOUNDARY_CONDITIONS_UNLOCKED /&

      BC_Hw_g, BC_Uw_g, BC_Vw_g, BC_Ww_g,                              &
      BC_Hw_s, BC_Uw_s, BC_Vw_s, BC_Ww_s, BC_JJ_PS,                    &
      BC_Hw_Theta_m, BC_Thetaw_m, BC_C_Theta_m,                        &
      BC_Hw_T_g, BC_Hw_T_s, BC_Tw_g, BC_Tw_s, BC_C_T_g, BC_C_T_s,      &
      BC_Hw_X_g, BC_Hw_X_s, BC_Xw_g, BC_Xw_s, BC_C_X_g, BC_C_X_s,      &
      BC_Hw_Scalar, BC_ScalarW, BC_C_Scalar,                           &

      BC_EP_g, BC_P_g, BC_ROP_s, BC_EP_s, BC_T_g, BC_T_s, BC_Theta_m,  &
      BC_X_g, BC_X_s,  BC_U_g, BC_U_s, BC_V_g, BC_V_s, BC_W_g, BC_W_s, &
      BC_VOLFLOW_g, BC_VOLFLOW_s, BC_MASSFLOW_g, BC_MASSFLOW_s,        &

      BC_DT_0, BC_Jet_g0, BC_DT_h, BC_Jet_gh, BC_DT_L, BC_Jet_GL,      &
      BC_Scalar, BC_K_Turb_G, BC_E_Turb_G, BC_VELMAG_G, BC_VELMAG_S,   &
      BC_JJ_M
