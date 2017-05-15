!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: INITIAL_CONDITIONS_LOCKED                                 !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / INITIAL_CONDITIONS_LOCKED /&

      IC_DES_FIT_TO_REGION, IC_PIC_CONST_NPC, IC_PIC_CONST_STATWT



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: INITIAL_CONDITIONS_UNLOCKED                               !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / INITIAL_CONDITIONS_UNLOCKED /&

      IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t,&
      IC_I_w, IC_I_e, IC_J_s, IC_J_n, IC_K_b, IC_K_t,&
      IC_EP_g, IC_P_g, IC_ROP_s, IC_T_g, IC_T_s,&
      IC_U_g, IC_U_s, IC_V_g, IC_V_s, IC_W_g, IC_W_s,&
      IC_X_g, IC_X_s, IC_TYPE, IC_P_star, IC_L_scale,&
      IC_GAMA_Rg, IC_T_Rg, IC_GAMA_Rs, IC_T_Rs,&
      IC_Theta_m, IC_Scalar, IC_K_Turb_G, IC_E_Turb_G, IC_EP_s

