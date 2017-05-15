!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: POINT_SOURCES_LOCKED                                      !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!      NAMELIST / POINT_SOURCES_LOCKED /&


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: POINT_SOURCES_UNLOCKED                                    !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / POINT_SOURCES_UNLOCKED /&

      PS_X_w, PS_X_e, PS_Y_s, PS_Y_n, PS_Z_b, PS_Z_t,&
      PS_I_w, PS_I_e, PS_J_s, PS_J_n, PS_K_b, PS_K_t,&
      PS_U_g, PS_V_g, PS_W_g, PS_MASSFLOW_g, &
      PS_U_s, PS_V_s, PS_W_s, PS_MASSFLOW_s, &
      PS_T_g, PS_T_s, PS_X_g, PS_X_s
