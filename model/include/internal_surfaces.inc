!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: INTERNAL_SURFACES_LOCKED                                  !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / INTERNAL_SURFACES_LOCKED /&

      IS_TYPE,                                                         &
      IS_X_w, IS_X_e, IS_Y_s, IS_Y_n, IS_Z_b, IS_Z_t,                  &
      IS_I_w, IS_I_e, IS_J_s, IS_J_n, IS_K_b, IS_K_t


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: INTERNAL_SURFACES_UNLOCKED                                !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / INTERNAL_SURFACES_UNLOCKED /&

      IS_PC, IS_VEL_s
