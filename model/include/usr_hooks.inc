!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: USER_HOOKS_LOCKED                                         !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!     NAMELIST / USER_HOOKS_LOCKED /&


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: USER_HOOKS_UNLOCKED                                       !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / USER_HOOKS_UNLOCKED /&

      CALL_USR, C, C_NAME, USR_DT, USR_TYPE, USR_VAR, USR_FORMAT,      &
      USR_EXT, USR_X_w, USR_X_e, USR_Y_s, USR_Y_n, USR_Z_b, USR_Z_t,   &
      USR_ROG, USR_CPG, USR_KG, USR_MUG, USR_DIFG,                     &
      USR_ROS, USR_CPS, USR_KS, USR_MUS, USR_DIFS,                     &
      USR_FGS, USR_FSS, USR_GAMA, CALL_USR_SOURCE
