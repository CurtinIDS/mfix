!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: NUMERICAL_PARAM_LOCKED                                    !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / NUMERICAL_PARAM_LOCKED /&

      FPFOI

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: NUMERICAL_PARAM_UNLOCKED                                  !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / NUMERICAL_PARAM_UNLOCKED /&

      MAX_NIT, DETECT_STALL, NORM_G, NORM_S, MAX_INLET_VEL_FAC, &
      TOL_RESID, TOL_RESID_T, TOL_RESID_X, TOL_RESID_TH, &
      TOL_RESID_SCALAR, TOL_DIVERGE, &
      LEQ_IT, LEQ_METHOD, LEQ_SWEEP, LEQ_TOL, LEQ_PC, DO_TRANSPOSE, &
      UR_FAC, UR_F_GS, UR_KTH_SML, &
      DISCRETIZE, DEF_COR, C_FAC, CHI_SCHEME, CN_ON

