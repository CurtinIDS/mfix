!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: TFM_SOLIDS_LOCKED                                         !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!     NAMELIST / TFM_SOLIDS_LOCKED /&


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: TFM_SOLIDS_UNLOCKED                                       !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / TFM_SOLIDS_UNLOCKED /&

      KT_TYPE, SIMONIN, AHMADI, FRICTION, SAVAGE, SCHAEFFER,  &
      BLENDING_STRESS, TANH_BLEND, SIGM_BLEND, YU_STANDISH,   &
      FEDORS_LANDEL, RDF_TYPE, ADDED_MASS, M_AM, E, C_e, R_p, &
      E_w, PHIP, PHIP0, C_f, PHI, PHI_w, EPS_f_min, EP_S_MAX, &
      SEGREGATION_SLOPE_COEFFICIENT, V_ex,  MU_s0, DIF_s0,    &
      EP_star, CLOSE_PACKED, M_MAX, JENKINS
