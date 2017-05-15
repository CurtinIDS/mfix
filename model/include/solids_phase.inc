!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: SOLIDS_PHASE_LOCKED                                       !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / SOLIDS_PHASE_LOCKED /&

      SOLIDS_MODEL, MMAX, NMAX_s, SPECIES_s, SPECIES_ALIAS_s

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: SOLIDS_PHASE_UNLOCKED                                     !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / SOLIDS_PHASE_UNLOCKED /&

      D_p0, RO_s0, X_s0, RO_Xs0, INERT_SPECIES, K_s0, C_ps0, MW_s, &
      DIL_FACTOR_VSD, DIL_INERT_X_VSD
