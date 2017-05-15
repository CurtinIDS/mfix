!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: OUTPUT_CONTROL_LOCKED                                     !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / OUTPUT_CONTROL_LOCKED /&

      SPX_DT, nRR, PHIP_OUT_JJ, bDist_IO , bStart_with_one_RES, &
      bWRITE_NetCDF, DBGPRN_LAYOUT, ENABLE_DMP_LOG, RES_BACKUPS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: OUTPUT_CONTROL_UNLOCKED                                   !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / OUTPUT_CONTROL_UNLOCKED /&

      RES_DT, OUT_DT, NLOG, FULL_LOG, RESID_STRING, GROUP_RESID, &
      REPORT_MASS_BALANCE_DT, REPORT_NEG_DENSITY, RES_BACKUP_DT, &
      REPORT_NEG_SPECIFICHEAT


