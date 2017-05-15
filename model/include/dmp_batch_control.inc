!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: DMP_BATCH_CONTROL_LOCKED                                  !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that CANNOT be edited during a run.     !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / DMP_BATCH_CONTROL_LOCKED /&

      NODESI , NODESJ , NODESK

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Namelist: DMP_BATCH_CONTROL_UNLOCKED                                !
!  Author: J.Musser                                   Date: 03-APR-15  !
!                                                                      !
!  Purpose: List MFIX keywords that may be edited during a run by      !
!  invoking a REFRESH.                                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      NAMELIST / DMP_BATCH_CONTROL_UNLOCKED /&

      CHK_BATCHQ_END, BATCH_WALLCLOCK, TERM_BUFFER, DEBUG_RESID, &
      IS_SERIAL, USE_DOLOOP, MINIMIZE_DOTPRODUCTS, ICHECK_BICGS, &
      SOLVER_STATISTICS, OPT_PARALLEL
