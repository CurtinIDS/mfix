! -*- f90 -*-
      MODULE STEP

      USE EXIT, ONLY: MFIX_EXIT
      USE MAIN, ONLY: NIT_TOTAL, DNCHECK, EXIT_SIGNAL, NCHECK
      USE RUN, ONLY: IER

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: TIME_STEP_INIT                                          !
!  Author: M. Syamlal                                 Date: 12-APR-96  !
!                                                                      !
!  Purpose: This module controls the iterations for solving equations  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE TIME_STEP_INIT

      !f2py threadsafe
      USE check, only: check_mass_balance
      USE compar, only: mype
      USE dashboard, only: run_status, write_dashboard
      USE discretelement, only: des_continuum_coupled, des_continuum_hybrid, discrete_element
      USE error_manager, only: err_msg
      USE error_manager, only: flush_err_msg
      USE output, only: res_dt
      USE output_man, only: output_manager
      USE param1, only: small_number, undefined
      USE qmom_kinetic_equation, only: qmomk
      USE run, only: auto_restart, automatic_restart, call_dqmom, call_usr, chk_batchq_end
      USE run, only: cn_on, dem_solids, dt, dt_min, dt_prev, ghd_2007, kt_type_enum
      USE run, only: nstep, nsteprst, odt, pic_solids, run_type, time, tstop, units, use_dt_prev
      USE stiff_chem, only: stiff_chemistry, stiff_chem_solver
      USE toleranc, only: max_inlet_vel
      USE utilities, only: max_vel_inlet

      IMPLICIT NONE

! Terminate MFIX normally before batch queue terminates.
      IF (CHK_BATCHQ_END) CALL CHECK_BATCH_QUEUE_END(EXIT_SIGNAL)

      IF (CALL_USR) CALL USR1

! Remove solids from cells containing very small quantities of solids
      IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
           DES_CONTINUUM_HYBRID) THEN
         IF(KT_TYPE_ENUM == GHD_2007) THEN
            CALL ADJUST_EPS_GHD
         ELSE
            CALL ADJUST_EPS
         ENDIF
      ENDIF

! sof modification: uncomment code below and modify MARK_PHASE_4_COR to
! use previous MFIX algorithm. Nov 22 2010.
! Mark the phase whose continuity will be solved and used to correct
! void/volume fraction in calc_vol_fr (see subroutine for details)
!      CALL MARK_PHASE_4_COR (PHASE_4_P_G, PHASE_4_P_S, DO_CONT, MCP,&
!           DO_P_S, SWITCH_4_P_G, SWITCH_4_P_S, IER)

! Set wall boundary conditions and transient flow b.c.'s
      CALL SET_BC1
! include here so they are set before calculations rely on value
! (e.g., calc_mu_g, calc_mu_s)
      CALL SET_EP_FACTORS


      CALL OUTPUT_MANAGER(EXIT_SIGNAL, .FALSE.)

! Update previous-time-step values of field variables
      CALL UPDATE_OLD

! Calculate coefficients
      CALL CALC_COEFF_ALL (0, IER)

! Calculate the stress tensor trace and cross terms for all phases.
      CALL CALC_TRD_AND_TAU()

! Calculate additional solids phase momentum source terms
      IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
         CALL CALC_EXPLICIT_MOM_SOURCE_S
      ENDIF

! Check rates and sums of mass fractions every NLOG time steps
      IF (NSTEP == NCHECK) THEN
         IF (DNCHECK < 256) DNCHECK = DNCHECK*2
         NCHECK = NCHECK + DNCHECK
! Upate the reaction rates for checking
         CALL CALC_RRATE(IER)
         CALL CHECK_DATA_30
      ENDIF

! Double the timestep for 2nd order accurate time implementation
      IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
           (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
         DT = 0.5d0*DT
         ODT = ODT * 2.0d0
      ENDIF

! Check for maximum velocity at inlet to avoid convergence problems
      MAX_INLET_VEL = MAX_VEL_INLET()

      END SUBROUTINE TIME_STEP_INIT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: TIME_STEP_INIT                                          !
!  Author: M. Syamlal                                 Date: 12-APR-96  !
!                                                                      !
!  Purpose: This module controls the iterations for solving equations  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE TIME_STEP_END

      USE check, only: check_mass_balance
      USE compar, only: mype
      USE dashboard, only: run_status, write_dashboard
      USE error_manager, only: err_msg
      USE error_manager, only: flush_err_msg
      USE iterate, only: nit
      USE leqsol, only: solver_statistics, report_solver_stats
      USE output, only: res_dt
      USE output_man, only: output_manager
      USE param1, only: small_number, undefined
      USE qmom_kinetic_equation, only: qmomk
      USE run, only: auto_restart, automatic_restart, call_dqmom, call_usr, chk_batchq_end
      USE run, only: cn_on, dem_solids, dt, dt_min, dt_prev, ghd_2007, kt_type_enum
      USE run, only: nstep, nsteprst, odt, pic_solids, run_type, time, tstop, units, use_dt_prev, steady_state
      USE stiff_chem, only: stiff_chemistry, stiff_chem_solver
      use discretelement, only: DES_EXPLICITLY_COUPLED
      IMPLICIT NONE

      IF(DT < DT_MIN) THEN

         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
1100 FORMAT('DT < DT_MIN.  Recovery not possible!')

         IF(WRITE_DASHBOARD) THEN
            RUN_STATUS = 'DT < DT_MIN.  Recovery not possible!'
            CALL UPDATE_DASHBOARD(NIT,0.0d0,'    ')
         ENDIF
         CALL MFIX_EXIT(MyPE)
      ENDIF

! Stiff Chemistry Solver.
      IF(STIFF_CHEMISTRY) THEN
         IF(DES_EXPLICITLY_COUPLED) CALL RXNS_GS_GAS1
         CALL STIFF_CHEM_SOLVER(DT, IER)
      ENDIF


! Check over mass and elemental balances.  This routine is not active by default.
! Edit the routine and specify a reporting interval to activate it.
      CALL CHECK_MASS_BALANCE (1)

! Other solids model implementations
      IF (DEM_SOLIDS) CALL DES_TIME_MARCH
      IF (PIC_SOLIDS) CALL PIC_TIME_MARCH
      IF (QMOMK) CALL QMOMK_TIME_MARCH
      IF (CALL_DQMOM) CALL USR_DQMOM

! Advance the time step and continue
      IF((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
           (CN_ON.AND.RUN_TYPE /= 'NEW' .AND. NSTEP >= (NSTEPRST+1))) THEN
! Double the timestep for 2nd order accurate time implementation
         DT = 2.d0*DT
         ODT = ODT * 0.5d0
! Perform the explicit extrapolation for CN implementation
         CALL CN_EXTRAPOL
      ENDIF

      IF (.NOT. STEADY_STATE) THEN
         IF(USE_DT_PREV) THEN
            TIME = TIME + DT_PREV
         ELSE
            TIME = TIME + DT
         ENDIF
         USE_DT_PREV = .FALSE.
         NSTEP = NSTEP + 1
      ENDIF

      NIT_TOTAL = NIT_TOTAL+NIT

      IF(SOLVER_STATISTICS) CALL REPORT_SOLVER_STATS(NIT_TOTAL, NSTEP)

! write (*,"('Compute the Courant number')")
! call get_stats(IER)

      FLUSH (6)

      CALL OUTPUT_MANAGER(EXIT_SIGNAL, .TRUE.)

      END SUBROUTINE TIME_STEP_END

      END MODULE STEP
