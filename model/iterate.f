! -*- f90 -*-
      MODULE ITERATE

      USE EXIT, ONLY: MFIX_EXIT
      USE RUN, ONLY: IER, TUNIT, TIME

! flag indicating convergence status with MUSTIT = 0,1,2 implying
! complete convergence, non-covergence and divergence respectively
      INTEGER :: MUSTIT

! Number of iterations completed for current timestep
      INTEGER :: NIT

! User defined maximum number of iterations
      INTEGER :: MAX_NIT

      LOGICAL :: CONVERGED, DIVERGED

! cpu time left
      DOUBLE PRECISION :: TLEFT
! Normalization factor for gas & solids pressure residual
      DOUBLE PRECISION :: NORMg, NORMs
! Set normalization factor for gas and solids pressure residual
      LOGICAL :: SETg, SETs
! gas & solids pressure residual
      DOUBLE PRECISION :: RESg, RESs
! Weight of solids in the reactor
      DOUBLE PRECISION :: SMASS
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: errorpercent
! Error Message
      CHARACTER(LEN=32) :: lMsg

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ITERATE                                                 !
!  Author: M. Syamlal                                 Date: 12-APR-96  !
!                                                                      !
!  Purpose: This module controls the iterations for solving equations  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ITERATE_INIT

      USE compar, only: mype, pe_io
      USE cutcell, only: cartesian_grid
      USE geometry, only: cyclic
      USE leqsol, only: leq_adjust
      USE output, only: full_log
      USE param1, only: one, small_number, undefined, zero
      USE run, only: dt, dt_prev, run_type, time, tstop, nstep, nsteprst, cn_on, get_tunit, steady_state
      USE time_cpu
      USE toleranc, only: norm_g, norm_s

      IMPLICIT NONE

! initializations
      DT_prev = DT
      NIT = 0
      MUSTIT = 0
      CONVERGED = .FALSE.
      DIVERGED  = .FALSE.
      RESG = ZERO
      RESS = ZERO

      IF(NORM_G == UNDEFINED) THEN
         NORMG = ONE
         SETG = .FALSE.
      ELSE
         NORMG = NORM_G
         SETG = .TRUE.
      ENDIF

      IF(NORM_S == UNDEFINED) THEN
         NORMS = ONE
         SETS = .FALSE.
      ELSE
         NORMS = NORM_S
         SETS = .TRUE.
      ENDIF

      LEQ_ADJUST = .FALSE.

! Initialize residuals
      CALL INIT_RESID ()

! Initialize the routine for holding gas mass flux constant with cyclic bc
      IF(CYCLIC) CALL GoalSeekMassFlux(NIT, MUSTIT, .false.)

! CPU time left
      IF (FULL_LOG) THEN
         TLEFT = (TSTOP - TIME)*CPUOS
         CALL GET_TUNIT (TLEFT, TUNIT)

         IF (STEADY_STATE) THEN
            CALL GET_SMASS (SMASS)
            IF(myPE.eq.PE_IO) THEN
               WRITE (*, '(/A,G12.5, A,F9.3,1X,A)') &
                  ' Starting solids mass = ', SMASS
            ENDIF
         ELSE
            IF(myPE.eq.PE_IO) THEN
               IF ((CN_ON.AND.NSTEP>1.AND.RUN_TYPE == 'NEW') .OR. &
                   (CN_ON.AND.RUN_TYPE /= 'NEW' .AND.&
                    NSTEP >= (NSTEPRST+1))) THEN
                  WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)')&
                     ' Time = ', TIME, '  Dt = ', 2.D0*DT
               ELSE
                  WRITE (*, '(/A,G12.5, A,G12.5, A,F9.3,1X,A)') &
                     ' Time = ', TIME, '  Dt = ', DT
               ENDIF
            ENDIF
         ENDIF   ! if/else(steady_state)
      ENDIF   ! if(full_log)

      CALL CALC_RESID_MB(0, errorpercent)

! Calculate the face values of densities and mass fluxes for the first
! solve_vel_star call.
      CALL CONV_ROP()
      CALL CALC_MFLUX ()
      CALL SET_BC1
      CALL SET_EP_FACTORS

! JFD: modification for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SET_OUTFLOW

! Default/Generic Error message
      lMsg = 'Run diverged/stalled'

      RETURN
      END SUBROUTINE ITERATE_INIT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ITERATE                                                 !
!  Author: M. Syamlal                                 Date: 12-APR-96  !
!                                                                      !
!  Purpose: This module controls the iterations for solving equations  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DO_ITERATION

      USE cont, only: solve_continuity
         USE cutcell, only: cartesian_grid
         USE discretelement, only: discrete_element, des_continuum_hybrid
         USE fldvar, only: ep_g, ro_g, rop_g, rop_s, p_star
         USE geometry, only: cyclic, cylindrical
         USE leqsol, only: leq_adjust
         USE mms, only: USE_MMS
         USE param1, only: small_number, zero, undefined, undefined_i
         USE physprop, only: mmax, ro_g0, smax
         USE pscor, only: k_cp, mcp
         USE qmom_kinetic_equation, only: qmomk
         USE residual, only: resid_p, resid
         USE run, only: auto_restart, automatic_restart, call_dqmom, call_usr, chk_batchq_end, friction, ghd_2007, granular_energy, k_epsilon, kt_type_enum, phip_out_iter, energy_eq, dt, ier, steady_state
         USE scalars, only: nscalar

      IMPLICIT NONE

      INTEGER :: M

      PHIP_OUT_ITER=NIT ! To record the output of phip
! mechanism to set the normalization factor for the correction
! after the first iteration to the corresponding residual found
! in the first iteration
      IF (.NOT.SETG) THEN
         IF (RESG > SMALL_NUMBER) THEN
            NORMG = RESG
            SETG = .TRUE.
         ENDIF
      ENDIF
      IF (.NOT.SETS) THEN
         IF (RESS > SMALL_NUMBER) THEN
            NORMS = RESS
            SETS = .TRUE.
         ENDIF
      ENDIF

! Call user-defined subroutine to set quantities that need to be updated
! every iteration
      IF (CALL_USR) CALL USR2

! Calculate coefficients, excluding density and reactions.
      CALL CALC_COEFF(IER, 1)
      IF (IER_MANAGER()) return

! Diffusion coefficient and source terms for user-defined scalars
      IF(NScalar /= 0) CALL SCALAR_PROP()

! Diffusion coefficient and source terms for K & Epsilon Eq.
      IF(K_Epsilon) CALL K_Epsilon_PROP()

! Update the stress tensor trace and cross terms each subiteration
! for MMS cases.
      IF(USE_MMS) CALL CALC_TRD_AND_TAU()

! Solve starred velocity components
      CALL SOLVE_VEL_STAR(IER)

! Correct the centerline velocity for cylindrical simulations.
      IF(CYLINDRICAL) CALL RADIAL_VEL_CORRECTION

! Calculate densities.
      CALL PHYSICAL_PROP(IER, 0)
      IF (IER_MANAGER()) return

! Calculate chemical reactions.
      CALL CALC_RRATE(IER)

! Solve solids volume fraction correction equation for close-packed
! solids phases
      IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
         DES_CONTINUUM_HYBRID) THEN
         IF (MMAX > 0) THEN
! MMS:  Solve gas continuity only.
            IF(USE_MMS) THEN
               CALL SOLVE_CONTINUITY(0,IER)
! Regular, non-MMS cases.
            ELSE
               IF(MMAX == 1 .AND. MCP /= UNDEFINED_I)THEN
! if second phase (m=1) can overpack (e.g., bubbles) then solve its
! continuity equation
                  CALL CALC_K_CP (K_CP)
                  CALL SOLVE_EPP (NORMS, RESS, IER)
                  CALL CORRECT_1 ()
               ELSE

! If one chooses to revert back to old mark_phase_4_cor wherein the
! continuity of the gas phase can get marked to be solved then this
! loop should start at 0.
                  DO M=1,SMAX ! mmax -> smax for GHD theory
! Volume fraction correction technique for one of the solids phase
! is not implemented.  This will only slow down convergence.
!                      IF (M .EQ. MCP) THEN
!                       CALL CALC_K_CP (K_CP, IER)
!                       CALL SOLVE_EPP (NORMS, RESS, IER)
!                       CALL CORRECT_1 (IER)
!                    ELSE
                        CALL SOLVE_CONTINUITY(M,IER)
!                    ENDIF
                  ENDDO
               ENDIF   ! end if/else (mmax==1 .and. mcp /= undefined)
            ENDIF ! end if/else (MMS)

            IF(KT_TYPE_ENUM == GHD_2007) CALL ADJUST_EPS_GHD

            CALL CALC_VOL_FR (P_STAR, RO_G, ROP_G, EP_G, ROP_S, IER)
            IF (IER_MANAGER()) return

         ENDIF  ! endif (mmax >0)

      ENDIF  ! end if (.not.discrete_element)


! Calculate P_star in cells where solids continuity equation is
! solved
      IF(.NOT.(DISCRETE_ELEMENT .OR. QMOMK) .OR. &
         DES_CONTINUUM_HYBRID) THEN
         IF (MMAX > 0 .AND. .NOT.FRICTION) &
            CALL CALC_P_STAR (EP_G, P_STAR)
      ENDIF

! Calculate the face values of densities.
      CALL CONV_ROP()

      IF (RO_G0 /= ZERO) THEN
! Solve fluid pressure correction equation
#ifndef FLAG_MMS
         CALL SOLVE_PP_G (NORMG, RESG, IER)
#endif
! Correct pressure, velocities, and density
         CALL CORRECT_0 ()
      ENDIF

! Recalculate densities.
      CALL PHYSICAL_PROP(IER, 0)
      IF (IER_MANAGER()) return

! Update wall velocities:
! modified by sof to force wall functions so even when NSW or FSW are
! declared, default wall BC will still be treated as NSW and no wall
! functions will be used
      IF(.NOT. K_EPSILON) CALL SET_WALL_BC ()

! Calculate the face values of mass fluxes
      CALL CALC_MFLUX ()
      CALL SET_BC1
      CALL SET_EP_FACTORS

! JFD: modification for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SET_OUTFLOW

! Solve energy equations
      IF (ENERGY_EQ) THEN
         CALL SOLVE_ENERGY_EQ (IER)
         IF (IER_MANAGER()) return
      ENDIF

! Solve granular energy equation
      IF (GRANULAR_ENERGY) THEN
         IF(.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
            CALL SOLVE_GRANULAR_ENERGY (IER)
            IF (IER_MANAGER()) return
         ENDIF
      ENDIF

! Solve species mass balance equations.
      CALL SOLVE_SPECIES_EQ (IER)
      IF (IER_MANAGER()) return

! Solve other scalar transport equations
      IF(NScalar /= 0) CALL SOLVE_Scalar_EQ (IER)

! Solve K & Epsilon transport equations
      IF(K_Epsilon) CALL SOLVE_K_Epsilon_EQ (IER)

! User-defined linear equation solver parameters may be adjusted after
! the first iteration
      IF (.NOT.CYCLIC) LEQ_ADJUST = .TRUE.

! Check for convergence
      CALL ACCUM_RESID ! Accumulating residuals from all the processors
      RESG = RESID(RESID_P,0)
      RESS = RESID(RESID_P,1)
      CALL CALC_RESID_MB(1, errorpercent)
      MUSTIT = 0
      CALL CHECK_CONVERGENCE (NIT, errorpercent(0), MUSTIT)

      IF(CYCLIC .AND. (MUSTIT==0 .OR. NIT >= MAX_NIT)) &
         CALL GoalSeekMassFlux(NIT, MUSTIT, .true.)

! Display residuals
      CALL DISPLAY_RESID (NIT)

      CALL END_ITERATION

      contains

      SUBROUTINE END_ITERATION
         IMPLICIT NONE

         ! Determine course of simulation: converge, non-converge, diverge?
         IF (MUSTIT == 0) THEN
            IF (STEADY_STATE .AND. NIT==1) RETURN   !Iterations converged
            CONVERGED = .TRUE.
            IER = 0
         ELSEIF (MUSTIT==2 .AND. .NOT.STEADY_STATE) THEN
            DIVERGED = .TRUE.
            IER = 1
         ENDIF

      END SUBROUTINE END_ITERATION

!----------------------------------------------------------------------!
! Function: IER_Manager                                                !
!                                                                      !
! Purpose: Identify and account for errors from called subroutines.    !
!          Returns .TRUE. for lErr >= 100, otherwise .FALSE.           !
!                                                                      !
! Reserved Error Blocks:                                               !
!                                                                      !
! [ 100,  109]: PHYSICAL_PROP                                          !
! [ 110,  119]: CALC_VOL_FR                                            !
! [ 120,  129]: SOLVE_ENERGY_EQ                                        !
! [ 130,  139]: SOLVE_SPECIES_EQ                                       !
! [ 140,  149]: SOLVE_GRANULAR_ENERGY                                  !
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION IER_MANAGER()

! Default case: do nothing.
      IF(IER < 100) THEN
         IER_MANAGER = .FALSE.
         return
      ENDIF

! Errors with an index greater than 100 will force an exit from iterate
! and in turn, reduce the step-size, and restart the time-step.
      IER_MANAGER = .TRUE.
      MUSTIT = 2

! Errors reported from PHYSICAL_PROP
!```````````````````````````````````````````````````````````````````````
      IF(IER <  110) THEN
         IF(IER ==  100) THEN
            lMsg = 'Negative gas density detected'
         ELSEIF(IER ==  101) THEN
            lMsg = 'Negative solids density detected'
         ELSE
            lMsg = 'UCE in PHYSICAL_PROP'
         ENDIF


! Errors reported from CALC_VOL_FR
!```````````````````````````````````````````````````````````````````````
      ELSEIF(IER <  120) THEN
         IF(IER ==  110) THEN
            lMsg = 'Negative void fraction detected'
         ELSE
            lMsg = 'UCE in CALC_VOL_FR'
         ENDIF


! Errors reported from SOLVE_ENERGY_EQ
!```````````````````````````````````````````````````````````````````````
      ELSEIF(IER <  130) THEN
         IF(IER ==  120) THEN
            lMsg = 'Energy Equation diverged'
         ELSE
            lMsg = 'UCE in SOLVE_ENERGY_EQ'
         ENDIF


! Errors reported from SOLVE_SPECIES_EQ
!```````````````````````````````````````````````````````````````````````
      ELSEIF(IER <  140) THEN
         IF(IER ==  130) THEN
            lMsg = 'Species Equation diverged'
         ELSE
            lMsg = 'UCE in SOLVE_SPECIES_EQ'
         ENDIF


! Errors reported from SOLVE_GRANULAR_ENERGY
!```````````````````````````````````````````````````````````````````````
      ELSEIF(IER <  150) THEN
         IF(IER ==  140) THEN
            lMsg = 'Granular Energy Eq diverged'
         ELSE
            lMsg = 'UCE in SOLVE_GRANULAR_ENERGY'
         ENDIF

! Unclassified Errors
!```````````````````````````````````````````````````````````````````````
      ELSE
         lMsg = 'Run diverged/stalled with UCE'
      ENDIF


      IF(STEADY_STATE) IER_MANAGER = .FALSE.

      CALL END_ITERATION

      return
      END FUNCTION IER_MANAGER

      END SUBROUTINE DO_ITERATION


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ITERATE                                                 !
!  Author: M. Syamlal                                 Date: 12-APR-96  !
!                                                                      !
!  Purpose: This module controls the iterations for solving equations  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE POST_ITERATE

      USE compar, only: mype, pe_io
      USE funits, only: dmp_log, unit_log, unit_out
      USE machine, only: start_log, end_log
      USE run, only: dt, time

      IMPLICIT NONE

      IF (CONVERGED) THEN
         CALL LOG_CONVERGED
      ELSEIF (DIVERGED) THEN
         CALL LOG_DIVERGED
      ELSE
         CALL GET_SMASS (SMASS)
         IF (myPE.eq.PE_IO) WRITE(UNIT_OUT, 5100) TIME, DT, NIT, SMASS
         CALL START_LOG
         IF(DMP_LOG) WRITE(UNIT_LOG, 5100) TIME, DT, NIT, SMASS
         CALL END_LOG
      ENDIF

      IF(.NOT.(CONVERGED .OR. DIVERGED)) THEN
         IER = 1
      ENDIF

5100  FORMAT(1X,'t=',F11.4,' Dt=',G11.4,' NIT>',I3,' Sm= ',G12.5, &
           'MbErr%=', G11.4)

      END SUBROUTINE POST_ITERATE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ITERATE                                                 !
!  Author: M. Syamlal                                 Date: 12-APR-96  !
!                                                                      !
!  Purpose: This module controls the iterations for solving equations  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE LOG_DIVERGED
      USE compar, only: mype, pe_io
      USE dashboard, only: f_dashboard, n_dashboard, run_status, write_dashboard
      USE funits, only: dmp_log, unit_log
      USE machine, only: start_log, end_log
      USE output, only: full_log
      USE run, only: dt, time, tstop, get_tunit
      USE time_cpu, only: cpuos

      IMPLICIT NONE

      CHARACTER(LEN=4) :: TUNIT

      IF (FULL_LOG) THEN
         CALL START_LOG
         CALL CALC_RESID_MB(1, errorpercent)

         IF(DMP_LOG) WRITE(UNIT_LOG,5200) TIME, DT, NIT, &
              errorpercent(0), trim(adjustl(lMsg))
         CALL END_LOG

         IF (myPE.EQ.PE_IO) WRITE(*,5200) TIME, DT, NIT, &
              errorpercent(0), trim(adjustl(lMsg))
      ENDIF

      ! JFD: modification for cartesian grid implementation
      IF(WRITE_DASHBOARD) THEN
         RUN_STATUS = 'Diverged/stalled...'
         N_DASHBOARD = N_DASHBOARD + 1
         IF(MOD(N_DASHBOARD,F_DASHBOARD)==0) THEN
            TLEFT = (TSTOP - TIME)*CPUOS
            CALL GET_TUNIT (TLEFT, TUNIT)
            CALL UPDATE_DASHBOARD(NIT,TLEFT,TUNIT)
         ENDIF
      ENDIF
5200  FORMAT(1X,'t=',F11.4,' Dt=',G11.4,' NIT=',&
                 I3,'MbErr%=', G11.4, ': ',A,' :-(')
      END SUBROUTINE LOG_DIVERGED

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: ITERATE                                                 !
!  Author: M. Syamlal                                 Date: 12-APR-96  !
!                                                                      !
!  Purpose: This module controls the iterations for solving equations  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE LOG_CONVERGED

      USE compar, only: mype, pe_io
      USE dashboard, only: f_dashboard, n_dashboard, run_status, write_dashboard
      USE error_manager, only: err_msg
      USE funits, only: dmp_log, unit_log
      USE geometry, only: cyclic_x, cyclic_y, cyclic_z
      USE geometry, only: do_i, do_j, do_k
      USE machine, only: start_log, end_log
      USE output, only: full_log, nlog
      USE physprop, only: mmax, smax
      USE run, only: dt, energy_eq, time, tstop, nstep, get_tunit
      USE time_cpu, only: cpu0, cpu_nlog, cpuos, time_nlog

      IMPLICIT NONE

      CHARACTER(LEN=4) :: TUNIT
! Perform checks and dump to screen every NLOG time steps
      IF (MOD(NSTEP,NLOG) == 0) CALL DUMP_TO_SCREEN

! JFD: modification for cartesian grid implementation
      IF(WRITE_DASHBOARD) THEN
         RUN_STATUS = 'In Progress...'
         N_DASHBOARD = N_DASHBOARD + 1
         IF(MOD(N_DASHBOARD,F_DASHBOARD)==0) THEN
            TLEFT = (TSTOP - TIME)*CPUOS
            CALL GET_TUNIT (TLEFT, TUNIT)
            CALL UPDATE_DASHBOARD(NIT,TLEFT,TUNIT)
         ENDIF
      ENDIF

      CONTAINS

      SUBROUTINE DUMP_TO_SCREEN
      USE error_manager, only: flush_err_msg
      IMPLICIT NONE

      ! phase index
      INTEGER :: M
! current cpu time used
      DOUBLE PRECISION :: CPU_NOW
! Heat loss from the reactor
      DOUBLE PRECISION :: HLOSS
! average velocity
      DOUBLE PRECISION :: Vavg

!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: VAVG_U_G, VAVG_V_G, VAVG_W_G, &
         VAVG_U_S, VAVG_V_S, VAVG_W_S

      CALL CPU_TIME (CPU_NOW)
      CPUOS = (CPU_NOW - CPU_NLOG)/(TIME - TIME_NLOG)
      CPU_NLOG = CPU_NOW
      TIME_NLOG = TIME
      CPU_NOW = CPU_NOW - CPU0

      CALL CALC_RESID_MB(1, errorpercent)
      CALL GET_SMASS (SMASS)
      IF (ENERGY_EQ) CALL GET_HLOSS (HLOSS)

      CALL START_LOG
      IF (ENERGY_EQ) THEN
         WRITE(ERR_MSG,5000)TIME, DT, NIT, SMASS, HLOSS, CPU_NOW
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ELSE
         WRITE(ERR_MSG,5001) TIME, DT, NIT, SMASS, CPU_NOW
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 5000 FORMAT(1X,'t=',F11.4,' Dt=',G11.4,' NIT=',I3,' Sm=',G12.5, &
         ' Hl=',G12.5,T84,'CPU=',F8.0,' s')

 5001 FORMAT(1X,'t=',F11.4,' Dt=',G11.4,' NIT=',I3,' Sm=',G12.5, &
         T84,'CPU=',F8.0,' s')

      IF(DMP_LOG)WRITE (UNIT_LOG, 5002) (errorpercent(M), M=0,MMAX)
      IF (FULL_LOG .and. myPE.eq.PE_IO) &
         WRITE (*, 5002) (errorpercent(M), M=0,MMAX)

 5002 FORMAT(3X,'MbError%(0,MMAX):', 5(1X,G11.4))

      IF (.NOT.FULL_LOG) THEN
         TLEFT = (TSTOP - TIME)*CPUOS
         CALL GET_TUNIT (TLEFT, TUNIT)
         IF(DMP_LOG)WRITE (UNIT_LOG, '(46X,A,F9.3,1X,A)')
      ENDIF

      IF (CYCLIC_X .OR. CYCLIC_Y .OR. CYCLIC_Z) THEN
         IF (DO_I) THEN
           Vavg = VAVG_U_G()
           IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'U_g = ', Vavg
         ENDIF
         IF (DO_J) THEN
           Vavg = VAVG_V_G()
           IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'V_g = ',  Vavg
         ENDIF
         IF (DO_K) THEN
           Vavg = VAVG_W_G()
           IF(DMP_LOG)WRITE (UNIT_LOG, 5050) 'W_g = ', Vavg
         ENDIF
         DO M = 1, SMAX
            IF (DO_I) Then
              Vavg = VAVG_U_S(M)
              IF(DMP_LOG)WRITE (UNIT_LOG, 5060) 'U_s(', M, ') = ', Vavg
            ENDIF
            IF (DO_J) Then
              Vavg = VAVG_V_S(M)
              IF(DMP_LOG)WRITE (UNIT_LOG, 5060) 'V_s(', M, ') = ', Vavg
            ENDIF
            IF (DO_K) Then
              Vavg = VAVG_W_S(M)
              IF(DMP_LOG)WRITE (UNIT_LOG, 5060) 'W_s(', M, ') = ', Vavg
            ENDIF
         ENDDO
      ENDIF   ! end if cyclic_x, cyclic_y or cyclic_z

      CALL END_LOG

5050  FORMAT(5X,'Average ',A,G12.5)
5060  FORMAT(5X,'Average ',A,I2,A,G12.5)

      END SUBROUTINE DUMP_TO_SCREEN

      END SUBROUTINE LOG_CONVERGED

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GoalSeekMassFlux                                        !
!  Author: M. Syamlal                                 Date: 12-APR-96  !
!                                                                      !
! Purpose:  In the following subroutine the mass flux across a periodic!
! domain with pressure drop is held constant at a user-specified value.!
! This module is activated only if the user specifies a value for the  !
! keyword flux_g in the mfix.dat file.                                 !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GoalSeekMassFlux(NIT, MUSTIT, doit)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE compar
      USE constant
      USE exit, only: mfix_exit
      USE geometry
      USE param1, only: one
      USE run
      USE time_cpu
      USE utilities, ONLY: mfix_isnan
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: NIT, MUSTIT
      LOGICAL, INTENT(IN) :: doit
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER, PARAMETER :: MAXOUTIT = 500
      DOUBLE PRECISION, PARAMETER          :: omega = 0.9
      DOUBLE PRECISION, PARAMETER          :: TOL = 1E-03
      INTEGER, SAVE :: OUTIT
      LOGICAL, SAVE :: firstPass = .true.

      DOUBLE PRECISION, SAVE  :: mdot_n, mdot_nm1, delp_n, delp_nm1, err
      DOUBLE PRECISION        :: mdot_0, delp_xyz

!-----------------------------------------------
! Functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: VAVG_Flux_U_G, VAVG_Flux_V_G, &
                                    VAVG_Flux_W_G
!-----------------------------------------------

      IF(CYCLIC_X_MF)THEN
         delp_n = delp_x
      ELSEIF(CYCLIC_Y_MF)THEN
         delp_n = delp_y
      ELSEIF(CYCLIC_Z_MF)THEN
         delp_n = delp_z
      ELSE
         RETURN
      ENDIF

      IF(.NOT.doit) THEN
         OUTIT = 0
         RETURN
      ENDIF

      OUTIT = OUTIT + 1
      IF(OUTIT > MAXOUTIT) THEN
         IF (myPE.EQ.PE_IO) write(*,5400) MAXOUTIT
         CALL mfix_exit(0)
      ENDIF

      mdot_0 = Flux_g


      ! calculate the average gas mass flux and error
      IF(CYCLIC_X_MF)THEN
        mdot_n = VAVG_Flux_U_G()
      ELSEIF(CYCLIC_Y_MF)THEN
        mdot_n = VAVG_Flux_V_G()
      ELSEIF(CYCLIC_Z_MF)THEN
        mdot_n = VAVG_Flux_W_G()
      ENDIF

      IF (mfix_isnan(mdot_n) .OR. mfix_isnan(delp_n)) THEN
         IF (myPE.eq.PE_IO) write(*,*) mdot_n, delp_n, &
            ' NaN being caught in GoalSeekMassFlux '
         RETURN
      ENDIF

      err = abs((mdot_n - mdot_0)/mdot_0)
      IF( err < TOL) THEN
         MUSTIT = 0
      ELSE
        MUSTIT = 1
        NIT = 1
      ENDIF

! correct delp
      if(.not.firstPass)then
!        delp_xyz = delp_n - omega * (delp_n - delp_nm1) * (mdot_n - mdot_0) &
!                          / (mdot_n - mdot_nm1)
! Fail-Safe Newton's method (below) works better than the regular
! Newton method (above)

         delp_xyz = delp_n - omega * (delp_n - delp_nm1) * &
                     ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0)) / &
                     ((mdot_n - mdot_0)/(mdot_nm1 - mdot_0) - ONE)
      else
         firstPass=.false.
         delp_xyz = delp_n*0.99
      endif

      IF(MUSTIT == 0) then
        IF(myPE.eq.PE_IO) Write(*,5500) TIME, OUTIT, delp_xyz, mdot_n
      ENDIF

      mdot_nm1 = mdot_n
      delp_nm1 = delp_n

      IF(CYCLIC_X_MF)THEN
        delp_x = delp_xyz
      ELSEIF(CYCLIC_Y_MF)THEN
        delp_y = delp_xyz
      ELSEIF(CYCLIC_Z_MF)THEN
        delp_z = delp_xyz
      ENDIF


      RETURN

5400 FORMAT(/1X,70('*')//' From: GoalSeekMassFlux',/&
      ' Message: Number of outer iterations exceeded ', I4,/1X,70('*')/)
5500  Format('  Time=', G12.5, ' MassFluxIterations=', I4, ' DelP=', &
      G12.5, ' Gas Flux=', G12.5)

      END SUBROUTINE GoalSeekMassFlux

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ADJUST_DT()                                            !
!  Author: M. Syamlal                                 Date: FEB-10-97  !
!                                                                      !
!  Purpose: Automatically adjust time step.                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      LOGICAL FUNCTION ADJUSTDT()

! Global Variables:
!---------------------------------------------------------------------//
! User defined type of run: new or restart
      use run, only: RUN_TYPE
! Integer flag: 0=Good, 100=initialize, otherwise bad.
      use run, only: IER
! User defined: min, max DT and adjustment factor
      use run, only: DT_MIN, DT_MAX, DT_FAC
! Flag: Use stored DT value for advancing TIME
      use run, only: USE_DT_PREV
! Flag: 2nd order time implementation
      use run, only: CN_ON
! Flag: Continue to run at DT_MIN
      use run, only: PERSISTENT_MODE
! The current number of time steps (value at restart).
      use run, only: NSTEP, NSTEPRST
! Current DT (1/DT) and direction of last change (+/-)
      use run, only: DT, oDT, DT_DIR, STEADY_STATE

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, ONE, UNDEFINED

! Module procedures:
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments:

! Local Variables:
!---------------------------------------------------------------------//
! Number of steps in between DT adjustments.
      INTEGER, PARAMETER :: STEPS_MIN = 5
! Number of time steps since last DT adjustment
      INTEGER, SAVE :: STEPS_TOT=0
! number of iterations since last DT adjustment
      INTEGER, SAVE :: NIT_TOT=0
! Iterations per second for last dt
      DOUBLE PRECISION, SAVE :: NIToS=0.0
! Current number of iterations per second
      DOUBLE PRECISION :: NITOS_NEW
! Flag to half/double the current time step
      LOGICAL :: CN_ADJUST_DT
!......................................................................!

! Initialize the function result.
      ADJUSTDT = .FALSE.
      USE_DT_PREV = .FALSE.

! Steady-state simulation.
      IF (STEADY_STATE .OR. DT<ZERO) RETURN

! Local flag for adjusting the time step for CN implementation.
      CN_ADJUST_DT = CN_ON .AND. ((RUN_TYPE=='NEW' .AND. NSTEP>1) .OR. &
         (RUN_TYPE/='NEW' .AND. NSTEP >= (NSTEPRST+1)))

! Iterate successfully converged.
!---------------------------------------------------------------------//
      IF(IER == 0) THEN

! Set back the timestep to original size which was halved previously for
! 2nd order accurate time implementation.
         IF(CN_ADJUST_DT) DT = 2.0D0*DT

! Calculate a new DT every STEPS_MIN time steps.
         IF(STEPS_TOT >= STEPS_MIN) THEN
            NITOS_NEW = DBLE(NIT_TOT)/(STEPS_TOT*DT)
            IF (NITOS_NEW > NITOS) DT_DIR = DT_DIR*(-1)
            STEPS_TOT = 0
            NITOS = NITOS_NEW
            NIT_TOT = 0
            IF (DT_DIR > 0) THEN
               IF(NIT < MAX_NIT) DT = MIN(DT_MAX,DT/DT_FAC)
            ELSE
               DT = DT*DT_FAC
               IF(PERSISTENT_MODE) DT = max(DT, DT_MIN)
            ENDIF

! DT was modified. Use the stored DT should be used to update TIME.
            USE_DT_PREV = .TRUE.

! Write the convergence stats to the screen/log file.
            WRITE(ERR_MSG,"('DT=',g11.4,3x,'NIT/s=',A)")  &
               DT, trim(iVal(nint(NITOS)))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., &
               FOOTER=.FALSE., LOG=.FALSE.)

         ELSE
            STEPS_TOT = STEPS_TOT + 1
            NIT_TOT = NIT_TOT + NIT
         ENDIF
! No need to iterate again
         ADJUSTDT = .FALSE.
! Cut the timestep into half for 2nd order accurate time implementation.
         IF(CN_ADJUST_DT) DT = 0.5d0*DT

! Iterate failed to converge.
!---------------------------------------------------------------------//
      ELSE

! Clear the error flag.
         IER = 0

! Reset the timestep to original size which was halved previously for
! 2nd order accurate time implementation.
         IF(CN_ADJUST_DT) DT = 2.0d0*DT

! Reset counters.
         DT_DIR = -1
         STEPS_TOT = 0
         NITOS = 0.
         NIT_TOT = 0

! Reduce the step size.
         DT = DT*DT_FAC

! The step size has decreased to the minimum.
         IF (DT_FAC >= ONE) THEN

            IF(PERSISTENT_MODE) THEN
               ADJUSTDT = .FALSE.
            ELSE
               WRITE(ERR_MSG,"(3X,A)") &
                  'DT_FAC >= 1. Recovery not possible!'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE., &
                  HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

         ELSEIF (DT > DT_MIN) THEN

            WRITE(ERR_MSG,"(3X,'Recovered: Dt=',G12.5,' :-)')") DT
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

            CALL RESET_NEW

! Iterate again with new dt
            ADJUSTDT = .TRUE.

! Cut the timestep for 2nd order accurate time implementation.
            IF(CN_ADJUST_DT) DT = 0.5d0*DT

! Set the return flag stop iterating.
         ELSE

! Prevent DT from dropping below DT_MIN.
            IF(PERSISTENT_MODE) DT = max(DT, DT_MIN)
            ADJUSTDT = .FALSE.
         ENDIF

      ENDIF

! Update ONE/DT variable.
      ODT = ONE/DT

      RETURN
      END FUNCTION ADJUSTDT

      END MODULE ITERATE
