!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_RUN_CONTROL                                       !
!  Purpose: Check the run control namelist section                     !
!                                                                      !
!  Author: P.Nicoletti                                Date: 27-NOV-91  !
!          J.Musser                                   Date: 31-JAN-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_RUN_CONTROL


! Global Variables:
!---------------------------------------------------------------------//
! New or restart
      USE run, only: RUN_TYPE
! Brief description of simulation.
      USE run, only: DESCRIPTION
! Simulation units: SI, CGS
      USE run, only: UNITS
! Simulation start/stop times.
      USE run, only: TIME, TSTOP
! Time step size, one over time step size.
      USE run, only: DT, ODT, STEADY_STATE
! Flag: Use K-Epsilon turbulence model.
      USE run, only: K_EPSILON
      USE run, only: ishii, jackson
! Turbulence lenghth scale and viscosity bound.
      use constant, only: L_SCALE0, MU_GMAX

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: UNDEFINED, UNDEFINED_C
      USE param1, only: ONE, ZERO

! Global Module procedures:
!---------------------------------------------------------------------//
      USE error_manager

      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//


!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_RUN_CONTROL")


! Clear out the run description if not specified.
      IF (DESCRIPTION == UNDEFINED_C) DESCRIPTION = ' '

! Verify UNITS input.
      IF(UNITS == UNDEFINED_C) THEN
         WRITE(ERR_MSG,1000) 'UNITS'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF((UNITS /= 'CGS') .AND. (UNITS /= 'SI')) THEN
         WRITE(ERR_MSG,1001) 'UNITS', UNITS
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Verify that DT is valid.
      IF (DT < ZERO) THEN
         WRITE(ERR_MSG,1002) 'DT', DT
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

! Steady-state simulation.
      ELSEIF(DT == UNDEFINED .OR. DT == ZERO) THEN
         STEADY_STATE = .TRUE.
         DT = ZERO
         ODT = ZERO
         TIME = ZERO

! Transient simulation.
      ELSE
         STEADY_STATE = .FALSE.
! Calculate reciprocal of initial timestep.
         ODT = ONE/DT
! Verify the remaining time settings.
         IF (TIME == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'TIME'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF (TSTOP == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) 'TSTOP'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF (TIME < ZERO) THEN
            WRITE(ERR_MSG,1002)'TIME', TIME
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF (TSTOP < ZERO) THEN
            WRITE(ERR_MSG,1002) 'TSTOP', TSTOP
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

! Verify the run type.
      IF(.NOT.(RUN_TYPE=='NEW' .OR. RUN_TYPE=='RESTART_1'              &
         .OR. RUN_TYPE=='RESTART_2')) THEN
         WRITE(ERR_MSG,1001) 'RUN_TYPE', RUN_TYPE
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Turbulence model:
      IF (K_Epsilon .AND. L_SCALE0 /= ZERO) THEN
         WRITE(ERR_MSG,2001)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 2001 FORMAT('Error 2001: Cannot set K_EPSILON = .T. and specify ',    &
         'L_SCALE0 /= ZERO')
      ENDIF

! Ishii and jackson form of governing equations cannot both be invoked
      IF (ISHII .AND. JACKSON) THEN
         WRITE(ERR_MSG,2002)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 2002 FORMAT('Error 2002: Cannot set both ISHII = .T. and JACKSON = ',&
             '.T.',/,'Please correct the mfix.dat file.')
      ENDIF

!  Check whether MU_gmax is specified for turbulence (sof)
      IF (K_Epsilon .AND. MU_GMAX==UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) 'MU_GMAX'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Clear the error manager
      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal or unknown input: ',A,' = ',G14.4,/  &
         'Please correct the mfix.dat file.')

 1003 FORMAT('Error 1003: Illegal or unknown input: ',A,' = ',I4,/     &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_RUN_CONTROL



