! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: run                                                    C
!  Purpose: Common block containing run control data                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE run

! Modules
!---------------------------------------------------------------------//
      use param, only: dim_M, dim_eqs
      use param1, only: UNDEFINED_I
      use derived_types
!---------------------------------------------------------------------//

! Main filename to be used for output files  Name must
! still be legal after extensions are added to it.
      CHARACTER(LEN=60) :: RUN_NAME

! Brief description of the problem.
      CHARACTER(LEN=60) :: DESCRIPTION

! Units for data input and output: CGS.
      CHARACTER(LEN=16) :: UNITS

! Type of run: NEW, RESTART
      CHARACTER(LEN=16) :: RUN_TYPE

! Variable which triggers automatic restart
      LOGICAL :: AUTOMATIC_RESTART

! counter to keep track of how many auto_retart were performed
      INTEGER :: ITER_RESTART

! version.release of software
      CHARACTER(LEN=10) :: ID_VERSION

! Start-time of the run.
      DOUBLE PRECISION :: TIME

! Stop-time of the run.
      DOUBLE PRECISION :: TSTOP

! Time step.
      DOUBLE PRECISION :: DT

! 1./Time step.
      DOUBLE PRECISION :: oDT

! Indicates whether simulation is steady-state
      LOGICAL :: STEADY_STATE

! Number of times steps completed.
      INTEGER :: NSTEP

! Declare a new variable to use on CN with RESTART cases
! Number of time steps when restart file was read
      INTEGER :: NSTEPRST

! Discretization scheme for different equations
      INTEGER :: DISCRETIZE(DIM_EQS)

! Use Chi scheme for discretizing certain equation sets
!  (species mass fractions)
      LOGICAL :: Chi_scheme

! If .TRUE. solve X momentum equations
      LOGICAL :: MOMENTUM_X_EQ(0:DIM_M)

! If .TRUE. solve Y momentum equations
      LOGICAL :: MOMENTUM_Y_EQ(0:DIM_M)

! If .TRUE. solve Z momentum equations
      LOGICAL :: MOMENTUM_Z_EQ(0:DIM_M)

! IF .TRUE. use Jackson form momentum equations
      LOGICAL :: JACKSON
! IF .TRUE. use Ishii form momentum equations
      LOGICAL :: ISHII

! If .TRUE. use Model-B momentum equations
      LOGICAL :: Model_B

! If .TRUE. include added (virtual) mass in momentum eq.
      LOGICAL :: Added_Mass

! phase number where added mass is applied.
      INTEGER :: M_AM

! If .TRUE. solve K_Epsilon turbulence eq.
      LOGICAL :: K_Epsilon

! If .TRUE. solve energy equations
      LOGICAL :: ENERGY_EQ

! If .TRUE. use the deferred correction method
      LOGICAL :: DEF_COR

! If .TRUE. use the fourth order interpolation
      LOGICAL :: FPFOI

! If .TRUE. activate 2nd order accurate time implementation
      LOGICAL :: CN_ON

! If .TRUE. solve granular energy equations
      LOGICAL :: GRANULAR_ENERGY

! If .TRUE. solve species balance equations
      LOGICAL :: SPECIES_EQ(0:DIM_M)

! If .TRUE. one of the species equations is being solved
      LOGICAL :: ANY_SPECIES_EQ

! If .TRUE. call user-defined subroutines
      LOGICAL :: CALL_USR

! If .TRUE. force time-step when NIT=MAX_NIT and DT=DT_MIN
      LOGICAL :: PERSISTENT_MODE

! If .TRUE. solve population balance  equations
      LOGICAL :: Call_DQMOM

! If .TRUE. incorporate the wall effects upon the calculation of the
! subgrid solids viscosity, solids pressure, and gas-solids drag
      LOGICAL :: SUBGRID_Wall
! the ratio of the FilterSize to the GridSize
      DOUBLE PRECISION :: filter_size_ratio

! Single particle drag correlation
      CHARACTER(64) :: CD_FUNCTION

! Parameter used to calculate lubrication interactions between
! different particles in HYS drag model
      DOUBLE PRECISION :: LAM_HYS

! If .TRUE. use Simonin model (k_epsilon must also be true)
      LOGICAL :: SIMONIN

! If .TRUE. use Ahmadi model (k_epsilon must also be true)
      LOGICAL :: AHMADI

! If .TRUE. calculate frictional stress terms
      LOGICAL :: FRICTION
! Form of friction model:
!             If 0: use S:S
!             If 1: use the form of Savage to compute S:S
!             If 2: use combination of both for frictional stress terms
      INTEGER :: SAVAGE

! If .TRUE. use Scheffer frictional stress (default set to .TRUE.)
      LOGICAL :: SCHAEFFER

! If .TRUE. use blending frictional/kinetic stresses
! (default set to .FALSE. do not blend)
      LOGICAL :: BLENDING_STRESS
      LOGICAL :: TANH_BLEND ! default set to true
      LOGICAL :: SIGM_BLEND ! default set to false

! If .TRUE. use Jenkins small friction BC
      LOGICAL :: JENKINS
! If .TRUE. use revised phip for JJ BC
      LOGICAL :: BC_JJ_M
! If .TRUE. output PHIP to JJ_PHIP.dat
      LOGICAL :: PHIP_OUT_JJ
! to write specularity
      INTEGER :: PHIP_OUT_ITER

! If .TRUE. treat system as if shearing
      LOGICAL :: SHEAR
! Shear Vel
      DOUBLE PRECISION :: V_sh

! If .TRUE. use Yu and Standish correlation to compute ep_star
      LOGICAL :: YU_STANDISH

! If .TRUE. use Fedors and Landel correlation to compute ep_star
      LOGICAL :: FEDORS_LANDEL

! STOP Trigger mechanism to terminate MFIX normally before batch
! queue terminates flag variable to check for end of batch queue when
! set to TRUE check performed at the beginning of each time step and
! termination of mfix triggered after saving all files if condition
! is met
      LOGICAL :: CHK_BATCHQ_END
! variable to store the total wall clock duration of the batch queue
! session wall clock time specified in seconds
! for jaguarcnl@NCCS max wall clock limit is 2.5 hr limit up to 512
! processors
      DOUBLE PRECISION :: BATCH_WALLCLOCK
! variable to set a buffer time before the batch queue session ends to
! make sure once MFIX is triggered to shutdown, there is sufficient
! time to save files, make copies to HPSS storage before batch queue
! time runs out. Current logic in MFIX checks for:
!    if CPU_TIME > (BATCH_WALLCLOCK - TERM_BUFFER) then
!    save all .RES .SP files and trigger shutdown
      DOUBLE PRECISION :: TERM_BUFFER

! If .TRUE. code will automatically restart for DT < DT_MIN
      LOGICAL :: AUTO_RESTART

! If .TRUE. code will automatically restart for DT < DT_MIN
      LOGICAL :: REINITIALIZING = .FALSE.

! Time-step failure rate:
! 1) Number of failed time steps
! 2) Observation window
      INTEGER :: TIMESTEP_FAIL_RATE(2)

! parameters for dynamically adjusting time step
! +1 -> increase dt; -1 decrease dt
      INTEGER :: DT_dir = -1

! Maximum Time step.
      DOUBLE PRECISION :: DT_MAX

! Minimum Time step.
      DOUBLE PRECISION :: DT_MIN

! Time step adjustment factor (<1.0)
      DOUBLE PRECISION :: DT_FAC

! The previous time step used in iterate (before it is
! changed by adjust_dt)
      DOUBLE PRECISION :: DT_prev

! in case iterations converged and DT modified, use old dt
! to advance time in time_march.
      LOGICAL :: use_DT_prev

! Slope limiter parameter (0 < C _FAC <= 1.0)
      DOUBLE PRECISION :: C_FAC

! If .TRUE. reduce time step when residuals do not decrease
      LOGICAL :: DETECT_STALL

! String which controls reduction of global sums for residual
! calculations
      LOGICAL :: DEBUG_RESID

       common /run_dp/ time      !for Linux


! Flags indicating variable solids density.
      LOGICAL :: SOLVE_ROs(DIM_M), ANY_SOLVE_ROs

! Specifies the type of solids: TFM, DEM, MPPIC
      CHARACTER(len=3), DIMENSION(DIM_M) :: SOLIDS_MODEL

! Flags for various solids phase models.
      LOGICAL :: TFM_SOLIDS
      LOGICAL :: DEM_SOLIDS
      LOGICAL :: PIC_SOLIDS
! The number of the various solids phases.
      INTEGER :: TFM_COUNT = 0
      INTEGER :: DEM_COUNT = 0
      INTEGER :: PIC_COUNT = 0

      ! Error index
      INTEGER :: IER

      ! CPU time unit.
      CHARACTER(LEN=4) :: TUNIT

      CONTAINS

         !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
         !  Purpose:  Given time in seconds, calculate time in days/hours/seconds
         !
         !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
         SUBROUTINE GET_TUNIT(TLEFT, TUNIT)

            !-----------------------------------------------
            ! Modules
            !-----------------------------------------------
            IMPLICIT NONE
            !-----------------------------------------------
            ! Dummy arguments
            !-----------------------------------------------
            DOUBLE PRECISION, INTENT(INOUT) :: TLEFT
            CHARACTER(LEN=4) :: TUNIT
            !-----------------------------------------------

            IF (TLEFT < 3600.0d0) THEN
               TUNIT = 's'
            ELSE
               TLEFT = TLEFT/3600.0d0
               TUNIT = 'h'
               IF (TLEFT >= 24.) THEN
                  TLEFT = TLEFT/24.0d0
                  TUNIT = 'days'
               ENDIF
            ENDIF

            RETURN
         END SUBROUTINE GET_TUNIT

      END MODULE RUN
