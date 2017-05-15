!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  MODULE: GET_DATA                                                    C
!  Purpose: read and verify input data, open files                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
MODULE read_input

  CONTAINS

      SUBROUTINE GET_DATA

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE check_data_cg, only: adjust_ijk_size, check_data_cartesian
      USE cut_cell_preproc, only: cut_cell_preprocessing
      USE compar
      USE constant, only: L_SCALE0
      USE cutcell
      USE dashboard
      USE des_allocate
      USE des_rxns
      USE des_thermo
      USE desgrid, only: DESGRID_INIT
      USE discretelement
      USE error_manager
      USE funits
      USE gridmap
      USE iterate, only: max_nit
      USE leqsol
      USE mfix_pic
      USE mpi_init_des, only: DESMPI_INIT
      USE parallel
      USE param
      USE param1
      USE qmom_kinetic_equation
      USE run
      USE stl_preproc_des, only: DES_STL_PREPROCESSING
      USE visc_g, only: L_SCALE

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! shift DX, DY and DZ values
      LOGICAL, PARAMETER :: SHIFT = .TRUE.

! This module call routines to initialize the namelist variables.
      CALL INIT_NAMELIST
! Read in the namelist variables from the ascii input file.
      CALL READ_NAMELIST(0,"mfix.dat")

! Initialize the error manager. This call occurs after the mfix.dat
! is read so that message verbosity can be set and the .LOG file
! can be opened.
      CALL INIT_ERROR_MANAGER

! Write header in the .LOG file and to screen.
! Not sure if the values are correct or useful
      CALL WRITE_HEADER

! Open files
      CALL OPEN_FILES(RUN_NAME, RUN_TYPE, N_SPX)

! These checks verify that sufficient information was provided
! to setup the domain indices and DMP gridmap.
      CALL CHECK_GEOMETRY_PREREQS
      CALL CHECK_DMP_PREREQS

! Set up the physical domain indicies (cell index max/min values).
      CALL SET_MAX2

! Set constants
      CALL SET_CONSTANTS

! Adjust partition for better load balance (done when RE_INDEXING is .TRUE.)
      CALL ADJUST_IJK_SIZE

! Partition the domain and set indices
      CALL GRIDMAP_INIT

! Check the minimum solids phase requirements.
      CALL CHECK_SOLIDS_MODEL_PREREQS

      CALL CHECK_RUN_CONTROL
      CALL CHECK_NUMERICS
      CALL CHECK_OUTPUT_CONTROL

      CALL CHECK_GAS_PHASE
      CALL CHECK_SOLIDS_PHASES
      CALL SET_PARAMETERS

! Basic geometry checks.
      CALL CHECK_GEOMETRY(SHIFT)
      IF(DISCRETE_ELEMENT) CALL CHECK_GEOMETRY_DES

! Set grid spacing variables.
      CALL SET_GEOMETRY
      IF(DISCRETE_ELEMENT) CALL SET_GEOMETRY_DES

      CALL CHECK_INITIAL_CONDITIONS
      CALL CHECK_BOUNDARY_CONDITIONS
      CALL CHECK_INTERNAL_SURFACES
      CALL CHECK_POINT_SOURCES

      CALL CHECK_CHEMICAL_RXNS
      CALL CHECK_ODEPACK_STIFF_CHEM



!----------------------  DOMAIN SPECIFIC CHECKS  --------------------!


! This call needs to occur before any of the IC/BC checks.
      CALL SET_ICBC_FLAG

! Compute area of boundary surfaces.
      CALL GET_BC_AREA

! Convert (mass, volume) flows to velocities.
      CALL SET_BC_FLOW

! Set the flags for identifying computational cells
      CALL SET_FLAGS
! Set arrays for computing indices
      CALL SET_INCREMENTS
      CALL SET_INCREMENTS3

! Cartesian grid implementation
      CALL CHECK_DATA_CARTESIAN
      IF(CARTESIAN_GRID) THEN
         CALL CUT_CELL_PREPROCESSING
      ELSE
         CALL ALLOCATE_DUMMY_CUT_CELL_ARRAYS
      ENDIF

      IF(DISCRETE_ELEMENT) THEN
         CALL DESGRID_INIT
         CALL DESMPI_INIT
         CALL DES_STL_PREPROCESSING
      ENDIF

!--------------------------  ARRAY ALLOCATION -----------------------!

! Allocate array storage.
      CALL ALLOCATE_ARRAYS
      IF(DISCRETE_ELEMENT) CALL DES_ALLOCATE_ARRAYS
      IF(QMOMK) CALL QMOMK_ALLOCATE_ARRAYS

! Initialize arrays.
      CALL INIT_FVARS
      IF(DISCRETE_ELEMENT) CALL DES_INIT_ARRAYS

! This is all that happens in SET_L_SCALE so it needs moved, maybe
! this should go in int_fluid_var.?
!     CALL SET_L_SCALE
      L_SCALE(:) = L_SCALE0

!======================================================================
! Data initialization for Dashboard
!======================================================================
      INIT_TIME = TIME
      SMMIN =  LARGE_NUMBER
      SMMAX = -LARGE_NUMBER

      DTMIN =  LARGE_NUMBER
      DTMAX = -LARGE_NUMBER

      NIT_MIN = MAX_NIT
      NIT_MAX = 0

      N_DASHBOARD = 0


      RETURN

    END SUBROUTINE GET_DATA

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_CONSTANTS                                           C
!  Purpose: Set various constants                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-JAN-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_CONSTANTS

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, one, undefined
      USE constant, only: gas_const
      USE constant, only: gravity, gravity_x, gravity_y, gravity_z
      USE constant, only: to_SI
      USE constant, only: k_scale
      USE run, only: LAM_HYS, UNITS
      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: flush_err_msg

      IMPLICIT NONE
!---------------------------------------------------------------------//

! Note that the cell flags are not set when this routine is called.

! Dimensionless constants
      K_SCALE = .08D0   ! this actually isn't used anywhere...

! Enter the value of all constants in various units (CGS or SI)
      IF (UNITS == 'SI') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 9.80665D0 ! m/s2
         GAS_CONST = 8314.56D0                !Pa.m3/kmol.K, or kg m2/s2 kmol K (Perry and Green, 1984)
         to_SI = 0.1D0                        !to convert dyne/cm2 to Pa. e.g. calc_mu_g.f
         IF (LAM_HYS == UNDEFINED) LAM_HYS = 0.000001d0    ! m
      ELSEIF (UNITS == 'CGS') THEN
         IF (GRAVITY == UNDEFINED) GRAVITY = 980.665D0 !cm/s2
         GAS_CONST = 8314.56D4                !g.cm2/s2.mol.K
         to_SI = ONE                          !does not do anything in CGS. e.g.: calc_mu_g.f
         IF (LAM_HYS == UNDEFINED) LAM_HYS = 0.0001d0    ! cm
      ELSE
! Initialize the error manager.
         CALL INIT_ERR_MSG("SET_CONSTANTS")
         WRITE(ERR_MSG,1005) UNITS
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         CALL FINL_ERR_MSG
 1005 FORMAT('Error 1005: Unknown UNITS = ',A,/ &
         'Please correct the mfix.dat file.')
      ENDIF

! If all components of the gravitational acceleration vector are
! undefined (zero), then use the default value for the negative
! y-direction. This ensures backward compatibility with the old
! (legacy) GRAVITY keyword. At this point GRAVITY is defined,
! either from mfix.dat or by default above
      IF(GRAVITY_X==ZERO.AND.GRAVITY_Y==ZERO.AND.GRAVITY_Z==ZERO) THEN
         GRAVITY_X = ZERO
         GRAVITY_Y = - GRAVITY
         GRAVITY_Z = ZERO
      ENDIF

      RETURN
      END SUBROUTINE SET_CONSTANTS
    END MODULE read_input
