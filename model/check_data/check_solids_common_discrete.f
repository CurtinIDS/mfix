!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_COMMON_DISCRETE                            !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE

! Modules
!---------------------------------------------------------------------//
! Runtime Flag: Generate initial particle configuration.
      USE discretelement, only: GENER_PART_CONFIG
! Runtime Flag: Store DES_*_OLD arrays.
      USE discretelement, only: DO_OLD
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! User specified integration method.
      USE discretelement, only: DES_INTG_METHOD
      USE discretelement, only: INTG_ADAMS_BASHFORTH
      USE discretelement, only: INTG_EULER
! User specified neighbor search method.
      USE discretelement, only: DES_NEIGHBOR_SEARCH
! User specified data out format (VTP, TecPlot)
      USE discretelement, only: DES_OUTPUT_TYPE
! Max/Min particle radii
      USE discretelement, only: MAX_RADIUS, MIN_RADIUS
! Runtime Flag: Periodic boundaries
      USE discretelement, only: DES_PERIODIC_WALLS
      USE discretelement, only: DES_PERIODIC_WALLS_X
      USE discretelement, only: DES_PERIODIC_WALLS_Y
      USE discretelement, only: DES_PERIODIC_WALLS_Z
! Use the error manager for posting error messages.
      use error_manager
! Runtime Flag: Invoke MPPIC model.
      USE mfix_pic, only: MPPIC
      USE mpi_utility

      use param1, only: undefined, undefined_c
      use param, only: dim_m
! number of continuous solids phases and
! solids 'phase' diameters and densities
      USE physprop, only: MMAX, D_p0, RO_s0
! Calculated baseline variable solids density.
      USE physprop, only: CLOSE_PACKED
! Runtime Flag: Solve energy equations
      USE run, only: ENERGY_EQ
! Runtime Flag: One or more species equations are solved.
      use run, only: ANY_SPECIES_EQ
! Flag: Solve variable solids density.
      use run, only: SOLVE_ROs
      use run, only: SOLIDS_MODEL
      USE run, only: MOMENTUM_X_EQ
      USE run, only: MOMENTUM_Y_EQ
      USE run, only: MOMENTUM_Z_EQ
      use run, only: RUN_TYPE
      implicit none

! Local Variables
!---------------------------------------------------------------------//
      INTEGER :: M 
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE")

! Determine the maximum particle size in the system (MAX_RADIUS), which
! in turn is used for various tasks
      MAX_RADIUS = -UNDEFINED
      MIN_RADIUS =  UNDEFINED
! For number of continuous solids phases (use MMAX rather than SMAX to
! accomodate GHD particularity)
      DO M = MMAX+1,DES_MMAX+MMAX
         MAX_RADIUS = MAX(MAX_RADIUS, 0.5d0*D_P0(M))
         MIN_RADIUS = MIN(MIN_RADIUS, 0.5d0*D_P0(M))
      ENDDO

! Set close_packed to true to prevent possible issues stemming from the
! pressure correction equation.  Specifically, if closed_packed is false
! then a mixture pressure correction equation is invoked and this is not
! correctly setup for DEM.  To do so would require ensuring that
! 1) the solids phase continuum quantities used in these equations are
!    correctly set based on their DEM counterparts and
! 2) the pressure correction coefficients for such solids phases are
!    also calculated (currently these calculations are turned off
!    when using DEM)
      CLOSE_PACKED((MMAX+1):DIM_M) = .TRUE.


! Turn off the 'continuum' equations for discrete solids if the user
! specified them.  We could make use of these flags.
      MOMENTUM_X_EQ((MMAX+1):DIM_M) = .FALSE.
      MOMENTUM_Y_EQ((MMAX+1):DIM_M) = .FALSE.
      MOMENTUM_Z_EQ((MMAX+1):DIM_M) = .FALSE.

! Derive periodicity from cyclic boundary flags.
      DES_PERIODIC_WALLS_X = CYCLIC_X .OR. CYCLIC_X_PD
      DES_PERIODIC_WALLS_Y = CYCLIC_Y .OR. CYCLIC_Y_PD
      DES_PERIODIC_WALLS_Z = CYCLIC_Z .OR. CYCLIC_Z_PD

      DES_PERIODIC_WALLS = (DES_PERIODIC_WALLS_X .OR.                  &
        DES_PERIODIC_WALLS_Y .OR. DES_PERIODIC_WALLS_Z)


! Overwrite for restart cases.
      IF(TRIM(RUN_TYPE) .NE. 'NEW') GENER_PART_CONFIG = .FALSE.

! Check for valid neighbor search option.
      SELECT CASE(DES_NEIGHBOR_SEARCH)
      CASE (1) ! N-Square
      CASE (2)
         WRITE(ERR_MSG,2001) 2, 'QUADTREE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      CASE (3)
         WRITE(ERR_MSG,2001) 3, 'OCTREE'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      CASE (4) ! Grid based
      CASE DEFAULT
         WRITE(ERR_MSG,2001) DES_NEIGHBOR_SEARCH,'UNKNOWN'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2001 FORMAT('Error 2001:Invalid DES_NEIGHBOR_SEARCH method: ',I2,1X,  &
         A,/'Please correct the mfix.dat file.')

      END SELECT


! Check the output file format
      IF(DES_OUTPUT_TYPE == UNDEFINED_C) DES_OUTPUT_TYPE = 'PARAVIEW'
      SELECT CASE(trim(DES_OUTPUT_TYPE))
      CASE ('PARAVIEW')
      CASE ('TECPLOT')
      CASE DEFAULT
         WRITE(ERR_MSG,2010) trim(DES_OUTPUT_TYPE)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2010 FORMAT('Error 2010:Invalid DES_OUTPUT_TYPE: ',A,/'Please ',       &
         'correct the mfix.dat file.')

      END SELECT


! Check for valid integration method
      SELECT CASE(trim(DES_INTG_METHOD))
      CASE ('EULER')
         INTG_EULER = .TRUE.
         INTG_ADAMS_BASHFORTH = .FALSE.
         !DES_INTG_METHOD_ENUM = 1
      CASE ('ADAMS_BASHFORTH')
         INTG_EULER = .FALSE.
         INTG_ADAMS_BASHFORTH = .TRUE.
         !DES_INTG_METHOD_ENUM = 2
      CASE DEFAULT
         WRITE(ERR_MSG,2020) trim(DES_INTG_METHOD)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 2020 FORMAT('Error 2020:Invalid DES_INGT_METHOD: ',A,/'Please ',      &
         'correct the mfix.dat file.')

      END SELECT

      DO_OLD = INTG_ADAMS_BASHFORTH .OR. MPPIC

! Check interpolation input.
      CALL CHECK_SOLIDS_COMMON_DISCRETE_INTERP

! Set flags for energy equations
      IF(ENERGY_EQ) CALL CHECK_SOLIDS_COMMON_DISCRETE_ENERGY

! Check thermodynamic properties of discrete solids.
      IF(ANY_SPECIES_EQ) &
         CALL CHECK_SOLIDS_COMMON_DISCRETE_THERMO

! Check geometry constrains.
      CALL CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY

      CALL FINL_ERR_MSG


      RETURN

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose: Check input parameters for solving discrete solids phase   !
!  energy equations.  Only DEM simulations (neither hybrid nor MPPIC)  !
!  can invoke particle-particle heat transfer. Therefore, checks for   !
!  those functions are reseved for later.                              !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY

! Modules
!---------------------------------------------------------------------//
      use param1, only: ZERO, UNDEFINED
      use des_thermo, only: DES_CONV_CORR
      use des_thermo, only: DES_CONV_CORR_ENUM
      use des_thermo, only: RANZ_1952
      use des_thermo, only: SB_CONST
      use des_thermo, only: DES_Em
      use des_thermo, only: CALC_CONV_DES ! Convection
      use des_thermo, only: CALC_COND_DES ! Conduction
      use des_thermo, only: CALC_RADT_DES ! Radiation
! Flag to explicitly couple source terms and DES
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use discretelement, only: DES_MMAX
      use discretelement, only: DES_CONTINUUM_COUPLED
! Use the error manager for posting error messages.
      use error_manager
! User input for DES interpolation scheme.
      use particle_filter, only: DES_INTERP_SCHEME
! Enumerated interpolation scheme for faster access
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE

      use physprop, only: MMAX
      use physprop, only: K_S0
      use run, only: UNITS
      use run, only: SOLIDS_MODEL

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter
      INTEGER :: M

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_ENERGY")


! Set runtime flags for which modes of heat transfer to calculate.
      CALC_CONV_DES = DES_CONTINUUM_COUPLED
      DO M = MMAX+1, MMAX+DES_MMAX
! Flag to calculate radiation.
         IF(DES_Em(M) > ZERO) CALC_RADT_DES(M) = .TRUE.
! Flag to calculate conduction.
         CALC_COND_DES(M) = (K_s0(M) > ZERO .AND. K_s0(M) /= UNDEFINED)
      ENDDO

! Gas/Solids convection:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Verify the selected convective heat transfer coefficient model
      SELECT CASE(TRIM(DES_CONV_CORR))
! Ranz, W.E. and Marshall, W.R., "Friction and transfer coefficients
! for single particles and packed beds,"  Chemical Engineering Science,
! Vol. 48, No. 5, pp 247-253, 1952.
      CASE ('RANZ_1952')
         DES_CONV_CORR_ENUM = RANZ_1952
! If the heat transfer coefficient correlation provided by the user does
! not match one of the models outlined above, flag the error and exit.
      CASE DEFAULT
         WRITE(ERR_MSG,1001)'DES_CONV_CORR', trim(DES_CONV_CORR)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      END SELECT


! Radiation Equation:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
! Verify that a emmisivity value is specified for each solids phase
      DO M = MMAX+1, MMAX+DES_MMAX
         IF(DES_Em(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_Em',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

! Set the value of the Stefan-Boltzman Constant based on the units
      IF(UNITS == 'SI')THEN
         SB_CONST = 5.6704d0*(10.0d0**(-8)) ! W/((m^2).K^4)
      ELSE
         SB_CONST = 1.355282d0*(10.0d0**(-12)) ! cal/((cm^2).sec.K^4)
      ENDIF


! Notify that interpolation is not support for thermo variables
      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_NONE)
      CASE DEFAULT
         WRITE(ERR_MSG,2000) trim(adjustl(DES_INTERP_SCHEME))
         CALL FLUSH_ERR_MSG()
      END SELECT

 2000 FORMAT('WARNING 2000: The selected interpolation scheme (',A,    &
         ') is not',/'supported by the DES energy equation implemen',  &
         'tation. All energy',/'equation variables will use the ',     &
         'centroid method for interphase',/'data exchange.')

      IF(DES_EXPLICITLY_COUPLED)THEN
         WRITE(ERR_MSG, 2100)
         CALL FLUSH_ERR_MSG!(ABORT=.TRUE.)
      ENDIF

 2100 FORMAT('Error 2100: The DES Energy equation implementation ',    &
         'does not',/'currently support explicit coupling (DES_',      &
         'EXPLICITLY_COUPLED).',/'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG


      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_ENERGY




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_THERMO                     !
!  Author: J.Musser                                   Date: 17-Jun-10  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_THERMO

! Modules
!---------------------------------------------------------------------//
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use error_manager
! User input for DES interpolation scheme.
      use particle_filter, only: DES_INTERP_SCHEME
! Enumerated interpolation scheme for faster access
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use stiff_chem, only: STIFF_CHEMISTRY
      IMPLICIT NONE

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_THERMO")

! Stiff chemistry solver is a TFM reaction model not for DES.
      IF(STIFF_CHEMISTRY) THEN
         WRITE(ERR_MSG,9003)
         CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
      ENDIF

 9003 FORMAT('Error 9003: The stiff chemistry solver is not ',         &
      'available in DES',/'simulations. Please correct the input file.')

! Notify that interpolation is not support for thermo variables
      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_NONE)
      CASE DEFAULT
         WRITE(ERR_MSG,2000) trim(adjustl(DES_INTERP_SCHEME))
         CALL FLUSH_ERR_MSG()
      END SELECT

 2000 FORMAT('WARNING 2000: The selected interpolation scheme (',A,    &
         ') is not',/'supported by the DES Species equation implemen', &
         'tation. All energy',/'equation variables will use the ',     &
         'centroid method for interphase',/'data exchange.')

      IF(DES_EXPLICITLY_COUPLED)THEN
         WRITE(ERR_MSG, 2100)
         CALL FLUSH_ERR_MSG!(ABORT=.TRUE.)
      ENDIF

 2100 FORMAT('Error 2100: The DES Species equation implementation ',   &
         'does not',/'currently support explicit coupling (DES_',      &
         'EXPLICITLY_COUPLED).',/'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_THERMO


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY                   !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Purpose: Check user input data                                      !
!                                                                      !
!  Comments: Geometry checks were moved here from CHECK_DES_DATA.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY

! Modules
!---------------------------------------------------------------------//
! Flag: Use Cartesian grid cut-cell implementation
      USE cutcell, only: CARTESIAN_GRID
! Flag: Use STL representation in CG
      USE cutcell, only: USE_STL
! Flag: Use DES E-L model
      USE discretelement, only: DES_CONTINUUM_COUPLED
      USE discretelement, only: MAX_RADIUS
      use error_manager
      USE geometry, only: COORDINATES
      USE geometry, only: NO_I, NO_J
      USE geometry, only: ZLENGTH
      IMPLICIT NONE

! Local Variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: MIN_DEPTH

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY")


! DEM/MPPIC is restricted to CARTESIAN coordinates.
      IF(COORDINATES == 'CYLINDRICAL') THEN
         WRITE (ERR_MSG, 1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error: 1100: DES and MPPIC models only support ',        &
         'CARTESIAN coordinates.')


! Check dimension. This is redundant with check_data_03.
      IF(NO_I .OR. NO_J) THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1200 FORMAT('Error 1200: Illegal geometry for DEM/MPPIC. 2D ',        &
         'simulations are',/'restricted to the XY plane. Please ',     &
         'correct the mfix.dat file.')


      IF(DES_CONTINUUM_COUPLED)THEN
! Check that the depth of the simulation exceeds the largest particle
! to ensure correct calculation of volume fraction. This is important
! for coupled simulations.
         MIN_DEPTH = 2.0d0*MAX_RADIUS
         IF(ZLENGTH < MIN_DEPTH)THEN
            WRITE(ERR_MSG, 1300)
            CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
         ENDIF
      ENDIF

 1300 FORMAT('Error 1300: The maximum particle diameter exceeds the ', &
         'simulation',/'depth (ZLENGTH). Please correct the mfix.dat ',&
         'file.')


      IF(CARTESIAN_GRID .AND. .NOT.USE_STL) THEN
         WRITE(ERR_MSG,1400)
         CALL FLUSH_ERR_MSG(ABORT =.TRUE.)
      ENDIF

 1400 FORMAT('Error 1400: Cartesian grid and discrete models (DEM or ',&
         'PIC) only',/'support STL wall representations. Quadrics ',   &
         'and polygons are not',/'supported.')


      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_GEOMETRY

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_COMMON_DISCRETE_INTERP                     !
!  Author: J.Musser                                   Date: 25-Nov-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_INTERP

! Modules
!---------------------------------------------------------------------//
! Runtime Flag: Invoke gas/solids coupled simulation.
      use discretelement, only: DES_CONTINUUM_COUPLED
      use error_manager
! Runtime FLag: 3D simulation
      use geometry, only: DO_K
! Runtime Flag: Invoke MPPIC model.
      USE mfix_pic, only: MPPIC
      use param1, only: UNDEFINED
! User input for DES interpolation scheme.
      use particle_filter, only: DES_INTERP_SCHEME
! Enumerated interpolation scheme for faster access
      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_NONE
      use particle_filter, only: DES_INTERP_GARG
      use particle_filter, only: DES_INTERP_DPVM
      use particle_filter, only: DES_INTERP_GAUSS
      use particle_filter, only: DES_INTERP_LHAT
! User specified filter width
      use particle_filter, only: DES_INTERP_WIDTH
! Flag: Diffuse DES field variables.
      use particle_filter, only: DES_DIFFUSE_MEAN_FIELDS
! Diffusion filter width
      use particle_filter, only: DES_DIFFUSE_WIDTH
! Flag: Interpolate continuum fields
      use particle_filter, only: DES_INTERP_MEAN_FIELDS
! Flag: Interplate variables for drag calculation.
      use particle_filter, only: DES_INTERP_ON
! Size of interpolation filter
      use particle_filter, only: FILTER_SIZE
      IMPLICIT NONE

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_COMMON_DISCRETE_INTERP")

! Set the runtime flag for diffusing mean fields
      DES_DIFFUSE_MEAN_FIELDS = (DES_DIFFUSE_WIDTH /= UNDEFINED)

! Set the interpolation ENUM value.
      SELECT CASE(trim(adjustl(DES_INTERP_SCHEME)))
      CASE ('NONE')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_NONE
! Cannot use interpolation when no scheme is selected.
         IF(DES_INTERP_ON)THEN
            WRITE(ERR_MSG,2001) 'DES_INTERP_ON'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(DES_INTERP_MEAN_FIELDS)THEN
            WRITE(ERR_MSG,2001) 'DES_INTERP_MEAN_FIELDS'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

         ELSEIF(DES_CONTINUUM_COUPLED) THEN
            IF(MPPIC) THEN
               WRITE(ERR_MSG,2002) 'MPPIC solids'
               CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
            ELSEIF(MPPIC) THEN
               WRITE(ERR_MSG,2002) 'Cartesian grid cut-cells'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

      CASE ('GARG_2012')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_GARG

      CASE ('SQUARE_DPVM')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_DPVM

      CASE ('GAUSS_DPVM')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_GAUSS

      CASE ('LINEAR_HAT')
         DES_INTERP_SCHEME_ENUM = DES_INTERP_LHAT

      CASE DEFAULT
         WRITE(ERR_MSG,2000) trim(adjustl(DES_INTERP_SCHEME))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      END SELECT

 2000 FORMAT('Error 2000: Invalid DES_INTERP_SCHEME: ',A,/'Please ',   &
         'correct the mfix.dat file.')

 2001 FORMAT('Error 2001: No interpolation scheme specified when ',A,/ &
         'is enabled. Please correct the mfix.dat file.')

 2002 FORMAT('Error 2002: DES simulations utilizing ',A,' require',/   &
         'interpolation (DES_INTERP_ON and DES_INTERP_MEANFIELDS). ',/ &
         'Please correct the mfix.dat file.')


      SELECT CASE(DES_INTERP_SCHEME_ENUM)

      CASE(DES_INTERP_NONE)

         IF(DES_INTERP_WIDTH /= UNDEFINED) THEN
            WRITE(ERR_MSG,2100) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2100 FORMAT('Error 2100: The selected interpolation scheme (',A,') ', &
         'does',/'not support an adjustable interpolation width.',/    &
         'Please correct the input file.')


      CASE(DES_INTERP_GARG)
         DES_INTERP_MEAN_FIELDS= .TRUE.

         IF(DES_INTERP_WIDTH /= UNDEFINED) THEN
            WRITE(ERR_MSG,2100) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(DES_DIFFUSE_MEAN_FIELDS) THEN
            WRITE(ERR_MSG,2110) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2110 FORMAT('Error 2110: The selected interpolation scheme (',A,') ', &
         'does not',/'support diffusive filtering of mean field ',     &
          'quantities. Please correct',/'the input file.')

      CASE(DES_INTERP_DPVM, DES_INTERP_GAUSS)

! Set the size of the interpolation filter.
         FILTER_SIZE = merge(27, 9, DO_K)

         IF(DES_INTERP_WIDTH == UNDEFINED) THEN
            WRITE(ERR_MSG,2120) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 2120 FORMAT('Error 2120: The selected interpolation scheme (',A,') ', &
         'requires',/'a DES_INTERP_WIDTH. Please correct the ',        &
         'input file.')


      CASE(DES_INTERP_LHAT)

! Set the size of the interpolation filter.
         FILTER_SIZE = merge(27, 9, DO_K)

         IF(DES_INTERP_WIDTH /= UNDEFINED) THEN
            WRITE(ERR_MSG,2100) trim(adjustl(DES_INTERP_SCHEME))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

      END SELECT

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_SOLIDS_COMMON_DISCRETE_INTERP
