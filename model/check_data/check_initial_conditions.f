!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_INITIAL_CONDITIONS                                !
!  Author: P. Nicoletti                               Date: 02-DEC-91  !
!  Author: J.Musser                                   Date: 01-MAR-14  !
!                                                                      !
!  Purpose: check the initial conditions input section                 !
!     - check geometry of any specified IC region                      !
!     - check specification of physical quantities                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_INITIAL_CONDITIONS

! Global Variables:
!---------------------------------------------------------------------//
! Flag: IC geometry was detected.
      use ic, only: IC_DEFINED
! Flag: DEM solids present.
      use run, only: DEM_SOLIDS
! Flag: New run or a restart.
      use run, only: RUN_TYPE
! Runtime flag specifying MPPIC solids
      use run, only: PIC_SOLIDS
! Flag: Adjusting runtime parameters
      use run, only: REINITIALIZING

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of IC.
      use param, only: DIMENSION_IC

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter for ICs
      INTEGER :: ICV
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_INITIAL_CONDITIONS")

! Determine which ICs are DEFINED
      CALL CHECK_IC_GEOMETRY

! Loop over all IC arrays.
      DO ICV=1, DIMENSION_IC

! Verify user input for defined defined IC.
         IF(IC_DEFINED(ICV)) THEN
! Gas phase checks.
            CALL CHECK_IC_GAS_PHASE(ICV)
! Generic solids phase checks.
            CALL CHECK_IC_SOLIDS_PHASES(ICV)

! Verify that no data was defined for unspecified IC. ICs are only
! defined for new runs, so these checks are restricted to new runs.
         ELSEIF(RUN_TYPE == 'NEW' .AND. .NOT.REINITIALIZING) THEN
            CALL CHECK_IC_OVERFLOW(ICV)
         ENDIF
      ENDDO



! Check the initial conditions for the DEM and MPPIC models as well
      IF(DEM_SOLIDS.OR.PIC_SOLIDS) &
      CALL CHECK_IC_COMMON_DISCRETE
      IF(DEM_SOLIDS) CALL CHECK_IC_DEM
      IF(PIC_SOLIDS) CALL CHECK_IC_MPPIC

! Finalize the error manager.
      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_INITIAL_CONDITIONS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GEOMETRY                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_IC_GEOMETRY


! Global Variables:
!---------------------------------------------------------------------//
! Flag: IC contains geometric data and/or specified type
      use ic, only: IC_DEFINED
! Flag: IC type.
      use ic, only: IC_TYPE
! User specified: IC geometry
      use ic, only: IC_X_e, IC_X_w, IC_I_e, IC_I_w
      use ic, only: IC_Y_n, IC_Y_s, IC_J_n, IC_J_s
      use ic, only: IC_Z_t, IC_Z_b, IC_K_t, IC_K_b
! User specified: System geometry
      use geometry, only: NO_I, IMAX, IMIN1, IMAX1, XLENGTH, DX, XMIN
      use geometry, only: NO_J, JMAX, JMIN1, JMAX1, YLENGTH, DY
      use geometry, only: NO_K, KMAX, KMIN1, KMAX1, ZLENGTH, DZ
! Flag: New run or a restart.
      use run, only: RUN_TYPE
      use run, only: REINITIALIZING

! Global Parameters:
!---------------------------------------------------------------------//
! The max number of ICs.
      use param, only: DIMENSION_IC
! Parameter constants
      use param1, only: ZERO, UNDEFINED, UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      implicit none


! Local Variables:
!---------------------------------------------------------------------//
! Loop/variable indices
      INTEGER :: ICV
! Local spatial indices.
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_GEOMETRY")

! Check geometry of any specified IC region
      DO ICV = 1, DIMENSION_IC

         IC_DEFINED(ICV) = .FALSE.
         IF (IC_X_W(ICV) /= UNDEFINED)   IC_DEFINED(ICV) = .TRUE.
         IF (IC_X_E(ICV) /= UNDEFINED)   IC_DEFINED(ICV) = .TRUE.
         IF (IC_Y_S(ICV) /= UNDEFINED)   IC_DEFINED(ICV) = .TRUE.
         IF (IC_Y_N(ICV) /= UNDEFINED)   IC_DEFINED(ICV) = .TRUE.
         IF (IC_Z_B(ICV) /= UNDEFINED)   IC_DEFINED(ICV) = .TRUE.
         IF (IC_Z_T(ICV) /= UNDEFINED)   IC_DEFINED(ICV) = .TRUE.
         IF (IC_I_W(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE.
         IF (IC_I_E(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE.
         IF (IC_J_S(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE.
         IF (IC_J_N(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE.
         IF (IC_K_B(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE.
         IF (IC_K_T(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE.

! An IC is defined for restart runs only if it is a 'PATCH'.
         IF(RUN_TYPE /= 'NEW' .AND. IC_TYPE(ICV) /= 'PATCH') &
            IC_DEFINED(ICV) = .FALSE.

! Ignore patched IC regions for new runs. It may be better to flag this as
! and error to avoid user confusion.
         IF(RUN_TYPE == 'NEW' .AND. IC_TYPE(ICV) == 'PATCH') &
            IC_DEFINED(ICV) = .FALSE.

! Enable only PATCH IC regions when initializing.
         IF(REINITIALIZING) IC_DEFINED(ICV)=(IC_TYPE(ICV)=='PATCH')

         IF(.NOT.IC_DEFINED(ICV)) CYCLE

         IF (IC_X_W(ICV)==UNDEFINED .AND. IC_I_W(ICV)==UNDEFINED_I) THEN
            IF (NO_I) THEN
               IC_X_W(ICV) = ZERO
            ELSE
               WRITE(ERR_MSG, 1100) ICV, 'IC_X_w and IC_I_w'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF (IC_X_E(ICV)==UNDEFINED .AND. IC_I_E(ICV)==UNDEFINED_I) THEN
            IF (NO_I) THEN
               IC_X_E(ICV) = XLENGTH
            ELSE
               WRITE(ERR_MSG, 1100) ICV, 'IC_X_e and IC_I_e'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF (IC_Y_S(ICV)==UNDEFINED .AND. IC_J_S(ICV)==UNDEFINED_I) THEN
            IF (NO_J) THEN
               IC_Y_S(ICV) = ZERO
            ELSE
               WRITE(ERR_MSG, 1100) ICV, 'IC_Y_s and IC_J_s'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF (IC_Y_N(ICV)==UNDEFINED .AND. IC_J_N(ICV)==UNDEFINED_I) THEN
            IF (NO_J) THEN
               IC_Y_N(ICV) = YLENGTH
            ELSE
               WRITE(ERR_MSG, 1100) ICV, 'IC_Y_n and IC_J_n'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF (IC_Z_B(ICV)==UNDEFINED .AND. IC_K_B(ICV)==UNDEFINED_I) THEN
            IF (NO_K) THEN
               IC_Z_B(ICV) = ZERO
            ELSE
               WRITE(ERR_MSG, 1100) ICV, 'IC_Z_b and IC_K_b'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF (IC_Z_T(ICV)==UNDEFINED .AND. IC_K_T(ICV)==UNDEFINED_I) THEN
            IF (NO_K) THEN
               IC_Z_T(ICV) = ZLENGTH
            ELSE
               WRITE(ERR_MSG, 1100) ICV, 'IC_Z_t and IC_K_t'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

      ENDDO   ! end loop over (icv = 1,dimension_ic)

 1100 FORMAT('Error 1100: Initial condition region ',I3,' is ill-',    &
         'defined.',/' > ',A,' are not specified.',/'Please correct ', &
         'the mfix.dat file.')


      DO ICV = 1, DIMENSION_IC

! Skip this check if the IC region is not specified.
         IF(.NOT.IC_DEFINED(ICV)) CYCLE

         IF (IC_X_W(ICV)/=UNDEFINED .AND. IC_X_E(ICV)/=UNDEFINED) THEN
            IF (NO_I) THEN
               I_W = 1
               I_E = 1
            ELSE
               CALL CALC_CELL (XMIN, IC_X_W(ICV), DX, IMAX, I_W)
               I_W = I_W + 1
               CALL CALC_CELL (XMIN, IC_X_E(ICV), DX, IMAX, I_E)
            ENDIF
            IF (IC_I_W(ICV)/=UNDEFINED_I .OR. IC_I_E(ICV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK (IC_I_W(ICV), I_W, ICV, 'IC - west')
               CALL LOCATION_CHECK (IC_I_E(ICV), I_E, ICV, 'IC - east')
            ELSE
               IC_I_W(ICV) = I_W
               IC_I_E(ICV) = I_E
            ENDIF
         ENDIF

! Report problems with calculated bounds.
         IF(IC_I_W(ICV) > IC_I_E(ICV)) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_W > IC_I_E'
             write(*,*)' dump:',IC_I_W(ICV),IC_I_E(ICV)
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_W(ICV) < IMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_W < IMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_W(ICV) > IMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_W > IMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_E(ICV) < IMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_I_E < IMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_I_E(ICV) > IMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_Z_t and IC_K_t'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (IC_Y_S(ICV)/=UNDEFINED .AND. IC_Y_N(ICV)/=UNDEFINED) THEN
            IF (NO_J) THEN
               J_S = 1
               J_N = 1
            ELSE
               CALL CALC_CELL (ZERO, IC_Y_S(ICV), DY, JMAX, J_S)
               J_S = J_S + 1
               CALL CALC_CELL (ZERO, IC_Y_N(ICV), DY, JMAX, J_N)
            ENDIF
            IF (IC_J_S(ICV)/=UNDEFINED_I .OR. IC_J_N(ICV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK (IC_J_S(ICV), J_S, ICV, 'IC - south')
               CALL LOCATION_CHECK (IC_J_N(ICV), J_N, ICV, 'IC - north')
            ELSE
               IC_J_S(ICV) = J_S
               IC_J_N(ICV) = J_N
            ENDIF
         ENDIF

         IF(IC_J_S(ICV) > IC_J_N(ICV)) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_S > IC_J_N'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_S(ICV)<JMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_S < JMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_S(ICV)>JMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_S >  JMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_N(ICV)<JMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_N < JMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_J_N(ICV)>JMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_J_N > JMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


         IF (IC_Z_B(ICV)/=UNDEFINED .AND. IC_Z_T(ICV)/=UNDEFINED) THEN
            IF (NO_K) THEN
               K_B = 1
               K_T = 1
            ELSE
               CALL CALC_CELL (ZERO, IC_Z_B(ICV), DZ, KMAX, K_B)
               K_B = K_B + 1
               CALL CALC_CELL (ZERO, IC_Z_T(ICV), DZ, KMAX, K_T)
            ENDIF
            IF (IC_K_B(ICV)/=UNDEFINED_I .OR. IC_K_T(ICV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK (IC_K_B(ICV), K_B, ICV, 'IC - bottom')
               CALL LOCATION_CHECK (IC_K_T(ICV), K_T, ICV, 'IC - top')
            ELSE
               IC_K_B(ICV) = K_B
               IC_K_T(ICV) = K_T
            ENDIF
         ENDIF

         IF(IC_K_B(ICV) > IC_K_T(ICV)) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_B > IC_K_T'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_B(ICV) < KMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_B < KMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_B(ICV) > KMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_B > KMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_T(ICV) < KMIN1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_T < KMIN1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_K_T(ICV) > KMAX1) THEN
             WRITE(ERR_MSG, 1101) ICV, 'IC_K_T > KMAX1'
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF


 1101 FORMAT('Error 1101: Initial condition region ',I2,' is ill-',    &
         'defined.',/3x,A,/'Please correct the mfix.dat file.')

      ENDDO   ! end loop over (icv=1,dimension_ic)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_IC_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_GAS_PHASE                                       !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Verify gas phase input variables in IC region.              !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_IC_GAS_PHASE(ICV)

! Global Variables:
!---------------------------------------------------------------------//
      use constant, only: l_scale0
! Gas phase volume fraction, pressure, temperature, species.
      use ic, only: IC_EP_g, IC_P_g, IC_T_g, IC_X_g
! Gas phase velocity components.
      use ic, only: IC_U_g, IC_V_g, IC_W_g
! Radiation model parameters.
      use ic, only: IC_GAMA_RG, IC_T_RG
! K-Epsilon model parameters.
      use ic, only: IC_K_TURB_G, IC_E_TURB_G
! L-scale model parameters
      use ic, only: IC_L_SCALE
! IC for user-defined scalar equation.
      use ic, only: IC_SCALAR
! IC Type: UNDEFINED or PATCH.
      use ic, only: IC_TYPE

! Flag. Solve Energy equations
      use run, only: ENERGY_EQ
! Flag. Solve Species equations
      use run, only: SPECIES_EQ
! Flag: Solve K-Epsilon; K-E parameters
      use run, only: K_Epsilon
! Specified constant gas density and viscosity.
      use physprop, only: RO_G0, MU_G0
! Specified average molecular weight
      use physprop, only: MW_AVG
! Number of gas phase species
      use physprop, only: NMAX
! Specified number of scalar equations.
      use scalars, only: NSCALAR
! Flag: Do not solve in specified direction.
      use geometry, only: NO_I, NO_J, NO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, ONE, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      use toleranc

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: ICV

! Local Variables:
!---------------------------------------------------------------------//
! Loop index
      INTEGER :: N
! Sum of mass fraction.
      DOUBLE PRECISION :: SUM
! Flag for regular IC regions
      LOGICAL :: BASIC_IC
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_GAS_PHASE")

! Patch ICs skip various checks.
      BASIC_IC = (IC_TYPE(ICV) /= 'PATCH')

! Check that gas phase velocity components are initialized.
      IF(BASIC_IC) THEN
         IF(IC_U_G(ICV) == UNDEFINED) THEN
            IF(NO_I) THEN
               IC_U_G(ICV) = ZERO
            ELSE
               WRITE(ERR_MSG, 1000) trim(iVar('IC_U_g',ICV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(IC_V_G(ICV) == UNDEFINED) THEN
            IF(NO_J) THEN
               IC_V_G(ICV) = ZERO
            ELSE
               WRITE(ERR_MSG, 1000) trim(iVar('IC_V_g',ICV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(IC_W_G(ICV) == UNDEFINED) THEN
            IF (NO_K) THEN
               IC_W_G(ICV) = ZERO
            ELSE
               WRITE(ERR_MSG, 1000) trim(iVar('IC_W_g',ICV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
           ENDIF
         ENDIF
      ENDIF

! Check that gas phase void fraction is initialized. Patched ICs may
! have an undefined volume fration. A second check is preformed on
! the solids.
      IF(IC_EP_G(ICV) == UNDEFINED .AND. BASIC_IC) THEN
         WRITE(ERR_MSG, 1000) trim(iVar('IC_EP_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check that if the gas phase pressure is initialized and the gas is
! compressible that the gas phase pressure is not zero or negative
      IF(IC_P_G(ICV) /= UNDEFINED) THEN
         IF(RO_G0==UNDEFINED .AND. IC_P_G(ICV)<=ZERO) THEN
            WRITE(ERR_MSG, 1100) trim(iVar('IC_P_g',ICV)),             &
               iVal(IC_P_G(ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

 1100 FORMAT('Error 1100: Pressure must be greater than 0.0 for ',     &
         'compressible flow',/'Illegal value: ',A,' = ',A,/'Please ',  &
         'correct the mfix.dat file.')

      IF(BASIC_IC) THEN
         IF(ENERGY_EQ .OR. RO_G0==UNDEFINED .OR. MU_G0==UNDEFINED) THEN
            IF(IC_T_G(ICV)==UNDEFINED) THEN
               WRITE(ERR_MSG, 1000) trim(iVar('IC_T_g',ICV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDIF

! Gas phase radiation values.
      IF (ENERGY_EQ) THEN
         IF (IC_GAMA_RG(ICV) < ZERO) THEN
            WRITE(ERR_MSG, 1001) trim(iVar('IC_GAMA_Rg',ICV)),         &
               iVal(IC_GAMA_RG(ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF (IC_GAMA_RG(ICV) > ZERO) THEN
            IF (IC_T_RG(ICV) == UNDEFINED) THEN
               WRITE(ERR_MSG, 1000) iVar('IC_T_Rg',ICV)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDIF


! First: sum together defiend gas phase species mass fractions.
      SUM = ZERO
      DO N = 1, NMAX(0)
         IF(IC_X_G(ICV,N) /= UNDEFINED) THEN
            SUM = SUM + IC_X_G(ICV,N)
         ELSEIF(BASIC_IC) THEN
            IC_X_G(ICV,N) = ZERO
         ENDIF
      ENDDO


! Enforce that the species mass fractions must sum to one.
      IF(.NOT.COMPARE(ONE,SUM)) THEN

! No error for ZERO sum PATCH IC data.
         IF(.NOT.BASIC_IC .AND. COMPARE(ZERO,SUM))THEN

! Error for regular IC regions when solving species equations.
         ELSEIF(SPECIES_EQ(0)) THEN
            WRITE(ERR_MSG, 1110) ICV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1110 FORMAT('Error 1110: IC_X_g(',I3,',:) do NOT sum to ONE and the ',&
         'gas phase',/'species equations are solved. Please correct ', &
         'the mfix.dat file.')

! Error for slightly compressible flows with MW_AVG specified.
         ELSEIF(RO_G0 == UNDEFINED .AND. MW_AVG == UNDEFINED) THEN
            WRITE(ERR_MSG, 1111) ICV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1111 FORMAT('Error 1111: IC_X_g(',I3,',:) do NOT sum to ONE and the ',&
         'gas phase',/'is compressible and MW_AVG is UNDEFINED.',/     &
         'Please correct the mfix.dat file.')

         ELSEIF(.NOT.COMPARE(SUM,ZERO)) THEN
            WRITE(ERR_MSG, 1112) ICV
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1112 FORMAT('Error 1112: IC_X_g(',I3,',:) do not sum to ONE or ZERO ',&
         'and they',/'are not needed. Please correct the mfix.dat ',   &
         'file.')

         ELSE
            IC_X_G(ICV,:) = ZERO
            IC_X_G(ICV,1) = ONE
         ENDIF
      ENDIF


      DO N = 1, NScalar
         IF(IC_Scalar(ICV,N) == UNDEFINED) IC_Scalar(ICV,N) = ZERO
      ENDDO


      IF(BASIC_IC) THEN
         IF (L_SCALE0 == ZERO) THEN
            IF (IC_L_SCALE(ICV) /= UNDEFINED) THEN
               WRITE(ERR_MSG, 1113) ICV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1113 FORMAT('Error 1113: IC_L_SCALE(',I3,',:) is defined but L_SCALE0 ',&
         'is equal',/,'to zero. A non-zero value must be specified to ',&
         'activate this model.',/,'Please correct the mfix.dat file.')
            ENDIF
         ELSEIF (L_SCALE0 < ZERO) THEN
            WRITE(ERR_MSG, 1001) 'L_SCALE0', iVal(L_scale0)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSE   ! l_scale0 is defined at not zero
            IF (IC_L_SCALE(ICV) < ZERO) THEN
               WRITE(ERR_MSG, 1001) iVar('IC_L_SCALE',ICV), &
                  iVal(IC_L_SCALE(ICV))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDIF


      IF(K_Epsilon .AND. BASIC_IC) THEN
         IF (IC_K_Turb_G(ICV) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) iVar('IC_K_Turb_G',ICV)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF (IC_E_Turb_G(ICV) == UNDEFINED) THEN
            WRITE(ERR_MSG, 1000) iVar('IC_E_Turb_G',ICV)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_IC_GAS_PHASE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_IC_SOLIDS_PHASES                                  !
!  Author: P. Nicoletti                               Date: 02-DEC-91  !
!  Author: J.Musser                                   Date: 01-MAR-14  !
!                                                                      !
!  Purpose: Verify solids phase(s) input variables in IC region.       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_IC_SOLIDS_PHASES(ICV)


! Global Variables:
!---------------------------------------------------------------------//
! Solids volume fraction, bulk density
      use ic, only: IC_EP_s, IC_ROP_s
! Solids velocity components.
      use ic, only: IC_U_s, IC_V_s, IC_W_s
! Solids temperature, mass fractions, granular energy
      use ic, only: IC_T_s, IC_X_s, IC_THETA_M
! Radiation model parameters.
      use ic, only: IC_GAMA_RS, IC_T_RS
! Gas phase volume fraction and temperature.
      use ic, only: IC_EP_g, IC_T_g
! IC Type: UNDEFINED or PATCH.
      use ic, only: IC_TYPE

! Flag. Solve Energy equations
      use run, only: ENERGY_EQ
! Flag. Solve Species equations
      use run, only: SPECIES_EQ
! Flag. Solve Granular Energy PDE
      use run, only: GRANULAR_ENERGY
! Flag. Solve variable solids density
      use run, only: SOLVE_ROs
! Baseline solids mass fraction, index of intert
      use physprop, only: X_S0, INERT_SPECIES
! Specified constant solids density.
      use physprop, only: RO_S0
! Number of gas phase species
      use physprop, only: NMAX
! Number of TFM solids phases.
      use physprop, only: SMAX
! Number of DEM or PIC solids.
      use discretelement, only: DES_MMAX
! Flag: Do not solve in specified direction.
      use geometry, only: NO_I, NO_J, NO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Maximum number of solids species.
      use param, only: DIM_M
! Parameter constants
      use param1, only: ZERO, ONE, UNDEFINED
! function to calculate solids density.
      USE eos, ONLY: EOSS

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      use toleranc

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Index of IC region.
      INTEGER, INTENT(IN) :: ICV
! Loop/variable index
      INTEGER :: M, N
! Various sums.
      DOUBLE PRECISION SUM, SUM_EP
! Solids phase density in IC region.
      DOUBLE PRECISION :: IC_ROs(1:DIM_M)
! Index of inert species
      INTEGER :: INERT
! Flag to skip checks on indexed solid phase.
      LOGICAL :: SKIP(1:DIM_M)
! Total number of solids phases (TFM + DEM + MPPIC)
      INTEGER :: MMAX_TOT
! Flag for PATCH IC regions
      LOGICAL :: BASIC_IC
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_SOLIDS_PHASES")

! Patch ICs skip various checks.
      BASIC_IC = (IC_TYPE(ICV) /= 'PATCH')

! The total number of solids phases (all models).
      MMAX_TOT = SMAX + DES_MMAX

! Calculate EP_s from EP_g if there is only one solids phase.
      IF(MMAX_TOT == 1 .AND. IC_EP_S(ICV,1) == UNDEFINED) THEN
         IF(IC_EP_g(ICV) /= UNDEFINED) IC_EP_S(ICV,1) = ONE-IC_EP_g(ICV)
      ENDIF

! Bulk density or solids volume fraction must be explicitly defined
! if there are more than one solids phase.
      IF(MMAX_TOT > 1 .AND. .NOT.COMPARE(IC_EP_g(ICV),ONE)) THEN
! IC_EP_g may be undefined for PATCH IC regions.
         IF(IC_EP_g(ICV) /= UNDEFINED) THEN
            DO M = 1, MMAX_TOT
               IF(IC_ROP_S(ICV,M) == UNDEFINED .AND. &
                  IC_EP_S(ICV,M) == UNDEFINED) THEN
                  WRITE(ERR_MSG, 1400) M, ICV, 'IC_ROP_s and IC_EP_s'
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO

! If IC_EP_G is undefined, then ROP_s and EP_s should be too.
         ELSE
            DO M = 1, MMAX_TOT
               IF(IC_ROP_S(ICV,M) /= UNDEFINED .AND. &
                  IC_EP_S(ICV,M) /= UNDEFINED) THEN
                  WRITE(ERR_MSG, 1400) M, ICV, 'IC_ROP_s and IC_EP_s'
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

 1400 FORMAT('Error 1400: Insufficient solids phase ',I2,' ',          &
         'information for IC',/'region ',I3,'. ',A,' not specified.',/ &
         'Please correct the mfix.dat file.')

! Determine which solids phases are present.
      DO M = 1, MMAX_TOT
         SKIP(M)=(IC_ROP_S(ICV,M)==UNDEFINED.OR.IC_ROP_S(ICV,M)==ZERO) &
            .AND.(IC_EP_S(ICV,M)==UNDEFINED .OR.IC_EP_S(ICV,M)==ZERO)
      ENDDO

      IF(MMAX_TOT == 1 .AND. IC_EP_g(ICV)/=ONE) SKIP(1) = .FALSE.

      DO M=1, MMAX_TOT

! check that solids phase m velocity components are initialized
         IF(BASIC_IC) THEN
            IF(IC_U_S(ICV,M) == UNDEFINED) THEN
               IF (SKIP(M) .OR. NO_I) THEN
                  IC_U_S(ICV,M) = ZERO
               ELSE
                  WRITE(ERR_MSG, 1000)trim(iVar('IC_U_s',ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF

            IF(IC_V_S(ICV,M) == UNDEFINED) THEN
               IF(SKIP(M) .OR. NO_J) THEN
                  IC_V_S(ICV,M) = ZERO
               ELSE
                  WRITE(ERR_MSG, 1000)trim(iVar('IC_V_s',ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF

            IF(IC_W_S(ICV,M) == UNDEFINED) THEN
               IF(SKIP(M) .OR. NO_K) THEN
                  IC_W_S(ICV,M) = ZERO
               ELSE
                  WRITE(ERR_MSG, 1000)trim(iVar('IC_W_s',ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF

            IF(ENERGY_EQ .AND. IC_T_S(ICV,M)==UNDEFINED) THEN
               IF(SKIP(M)) THEN
                  IC_T_S(ICV,M) = IC_T_G(ICV)
               ELSE
                  WRITE(ERR_MSG, 1000)trim(iVar('IC_T_s',ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
            ENDIF

            IF(GRANULAR_ENERGY .AND. IC_THETA_M(ICV,M)==UNDEFINED) THEN
               IF(SKIP(M)) THEN
                  IC_THETA_M(ICV,M) = ZERO
               ELSE
                  WRITE(ERR_MSG, 1000)trim(iVar('IC_Theta_M',ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
           ENDIF

            IF(ENERGY_EQ) THEN
               IF(IC_GAMA_RS(ICV,M) < ZERO) THEN
                  WRITE(ERR_MSG, 1001)trim(iVar('IC_GAMA_Rs',ICV,M)),  &
                     iVal(IC_GAMA_RS(ICV,M))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ELSEIF (IC_GAMA_RS(ICV,M) > ZERO) THEN
                  IF(IC_T_RS(ICV,M) == UNDEFINED) THEN
                     WRITE(ERR_MSG, 1001)trim(iVar('IC_T_Rs',ICV,M))
                     CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF


! First: sum together defiend species mass fractions.
         SUM = ZERO
         DO N = 1, NMAX(M)
            IF(IC_X_S(ICV,M,N) /= UNDEFINED) THEN
               SUM = SUM + IC_X_S(ICV,M,N)
            ELSEIF(BASIC_IC) THEN
               IC_X_S(ICV,M,N) = ZERO
            ENDIF
         ENDDO

! Enforce that the species mass fractions must sum to one.
         IF(.NOT.COMPARE(ONE,SUM)) THEN

! Summing to zero is not an error for PATCH IC regions.
            IF(.NOT.BASIC_IC .AND. COMPARE(ZERO,SUM)) THEN

            ELSEIF(SPECIES_EQ(M) .AND. .NOT.SKIP(M)) THEN
               WRITE(ERR_MSG, 1402) ICV, M
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1402 FORMAT('Error 1402: IC_X_s(',I3,',',I2,',:) do NOT sum to ONE ', &
         'and the solids phase',/'species equations are solved. ',     &
         'Please correct the mfix.dat file.')

            ELSEIF(SOLVE_ROS(M) .AND. .NOT.SKIP(M)) THEN
               WRITE(ERR_MSG, 1403) ICV, M
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1403 FORMAT('Error 1403: IC_X_s(',I3,',',I2,':) do NOT sum to ONE ',  &
         'and the solids phase',/'density is calculated. Please ',     &
         'correct the mfix.dat file.')

            ELSEIF(.NOT.COMPARE(SUM,ZERO)) THEN
                WRITE(ERR_MSG, 1404) ICV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1404 FORMAT('Error 1404: IC_X_s(',I3,',',I2,':) do NOT sum to ONE ',  &
         'or ZERO and',/'they are not needed. Please correct the ',    &
         'mfix.dat file.')

            ELSE
               IC_X_S(ICV,M,:) = ZERO
               IC_X_S(ICV,M,1) = ONE
            ENDIF
         ENDIF

! Set the solids density for the IC region.
         IF(SKIP(M)) THEN
            IC_ROs(M) = merge(RO_s0(M), RO_s0(M), SOLVE_ROs(M))

         ELSEIF(.NOT.SOLVE_ROs(M)) THEN
            IC_ROs(M) = RO_s0(M)

         ELSE
! Verify that the species mass fraction for the inert material is not
! zero in the IC region when the solids is present.
            INERT = INERT_SPECIES(M)
            IF(IC_X_S(ICV,M,INERT) == ZERO) THEN
               WRITE(ERR_MSG,1405) M, ICV
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

 1405 FORMAT('Error 1405: No inert species for phase ',I2,' in IC ',   &
         'region ',I3,'.',/'Unable to calculate solids phase density. ',&
         'Please refer to the Readme',/' file for required variable ', &
         'solids density  model input parameters and',/' make the ',   &
         'necessary corrections to the data file.')

! Calculate the solids density.
            IC_ROs(M) = EOSS(RO_s0(M), X_s0(M,INERT),                  &
               IC_X_S(ICV,M,INERT))
         ENDIF

      ENDDO   ! end loop over (m=1,smax)


! Initialize the sum of the total volume fraction.
      SUM_EP = IC_EP_G(ICV)

      DO M=1, MMAX_TOT

! Clear out both varaibles if this phase is skipped.
         IF(BASIC_IC .AND. SKIP(M)) THEN
            IC_EP_S(ICV,M)  = ZERO
            IC_ROP_S(ICV,M) = ZERO

! Leave everything undefined for PATCH ICs that are not specifed.
         ELSEIF(.NOT.BASIC_IC .AND. (IC_ROP_S(ICV,M) == UNDEFINED      &
            .AND. IC_EP_S(ICV,M) == UNDEFINED)) THEN

! If both input parameters are defined. Make sure they are equivalent.
         ELSEIF(IC_ROP_S(ICV,M) /= UNDEFINED .AND.                     &
            IC_EP_S(ICV,M) /= UNDEFINED) THEN

            IF(.NOT.COMPARE(IC_EP_S(ICV,M)*IC_ROs(M),                  &
               IC_ROP_S(ICV,M))) THEN

! BASIC_IC regions require that the IC_ROP_s and IC_EP_s specifications
! match although it is unlikely that anyone would specify both.
               IF(BASIC_IC) THEN

                  WRITE(ERR_MSG,1406) M, ICV
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

1406 FORMAT('Error 1406: IC_EP_s and IC_ROP_s are inconsistent for ',&
         'phase ',I2,/,'in IC region ', I3,'. Please correct the ',&
         'mfix.dat file.')


! PachedeIC regions defer to IC_EP_s if the values do not match. This
! prevents a dead lock or the need to define both. This case is rather
! common as a defined IC_EP_s is converted to IC_ROP_s. Therefore, if
! a patch region is used more than once, these values may not match.
               ELSE

                  WRITE(ERR_MSG,1407) trim(iVar('IC_ROP_s',ICV,M)), &
                     trim(iVAL(IC_ROP_S(ICV,M))), trim(iVar('IC_EP_s',&
                     ICV,M)), trim(iVAL(IC_EP_S(ICV,M)))
                  CALL FLUSH_ERR_MSG()

1407 FORMAT('Warning 1407: IC_EP_s and IC_ROP_s are inconsistent:',    &
         2(/3x,A,' = ',A),/'Deferring to IC_EP_s to overcome conflict.')

                  IC_ROP_S(ICV,M) = IC_EP_S(ICV,M)*IC_ROs(M)

               ENDIF
            ENDIF


! Compute IC_EP_s from IC_ROP_s
         ELSEIF(IC_EP_S(ICV,M) == UNDEFINED)THEN
            IC_EP_S(ICV,M) = IC_ROP_S(ICV,M) / IC_ROs(M)

! Compute IC_ROP_s from IC_EP_s and IC_ROs
         ELSEIF(IC_ROP_S(ICV,M) == UNDEFINED) THEN
            IC_ROP_S(ICV,M) = IC_EP_S(ICV,M) * IC_ROs(M)
! This is a sanity check.
         ELSE

         ENDIF
! Add this phase to the total volume fraction.
         SUM_EP = SUM_EP + IC_EP_S(ICV,M)
      ENDDO

! Verify that the volume fractions sum to one.
      IF(BASIC_IC .AND. .NOT.COMPARE(SUM_EP,ONE)) THEN
         WRITE(ERR_MSG,1410) ICV
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1410 FORMAT('Error 1410: Illegal initial condition region : ',I3,/    &
         'Sum of volume fractions does NOT equal ONE. Please correct',/&
         'the mfix.dat file.')


      CALL FINL_ERR_MSG

      RETURN


 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Illegal data in initial condition region ', &
         I3,/A,' defined for unknown solids phase ',I2,'.',/          &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_IC_SOLIDS_PHASES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_IC_OVERFLOW                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Verify that no data was defined for unspecified IC.         !
!                                                                      !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!

      SUBROUTINE CHECK_IC_OVERFLOW(ICV)

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase volume fraction, pressure, temperature, species.
      use ic, only: IC_EP_g, IC_T_g, IC_X_g
! Gas phase velocity components.
      use ic, only: IC_U_g, IC_V_g, IC_W_g
! Radiation model parameters.
      use ic, only: IC_T_RG
! K-Epsilon model parameters.
      use ic, only: IC_K_TURB_G, IC_E_TURB_G
! IC for user-defined scalar equation.
      use ic, only: IC_SCALAR
! Solids volume fraction, bulk density
      use ic, only: IC_ROP_s
! Solids velocity components.
      use ic, only: IC_U_s, IC_V_s, IC_W_s
! Solids temperature, mass fractions, granular energy
      use ic, only: IC_T_s, IC_X_s
! Radiation model parameters.
      use ic, only: IC_T_RS
! IC Type: UNDEFINED or PATCH.
      use ic, only: IC_TYPE

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constant
      use param1, only: UNDEFINED
! Maximum number of solids phases and user-defined scalars.
      use param, only: DIM_M, DIM_SCALAR
! Maximum number of gas and solids phase species
      use param, only: DIM_N_g, DIM_N_s

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Dummy Arguments:
!---------------------------------------------------------------------//
! IC region index.
      INTEGER, INTENT(IN) :: ICV

! Local Variables:
!---------------------------------------------------------------------//
! Loop counters
      INTEGER :: M, N
!......................................................................!

      IF (IC_TYPE(ICV) == 'PATCH') RETURN

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_IC_OVERFLOW")

! GAS PHASE quantities
! -------------------------------------------->>>
      IF(IC_U_G(ICV) /= UNDEFINED) THEN
          WRITE(ERR_MSG, 1010) trim(iVar('IC_U_g',ICV))
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IC_V_G(ICV) /= UNDEFINED) THEN
         WRITE(ERR_MSG, 1010) trim(iVar('IC_V_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IC_W_G(ICV) /= UNDEFINED) THEN
         WRITE(ERR_MSG, 1010) trim(iVar('IC_W_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IC_EP_G(ICV) /= UNDEFINED) THEN
         WRITE(ERR_MSG, 1010) trim(iVar('IC_EP_g',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IC_T_G(ICV) /= UNDEFINED) THEN
          WRITE(ERR_MSG, 1010) trim(iVar('IC_T_g',ICV))
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IC_T_RG(ICV) /= UNDEFINED) THEN
          WRITE(ERR_MSG, 1010) trim(iVar('IC_T_Rg',ICV))
          CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IC_K_Turb_G(ICV) /= UNDEFINED) THEN
         WRITE(ERR_MSG, 1010) trim(iVar('IC_K_Turb_G',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(IC_E_Turb_G(ICV) /= UNDEFINED) THEN
         WRITE(ERR_MSG, 1010) trim(iVar('IC_E_Turb_G',ICV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      DO N = 1, DIM_N_G
         IF(IC_X_G(ICV,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_X_g',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO
      DO N = 1, DIM_SCALAR
         IF(IC_Scalar(ICV,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_Scalar',ICV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO
! --------------------------------------------<<<


! SOLIDS PHASE quantities
! -------------------------------------------->>>
      DO M=1, DIM_M
         IF(IC_ROP_S(ICV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_ROP_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_U_S(ICV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_U_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_V_S(ICV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_V_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_W_S(ICV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_W_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_T_S(ICV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_T_s',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(IC_T_RS(ICV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG, 1010) trim(iVar('IC_T_Rs',ICV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         DO N = 1, DIM_N_S
            IF(IC_X_S(ICV,M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG, 1010) trim(iVar('IC_X_s',ICV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO
      ENDDO
! --------------------------------------------<<<

      CALL FINL_ERR_MSG
      RETURN

 1010 FORMAT('Error 1010: ',A,' specified in an undefined IC region')

      END SUBROUTINE CHECK_IC_OVERFLOW
