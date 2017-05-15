!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_DES_SOLIDS                                        !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DEM")

! Particle-particle collision parameters.
      CALL CHECK_SOLIDS_DEM_COLLISION
! DES cohesion model parameters.
      CALL CHECK_SOLIDS_DEM_COHESION
! Particle-particle conduction model parameters.
      CALL CHECK_SOLIDS_DEM_ENERGY

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_DEM

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE CHECK_SOLIDS_DEM_ENERGY                                  !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_ENERGY

! Global Variables
!---------------------------------------------------------------------//

      use des_thermo, only: DES_MIN_COND_DIST
      use run, only: UNITS

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED
      use error_manager
      use des_thermo, only : CALC_COND_DES
      use des_thermo_cond, only: DO_AREA_CORRECTION
      use discretelement
      use param1, only : one
      use physprop, only : mmax, ro_s0, d_p0
      use run, only : solids_model
      use constant, only : pi

      IMPLICIT NONE
!......................................................................!
      INTEGER :: M, L  ! Loop indices for DEM solids
! Calculate phase index offset for certain inputs until it can be
! addressed in other ways. should not matter unlesss hybrid
      CHARACTER(len=64) :: MSG
      INTEGER :: lent, lend, lenc, LC
      INTEGER :: MMAX_TOT
      LOGICAL :: ANY_CONDUCTION = .FALSE.

      DOUBLE PRECISION :: MASS_L, MASS_M, MASS_EFF
      DOUBLE PRECISION :: R_EFF, E_EFF
      
! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DEM_ENERGY")

! Conduction Equations:
!'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

! Set the default value for the minimum distance separating particles'
! surfaces.
      IF(DES_MIN_COND_DIST == UNDEFINED)THEN
         DES_MIN_COND_DIST = 1.0D-04 ! cm
         IF (UNITS == 'SI') DES_MIN_COND_DIST = &
            DES_MIN_COND_DIST/100.0  ! m
      ENDIF

! Setup code for conduction correction terms for artificial softening
      ! Calculate masses used for collision calculations.

! Shift the phase index for certain inputs to match the global phase
! index until this matter can be addressed otherwise (i.e., require
! the user specify correct indexing in mfix.dat). This should have no
! impact if not running a hybrid case
      MMAX_TOT = DES_MMAX+MMAX
      e_young_actual((MMAX+1):MMAX_TOT) = e_young_actual(1:DES_MMAX)
      v_poisson_actual((MMAX+1):MMAX_TOT) = v_poisson_actual(1:DES_MMAX)
      lent = MMAX_TOT+MMAX_TOT*(MMAX_TOT-1)/2
      lend = DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
      lenc = lent-lend
      DO_AREA_CORRECTION = .TRUE.
      DO M=MMAX+1,MMAX_TOT
         IF (SOLIDS_MODEL(M) /= 'DEM') CYCLE
         IF (.NOT.CALC_COND_DES(M))CYCLE
         ANY_CONDUCTION = .TRUE.
         IF(E_YOUNG_ACTUAL(M) == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"('Phase ',I2,' Actual EYoungs')") M
            WRITE(ERR_MSG,2002) 'E_YOUNG_ACTUAL', MSG
            DO_AREA_CORRECTION = .FALSE.
            CALL FLUSH_ERR_MSG
         ENDIF
         IF(V_POISSON_ACTUAL(M) == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"('Phase ',I2,' Actual poissons ratio')") M
            WRITE(ERR_MSG,2002) 'V_POISSON_ACTUAL', MSG
            DO_AREA_CORRECTION = .FALSE.
            CALL FLUSH_ERR_MSG
         ELSEIF(V_POISSON_ACTUAL(M) > 0.5d0 .OR. &
                V_POISSON_ACTUAL(M) <= -ONE) THEN
            WRITE(ERR_MSG,1001) trim(iVar('V_POISSON_ACTUAL',M)),  &
               iVal(V_POISSON_ACTUAL(M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO
      
      IF(ANY_CONDUCTION)THEN
         IF(EW_YOUNG_ACTUAL == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"('Actual EYoungs (wall)')")
            WRITE(ERR_MSG,2002) 'EW_YOUNG_ACTUAL', MSG
            DO_AREA_CORRECTION = .FALSE.
            CALL FLUSH_ERR_MSG
         ENDIF
         IF(VW_POISSON_ACTUAL == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"(' Actual poisson ratio (wall)')")
            WRITE(ERR_MSG,2002) 'VW_POISSON_ACTUAL', MSG
            DO_AREA_CORRECTION = .FALSE.
            CALL FLUSH_ERR_MSG
         ENDIF
      ENDIF

! see above index shift
      LC = LENC
      DO M=MMAX+1,MMAX_TOT
         IF(SOLIDS_MODEL(M) /='DEM') CYCLE
         IF(.NOT.DO_AREA_CORRECTION)CYCLE
         IF(.NOT.CALC_COND_DES(M)) CYCLE
! Calculate the mass of a phase M particle.
         MASS_M = (PI/6.d0)*(D_P0(M)**3)*RO_S0(M)

! Particle-Particle Collision Parameters ------------------------------>
         DO L=M,MMAX_TOT
            LC = LC + 1
            IF(.NOT.CALC_COND_DES(L)) CYCLE

            MASS_L = (PI/6.d0)*(D_P0(L)**3)*RO_S0(L)
            MASS_EFF = (MASS_M*MASS_L)/(MASS_M+MASS_L)

! Calculate the effective radius, Youngs modulus, and shear modulus.
            R_EFF = 0.5d0*(D_P0(M)*D_P0(L)/(D_P0(M) + D_P0(L)))
            E_EFF = E_YOUNG_ACTUAL(M)*E_YOUNG_ACTUAL(L) /  &
            &  (E_YOUNG_ACTUAL(M)*(1.d0 - V_POISSON_ACTUAL(L)**2) + &
            &   E_YOUNG_ACTUAL(L)*(1.d0 - V_POISSON_ACTUAL(M)**2))
 

! Calculate the spring properties and store in symmetric matrix format.
            HERT_KN_ACTUAL(M,L)=(4.d0/3.d0)*SQRT(R_EFF)*E_EFF

! Compute baseline for Hertzian collision time
            TAU_C_BASE_ACTUAL(M,L)=3.21D0*(MASS_Eff/HERT_KN_ACTUAL(M,L))**0.4 
            ! Can compute actual collision time via:
            !    TAU_C_ACTUAL = TAU_C_BASE_ACTUAL * (1/ImpactVel)^0.2

! Compute base for simulated collision time.  If Hertzian, only include base (without impact vel component).
            IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN)THEN
               TAU_C_BASE_SIM(M,L)=3.21D0*(MASS_Eff/HERT_KN(M,L))**0.4
            ELSE
               TAU_C_BASE_SIM(M,L)=PI/SQRT(KN/MASS_EFF - &
               ((DES_ETAN(M,L)/MASS_EFF)**2)/4.d0)
            ENDIF
         ENDDO

! Do particle-wall calculations
         MASS_EFF = MASS_M
         R_EFF = 0.5d0*D_P0(M)
         E_EFF = E_YOUNG_ACTUAL(M)*EW_YOUNG_ACTUAL /  &
         &  (E_YOUNG_ACTUAL(M)*(1.d0 - VW_POISSON_ACTUAL**2) + &
         &   EW_YOUNG_ACTUAL*(1.d0 - V_POISSON_ACTUAL(M)**2))
         
         HERT_KWN_ACTUAL(M) = (4.d0/3.d0)*SQRT(R_EFF)*E_EFF
         TAUW_C_BASE_ACTUAL(M) = 3.21D0 * (MASS_Eff/HERT_KWN_ACTUAL(M))**0.4 

         IF (DES_COLL_MODEL_ENUM .EQ. HERTZIAN)THEN
            TAUW_C_BASE_SIM(M)=3.21D0*(MASS_Eff/HERT_KWN(M))**0.4
         ELSE
            TAUW_C_BASE_SIM(M)=PI/SQRT(KN_w/MASS_EFF - &
            ((DES_ETAN_WALL(M)/MASS_EFF)**2)/4.d0)
         ENDIF

      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 1001 FORMAT('Warning 2001: Illegal or unknown input: ',A,' = ',A,/      &
      'Please correct the mfix.dat file.')
 2002 FORMAT('Warning 2002: Recommended input not specified: ',A,/          &
      'Description:',A,/'Not correcting contact area.')


      END SUBROUTINE CHECK_SOLIDS_DEM_ENERGY



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_SOLIDS_DEM_COHESION                              !
!  Author: J.Musser                                   Date: 11-DEC-13  !
!                                                                      !
!  Purpose: Check/set parameters for DES cohesion models.              !
!                                                                      !
!  Comments: Original code moved from CHECK_DES_DATA                   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_COHESION

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Invoke a cohesion model for DES simulation
      use discretelement, only: USE_COHESION
! Largest discrete particle diameter.
      use discretelement, only: MAX_RADIUS
! Runtime Flag: Invoke Square Well
      use discretelement, only: SQUARE_WELL
! Runtime Flag: Invoke Van der Waals model.
      use discretelement, only: VAN_DER_WAALS
! User specified parameter for increase neighbor search area.
      use discretelement, only: FACTOR_RLM
! User specified parameters for Van der Waals model.
      use discretelement, only: VDW_INNER_CUTOFF
      use discretelement, only: VDW_OUTER_CUTOFF
      use discretelement, only: HAMAKER_CONSTANT
      use discretelement, only: WALL_VDW_INNER_CUTOFF
      use discretelement, only: WALL_VDW_OUTER_CUTOFF
      use discretelement, only: WALL_HAMAKER_CONSTANT
      use discretelement, only: ASPERITIES
      use discretelement, only: SURFACE_ENERGY
      use discretelement, only: WALL_SURFACE_ENERGY
      use error_manager

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED
      use param1, only: ZERO
      use constant, only: Pi
      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Neighborhood size for Van der Waals force.
      DOUBLE PRECISION :: VDW_NEIGHBORHOOD
!......................................................................!


! Override the following settings if cohesion not used.
      IF(.NOT.USE_COHESION) THEN
!No more square well in the code
         SQUARE_WELL = .FALSE.
         VAN_DER_WAALS = .FALSE.
         WALL_VDW_OUTER_CUTOFF = ZERO
         RETURN
      ENDIF


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DEM_COHESION")


! Verify that only one cohesion model is specified.
      IF (SQUARE_WELL .AND. VAN_DER_WAALS) THEN
         WRITE(ERR_MSG,1002)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

! Verify that at a cohesion model is specified.
      ELSEIF(.NOT.SQUARE_WELL .AND. .NOT.VAN_DER_WAALS) THEN
         WRITE(ERR_MSG,1003)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


! Van der Waals model checks.
      IF (VAN_DER_WAALS) THEN

         IF (VDW_INNER_CUTOFF .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'VDW_INNER_CUTOFF'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(VDW_OUTER_CUTOFF .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'VDW_OUTER_CUTOFF'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'HAMAKER_CONSTANT'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (WALL_VDW_INNER_CUTOFF .EQ. UNDEFINED)THEN
            WRITE(ERR_MSG,1201) 'WALL_VDW_INNER_CUTOFF'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (WALL_VDW_OUTER_CUTOFF .EQ. UNDEFINED)THEN
            WRITE(ERR_MSG,1201) 'WALL_VDW_OUTER_CUTOFF'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(WALL_HAMAKER_CONSTANT .EQ. UNDEFINED) THEN
            WRITE(ERR_MSG,1201) 'WALL_HAMAKER_CONSTANT'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         VDW_NEIGHBORHOOD = 1.0d0 + (VDW_OUTER_CUTOFF/(2.d0*MAX_RADIUS))
         IF (FACTOR_RLM < VDW_NEIGHBORHOOD) THEN
            WRITE(ERR_MSG,1202)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF (ASPERITIES < ZERO) THEN
            WRITE(ERR_MSG,1001) 'ASPERITIES', trim(iVal(ASPERITIES))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         SURFACE_ENERGY=HAMAKER_CONSTANT/&
            (24.d0*Pi*VDW_INNER_CUTOFF**2)

         WALL_SURFACE_ENERGY=WALL_HAMAKER_CONSTANT/&
            (24.d0*Pi*WALL_VDW_INNER_CUTOFF**2)

      ENDIF


      CALL FINL_ERR_MSG

      RETURN

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A, ' = ',A,/     &
             'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1000: Cannot use SQUARE_WELL and VAN_DER_WAALS ',  &
         'cohesion',/'models simultaneously.')

 1003 FORMAT('Error 1001: A cohesion model was not selected. Specify ',&
         'one of the available models in the mfix.dat file.')


!<------------------- Van der Waals model messages. ----------------->!

 1201 FORMAT('Error 1201: Missing input data for Van der Waals ',      &
         'cohesion model.',/'Input parameter ',A,' is UNDEFINED.')

 1202 FORMAT('Error 1202: VDW_OUTER_CUTOFF outside of the neighbor ',  &
         'search distance.',/'Increase FACTOR_RLM to increase the ',   &
         'search distance.')

      END SUBROUTINE CHECK_SOLIDS_DEM_COHESION


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_DEM_COLLISION                              !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input data for DES collision calculations.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_COLLISION


! Global Variables:
!---------------------------------------------------------------------//

! User specified collision model
      USE discretelement, only: DES_COLL_MODEL
      USE discretelement, only: DES_COLL_MODEL_ENUM
      USE discretelement, only: LSD
      USE discretelement, only: HERTZIAN
! Particle and wall friction coeff.
      USE discretelement, only: MEW, MEW_W
! Parameter constatns.
      USE param1, only: ONE, ZERO, UNDEFINED

!      USE mpi_utility
      use error_manager

      IMPLICIT NONE

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DEM_COLLISION")

! Check coefficient friction
      IF(MEW == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'MEW'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF (MEW < ZERO .OR. MEW_W > ONE) THEN
         WRITE(ERR_MSG,1001) 'MEW', trim(iVal(MEW))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(MEW_W == UNDEFINED) THEN
         WRITE(ERR_MSG,1000) 'MEW_W'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF(MEW_w < ZERO .OR. MEW_W > ONE) THEN
         WRITE(ERR_MSG,1001) 'MEW_W', trim(iVal(MEW_W))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check collision model specific parameters.
      SELECT CASE (trim(DES_COLL_MODEL))
! Linear spring-dashpot model.
      CASE('LSD')
        DES_COLL_MODEL_ENUM = LSD
        CALL CHECK_SOLIDS_DEM_COLL_LSD
! Hertzian collision model.
      CASE('HERTZIAN')
         DES_COLL_MODEL_ENUM = HERTZIAN
         CALL CHECK_SOLIDS_DEM_COLL_HERTZ
! Unknown collision model.
      CASE DEFAULT
         WRITE(ERR_MSG,2000) TRIM(DES_COLL_MODEL)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      END SELECT

 2000 FORMAT('Error 2000: Invalid particle-particle collision model:',&
         A,/'Please correct the mfix.dat file.')


      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_DEM_COLLISION

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_DEM_COLL_LSD                               !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input data for DES collision calculations.      !
!                                                                      !
!  References:                                                         !
!   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
!   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
!      65-149 (pages 94-95)                                            !
!   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_COLL_LSD

! Modules
!---------------------------------------------------------------------//
      use constant, only: PI
! Number of discrete solids phases
      USE discretelement, only: DES_MMAX
! Particle and wall normal and tangential spring constants
      USE discretelement, only: KN, KN_W
      USE discretelement, only: KT, KT_W
! Particle and wall tangential spring factor := KN/KT
      USE discretelement, only: KT_FAC, KT_W_FAC

      use discretelement, only: DES_ETAN, DES_ETAN_WALL
      use discretelement, only: DES_ETAT, DES_ETAT_WALL

! Coefficients of restitution: Normal and Tangential
      USE discretelement, only: DES_EN_INPUT, DES_EN_WALL_INPUT
      USE discretelement, only: DES_ET_INPUT, DES_ET_WALL_INPUT

! Tangential damping factors := ET/EN
      USE discretelement, only: DES_ETAT_FAC, DES_ETAT_W_FAC

      use discretelement, only: DTSOLID

! Parameter constants.
      USE param1, only: ZERO, HALF, ONE, UNDEFINED
!
      use physprop, only: mmax, d_p0, ro_s0
      use run, only: solids_model
! USE mpi_utility
      use error_manager
      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M, L, LC, MMAX_TOT
! Calculate phase index offset for certain inputs until it can be
! addressed in other ways. should not matter unlesss hybrid
      INTEGER :: lent, lend, lenc
! Flag to warn user.
      LOGICAL :: FLAG_WARN
! Collision length scale.
      DOUBLE PRECISION :: TCOLL, TCOLL_TMP
! Collision length scale.
      DOUBLE PRECISION :: MASS_M, MASS_L, MASS_EFF
! Alias for coefficient restitution
      DOUBLE PRECISION :: EN
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DEM_COLL_LSD")

! Initialize.
      TCOLL = UNDEFINED

! Check for particle-particle normal spring constants.
      IF(KN == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) 'KN'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check for particle-particle tangential spring constant factors.
      IF(KT_FAC == UNDEFINED) THEN
         WRITE (ERR_MSG, 2100) 'KT_FAC'
         CALL FLUSH_ERR_MSG()
         KT_FAC = 2.0d0/7.0d0
      ELSEIF(KT_FAC > ONE .OR. KT_FAC < ZERO) THEN
         WRITE(ERR_MSG,1001) 'KT_FAC', trim(iVal(KT_FAC))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
! Calculate the particle-particle tangential spring factor.
      KT = KT_FAC*KN

! Check for particle-wall normal spring constants.
      IF(KN_W == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) 'KN_W'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check for particle-wall tangential spring constant factors.
      IF(KT_W_FAC == UNDEFINED) THEN
         WRITE (ERR_MSG, 2100) 'KT_W_FAC'
         CALL FLUSH_ERR_MSG()
         KT_W_FAC = 2.0d0/7.0d0
      ELSEIF(KT_W_FAC > ONE .OR. KT_W_FAC < ZERO) THEN
         WRITE(ERR_MSG,1001) 'KT_W_FAC', trim(iVal(KT_W_FAC))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
! Calculate the particle-wall tangential spring factor.
      KT_W = KT_W_FAC*KN_W

 2100 FORMAT('Warning 2100: Tangential spring factor ',A,' not ',      &
         'specified in mfix.dat.',/'Setting to default: (2/7).')

! Check for particle-particle tangential damping coefficients
      IF(DES_ETAT_FAC == UNDEFINED) THEN
         WRITE (ERR_MSG, 2101) 'DES_ETAT_FAC'
         CALL FLUSH_ERR_MSG
         DES_ETAT_FAC = HALF
      ELSEIF(DES_ETAT_FAC > ONE .OR. DES_ETAT_FAC < ZERO) THEN
         WRITE(ERR_MSG,1001) 'DES_ETAT_FAC', iVal(DES_ETAT_FAC)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Check for particle-wall tangential damping coefficients
      IF(DES_ETAT_W_FAC == UNDEFINED) THEN
         WRITE (ERR_MSG, 2101) 'DES_ETAT_W_FAC'
         CALL FLUSH_ERR_MSG
         DES_ETAT_W_FAC = HALF
      ELSEIF(DES_ETAT_W_FAC > ONE .OR. DES_ETAT_W_FAC < ZERO) THEN
         WRITE(ERR_MSG,1001) 'DES_ETAT_W_FAC', iVal(DES_ETAT_W_FAC)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 2101 FORMAT('Warning 2101: Tangential damping factor ',A,' not ', &
         'specified',/'in mfix.dat. Setting to default: (1/2).')

! Shift the phase index for certain inputs to match the global phase
! index until this matter can be addressed otherwise (i.e., require
! the user specify correct indexing in mfix.dat). This should have no
! impact if not running a hybrid case
      MMAX_TOT = DES_MMAX+MMAX
      des_en_wall_input((MMAX+1):MMAX_TOT) = des_en_wall_input(1:DES_MMAX)
      des_et_wall_input((MMAX+1):MMAX_TOT) = des_et_wall_input(1:DES_MMAX)
      lent = MMAX_TOT+MMAX_TOT*(MMAX_TOT-1)/2
      lend = DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
      lenc = lent-lend
      des_en_input((lenc+1):lent) = des_en_input(1:lend)
      des_et_input((lenc+1):lent) = des_et_input(1:lend)
      LC = lenc

      DO M = MMAX+1, MMAX_TOT
         IF (SOLIDS_MODEL(M) /= 'DEM') CYCLE

! Calculate the mass of a phase M particle.
         MASS_M = (PI/6.d0)*(D_P0(M)**3)*RO_S0(M)

! Particle-Particle Collision Parameters ------------------------------>
         DO L = M, MMAX_TOT
! not necessary to cycle here given requirements on ordering of solids input
            LC = LC+1

! Check particle-particle normal restitution coefficient
            IF(DES_EN_INPUT(LC) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('DES_EN_INPUT',LC))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(DES_EN_INPUT(LC) > ONE .OR. &
                   DES_EN_INPUT(LC) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('DES_EN_INPUT',LC)), &
                  trim(iVal(DES_EN_INPUT(LC)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            EN = DES_EN_INPUT(LC)

! Calculate masses used for collision calculations.
            MASS_L = (PI/6.d0)*(D_P0(L)**3)*RO_S0(L)
            MASS_EFF = MASS_M*MASS_L/(MASS_M+MASS_L)

! Calculate the M-L normal and tangential damping coefficients.
            IF(EN .NE. ZERO) THEN
               DES_ETAN(M,L) = 2.0D0*SQRT(KN*MASS_EFF) * ABS(LOG(EN))
               DES_ETAN(M,L) = DES_ETAN(M,L)/SQRT(PI*PI + (LOG(EN)**2))
            ELSE
               DES_ETAN(M,L) = 2.0D0*SQRT(KN*MASS_EFF)
            ENDIF
            DES_ETAT(M,L) = DES_ETAT_FAC*DES_ETAN(M,L)

! Store the entries in the symmetric matrix.
            DES_ETAN(L,M) = DES_ETAN(M,L)
            DES_ETAT(L,M) = DES_ETAT(M,L)

! Calculate the collision time scale.
            TCOLL_TMP = PI/SQRT(KN/MASS_EFF - &
               ((DES_ETAN(M,L)/MASS_EFF)**2)/4.d0)
            TCOLL = MIN(TCOLL_TMP, TCOLL)
         ENDDO


! Particle-Wall Collision Parameters ---------------------------------->
! Check particle-wall normal restitution coefficient.
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_EN_WALL_INPUT',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(DES_EN_WALL_INPUT(M) > ONE .OR. &
            DES_EN_WALL_INPUT(M) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('DES_EN_WALL_INPUT',M)), &
               trim(iVal(DES_EN_WALL_INPUT(M)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         EN = DES_EN_WALL_INPUT(M)

! Calculate masses used for collision calculations.
         MASS_EFF = MASS_M

! Calculate the M-Wall normal and tangential damping coefficients.
         IF(EN .NE. ZERO) THEN
            DES_ETAN_WALL(M) = 2.d0*SQRT(KN_W*MASS_EFF)*ABS(LOG(EN))
            DES_ETAN_WALL(M) = DES_ETAN_WALL(M)/SQRT(PI*PI+(LOG(EN))**2)
         ELSE
            DES_ETAN_WALL(M) = 2.D0*SQRT(KN_W*MASS_EFF)
         ENDIF
         DES_ETAT_WALL(M) = DES_ETAT_W_FAC*DES_ETAN_WALL(M)

! Calculate the collision time scale.
         TCOLL_TMP = PI/SQRT(KN_W/MASS_EFF -                           &
            ((DES_ETAN_WALL(M)/MASS_EFF)**2.d0)/4.d0)
!         TCOLL = MIN(TCOLL_TMP, TCOLL)
      ENDDO

! if following are assigned warn user they are discarded
      FLAG_WARN = .FALSE.
      LC = lenc
      DO M = MMAX+1, MMAX_TOT
         IF(SOLIDS_MODEL(M) /= 'DEM') CYCLE
         DO L = M, MMAX_TOT
            LC = LC + 1
            IF(DES_ET_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
         ENDDO
      ENDDO
      IF (FLAG_WARN) THEN
         WRITE(ERR_MSG,2102) 'DES_ET_INPUT'
         CALL FLUSH_ERR_MSG
      ENDIF

      FLAG_WARN = .FALSE.
      DO M = MMAX+1, MMAX_TOT
         IF (SOLIDS_MODEL(M) /= 'DEM') CYCLE
         IF(DES_ET_WALL_INPUT(M) .NE. UNDEFINED) FLAG_WARN = .TRUE.
      ENDDO
      IF (FLAG_WARN)THEN
         WRITE(ERR_MSG,2102) 'DES_ET_WALL_INPUT'
         CALL FLUSH_ERR_MSG
      ENDIF

 2102 FORMAT('Warning 2102: ',A,' values are not used ',/' with the',  &
         ' linear spring-dashpot collision model.')

! Store the smalled calculated collision time scale. This value is used
! in time-marching the DEM solids.
      DTSOLID = TCOLL/50.d0


      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_DEM_COLL_LSD

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_SOLIDS_DEM_COLL_HERTZ                             !
!  Author: J.Musser                                   Date: 11-Dec-13  !
!                                                                      !
!  Purpose: Check user input data for Hertzian collisions.             !
!                                                                      !
!  References:                                                         !
!   - Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13) !
!   -  Van der Hoef et al., Advances in Chemical Engineering, 2006, 31,!
!      65-149 (pages 94-95)                                            !
!   - Silbert et al., Physical Review E, 2001, 64, 051302 1-14 (page 5)!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_DEM_COLL_HERTZ

! Modules
!---------------------------------------------------------------------//
! The constant PI
      use constant, only: PI
! Number of discrete solids phases, diameters and densities
      USE discretelement, only: DES_MMAX
! User defined coefficients of restitution: Normal and Tangential
      USE discretelement, only: DES_EN_INPUT, DES_EN_WALL_INPUT
      USE discretelement, only: DES_ET_INPUT, DES_ET_WALL_INPUT
! Damping coefficients: Normal and Tangential
      use discretelement, only: DES_ETAN, DES_ETAN_WALL
      use discretelement, only: DES_ETAT, DES_ETAT_WALL
! Hertzian spring constants: Normal and Tangential
      use discretelement, only: HERT_KN, HERT_Kwn
      use discretelement, only: HERT_KT, HERT_Kwt
! Tangential damping factors := ET/EN
      USE discretelement, only: DES_ETAT_FAC, DES_ETAT_W_FAC
! Particle and wall Young's modulus and Shear modulus
      USE discretelement, only: E_YOUNG, Ew_YOUNG
! Particle and wall Poisson ratio
      USE discretelement, only: V_POISSON, Vw_POISSON
! Solids time step-size.
      use discretelement, only: DTSOLID
! Particle and wall normal spring constants
      USE discretelement, only: KN, KN_W
! Particle and wall tangential spring factor := KN/KT
      USE discretelement, only: KT_FAC, KT_W_FAC

      use param, only: DIM_M
! Parameter constatns.
      USE param1, only: ZERO, ONE, UNDEFINED

      USE physprop, only: mmax, ro_s0, d_p0
      USE run, only: solids_model
! USE mpi_utility
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Loop index.
      INTEGER :: M, L, LC, MMAX_TOT
! Calculate phase index offset for certain inputs until it can be
! addressed in other ways. should not matter unlesss hybrid
      INTEGER :: lent, lend, lenc
! Message for formatted output.
      CHARACTER(len=64) :: MSG
! Collision length scale.
      DOUBLE PRECISION :: TCOLL, TCOLL_TMP
! Particle and effective mass.
      DOUBLE PRECISION :: MASS_M, MASS_L, MASS_EFF
! Effective physical quantities. Radius, Youngs, Shear
      DOUBLE PRECISION :: R_EFF, E_EFF, G_MOD_EFF, RED_MASS_EFF
! Alias for coefficient restitution
      DOUBLE PRECISION :: EN, ET
! shear modulus
      DOUBLE PRECISION :: G_MOD(DIM_M)
! Shear modules for wall
      DOUBLE PRECISION :: G_MOD_WALL
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_DEM_COLL_HERTZ")

! Initialize.
      TCOLL = UNDEFINED

! check young's modulus and poisson ratio
      IF(Ew_YOUNG == UNDEFINED ) THEN
         MSG='Wall value for Youngs modulus'
         WRITE(ERR_MSG,1002) 'Ew_YOUNG', MSG
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(Vw_POISSON == UNDEFINED) THEN
         MSG='Wall value for Poissons ratio'
         WRITE(ERR_MSG,1002) 'Vw_POISSON', MSG
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ELSEIF (Vw_POISSON > 0.5d0 .OR. Vw_POISSON <= -ONE) THEN
         WRITE(ERR_MSG,1001) 'Vw_POISSON',iVal(Vw_POISSON)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      G_MOD_WALL = 0.5d0*Ew_YOUNG/(1.d0+Vw_POISSON)

! Shift the phase index for certain inputs to match the global phase
! index until this matter can be addressed otherwise (i.e., require
! the user specify correct indexing in mfix.dat). This should have no
! impact if not running a hybrid case
      MMAX_TOT = DES_MMAX+MMAX
      e_young((MMAX+1):MMAX_TOT) = e_young(1:DES_MMAX)
      v_poisson((MMAX+1):MMAX_TOT) = v_poisson(1:DES_MMAX)
      des_en_wall_input((MMAX+1):MMAX_TOT) = des_en_wall_input(1:DES_MMAX)
      des_et_wall_input((MMAX+1):MMAX_TOT) = des_et_wall_input(1:DES_MMAX)
      lent = MMAX_TOT+MMAX_TOT*(MMAX_TOT-1)/2
      lend = DES_MMAX+DES_MMAX*(DES_MMAX-1)/2
      lenc = lent-lend
      des_en_input((lenc+1):lent) = des_en_input(1:lend)
      des_et_input((lenc+1):lent) = des_et_input(1:lend)

      DO M=MMAX+1,MMAX_TOT
         IF (SOLIDS_MODEL(M) /= 'DEM') CYCLE

         IF(E_YOUNG(M) == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"('Phase ',I2,' Youngs modulus')") M
            WRITE(ERR_MSG,1002) 'E_YOUNG', MSG
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(V_POISSON(M) == UNDEFINED) THEN
            MSG=''; WRITE(MSG,"('Phase ',I2,' Poissons ratio')") M
            WRITE(ERR_MSG,1002) 'V_POISSON', MSG
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(V_POISSON(M) > 0.5d0 .OR. &
                V_POISSON(M) <= -ONE) THEN
            WRITE(ERR_MSG,1001) trim(iVar('V_POISSON',M)),  &
               iVal(V_POISSON(M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
! Calculate the shear modulus for phase M.
         G_MOD(M) = 0.5d0*E_YOUNG(M)/(1.d0+V_POISSON(M))
      ENDDO

! see above index shift
      LC = LENC
      DO M=MMAX+1,MMAX_TOT
         IF(SOLIDS_MODEL(M) /='DEM') CYCLE

! Calculate the mass of a phase M particle.
         MASS_M = (PI/6.d0)*(D_P0(M)**3)*RO_S0(M)

! Particle-Particle Collision Parameters ------------------------------>
         DO L=M,MMAX_TOT
! not necessary to cycle here given requirements on ordering of solids input
            LC = LC+1

! Check particle-particle normal restitution coefficient
            IF(DES_EN_INPUT(LC) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('DES_EN_INPUT',LC))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(DES_EN_INPUT(LC) > ONE .OR. &
               DES_EN_INPUT(LC) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('DES_EN_INPUT',LC)), &
                  trim(iVal(DES_EN_INPUT(LC)))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            EN = DES_EN_INPUT(LC)

! Check particle-particle tangential restitution coefficient
            IF(DES_ET_INPUT(LC) == UNDEFINED) THEN
               WRITE(ERR_MSG,1000) trim(iVar('DES_ET_INPUT',LC))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ELSEIF(DES_ET_INPUT(LC) > ONE .OR. &
                   DES_ET_INPUT(LC) < ZERO) THEN
               WRITE(ERR_MSG,1001) trim(iVar('DES_ET_INPUT',LC)), &
                  iVal(DES_ET_INPUT(LC))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
            ET = DES_ET_INPUT(LC)


! Calculate masses used for collision calculations.
            MASS_L = (PI/6.d0)*(D_P0(L)**3)*RO_S0(L)
            MASS_EFF = (MASS_M*MASS_L)/(MASS_M+MASS_L)
            RED_MASS_EFF = (2.d0/7.d0)*MASS_EFF
! Calculate the effective radius, Youngs modulus, and shear modulus.
            R_EFF = 0.5d0*(D_P0(M)*D_P0(L)/(D_P0(M) + D_P0(L)))
            E_EFF = E_YOUNG(M)*E_YOUNG(L) /  &
               (E_YOUNG(M)*(1.d0 - V_POISSON(L)**2) + &
                E_YOUNG(L)*(1.d0 - V_POISSON(M)**2))
            G_MOD_EFF = G_MOD(M)*G_MOD(L)/ &
               (G_MOD(M)*(2.d0 - V_POISSON(L)) + &
                G_MOD(L)*(2.d0 - V_POISSON(M)))

! Calculate the spring properties and store in symmetric matrix format.
            HERT_KN(M,L)=(4.d0/3.d0)*SQRT(R_EFF)*E_EFF
            HERT_KT(M,L)= 8.d0*SQRT(R_EFF)*G_MOD_EFF

            HERT_KN(L,M) = HERT_KN(M,L)
            HERT_KT(L,M) = HERT_KT(M,L)

! Calculate the normal damping coefficient.
            IF(EN .NE. ZERO) THEN
               DES_ETAN(M,L) = 2.d0*SQRT(HERT_KN(M,L)*MASS_EFF)* &
                  ABS(LOG(EN))
               DES_ETAN(M,L) = DES_ETAN(M,L)/ &
                  SQRT(PI*PI + (LOG(EN))**2)
            ELSE
               DES_ETAN(M,L) = 2.d0*SQRT(HERT_KN(M,L)*MASS_EFF)
            ENDIF
            DES_ETAN(L,M) = DES_ETAN(M,L)

! Calculate the tangential coefficients.
            IF(ET .NE. ZERO) THEN
               DES_ETAT(M,L) = 2.d0*SQRT(HERT_KT(M,L)*RED_MASS_EFF)* &
                  ABS(LOG(ET))
               DES_ETAT(M,L) = DES_ETAT(M,L)/ SQRT(PI*PI+(LOG(ET))**2)
            ELSE
               DES_ETAT(M,L) = 2.d0*SQRT(HERT_KT(M,L)*RED_MASS_EFF)
            ENDIF
            DES_ETAT(L,M) = DES_ETAT(M,L)

            TCOLL_TMP = PI/SQRT(HERT_KN(M,L)/MASS_EFF - &
               ((DES_ETAN(M,L)/MASS_EFF)**2)/4.d0)
            TCOLL = MIN(TCOLL_TMP, TCOLL)
         ENDDO

! Particle-Wall Collision Parameters ---------------------------------->
! Check particle-wall normal restitution coefficient.
         IF(DES_EN_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_EN_WALL_INPUT',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(DES_EN_WALL_INPUT(M) > ONE .OR. &
                DES_EN_WALL_INPUT(M) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('DES_EN_WALL_INPUT',M)), &
               trim(iVal(DES_EN_WALL_INPUT(M)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         EN = DES_EN_WALL_INPUT(M)

! Check particle-wall tangential restitution coefficient
         IF(DES_ET_WALL_INPUT(M) == UNDEFINED) THEN
            WRITE(ERR_MSG,1000) trim(iVar('DES_ET_WALL_INPUT',M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ELSEIF(DES_ET_WALL_INPUT(M) > ONE .OR. &
                DES_ET_WALL_INPUT(M) < ZERO) THEN
            WRITE(ERR_MSG,1001) trim(iVar('DES_ET_WALL_INPUT',M)), &
               trim(iVal(DES_ET_WALL_INPUT(M)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         ET = DES_ET_WALL_INPUT(M)

! Calculate masses used for collision calculations.
         MASS_EFF = MASS_M
         RED_MASS_EFF = (2.d0/7.d0)*MASS_EFF
! Calculate the effective radius, Youngs modulus, and shear modulus.
         R_EFF = 0.5d0*D_P0(M)
         E_EFF = E_YOUNG(M)*Ew_YOUNG /                                 &
            (E_YOUNG(M)*(1.d0-Vw_POISSON**2) +                         &
             Ew_YOUNG  *(1.d0-V_POISSON(M)**2))
         G_MOD_EFF = G_MOD(M)*G_MOD_WALL /                             &
            (G_MOD(M)*(2.d0 - Vw_POISSON) +                            &
             G_MOD_WALL*(2.d0 - V_POISSON(M)))

! Calculate the spring properties.
         HERT_Kwn(M) = (4.d0/3.d0)*SQRT(R_EFF)*E_EFF
         HERT_Kwt(M) = 8.0*SQRT(R_EFF)*G_MOD_EFF

! Calculate the tangential coefficients.
         IF(EN /= ZERO) THEN
            DES_ETAN_WALL(M) = 2.d0*SQRT(HERT_Kwn(M)*MASS_EFF)*&
               ABS(LOG(EN))
            DES_ETAN_WALL(M) = DES_ETAN_WALL(M)/&
               SQRT(PI*PI + (LOG(EN))**2)
         ELSE
            DES_ETAN_WALL(M) = 2.d0*SQRT(HERT_Kwn(M)*MASS_EFF)
         ENDIF

         IF(ET /= ZERO) THEN
            DES_ETAT_WALL(M) = 2.d0*SQRT(HERT_Kwt(M)*RED_MASS_EFF)*    &
                ABS(LOG(ET))
            DES_ETAT_WALL(M) = DES_ETAT_WALL(M)/SQRT(PI*PI+(LOG(ET))**2)
         ELSE
            DES_ETAT_WALL(M) = 2.d0*SQRT(HERT_Kwt(M)*RED_MASS_EFF)
         ENDIF

! Calculate the collision time scale.
         TCOLL_TMP = PI/SQRT(HERT_Kwn(M)/MASS_EFF -                    &
            ((DES_ETAN_WALL(M)/MASS_EFF)**2)/4.d0)
      ENDDO


! If following are assigned warn user they are discarded.
       IF(KN .NE. UNDEFINED) THEN
          WRITE(ERR_MSG, 2200) 'KN'
          CALL FLUSH_ERR_MSG
       ENDIF
       IF(KN_W .NE. UNDEFINED) THEN
          WRITE(ERR_MSG, 2200) 'KN_W'
          CALL FLUSH_ERR_MSG
       ENDIF
       IF(KT_FAC .NE. UNDEFINED) THEN
          WRITE(ERR_MSG, 2200) 'KT_FAC'
          CALL FLUSH_ERR_MSG
       ENDIF
       IF(KT_W_FAC .NE. UNDEFINED) THEN
          WRITE(ERR_MSG, 2200) 'KT_W_FAC'
          CALL FLUSH_ERR_MSG
       ENDIF
       IF(DES_ETAT_FAC .NE. UNDEFINED) THEN
          WRITE(ERR_MSG, 2200) 'DES_ETAT_FAC'
          CALL FLUSH_ERR_MSG
       ENDIF
       IF(DES_ETAT_W_FAC .NE. UNDEFINED) THEN
          WRITE(ERR_MSG, 2200) 'DES_ETAT_W_FAC'
          CALL FLUSH_ERR_MSG
       ENDIF

 2200 FORMAT('Warning 2200: ',A,' values are not used ',/' with the',  &
         ' linear spring-dashpot collision model.')


! Store the smalled calculated collision time scale. This value is used
! in time-marching the DEM solids.
      DTSOLID = TCOLL/50.d0


      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/      &
         'Please correct the mfix.dat file.')

 1002 FORMAT('Error 1002: Required input not specified: ',A,/          &
         'Description:',A,/'Please correct the mfix.dat file.')

      END SUBROUTINE CHECK_SOLIDS_DEM_COLL_HERTZ










