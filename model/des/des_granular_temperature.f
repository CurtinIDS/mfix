!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DES_GRANULAR_TEMPERATURE                                C
!  Purpose: Calculate the DES granular temperature                     C
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!  Revision: For parallel processing modifications are made to         C
!            accomodate ghost cells  (Pradeep G.)                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_GRANULAR_TEMPERATURE

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE discretelement
      USE fldvar, only: u_s, v_s, w_s, theta_m
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE mpi_utility
      USE param
      USE param1
      USE physprop, only: mmax
      USE run
      USE sendrecv
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! indices
      INTEGER :: I, J, K, IJK
!
      INTEGER :: M, LL
! counter for no. of particles in phase m in cell ijk
      INTEGER :: NP_PHASE(DIMENSION_3, DIMENSION_M)
! temporary variable for mth phase granular temperature in cell ijk
      DOUBLE PRECISION :: TEMP(DIMENSION_3, DIMENSION_M)
! accounted for particles
      INTEGER :: PC
! squared particle velocity v.v
      DOUBLE PRECISION :: SQR_VEL, SQR_ROT_VEL
!-----------------------------------------------

! Calculate a local species granular temperature for current instant of
! time.  Note that the following calculation of species granular
! temperature employs a fluctuating particle velocity that is defined
! as the difference between a particles velocity and the corresponding
! local mean velocity of that particles species evaluated at the same
! instant of time.
!----------------------------------------------------------------->>>
! The following calculations are performed on the 'fluid' grid
      TEMP(:,:) = ZERO
      NP_PHASE(:,:) = 0
      PC = 0
      DO LL = 1, MAX_PIP
! skipping particles that do not exist
         IF(IS_NONEXISTENT(LL)) CYCLE
         PC = PC + 1
! skipping ghost particles
         IF(IS_GHOST(LL) .or. IS_ENTERING_GHOST(LL) .or. IS_EXITING_GHOST(LL)) CYCLE

         I = PIJK(LL,1)
         J = PIJK(LL,2)
         K = PIJK(LL,3)
         IJK = PIJK(LL,4)

         M = PIJK(LL,5)
         NP_PHASE(IJK,M) = NP_PHASE(IJK,M) + 1

         TEMP(IJK,M) = TEMP(IJK,M) + &
            (DES_VEL_NEW(LL,1)-U_s(IJK,M))**2
         TEMP(IJK,M) = TEMP(IJK,M) + &
            (DES_VEL_NEW(LL,2)-V_s(IJK,M))**2
         IF(DO_K) THEN
            TEMP(IJK,M) = TEMP(IJK,M) + &
               (DES_VEL_NEW(LL,3)-W_s(IJK,M))**2
         ENDIF

         IF(PC .EQ. PIP) EXIT
      ENDDO

! loop over all fluid cells
      DO IJK = IJKSTART3, IJKEND3
         IF(FLUID_AT(IJK)) THEN

            DO M = MMAX+1,DES_MMAX+MMAX
               IF (NP_PHASE(IJK,M) > 0 ) THEN
                  THETA_M(IJK,M) = TEMP(IJK,M)/&
                     DBLE(DIMN*NP_PHASE(IJK,M))
               ELSE
                  THETA_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDIF
      ENDDO


! Calculate global quantities: granular temperature, kinetic energy,
! potential energy and average velocity at the current instant of time
!----------------------------------------------------------------->>>

! initialization for calculations
      DES_KE = ZERO
      DES_ROTE = ZERO
      DES_PE = ZERO
      DES_VEL_AVG(:) = ZERO

! Calculate global average velocity, kinetic energy &
! potential energy
      PC = 0
      DO LL = 1, MAX_PIP
! skipping ghost particles and particles that don't exist
         IF(IS_NONEXISTENT(LL)) CYCLE
         PC = PC + 1
         IF(IS_GHOST(LL) .OR. IS_ENTERING_GHOST(LL) .OR. IS_EXITING_GHOST(LL)) CYCLE

         SQR_VEL = ZERO
         SQR_ROT_VEL = ZERO
         DO I = 1, DIMN
            SQR_VEL = SQR_VEL + DES_VEL_NEW(LL,I)**2
         ENDDO

         DO I = 1, merge(1,3,NO_K)
            SQR_ROT_VEL = SQR_ROT_VEL + OMEGA_NEW(LL,I)**2
         ENDDO


         DES_KE = DES_KE + PMASS(LL)/2.d0 * SQR_VEL
! Calculation of rotational kinetic energy
         DES_ROTE = DES_ROTE +                                         &
            (0.4D0*PMASS(LL)*DES_RADIUS(LL)**2)/2.d0 * SQR_ROT_VEL
         DES_PE = DES_PE + PMASS(LL)*DBLE(ABS(GRAV(2)))*&
            DES_POS_NEW(LL,2)
         DES_VEL_AVG(:) =  DES_VEL_AVG(:) + DES_VEL_NEW(LL,:)

         IF(PC .EQ. PIP) EXIT
      ENDDO

!Calculating total number of particles in the entire domain
!Needed for correct average velocities and granular temp.
      CALL GLOBAL_ALL_SUM(PIP-iGHOST_CNT,TOT_PAR)
      CALL GLOBAL_ALL_SUM(DES_VEL_AVG(1:DIMN))

!J.Musser changed PARTICLES TO PIP
      IF(TOT_PAR > 0) DES_VEL_AVG(:) = DES_VEL_AVG(:)/DBLE(TOT_PAR)

! The following quantities are primarily used for debugging/developing
! and allow a quick check of the energy conservation in the system.
! In their current form they are best applied to monodisperse cases.
! Calculate x,y,z components of global energy & granular temperature
      GLOBAL_GRAN_ENERGY(:) = ZERO
      GLOBAL_GRAN_TEMP(:)  = ZERO
      PC = 0
      DO LL = 1, MAX_PIP

! skipping ghost particles and particles that don't exist
         IF(IS_NONEXISTENT(LL)) CYCLE
         PC = PC + 1
         IF(IS_GHOST(LL) .OR. IS_ENTERING_GHOST(LL) .OR. IS_EXITING_GHOST(LL)) CYCLE

         GLOBAL_GRAN_ENERGY(:) = GLOBAL_GRAN_ENERGY(:) + &
            0.5d0*PMASS(LL)*(DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2
         GLOBAL_GRAN_TEMP(:) = GLOBAL_GRAN_TEMP(:) + &
            (DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2

         IF(PC .EQ. PIP) EXIT
      ENDDO

      CALL GLOBAL_ALL_SUM(GLOBAL_GRAN_TEMP)
      CALL GLOBAL_ALL_SUM(GLOBAL_GRAN_ENERGY)
      CALL GLOBAL_ALL_SUM(DES_KE)
      CALL GLOBAL_ALL_SUM(DES_PE)
      CALL GLOBAL_ALL_SUM(DES_ROTE)

      IF(TOT_PAR > 0) GLOBAL_GRAN_ENERGY(:) =                          &
         GLOBAL_GRAN_ENERGY(:)/DBLE(TOT_PAR)
      IF(TOT_PAR > 0) GLOBAL_GRAN_TEMP(:) =                            &
         GLOBAL_GRAN_TEMP(:)/DBLE(TOT_PAR)

      RETURN

      END SUBROUTINE DES_GRANULAR_TEMPERATURE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DES_BEDHEIGHT                                      !
!  Purpose: Calculate an average bed height for each solids phase      !
!           present                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE CALC_DES_BEDHEIGHT

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE des_bc
      USE discretelement
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE param, only: dimension_m
      USE param1, only: zero
      USE physprop
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle index
      INTEGER :: L
! accounted for particles
      INTEGER :: PC
! solids phase no.
      INTEGER :: M
! ijk indices
      INTEGER :: J, IJK
! average height of fluid cell
      DOUBLE PRECISION :: hcell
! height of particle (y-position)
      DOUBLE PRECISION :: hpart
! volume fraction of phase M in fluid cell
      DOUBLE PRECISION :: EP_SM
! tmp variables for calculations
      DOUBLE PRECISION :: tmp_num(DIMENSION_M), tmp_den(DIMENSION_M)
!-----------------------------------------------

! calculation of bed height following the formulation of Goldschmidt et
! al., Powder Technology 138, 2003, 135-159
      tmp_num(:) = ZERO
      tmp_den(:) = ZERO
      DO IJK = IJKSTART3, IJKEND3
         J = J_OF(IJK)
         DO M = MMAX+1, DES_MMAX+MMAX
            IF(ROP_S(IJK,M) > ZERO) THEN
               hcell = 0.5d0*(YN(J)+YN(J-1))
               EP_SM = EP_S(IJK,M)
               tmp_num(M) = tmp_num(M) + EP_SM*hcell*VOL(IJK)
               tmp_den(M) = tmp_den(M) + EP_SM*VOL(IJK)
            ENDIF
         ENDDO
      ENDDO
! calculate avg height for each phase
      IF (PIP >0) bed_height(:) = tmp_num(:)/tmp_den(:)

! alternative method to calculating bed height (turned off atm)
      IF(.FALSE.) THEN
      tmp_num(:) = ZERO
      tmp_den(:) = ZERO

      PC = 0
      DO L = 1, MAX_PIP
! skipping ghost particles and particles that don't exist
         IF(IS_NONEXISTENT(L)) CYCLE
         PC = PC + 1
         IF(IS_GHOST(L) .OR. IS_ENTERING_GHOST(L) .OR. IS_EXITING_GHOST(L)) CYCLE

         M = PIJK(L,5)
         hpart = DES_POS_NEW(L,2)
         tmp_num(M) = tmp_num(M) + hpart
         tmp_den(M) = tmp_den(M) + 1
         IF(PC .EQ. PIP) EXIT
      ENDDO
      ! calculate avg height for each phase
      IF (PIP >0) bed_height(:) = tmp_num(:)/tmp_den(:)
      ENDIF

      RETURN

      END SUBROUTINE CALC_DES_BEDHEIGHT
