!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_PHASE_INDEX                                         !
!                                                                      !
!  Purpose: Set the index of all particles based on their diameter and !
!  density.                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_PHASE_INDEX

! Modules
!---------------------------------------------------------------------//
      use discretelement, only: PIJK
      USE discretelement, only: DES_RADIUS, RO_SOL
      USE discretelement, only: DES_MMAX
      USE discretelement, only: MAX_PIP
      USE functions, only: IS_NONEXISTENT, IS_GHOST, IS_ENTERING_GHOST
      USE functions, only: IS_EXITING_GHOST
      USE error_manager
      use mpi_funs_des, only: des_par_exchange
      use mpi_utility
      use param1, only: small_number
      USE physprop, only: MMAX, D_p0, RO_s0
      USE run, only: RUN_TYPE, solids_model
      USE run, only: ANY_SPECIES_EQ
      use sendrecv
      IMPLICIT NONE

! Local Variables
!---------------------------------------------------------------------//
! particle no.
      INTEGER :: L
! solids phase no.
      INTEGER :: M
! May need to offset index when using d_p0 and ro_s
      INTEGER :: DM
! IER for error reporting
      INTEGER :: IER
! Difference between a particles diameter (density) and the diameter
! (density) of a phase specified in the data file.
      DOUBLE PRECISION dDp, dRho
!......................................................................!


! The restart file contains the phase index for reacting cases as the
! diameter and/or density of the particle may have changed.
      IF(RUN_TYPE /= 'NEW' .AND. ANY_SPECIES_EQ) RETURN

! Initialize the error flag.
      IER = 0

! solids phase index of particle.
! ---------------------------------------------------------------->>>
      DO L = 1, MAX_PIP
         IF(IS_NONEXISTENT(L)) CYCLE
         IF(IS_GHOST(L) .OR. IS_ENTERING_GHOST(L) .OR. IS_EXITING_GHOST(L)) CYCLE

! Determining the solids phase of each particle by matching the diameter
! and density to those specified in the data file.
         M_LP: DO M = MMAX+1, MMAX+DES_MMAX
            dDp  = ABS(2.0d0*DES_RADIUS(L)-D_P0(M))
            dRho = ABS( RO_Sol(L)-RO_S0(M))
            IF( dDp < SMALL_NUMBER .AND. dRho < SMALL_NUMBER) THEN
               PIJK(L,5) = M
               EXIT M_LP
            ENDIF
         ENDDO M_LP
! Flag error if no match is found.
         IF(PIJK(L,5).EQ.0) IER = 1
      ENDDO

! Sync up the error flag across all processes.
      CALL GLOBAL_ALL_SUM(IER)
      IF(IER == 0) RETURN

! Point of no return: Report errors and abort
!----------------------------------------------------------------------
      CALL INIT_ERR_MSG("SET_PHASE_INDEX")

      CALL OPEN_PE_LOG(IER)

      WRITE(ERR_MSG, 1100)
      CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)

 1100 FORMAT('Error 1100: Unable to determine the phase of one or ',&
         'more particles.',/8x,'ID',4X,'Diameter',6x,'Density',/)

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(IS_NONEXISTENT(L)) CYCLE
         IF(IS_GHOST(L) .OR. IS_ENTERING_GHOST(L) .OR. IS_EXITING_GHOST(L)) CYCLE

! Flag as an error if no match is found.
         IF(PIJK(L,5).EQ.0) THEN
            WRITE(ERR_MSG,9000) L,  2.0*DES_RADIUS(L), Ro_Sol(L)
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ENDDO

      WRITE(ERR_MSG, 1101)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1101 FORMAT(' ',/'Defined phase parameters from mfix.dat:',/3x,'ID',&
         5X,'Diameter',5x,'Density')

      DO M = MMAX+1, DES_MMAX+MMAX
         WRITE(ERR_MSG, 9000) M, D_P0(M), RO_S0(M)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDDO

      WRITE(ERR_MSG, 1102)
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., ABORT=.TRUE.)

 1102 FORMAT('Please correct the mfix.dat or particle_input.dat files.')

 9000 FORMAT(I10,2(2x,g12.5))

      END SUBROUTINE SET_PHASE_INDEX
