!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_THERMO_NEWVALUES                                   !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_THERMO_NEWVALUES

      USE compar
      Use des_thermo
      Use des_rxns
      Use discretelement
      USE geometry
      USE indices
      Use param1
      Use physprop
      use run, only: ENERGY_EQ
      use functions
      use funits, only: dmp_log
      use run, only: ANY_SPECIES_EQ
      USE des_thermo_cond, only: DES_QW_cond
      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! NONE

! Local variables
!---------------------------------------------------------------------//
! Index of neighbor particle of particle I such that I < J
      INTEGER IJK
! Loop index for particles.
      INTEGER NP, lNP
! Logical for Adams-Bashfort integration.
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
!---------------------------------------------------------------------//

      IF(.NOT.ENERGY_EQ) RETURN

! Second-order Adams-Bashforth scheme defaults to Euler on first pass.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE) &
            Q_Source0(:MAX_PIP) = Q_Source(:MAX_PIP)/       &
            (PMASS(:MAX_PIP)*DES_C_ps(:MAX_PIP))
      ENDIF
      FIRST_PASS = .FALSE.

! First-order method
      IF (INTG_EULER) THEN
         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)   &
            DES_T_s(:MAX_PIP) = DES_T_s(:MAX_PIP) +   &
            DTSOLID*(Q_Source(:MAX_PIP)/(PMASS(:MAX_PIP)*     &
            DES_C_ps(:MAX_PIP)))

! Second-order Adams-Bashforth scheme
      ELSE
         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)
            DES_T_s(:MAX_PIP) = DES_T_s(:MAX_PIP) + DTSOLID *         &
               (1.5d0*Q_Source(:MAX_PIP) -0.5d0*Q_Source0(:MAX_PIP))/ &
               (PMASS(:MAX_PIP)*DES_C_ps(:MAX_PIP))
            Q_Source0(:MAX_PIP) = Q_Source(:MAX_PIP)
         ENDWHERE
      ENDIF


      Q_Source(:) = ZERO
      IF(ALLOCATED(DES_QW_Cond)) &
         DES_QW_Cond(:,:) = ZERO

! Update particle from reactive chemistry process.
      IF(ANY_SPECIES_EQ .AND. .NOT.DES_EXPLICITLY_COUPLED)&
         CALL DES_REACTION_MODEL

      RETURN

      END SUBROUTINE DES_THERMO_NEWVALUES
