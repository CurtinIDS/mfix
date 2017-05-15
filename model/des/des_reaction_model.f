!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_REACTION_MODEL                                     !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_REACTION_MODEL

      USE compar
      USE constant
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE geometry
      USE indices
      USE param, only: dimension_n_s
      USE param1, only: zero
      USE run, only: ANY_SPECIES_EQ, SPECIES_EQ
      USE physprop, only: NMAX
      USE run, only: DT
      USE run, only: SOLVE_ROs
      USE functions

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! None

! Local variables
!-----------------------------------------------
! Loop counter
      INTEGER :: NN
! total rate of consumption/production of species (g/sec)
      DOUBLE PRECISION  :: SUM_DES_Rs(1:MAX_PIP)

      DOUBLE PRECISION :: PIx4o3
      DOUBLE PRECISION :: o3 = 1.0d0/3.0d0

      DOUBLE PRECISION :: lDT, lOoDT
! Logical for Adams-Bashfort integration.
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

!---------------------------------------------------------------------//

      IF(.NOT.ANY_SPECIES_EQ) RETURN

      PIx4o3 = Pi*4.0d0/3.0d0

      lDT = merge(DT, DTSOLID, DES_EXPLICITLY_COUPLED)
      lOoDT = -1.0d0/lDT

! Bound the amount of mass loss.
      FORALL(NN=1:DIMENSION_N_S)
         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)         &
            DES_R_s(:MAX_PIP,NN) = max(DES_R_s(:MAX_PIP,NN),        &
            DES_X_s(:MAX_PIP,NN)*PMASS(:MAX_PIP)*lOoDT)
      END FORALL

! First-order method: Euler
      IF(INTG_EULER) THEN
         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)            &
            SUM_DES_Rs(:MAX_PIP) = sum(DES_R_s(:MAX_PIP,:),DIM=2)

         WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)            &
            PMASS(:MAX_PIP) = PMASS(:MAX_PIP) + lDT*SUM_DES_Rs(:MAX_PIP)

         FORALL(NN=1:DIMENSION_N_S)
            WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE)         &
               DES_X_s(:MAX_PIP,NN) = max(DES_X_s(:MAX_PIP,NN) + lDT*  &
               (DES_R_s(:MAX_PIP,NN) - DES_X_s(:MAX_PIP,NN)*           &
               SUM_DES_Rs(:MAX_PIP))/PMASS(:MAX_PIP), ZERO)
         END FORALL

      ELSE
         IF(FIRST_PASS) THEN
         ENDIF
      ENDIF

      DO NN=1,MAX_PIP
         IF(IS_NORMAL(NN)) THEN
            IF(SOLVE_ROs(PIJK(NN,5))) THEN
               RO_Sol(NN)= PMASS(NN)/PVOL(NN)
            ELSE
               DES_RADIUS(NN) = (PMASS(NN)/(Pix4o3*RO_SOL(NN)))**o3
               PVOL(NN) = PMASS(NN)/RO_SOL(NN)
            ENDIF
         ENDIF
      ENDDO

! Clear the necessary variables.
      DES_R_s = ZERO

! Flag that the first pass is over
      FIRST_PASS = .FALSE.

      RETURN
      END SUBROUTINE DES_REACTION_MODEL
