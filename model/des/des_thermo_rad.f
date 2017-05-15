!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_RADIATION                                          !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 25-Jun-10  !
!                                                                      !
!  Commen:                                                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_RADIATION

      USE constant
      USE des_thermo
      USE discretelement
      USE fldvar
      USE param1
      USE physprop, only: SMAX
      USE toleranc
      use functions, only: FLUID_AT
      use functions, only: IS_NORMAL
      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! Global index of particle
      INTEGER :: NP
! Solid phase of particle I
      INTEGER :: lM
! Fluid cell index containing particle I
      INTEGER :: IJK

! Local variables
!---------------------------------------------------------------------//
! Environment temperature
      DOUBLE PRECISION :: Tenv
! SB constant TIMES particle surface area
      DOUBLE PRECISION :: SBx4Pi
      INTEGER, PARAMETER :: lUpdateFreq=5
      INTEGER, SAVE :: lUpdate_avgTs=0
!......................................................................!


! Coupled simulations update the average solids temperature at the start
! of the DEM time march. Granular flow (no gas) simulations update the
! average solids temperature every "lUpdateFreq" time steps.
      IF(.NOT.DES_CONTINUUM_COUPLED) THEN
         IF(MOD(lUpdate_avgTs,lUpdateFreq) == 0) THEN
            CALL PARTICLES_IN_CELL
            CALL CALC_avgTs
            lUpdate_avgTs = 0
         ELSE
            lUpdate_avgTs = lUpdate_avgTs + 1
         ENDIF
      ENDIF

      SBx4Pi = SB_CONST*4.0d0*Pi

      DO NP=1,MAX_PIP
         IF(IS_NORMAL(NP)) THEN
            lM = PIJK(NP,5)
            IF(CALC_RADT_DES(lM)) THEN

               IJK = PIJK(NP,4)
               IF(FLUID_AT(IJK)) THEN
                  Tenv = EP_g(IJK)*T_g(IJK) + &
                     (ONE-EP_g(IJK))*avgDES_T_s(IJK)
               ELSE
                  Tenv = avgDES_T_s(IJK)
               ENDIF
! Calculate the heat source.
               Q_Source(NP) = Q_Source(NP) + DES_Em(lM)* SBx4Pi * &
                  (DES_RADIUS(NP)**2)*(Tenv**4 - DES_T_s(NP)**4)
! Update the thermal source term.
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE DES_RADIATION

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_avgTs                                             !
!  Author: J.Musser                                   Date: 06-NOV-12  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_avgTs

      USE compar
      USE derived_types, only: PIC
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE functions
      USE geometry
      USE indices
      USE param1
      USE physprop
      USE run, only: ENERGY_EQ

      IMPLICIT NONE

! Passed variables
!-----------------------------------------------
! NONE

! Local variables
!---------------------------------------------------------------------//
! Index of neighbor particle of particle I such that I < J
      INTEGER :: IJK
! Loop index for particles.
      INTEGER :: NP, lNP
! Sum of particle temperatures in fluid cell.
      DOUBLE PRECISION :: SUM_T_s
      INTEGER, SAVE :: PASS=0
!---------------------------------------------------------------------//

      IF(.NOT.ENERGY_EQ) RETURN

! Loop over fluid cells.
      IJK_LP: DO IJK = IJKSTART3, IJKEND3

         avgDES_T_s(IJK) = ZERO
         IF(.NOT.FLUID_AT(IJK)) CYCLE IJK_LP

         IF(PINC(IJK) > 0) THEN
! Initialize local solids temperature.
            SUM_T_s = ZERO
! Loop over all particles in cell IJK.
            lNP_LP: DO lNP = 1, PINC(IJK)
               NP = PIC(IJK)%p(lNP)
! Update the sum of particle temperatures in fluid cell IJK.
               IF(IS_NORMAL(NP)) SUM_T_s = SUM_T_s + DES_T_s(NP)
            ENDDO lNP_LP

! Average solids temperature in fluid cell IJK. The average method
! (over particles) will need changed for Hybrid model (mass? volume?).
            avgDES_T_s(IJK) = SUM_T_s/PINC(IJK)
         ENDIF
      ENDDO IJK_LP

      RETURN
      END SUBROUTINE CALC_avgTs
