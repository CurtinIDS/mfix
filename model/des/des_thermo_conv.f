#include "version.inc"
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CONV_GS_DES1                                           !
!  Author: J.Musser: 16-Jun-10                                         !
!                                                                      !
!  Purpose: This routine is called from the DISCRETE side to calculate !
!  the gas-particle convective heat transfer.                          !
!                                                                      !
!  Comments: Explicitly coupled simulations use a stored convective    !
!  heat transfer coefficient. Otherwise, the convective heat transfer  !
!  coeff is calculated every time step and the total interphase energy !
!  transfered is 'stored' and used explictly in the gas phase. The     !
!  latter conserves all energy
!                                                                      !
!  REF: Zhou, Yu, and Zulli, "Particle scale study of heat transfer in !
!       packed and bubbling fluidized beds," AIChE Journal, Vol. 55,   !
!       no 4, pp 868-884, 2009.                                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CONV_GS_DES1

      use constant, only: Pi
      Use des_thermo
      Use discretelement
      Use fldvar
      Use interpolation
      Use param1
      use des_thermo, only: GAMMAxSA
      use geometry, only: NO_K
      use particle_filter, only: DES_INTERP_ON
      use particle_filter, only: FILTER_CELL
      use particle_filter, only: FILTER_WEIGHT
      use particle_filter, only: FILTER_SIZE
      use functions, only: FLUID_AT
      use functions, only: IS_NORMAL
      IMPLICIT NONE

      DOUBLE PRECISION :: lTg, GAMMA
      DOUBLE PRECISION :: Qcv_DT, Qcv
      DOUBLE PRECISION :: l4Pi
      INTEGER :: IJK, LC, NP

      l4Pi = 4.0d0*Pi

      DO NP=1,MAX_PIP
         IF(.NOT.IS_NORMAL(NP)) CYCLE

! Calculate the gas temperature.
         IF(DES_INTERP_ON) THEN
            lTg = ZERO
            DO LC=1,FILTER_SIZE
               IJK = FILTER_CELL(LC,NP)
               lTg = lTg + T_G(IJK)*FILTER_WEIGHT(LC,NP)
            ENDDO
         ELSE
            IJK = PIJK(NP,4)
            lTg = T_G(IJK)
         ENDIF

         IJK = PIJK(NP,4)

! Avoid convection calculations in cells without fluid (cut-cell)
         IF(.NOT.FLUID_AT(IJK)) THEN
            GAMMAxSA(NP) = ZERO
            CONV_Qs(NP) = ZERO

! For explicit coupling, use the heat transfer coefficient calculated
! for the gas phase heat transfer calculations.
         ELSEIF(DES_EXPLICITLY_COUPLED) THEN
            CONV_Qs(NP) = GAMMAxSA(NP)*(lTg - DES_T_s(NP))

         ELSE

! Calculate the heat transfer coefficient.
            CALL CALC_GAMMA_DES(NP, GAMMA)
            GAMMAxSA(NP) = GAMMA* l4Pi*DES_RADIUS(NP)*DES_RADIUS(NP)

! Calculate the rate of heat transfer to the particle
            Qcv = GAMMAxSA(NP)*(lTg - DES_T_s(NP))
! Store convection source in global energy source array.
            Q_Source(NP) = Q_Source(NP) + Qcv

! Calculate the gas phase source term components.
            Qcv_DT = Qcv*DTSOLID
            IF(DES_INTERP_ON) THEN
               DO LC=1,FILTER_SIZE
                  CONV_SC(IJK)=CONV_Sc(IJK)-Qcv_DT*FILTER_WEIGHT(LC,NP)
               ENDDO
            ELSE
               CONV_SC(IJK) = CONV_Sc(IJK) - Qcv_DT
            ENDIF
         ENDIF
      ENDDO

! Note that MPI sync is managed at the end of des_time_march for
! non-explicitly coupled cases that use interpolation.

      RETURN
      END SUBROUTINE CONV_GS_DES1


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CONV_GS_GAS1                                            !
!  Author: J.Musser                                   Date: 21-NOV-14  !
!                                                                      !
!                                                                      !
!  Purpose: This routine is called from the CONTINUUM. It calculates   !
!  the scalar cell center drag force acting on the fluid using         !
!  interpolated values for the gas velocity and volume fraction. The   !
!  The resulting sources are interpolated back to the fluid grid.      !
!                                                                      !
!  NOTE: The loop over particles includes ghost particles so that MPI  !
!  communications are needed to distribute overlapping force between   !
!  neighboring grid cells. This is possible because only cells "owned" !
!  by the current process will have non-zero weights.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CONV_GS_GAS1

! Flag: The fluid and discrete solids are explicitly coupled.
      use discretelement, only: DES_EXPLICITLY_COUPLED
! Size of particle array on this process.
      use discretelement, only: MAX_PIP
! Flag to use interpolation
      use particle_filter, only: DES_INTERP_ON
! Interpolation cells and weights
      use particle_filter, only: FILTER_CELL, FILTER_WEIGHT, FILTER_SIZE
! IJK of fluid cell containing particles center
      use discretelement, only: PIJK
! Particle temperature
      use des_thermo, only: DES_T_s
! Gas phase energy equation sources
      use des_thermo, only: CONV_Sp, CONV_Sc
      Use discretelement, only: DES_RADIUS
! Heat transfer coefficint (GAMMA) multiplied by sufrace area
      use des_thermo, only: GAMMAxSA
! Funtion for identifying fluid cells and normal particles.
      use functions, only: FLUID_AT
      use functions, only: IS_NORMAL
! MPI function for collecting interpolated data from ghost cells.
      use sendrecvnode, only: DES_COLLECT_gDATA
! MPI wrapper for halo exchange.
      use sendrecv, only: SEND_RECV

! Global Parameters:
!---------------------------------------------------------------------//
! Double precision values.
      use param1, only: ZERO, ONE
      use constant, only: Pi

      IMPLICIT NONE

! Loop counters: Particle, fluid cell, neighbor cells
      INTEGER :: NP, IJK, LC
! Interpolation weight
      DOUBLE PRECISION :: WEIGHT
      DOUBLE PRECISION :: GAMMAxSAxTp, GAMMA
      DOUBLE PRECISION :: l4Pi

      l4Pi = 4.0d0*Pi

! Initialize fluid cell values.
      CONV_Sc = ZERO
      CONV_Sp = ZERO

! Calculate the gas phase forces acting on each particle.
      DO NP=1,MAX_PIP

         IF(.NOT.IS_NORMAL(NP)) CYCLE
         IF(.NOT.FLUID_AT(PIJK(NP,4))) CYCLE

! Calculate the heat transfer coefficient.
         CALL CALC_GAMMA_DES(NP, GAMMA)

! Calculate the surface area of the particle

         GAMMAxSA(NP) = GAMMA*l4Pi*DES_RADIUS(NP)*DES_RADIUS(NP)
         GAMMAxSAxTp = GAMMAxSA(NP)*DES_T_s(NP)

         IF(DES_INTERP_ON) THEN
            DO LC=1,FILTER_SIZE
               IJK = FILTER_CELL(LC,NP)
               WEIGHT = FILTER_WEIGHT(LC,NP)

               CONV_Sc(IJK) = CONV_Sc(IJK) + WEIGHT*GAMMAxSAxTp
               CONV_Sp(IJK) = CONV_Sp(IJK) + WEIGHT*GAMMAxSA(NP)
            ENDDO
         ELSE
            IJK = PIJK(NP,4)

            CONV_Sc(IJK) = CONV_Sc(IJK) + GAMMAxSAxTp
            CONV_Sp(IJK) = CONV_Sp(IJK) + GAMMAxSA(NP)
         ENDIF

      ENDDO

! Add in data stored in ghost cells from interpolation. This call must
! preceed the SEND_RECV to avoid overwriting ghost cell data.
      IF(DES_INTERP_ON) THEN
         CALL DES_COLLECT_gDATA(CONV_SC)
         CALL DES_COLLECT_gDATA(CONV_SP)
      ENDIF

! Update the drag force and sources in ghost layers.
      CALL SEND_RECV(CONV_SC, 2)
      CALL SEND_RECV(CONV_SP, 2)

      RETURN
      END SUBROUTINE CONV_GS_GAS1




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: ZERO_ENERGY_SOURCE                                      !
!                                                                      !
!  Purpose: ZERO out the array that passes energy source terms back to !
!  the continuum model. Additional entries may be needed to include    !
!  heat transfer to the hybrid mode.                                   !
!                                                                      !
!  Author: J.Musser                                   Date: 15-Jan-11  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ZERO_ENERGY_SOURCE

      Use des_thermo
      Use param1

      IMPLICIT NONE

      CONV_Sc = ZERO
      CONV_Sp = ZERO

      RETURN
      END SUBROUTINE ZERO_ENERGY_SOURCE

