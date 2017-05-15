!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose: Hook for user defined reaction rates.                      !
!                                                                      !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Comments: Write reaction rates in units of moles/sec (cgs and SI).  !
!                                                                      !
!  WARNING: Only discrete phase reactions should be specifed here.     !
!  Homogeneous gas phase reactions in DEM simulations must be given    !
!  in the continuum reaction hook (usr_rates.f).                       !
!                                                                      !
!  The call to usr_rates_des is made from inside a particle loop which !
!  is nested inside an IJK loop. Fluid grid calculations independent   !
!  of particle properties can be carried out in des/usr4_des.f to      !
!  reduce redundant calculations.                                      !
!                                                                      !
!  Example: Evaporation                                                !
!                                                                      !
!  mfix.dat input:                                                     !
!``````````````````````````````````````````````````````````````````````!
!    @(DES_RXNS)                                                       !
!      Evap { chem_eq = "Liquid --> Vapor" }                           !
!    @(DES_END)                                                        !
!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!  des/usr_rates_des.f input:                                          !
!``````````````````````````````````````````````````````````````````````!
!    Sa = Pi*(2.0d0*DES_RADIUS(NP))**2    ! Particle surface area      !
!    H2O_xfr = ...  ! An expression for mass transfer coefficient      !
!    Cmg_H2O = ...  ! Molar concentration grad of water vapor          !
!    DES_RATES(EVAP) = Sa * H2O_xfr * Cmg_H2O                          !
!``````````````````````````````````````````````````````````````````````!
!  * Species alias and reaction names given in the data file can be    !
!    used in reference to the reaction index in DES_RATES and a        !
!    species index in gas/solids phase variables.                      !
!                                                                      !
!  * Additional information is provided in section 4.11 of the code    !
!    Readme.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES_DES(NP, pM, IJK, DES_RATES)

      USE compar
      USE constant
      USE des_thermo
      USE des_rxns
      USE discretelement
      USE energy
      USE fldvar
      USE funits
      USE geometry
      USE indices
      USE param
      USE param1
      USE physprop
      USE rxns
      USE run
      USE usr
      USE fun_avg
      USE functions

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NP  ! Global index of particle
      INTEGER, INTENT(IN) :: pM  ! Solid phase index of particle NP
      INTEGER, INTENT(IN) :: IJK ! Fluid cell index containing NP

! Calculated reaction rates. (reacted moles per sec)
      DOUBLE PRECISION, INTENT(OUT) :: DES_RATES(NO_OF_DES_RXNS)


! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//

      INCLUDE 'species.inc'

! Reaction rates:
!`````````````````````````````````````````````````````````````````````//
! Include reaction rates here. Reaction rates should be stored in the
! variable DES_RATES. The reaction name given in the data file can be
! used to store the rate in the appropriate array location. Additional
! input format parameters are given in Section 4.11 of the code Readme.

      DES_RATES(:) = ZERO

      RETURN

      END SUBROUTINE USR_RATES_DES
