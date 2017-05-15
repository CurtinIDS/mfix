!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose: Hook for user defined reaction rates.                      !
!                                                                      !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Comments: Write reaction rates in units of moles/sec.cm^3 (cgs) or  !
!  kmoles/sec.m^3 (SI). Units should match those specified in the data !
!  file.
!                                                                      !
!  Example reaction: Methane combustion                                !
!                                                                      !
!  mfix.dat input:                                                     !
!``````````````````````````````````````````````````````````````````````!
!    @(RXNS)                                                           !
!      CH4_Comb { chem_eq = "CH4 + 2.0*O2 --> CO2 + 2.0*H2O" }         !
!    @(END)                                                            !
!``````````````````````````````````````````````````````````````````````!
!                                                                      !
!  usr_rates.f input:                                                  !
!``````````````````````````````````````````````````````````````````````!
!    c_O2  = (RO_g(IJK)*X_g(IJK,O2)/MW_g(O2))                          !
!    c_CH4 = (RO_g(IJK)*X_g(IJK,CH4)/MW_g(CH4))                        !
!    RATES(CH4_Comb) = 2.0d5 * EP_g(IJK) * c_O2 * c_CH4                !
!``````````````````````````````````````````````````````````````````````!
!  * Species alias and reaction names given in the data file can be    !
!    used in reference to the reaction index in RATES and a species    !
!    index in gas/solids phase variables.                              !
!                                                                      !
!  * Additional information is provided in section 5.11 of the code    !
!    Readme.                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES(IJK, RATES)

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits
      USE compar
      USE sendrecv
      USE toleranc
      USE usr
      USE fun_avg
      USE functions

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJK

      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

      INCLUDE 'species.inc'
      INCLUDE 'usrnlst.inc'

! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
! Bounded phase temperatures (K)
      DOUBLE PRECISION :: TgX   ! gas phase
      DOUBLE PRECISION :: TsX   ! solids phase 1
      DOUBLE PRECISION :: TgsX  ! Average of gas/solids 1

! Molar concentration of solids species (mol/cm^3)
      DOUBLE PRECISION :: c_Biomass
      DOUBLE PRECISION :: c_Moisture

! Bound the gas and solids phase temperatures.
      TgX  = min(TMAX, T_g(IJK))
      TsX  = min(TMAX, T_s(IJK,1))

! Compute the gas/solids average bounded temperature.
      TgsX = HALF * (TgX + TsX)

! Compute concentration of moisture (gmole/cm^3)
      c_Moisture = ROP_s(IJK,1)*X_s(IJK,1,Moisture)/MW_s(1,Moisture)

! Compute concentration of carbon (gmole/cm^3)
      c_Biomass = ROP_s(IJK,1)*X_s(IJK,1,Biomass)/MW_s(1,Biomass)


! Reaction rates:
!`````````````````````````````````````````````````````````````````````//

      RATES(Drying)    = 5.13d6 * exp(-1.058d4/TsX) * c_Moisture

      RATES(Pyrolysis) = 1.30d8 * exp(-1.688d4/TsX) * c_Biomass

      RATES(Tarring)   = 2.00d8 * exp(-1.602d4/TsX) * c_Biomass
      RATES(Charring)  = 1.08d7 * exp(-1.460d4/TsX) * c_Biomass


!      RATES(:) = ZERO

      RETURN

      END SUBROUTINE USR_RATES
