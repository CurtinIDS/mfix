!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose: Hook for user defined reaction rates.                      !
!                                                                      !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Comments: Write reaction rates in units of moles/sec.cm^3 (cgs) or  !
!  modles/sec.m^3 (SI). Units should match those specified in the data !
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

      INTEGER, INTENT(IN) :: NP, pM, IJK
      DOUBLE PRECISION, INTENT(OUT) :: DES_RATES(NO_OF_DES_RXNS)


! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
! Concentration of water vapor at droplet surface (mol-gH2O/cm^3)
      DOUBLE PRECISION Cm_H2O

! Mole fraction of water vapor in gas phase.  (mol-gH2O/mol-MIX)
      DOUBLE PRECISION Y_H2O
! Concentration of water vapor in bluk gas phase (mol-gH2O/cm^3)
      DOUBLE PRECISION Cg_H2O

! Difference in water vapor molar concentration between saturated vapor
! at the solids surface and bulk gas phase.
      DOUBLE PRECISION Cmg_H2O

! Diffusion coefficient of water vapor in Air (cm^2/sec)
      DOUBLE PRECISION Diff_H2OinAIR

! Dimensionless numbers:
      DOUBLE PRECISION N_Re ! Reynolds
      DOUBLE PRECISION N_Sc ! Schmidt
      DOUBLE PRECISION N_Sh ! Sherwood

! Water vapor mass transfer coefficient (cm/sec)
      DOUBLE PRECISION H2O_xfr

! Total surface area in CV  (cm^2/cm^3  --> 1/cm)
      DOUBLE PRECISION Sa

      DOUBLE PRECISION Tp, Dp0

      INCLUDE 'species.inc'
      INCLUDE 'usrnlst.inc'

! Alias particle temperature.
      Tp = DES_T_s(NP)
! Calculate particle diameter.
      Dp0 = 2.0d0 * DES_RADIUS(NP)

! Reaction rates:
!`````````````````````````````````````````````````````````````````````//

! Concentration of water vapor at droplet surface (mol-Vapor/cm^3)
      Cm_H2O = Psat_H2O(Tp) / (8314.56d4 * T_G(IJK))
! Mole fraction of water vapor in gas phase.  (mol-Vapor/mol-MIX)
      Y_H2O = X_g(IJK,Vapor) * (MW_MIX_g(IJK)/MW_g(Vapor))
! Concentration of water vapor in bluk gas phase (mol-Vapor/cm^3)
      Cg_H2O = (Y_H2O * P_G(IJK)) / (8314.56d4 * T_G(IJK))
! Water vapor concentration gradient (mol-Vapor/cm^3)
      Cmg_H2O = Cm_H2O - Cg_H2O

      IF(Cmg_H2O > ZERO) THEN
! Diffusion coefficient of water vapor in air. (cm^2/sec)
         Diff_H2OinAIR = Diff_H2O_Air(T_G(IJK), P_G(IJK))
! Reynolds Number
         N_Re = cal_NRe(1)
! Schmidt Number
         N_Sc = cal_NSc(Diff_H2OinAIR)
! Sherwood Number (Ranz and Marshal, 1952)
         N_Sh = cal_NSh(N_Re, N_Sc)
! Mass transfer coefficient (cm/sec)
         H2O_xfr = (N_Sh * Diff_H2OinAIR) / Dp0
! Solids phase surface area (cm^2)
         Sa = Pi * Dp0 * Dp0
! Mass transfer rate. (moles/sec)
         DES_RATES(Evaporation) = Sa * H2O_xfr * Cmg_H2O

      ELSE
         DES_RATES(Evaporation) = ZERO
      ENDIF

      RETURN

      CONTAINS
!----------------------------------------------------------------------!
! Function: Psat_H2O                                                   !
!                                                                      !
! Purpose: Calculate the saturation pressure of water vapor as a       !
! function of temperature. The saturation pressure is returned in      !
! Pascals for SI units and Baryes for CGS units.                       !
!                                                                      !
! REF: Reynolds, W.C., "Thermodynamic properties in SI," Dept. of      !
! Mechanical Engineer, Stanford University. (1979)                     !
!                                                                      !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION Psat_H2O(sat_T)

! Temperature at which to calculate the saturation pressure.
      DOUBLE PRECISION, INTENT(IN) :: sat_T
! Algebraic coefficients for equation
      DOUBLE PRECISION, PARAMETER, DIMENSION(8) :: COEFFS = &
         (/-7.419242d0, 2.97221d-1,-1.155286d-1, 8.68563d-3, &
            1.094098d-3,-4.39993d-3, 2.520658d-3, -5.218684d-4/)
! Reference temperature (K)
      DOUBLE PRECISION, PARAMETER :: ref_T = 338.15d0
! Reference pressure (Pascals)
      DOUBLE PRECISION, PARAMETER :: ref_P = 0.01d0
! Critical point temperature. (K)
      DOUBLE PRECISION, PARAMETER :: CP_T = 647.286d0
! Critical point pressure. (Pascals)
      DOUBLE PRECISION, PARAMETER :: CP_P = 22.089d6

! Local temporary variables.
      DOUBLE PRECISION val, sum1
! Generic loop index
      INTEGER I

! Initialize variables
      val = ref_P * (sat_T - ref_T)
      sum1 = ZERO
      DO I = 1, 8
         sum1 = sum1 + COEFFS(I)*((val)**(I-1))
      ENDDO
      val = sum1 * ((CP_T/sat_T) - 1.0d0)
! Calculate the saturation pressure (Pascals)
      Psat_H2O = CP_P * exp(val)
! Convert to the necessary units.
      Psat_H2O = 10.0d0 * Psat_H2O ! (barye)

      RETURN
      END FUNCTION Psat_H2O


!----------------------------------------------------------------------!
! Function: Diff_H2O_Air                                               !
!                                                                      !
! Purpose: Calculate the diffusion coefficient of water vapor in air.  !
!                                                                      !
! The diffusion coefficient is returned in m^2/sec for SI units and    !
! cm^2/sec for CGS units.                                              !
!                                                                      !
! REF: Mills, A.F., "Basic Heat and Mass Transfer," 2nd Edition,       !
! Prentice Hall, pg. 947 (1998)                                        !
!                                                                      !
! Original Source:                                                     !
! Marreo, T.R., and Masion, E.A., "Gaseous diffusion coefficients,"    !
! Journal of Physical and Chemical Reference Data, 1, 3-118 (1972)     !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION Diff_H2O_Air(ref_T, ref_P)

! Reference temperature and pressure
      DOUBLE PRECISION, INTENT(IN) :: ref_T, ref_P
! Reference pressure in atm
      DOUBLE PRECISION ref_Patm


      ref_Patm = ref_P/1.01325d6

      IF( ref_T < 450.0d0) THEN
         Diff_H2O_Air = (1.87d-6/ref_Patm) * (ref_T**2.072d0)
      ELSE
         Diff_H2O_Air = (27.5d-6/ref_Patm) * (ref_T**1.632d0)
      ENDIF

      RETURN
      END FUNCTION Diff_H2O_Air

!----------------------------------------------------------------------!
! Function: calc_NRe(M)                                                !
!                                                                      !
! Purpose: Calculate the Reynolds number.                              !
!                                                                      !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION cal_NRe(M)

      INTEGER M

! Various fluid cell indicies
      INTEGER I, IMJK, IJMK, IJKM
! Gas Velocity - cell centered
      DOUBLE PRECISION UGC, VGC, WGC
! Solids Velocity - cell centered
      DOUBLE PRECISION USCM, VSCM, WSCM
! Relative velocity.
      DOUBLE PRECISION VREL

! Initialize fluid cell variables
      I =  I_OF(IJK)
      IMJK  = IM_OF(IJK)
      IJMK  = JM_OF(IJK)
      IJKM  = KM_OF(IJK)

! Calculate velocity components at i, j, k
! Gas
      UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
      VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
      WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))
! Solids
      USCM = DES_VEL_NEW(NP,1)
      VSCM = DES_VEL_NEW(NP,2)
      WSCM = DES_VEL_NEW(NP,3)

! magnitude of gas-solids relative velocity
      VREL = SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + (WGC - WSCM)**2)

! Reynods Number
      IF(MU_g(IJK) > ZERO) THEN
         cal_NRe = 2.0d0 * DES_RADIUS(NP) * VREL * RO_g(IJK) / MU_g(IJK)
      ELSE
         cal_NRe = LARGE_NUMBER
      ENDIF

      RETURN
      END FUNCTION cal_NRe

!----------------------------------------------------------------------!
! Function: cal_NSc                                                    !
!                                                                      !
! Purpose: Calculate the Schmidt Number.                               !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION cal_NSc(Diff_Coeff)
! Diffustion coefficient
      DOUBLE PRECISION, intent(IN) :: Diff_Coeff

! Schmidt Number
      cal_NSc = MU_g(IJK)/(RO_g(IJK)*Diff_Coeff)

      RETURN
      END FUNCTION cal_NSc


!----------------------------------------------------------------------!
! Function: cal_NSh                                                    !
!                                                                      !
! Purpose: Calculate the Sherwood Number.                              !
!''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''!
      DOUBLE PRECISION FUNCTION cal_NSh(lRe, lSc)

! Reynolds and Schmidt Numbers.
      DOUBLE PRECISION, intent(IN) :: lRe, lSc

! Sherwood Number: Ranz and Marshall, 1952
      cal_NSh = 2.0d0 + 0.60d0*(lRe**(1.0d0/2.0d0))*(lSc**(1.0d0/3.0d0))

      RETURN

      END FUNCTION cal_NSh

      END SUBROUTINE USR_RATES_DES
