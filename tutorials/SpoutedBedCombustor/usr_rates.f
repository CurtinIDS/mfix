!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES(IJK, RATES)

      USE compar
      USE constant
      USE energy
      USE fldvar
      USE fun_avg
      USE functions
      USE funits
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run
      USE rxns
      USE sendrecv
      USE toleranc
      USE usr

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: IJK

      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

!-----------------------------------------------
      INCLUDE 'species.inc'
      INCLUDE 'usrnlst.inc'

! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
! Bounded phase temperatures (K)
      DOUBLE PRECISION TgX   ! gas phase
      DOUBLE PRECISION Ts1X  ! solids phase 1
      DOUBLE PRECISION Ts2X  ! solids phase 2

! Average of gas and solids phase (bounded) temperatures (K)
      DOUBLE PRECISION Tgs1X ! Average of gas/solids 1
      DOUBLE PRECISION Tgs2X ! Average of gas/solids 2

      DOUBLE PRECISION Pg_atm      ! Gas pressure (atm)
      DOUBLE PRECISION Pg_atmMW    ! Gas pressure * MW ()

! Gas phase species partial pressures (atm)
      DOUBLE PRECISION Pg_O2       ! Oxygen
      DOUBLE PRECISION Pg_CO2      ! Carbon dioxide
      DOUBLE PRECISION Pg_CO       ! Carbon monoxide
      DOUBLE PRECISION Pg_CO2_star ! (reverse reaction)
      DOUBLE PRECISION dPg_CO2     ! Pg_CO2 - Pg_CO2_star

! Gas phase concentrations (mol/cm^3)
      DOUBLE PRECISION c_O2        ! Oxygen
      DOUBLE PRECISION c_CO        ! Carbon monoxide
      DOUBLE PRECISION c_H2O       ! Water vapor (artificial)

! Solids volume fraction (aliases)
      DOUBLE PRECISION EP_s1 ! Solids phase 1
      DOUBLE PRECISION EP_s2 ! Solids phase 2

! Concentration of carbon (mol/cm^3)
      DOUBLE PRECISION c_FC1 ! Solids phase 1
      DOUBLE PRECISION c_FC2 ! Solids phase 2

! Reaction resistances (char combustion)
      DOUBLE PRECISION K_a   ! Ash layer diffusion
      DOUBLE PRECISION K_f   ! Gas film diffusion
      DOUBLE PRECISION K_r   ! Chemical reaction
      DOUBLE PRECISION K_eff ! Effective resistance

! Radio of particle (outer) diameter to unreacted core diameter.
      DOUBLE PRECISION R_D1  ! Solids phase 1
      DOUBLE PRECISION R_D2  ! Solids phase 2

      DOUBLE PRECISION R_O2, DIFF

      R_O2 = GAS_CONST*9.86923E-7/32.0  !cm^3.atm/g.K

! Bound the gas and solids phase temperatures.
      TGX   = min(TMAX, T_g(IJK))
      TS1X  = min(TMAX, T_s(IJK,1))
      TS2X  = min(TMAX, T_s(IJK,2))

! Compute the gas/solids average bounded temperature.
      TGS1X = HALF * (TGX + TS1X)
      TGS2X = HALF * (TGX + TS2X)

! Compute gas phase partial pressures (atm)
      Pg_atm    = P_g(IJK) / 1.013d6
      Pg_atmMW =  Pg_atm * MW_MIX_g(IJK)
      Pg_O2     = Pg_atmMW * X_g(IJK, O2) / MW_g(O2)
      Pg_CO     = Pg_atmMW * X_g(IJK, CO) / MW_g(CO)
      Pg_CO2    = Pg_atmMW * X_g(IJK, CO2) / MW_g(CO2)

! Alias the solids phase volume fractions
      EP_s1  = EP_s(IJK,1)
      EP_s2  = EP_s(IJK,2)

! Compute concentration of carbon (gmole/cm^3)
      c_FC1 = ROP_s(IJK,1) * X_s(IJK,1,FC1) / MW_s(1,FC1)
      c_FC2 = ROP_s(IJK,2) * X_s(IJK,2,FC2) / MW_s(2,FC2)


! Combustion: 2C + O2 --> 2CO  (g-mol/cm^3.sec)
!---------------------------------------------------------------------//
! Ref: Wen at al. (1982), Syamlal et al. (1993), Desai and Wen (1978)

      IF(Pg_O2 .GT. ZERO .AND. .NOT.COMPARE(EP_g(IJK), ONE)) THEN

! Calculate the diffusion coefficient for O2 in N2. Field, 1967.
         DIFF = 4.26d0 * ((TGX/1.8d3)**1.75d0) / Pg_atm

! Combustion_s1: 2FC1 + O2 --> 2CO  (g-mol/cm^3.sec)
!```````````````````````````````````````````````````````````````````````
! Calculate the ratio of outer diamter and unreacted core diameter.
         IF(X_s(IJK,1,Ash1) .GT. ZERO) THEN
            R_D1 = (X_s(IJK,1,FC1) * PAA / &
               (X_s(IJK,1,Ash1) * PAFC))**(1.0d0/3.0d0)
            R_D1 = MIN(ONE, R_D1)
         ELSE
            R_D1 = ZERO
         ENDIF
! If all of the carbon is consumed in solids phase 1, set the reaction
! rate to zero.

         IF(R_D1 .EQ. ZERO .OR. EP_s1 .EQ. ZERO) THEN
            RATES(Combustion_s1) = ZERO
         ELSE
! Calculate film layer resistance.
            K_f = DIFF * N_sh(IJK,1) / (D_p0(1) * R_O2 * TGX)
! Calculate kinetic resistance.
            K_r = 8.71d3 * EXP(-1.359d4/TS1X) * R_D1*R_D1
! Calculate ash layer resistance.
            IF(R_D1 .GE. ONE) THEN
! No ash layer resistance.
               K_eff = ONE/(ONE/K_f + ONE/K_r)
            ELSE
               K_a = 2.0d0 * DIFF * f_EP_A * R_D1 / &
                   (D_p0(1) * (ONE - R_D1) * R_O2 * TS1X)
               k_eff = ONE/(ONE/K_f + ONE/K_a + ONE/ K_r)
            ENDIF
! Combustion rate.
            RATES(Combustion_s1) = K_eff * Pg_O2 * 6.0 * EP_s1 / &
               (D_p0(1) * MW_g(O2))
         ENDIF

! Combustion_s2: 2FC2 + O2 --> 2CO  (g-mol/cm^3.sec)
!```````````````````````````````````````````````````````````````````````
! Calculate the ratio of outer diamter and unreacted core diameter.
         IF(X_s(IJK,2,Ash2) .GT. ZERO) THEN
            R_D2 = (X_s(IJK,2,FC2) * PAA / &
               (X_s(IJK,2,Ash2) * PAFC))**(1.0d0/3.0d0)
            R_D2 = MIN(ONE, R_D2)
         ELSE
            R_D2 = ZERO
         ENDIF
! If all of the carbon is consumed in solids phase 2, set the reaction
! rate to zero.
         IF(R_D2 .EQ. ZERO .OR. EP_s2 .EQ. ZERO) THEN
            RATES(Combustion_s2) = ZERO
         ELSE
! Calculate film layer resistance.
            K_f = DIFF * N_sh(IJK,2) / (D_p0(2) * R_O2 * TGX)
! Calculate kinetic resistance.
            K_r = 8.71d3 * EXP(-1.359d4/TS2X) * R_D2*R_D2
! Calculate ash layer resistance.
            IF(R_D2 .GE. ONE) THEN
! No ash layer resistance.
               K_eff = ONE/(ONE/K_f + ONE/K_r)
            ELSE
               K_a = 2.0d0 * DIFF * f_EP_A * R_D2 / &
                   (D_p0(2) * (ONE - R_D2) * R_O2 * TS2X)
               k_eff = ONE/(ONE/K_f + ONE/K_a + ONE/K_r)
            ENDIF
! Combustion rate.
            RATES(Combustion_s2) = K_eff * Pg_O2 * 6.0 * EP_s2 / &
               (D_p0(2) * MW_g(O2))
         ENDIF
      ENDIF



! CO2 Gasification: C + CO2 <--> 2CO  (g-mol/cm^3.sec)
!---------------------------------------------------------------------//
! Ref: Wen et al. (1982) rate was increased by a factor of 1000.

! Char_CO2_s1:  FC1 + CO2 --> 2CO
! Char_CO2_R: 2CO --> Soot + CO2 (catalyzed by FC1)
!```````````````````````````````````````````````````````````````````````
      IF(EP_s1 .GT. ZERO) THEN
         Pg_CO2_star = (Pg_CO * Pg_CO)/EXP (2.09238d1 - 2.02818d4/TGS1X)
         dPg_CO2 = Pg_CO2 - Pg_CO2_star
! Net Forward reaction.
         IF(dPg_CO2 .GT. ZERO) THEN
            RATES(Char_CO2_s1) = 9.3d2 * dPg_CO2 * &
               EXP(-2.265d4/TGS1X)*c_FC1
            RATES(Char_CO2_s1r) = ZERO
! Net Reverse reaction
         ELSE
            RATES(Char_CO2_s1) = ZERO
            RATES(Char_CO2_s1r) = 9.3d2 * (-dPg_CO2) * &
               EXP(-2.265d4/TGS1X)*c_FC1
         ENDIF
      ELSE
         RATES(Char_CO2_s1)  = ZERO
         RATES(Char_CO2_s1r) = ZERO
      ENDIF



! Char_CO2_s2:  2FC2 + O2 --> 2CO
! Char_CO2_s2r: 2CO --> Soot + CO2 (catalyzed by FC2)
!```````````````````````````````````````````````````````````````````````
      IF(EP_s2 .GT. ZERO) THEN
         Pg_CO2_star = (Pg_CO * Pg_CO)/EXP (2.09238d1 - 2.02818d4/TGS2X)
         dPg_CO2 = Pg_CO2 - Pg_CO2_star
! Net Forward reaction.
         IF(dPg_CO2 .GT. ZERO) THEN
            RATES(Char_CO2_s2) = 9.3d2 * dPg_CO2 * &
               EXP(-2.265d4/TGS2X)*c_FC2
            RATES(Char_CO2_s2r) = ZERO
! Net Reverse reaction
         ELSE
            RATES(Char_CO2_s2) = ZERO
            RATES(Char_CO2_s2r) = 9.3d2 * (-dPg_CO2) * &
               EXP(-2.265d4/TGS2X)*c_FC2
         ENDIF
      ELSE
         RATES(Char_CO2_s2)  = ZERO
         RATES(Char_CO2_s2r) = ZERO
      ENDIF



! CO Combustion: C0 + 0.5O2 --> CO2 (g-mol/cm^3.sec)
!---------------------------------------------------------------------//
! Ref: Westbrook and Dryer (1981)

      IF(Pg_CO .GT. ZERO .AND. Pg_O2 .GT. ZERO) THEN

! Combute gas phase concentrations of O2, CO, and water vapor.
         c_O2  = RO_g(IJK)*X_g(IJK,O2)/MW_g(O2)
         c_CO  = RO_g(IJK)*X_g(IJK,CO)/MW_g(CO)
!        c_H2O = RO_g(IJK)*X_g(IJK,H2O)/MW_g(H2O)

! Assume that the mass fraction of water vapor is 0.1 although it is
! not included as a gas phase species.
         c_H2O = RO_g(IJK)*0.1d0/1.8d1

         RATES(CO_Combustion) =  3.98d14 * EXP(-2.013d4/TGX) * &
            EP_g(IJK) * (c_O2)**0.25d0 * (c_CO) * (c_H2O)**0.5d0
      ELSE
         RATES(CO_Combustion) =  ZERO
      ENDIF


! Phase Change: C + Ash (phase 2) --> C + Ash (phase 2)
!---------------------------------------------------------------------//
! Note: Psuedo-reactions that convert coal to char when ash fraction
!       exceeds a threshold;  g/(cm^3.s)

       IF(EP_s2 .GT. ZERO .AND. X_s(IJK,2,Ash2) .GT. C(1) )THEN

! Char_to_Char: FC2 --> FC1
!```````````````````````````````````````````````````````````````````````
         RATES(Char_to_Char) = C(2) * c_FC2

! Ash_to_Ash: Ash2 --> Ash1
!```````````````````````````````````````````````````````````````````````
         RATES(Ash_to_Ash) = C(2) * ROP_s(IJK,2) * X_s(IJK,2,Ash2) /   &
            MW_s(2,Ash2)
      ELSE
         RATES(Char_to_Char) = ZERO
         RATES(Ash_to_Ash)   = ZERO
      ENDIF


!  Write the reaction rate data into .SPA file and for visualizing.
!---------------------------------------------------------------------//
      IF(nRR >= Combustion_s1) &!   2FC1 + O2 --> 2CO
         ReactionRates(IJK,Combustion_s1) = RATES(Combustion_s1)
      IF(nRR >= Combustion_s2) &!   2FC2 + O2 --> 2CO
         ReactionRates(IJK,Combustion_s2) = RATES(Combustion_s2)
      IF(nRR >= Char_CO2_s1)   &!   FC1 + CO2 --> 2CO
         ReactionRates(IJK,Char_CO2_s1) = RATES(Char_CO2_s1)
      IF(nRR >= Char_CO2_s1r)  &!   2CO + 0.FC1 --> Soot + CO2
         ReactionRates(IJK,Char_CO2_s1r) = RATES(Char_CO2_s1r)
      IF(nRR >= Char_CO2_s2)   &!   FC2 + CO2 --> 2CO
         ReactionRates(IJK,Char_CO2_s2) = RATES(Char_CO2_s2)
      IF(nRR >= Char_CO2_s2r)  &!   2CO + 0.FC2 --> Soot + CO2
         ReactionRates(IJK,Char_CO2_s2r) = RATES(Char_CO2_s2r)
      IF(nRR >= CO_Combustion) &!   CO + 0.5O2 --> CO2
         ReactionRates(IJK,CO_Combustion) = RATES(CO_Combustion)
      IF(nRR >= Char_to_Char)  &!   FC2 --> FC1
         ReactionRates(IJK,Char_to_Char) = RATES(Char_to_Char)
      IF(nRR >= Ash_to_Ash)    &!   Ash2 --> Ash1
         ReactionRates(IJK,Ash_to_Ash) = RATES(Ash_to_Ash)


      RETURN

      END SUBROUTINE USR_RATES
