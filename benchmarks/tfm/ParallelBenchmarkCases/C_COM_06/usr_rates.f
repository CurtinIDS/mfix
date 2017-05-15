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
      DOUBLE PRECISION :: TgX  ! gas phase
      DOUBLE PRECISION :: TSX  ! solids phase 1

! Average of gas and solids phase (bounded) temperatures (K)
      DOUBLE PRECISION :: TgsX ! Average of gas/solids 1

      DOUBLE PRECISION :: Pg_atm      ! Gas pressure (atm)
      DOUBLE PRECISION :: Pg_atmMW    ! Gas pressure * MW ()

! Gas phase species partial pressures (atm)
      DOUBLE PRECISION :: Pg_O2       ! Oxygen
      DOUBLE PRECISION :: Pg_CO2      ! Carbon dioxide
      DOUBLE PRECISION :: Pg_CO       ! Carbon monoxide
      DOUBLE PRECISION :: Pg_CO2_star ! (reverse reaction)
      DOUBLE PRECISION :: dPg_CO2     ! Pg_CO2 - Pg_CO2_star

! Gas phase concentrations (mol/cm^3)
      DOUBLE PRECISION :: c_O2        ! Oxygen
      DOUBLE PRECISION :: c_CO        ! Carbon monoxide
      DOUBLE PRECISION :: c_H2O       ! Water vapor (artificial)

! Solids volume fraction (aliases)
      DOUBLE PRECISION :: EPs ! Solids phase 1

! Concentration of carbon (mol/cm^3)
      DOUBLE PRECISION :: c_FC1 ! Solids phase 1

! Reaction resistances (char combustion)
      DOUBLE PRECISION :: K_a   ! Ash layer diffusion
      DOUBLE PRECISION :: K_f   ! Gas film diffusion
      DOUBLE PRECISION :: K_r   ! Chemical reaction
      DOUBLE PRECISION :: K_eff ! Effective resistance

! Radio of particle (outer) diameter to unreacted core diameter.
      DOUBLE PRECISION :: R_D1  ! Solids phase 1
      DOUBLE PRECISION :: R_D2  ! Solids phase 2

      DOUBLE PRECISION :: R_O2, DIFF

! Generic loop variable
      INTEGER :: L

      R_O2 = GAS_CONST*9.86923E-7/32.0  !cm^3.atm/g.K

! Bound the gas and solids phase temperatures.
      TGX   = min(TMAX, T_g(IJK))
      TSX  = min(TMAX, T_s(IJK,1))

! Compute the gas/solids average bounded temperature.
      TgsX = HALF * (TGX + TSX)

! Compute gas phase partial pressures (atm)
      Pg_atm    = P_g(IJK) / 1.013d6
      Pg_atmMW =  Pg_atm * MW_MIX_g(IJK)
      Pg_O2     = Pg_atmMW * X_g(IJK, O2) / MW_g(O2)
      Pg_CO     = Pg_atmMW * X_g(IJK, CO) / MW_g(CO)
      Pg_CO2    = Pg_atmMW * X_g(IJK, CO2) / MW_g(CO2)

! Alias the solids phase volume fractions
      EPs  = EP_s(IJK,1)

! Compute concentration of carbon (gmole/cm^3)
      c_FC1 = ROP_s(IJK,1) * X_s(IJK,1,FC1) / MW_s(1,FC1)


! Combustion: 2C + O2 --> 2CO  (g-mol/cm^3.sec)
!---------------------------------------------------------------------//
! Ref: Wen at al. (1982), Syamlal et al. (1993), Desai and Wen (1978)

      IF(Pg_O2 .GT. ZERO .AND. .NOT.COMPARE(EP_g(IJK), ONE)) THEN

! Calculate the diffusion coefficient for O2 in N2. Field, 1967.
         DIFF = 4.26d0 * ((TGX/1.8d3)**1.75d0) / Pg_atm

! Combustion: 2FC1 + O2 --> 2CO  (g-mol/cm^3.sec)
!```````````````````````````````````````````````````````````````````````
! Calculate the ratio of outer diamter and unreacted core diameter.
         IF(X_s(IJK,1,ASH) .GT. ZERO) THEN
            R_D1 = (X_s(IJK,1,FC1) * PAA / &
               (X_s(IJK,1,ASH) * PAFC))**(1.0d0/3.0d0)
            R_D1 = MIN(ONE, R_D1)
         ELSE
            R_D1 = ZERO
         ENDIF
! If all of the carbon is consumed in solids phase 1, set the reaction
! rate to zero.

         IF(R_D1 .EQ. ZERO .OR. EPs .EQ. ZERO) THEN
            RATES(CHAR_COMBUSTION) = ZERO
         ELSE
! Calculate film layer resistance.
            K_f = DIFF * N_sh(IJK,1) / (D_p0(1) * R_O2 * TGX)
! Calculate kinetic resistance.
            K_r = 8.71d3 * EXP(-1.359d4/TSX) * R_D1*R_D1
! Calculate ash layer resistance.
            IF(R_D1 .GE. ONE) THEN
! No ash layer resistance.
               K_eff = ONE/(ONE/K_f + ONE/K_r)
            ELSE
               K_a = 2.0d0 * DIFF * f_EP_A * R_D1 / &
                   (D_p0(1) * (ONE - R_D1) * R_O2 * TSX)
               k_eff = ONE/(ONE/K_f + ONE/K_a + ONE/ K_r)
            ENDIF
! Combustion rate.
            RATES(CHAR_COMBUSTION) = K_eff * Pg_O2 * 6.0 * EPs /     &
               (D_p0(1) * MW_g(O2))
         ENDIF
      ENDIF

! CO2 Gasification: C + CO2 <--> 2CO  (g-mol/cm^3.sec)
!---------------------------------------------------------------------//
! Ref: Wen et al. (1982) rate was increased by a factor of 1000.

! CHAR_CO2:  FC1 + CO2 --> 2CO
! Char_CO2r: 2CO --> Soot + CO2 (catalyzed by FC1)
!```````````````````````````````````````````````````````````````````````
      IF(EPs .GT. ZERO) THEN
         Pg_CO2_star = (Pg_CO * Pg_CO)/EXP (2.09238d1 - 2.02818d4/TgsX)
         dPg_CO2 = Pg_CO2 - Pg_CO2_star
! Net Forward reaction.
         IF(dPg_CO2 .GT. ZERO) THEN
            RATES(CHAR_CO2) = 9.3d2 * dPg_CO2 * &
               EXP(-2.265d4/TgsX)*c_FC1
            RATES(CHAR_CO2r) = ZERO
! Net Reverse reaction
         ELSE
            RATES(CHAR_CO2) = ZERO
            RATES(CHAR_CO2r) = 9.3d2 * (-dPg_CO2) * &
               EXP(-2.265d4/TgsX)*c_FC1
         ENDIF
      ELSE
         RATES(CHAR_CO2)  = ZERO
         RATES(CHAR_CO2r) = ZERO
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

         RATES(CO_COMBUSTION) =  3.98d14 * EXP(-2.013d4/TGX) * &
            EP_g(IJK) * (c_O2)**0.25d0 * (c_CO) * (c_H2O)**0.5d0
      ELSE
         RATES(CO_Combustion) =  ZERO
      ENDIF


!  Write the reaction rate data into .SPA file and for visualizing.
!---------------------------------------------------------------------//
      DO L=1,min(nRR, NO_OF_RXNS)
         ReactionRates(IJK,L) = RATES(L)
      ENDDO

      RETURN

      END SUBROUTINE USR_RATES
