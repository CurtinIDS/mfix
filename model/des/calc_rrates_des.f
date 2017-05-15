!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CALC_RRATES_DES                                        !
!  Author: J. Musser                                  Date: 10-Oct-12  !
!                                                                      !
!  Purpose: Calculate reaction rates for various reactions for DES.    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_RRATES_DES(NP, pRgp, pRgc, pRPhase, pHoRg, pSUMRg)

      USE compar
      USE constant
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE energy
      USE fldvar
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
      USE usr
      Use parse
      use functions
      use toleranc, only: ZERO_X_gs, COMPARE
      use physprop, only: NMAX
      use param1, only: DIMENSION_LM
      use param, only: DIMENSION_M

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
      INTEGER, INTENT(IN) :: NP   ! particle index

! Local gas phase values.
      DOUBLE PRECISION, INTENT(OUT) :: pRgp(NMAX(0)) ! Rate of production
      DOUBLE PRECISION, INTENT(OUT) :: pRgc(NMAX(0)) ! Rate of consumption
      DOUBLE PRECISION, INTENT(OUT) :: pHoRg   ! Heat of reaction
      DOUBLE PRECISION, INTENT(OUT) :: pSUMRg
      DOUBLE PRECISION, INTENT(OUT) :: pRPhase(DIMENSION_LM+DIMENSION_M-1)

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: pM   ! Global Solids Phase index
      INTEGER :: IJK  ! fluid cell index

      INTEGER :: H    ! Reaction loop counter
      INTEGER :: M    ! Global Phase index loop counter
      INTEGER :: NN   ! Global species index
      INTEGER :: lN   ! Local reaction speices index/loop counter
      INTEGER :: LM   !

      INTEGER :: mXfr ! Global phase index for mass transfer

! User-defined reaction rates returned from USR_RATES
      DOUBLE PRECISION :: DES_RATES(NO_OF_DES_RXNS)

      DOUBLE PRECISION :: lRate
      DOUBLE PRECISION :: lTp
      DOUBLE PRECISION :: lHoRs, llHoRg

      DOUBLE PRECISION :: RxH

! Reaction limiters. If a species mass fraction is less than this
! value, then the reaction is suppressed.
      DOUBLE PRECISION :: speciesLimiter

! External functions
!---------------------------------------------------------------------//
! Enthalpy calculations (cal/gram)
      DOUBLE PRECISION, EXTERNAL :: CALC_H

! Alias particle temperature.
      lTp = DES_T_s(NP)
      IJK = PIJK(NP,4)
      pM = PIJK(NP,5)

! Initialize storage arrays
      pRgp(:) = ZERO
      pRgc(:) = ZERO
      pHoRg = ZERO

      pRPhase = ZERO
      pSUMRg = ZERO

! Set the species limiter:
      speciesLimiter = ZERO_X_gs

! Calculate user defined reaction rates.
      DES_RATES(:) = ZERO
      CALL USR_RATES_DES(NP, pM, IJK, DES_RATES)

! Loop over reactions.
      RXN_LP: DO H = 1, NO_OF_DES_RXNS

! Skip empty reactions
         IF(DES_Reaction(H)%nSpecies == 0) CYCLE RXN_LP
         IF(COMPARE(DES_RATES(H),ZERO)) CYCLE RXN_LP

! Initialize local loop arrays
         llHoRg = ZERO
         lHoRs = ZERO
         RxH = ZERO

! Calculate the rate of formation/consumption for each species.
!---------------------------------------------------------------------//
         DO lN = 1, DES_Reaction(H)%nSpecies
! Global phase index.
            M = DES_Reaction(H)%Species(lN)%pMap
! Global species index.
            NN = DES_Reaction(H)%Species(lN)%sMap
! Index for interphase mass transfer. For a gas/solid reaction, the
! index is stored with the gas phase.
            mXfr = DES_Reaction(H)%Species(lN)%mXfr
            lRate = DES_RATES(H) * DES_Reaction(H)%Species(lN)%MWxStoich
! Gas Phase:
            IF(M == 0) THEN
! Consumption of gas phase species.
               IF(lRate < ZERO) THEN
                  IF(X_g(IJK,NN) > speciesLimiter) THEN
                     pRgc(NN) = pRgc(NN) - lRate
! Enthalpy transfer associated with mass transfer. (gas/solid)
                     IF(M /= mXfr) RxH = RxH +                         &
                        lRate*CALC_H(T_g(IJK),0,NN)
                  ELSE
! There is an insignificant amount of reactant. Skip this reaction.
                     DES_RATES(H) = ZERO
                     CYCLE RXN_LP
                  ENDIF
               ELSE
! Formation of gas phase species.
                  pRgp(nn) = pRgp(nn) + lRate
! Enthalpy transfer associated with mass transfer. (gas/solid)
                  IF(M /= mXfr) RxH = RxH + lRate*CALC_H(lTp,0,NN)
               ENDIF
! Discrete Solids Phase:
            ELSE
! Formation/consumption of solids phase species.
               DES_R_s(NP,nn) = DES_R_s(NP,nn) + lRate
            ENDIF
         ENDDO ! Loop of species


! Calculate and store the heat of reaction.
!---------------------------------------------------------------------//
         IF(ENERGY_EQ) THEN
! Automated heat of reaction calculations
            IF(DES_Reaction(H)%Calc_DH) THEN
! Loop over reaction species.
               DO lN = 1, DES_Reaction(H)%nSpecies
! Global phase index.
                  M = DES_Reaction(H)%Species(lN)%pMap
! Global species index.
                  NN = DES_Reaction(H)%Species(lN)%sMap
! Rate of formation/consumption for speices N
                  lRate = DES_RATES(H) * &
                     DES_Reaction(H)%Species(lN)%MWxStoich
! Gas phase enthalpy chnage from energy equation derivation.
                  IF(M == 0) THEN
                     llHORg = llHORg + CALC_H(T_g(IJK),0,nn) * lRate
! Solid phase enthalpy change from energy equation derivation.
                  ELSE
                     lHORs = lHORs + CALC_H(lTp,M,nn) * lRate
                  ENDIF
               ENDDO

! Apply enthalpy transfer associated with mass transfer to get the
! complete heat of reaction for Reaction H.
               llHORg = llHORg - RxH
               lHORs = lHORs + RxH

! Convert the heat of reaction to the appropriate units (if SI), and
! store in the global array.
               IF(UNITS == 'SI') THEN
                  pHoRg = pHoRg + 4.183925d3*llHORg
                  Q_Source(NP) = Q_Source(NP) - 4.183925d3*lHORs
               ELSE
                  pHORg = pHORg + llHORg
                  Q_Source(NP) = Q_Source(NP) - lHORs
               ENDIF
            ELSE
! User-defined heat of reaction.
               HOR_g(IJK) = HOR_g(IJK) + &
                  DES_Reaction(H)%HoR(0) * DES_RATES(H)
               Q_Source(NP) = Q_Source(NP) - &
                  DES_Reaction(H)%HoR(pM) * DES_RATES(H)
            ENDIF
         ENDIF

! Update rate of interphase mass transfer.
!---------------------------------------------------------------------//
         LM = 1 + (pM - 1)*pM/2
         pRPhase(LM) = pRPhase(LM) + &
            DES_RATES(H) * DES_Reaction(H)%rPHASE(LM)

      ENDDO RXN_LP ! Loop over reactions.

! Calculate the toal rate of formation and consumption for each species.
!---------------------------------------------------------------------//
      IF(SPECIES_EQ(0)) THEN
         pSUMRg = SUM(pRgp(:NMAX(0)) - pRgc(:NMAX(0)))
      ELSE
         DO H=1, NO_OF_DES_RXNS
            IF(DES_Reaction(H)%nPhases <= 0) CYCLE
            LM = 1 + ((pM-1)*pM)/2
            pSUMRg = pSUMRg + &
               DES_RATES(H) * DES_Reaction(H)%rPHASE(LM)
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CALC_RRATES_DES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine name: CALC_RRATE_DES                                     !
!                                                                      !
!  Purpose: This routine manages gas-solid reactions for the continuum !
!  phase.                                                              !
!                                                                      !
!  Author: J.Musser                                   Date: 16-May-11  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ZERO_RRATE_DES

      USE des_rxns
      USE param1, only: zero

      IMPLICIT NONE

      DES_R_gp(:,:) = ZERO
      DES_R_gc(:,:) = ZERO
      DES_R_PHASE(:,:) = ZERO
      DES_HOR_G(:) = ZERO
      DES_SUM_R_g(:) = ZERO

      RETURN
      END SUBROUTINE ZERO_RRATE_DES
