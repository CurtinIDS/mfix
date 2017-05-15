!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES0(IER)                                           C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!           using information from the data file                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 3-10-98    C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose:Replaced routines with new proceedures for automated        C
!          reaction rate calculations.                                 C
!  Author: J. Musser                                  Date: 10-Oct-12  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX, IJK, T_g, T_s1, D_p, X_g, X_s, EP_g,    C
!            P_g, HOR_g, HOR_s                                         C
!                                                                      C
!                                                                      C
!  Variables modified: M, N, R_gp, R_sp, RoX_gc, RoX_sc, SUM_R_g,      C
!                      SUM_R_s                                         C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE RRATES0()

! Global domain parameters.
!`````````````````````````````````````````````````````````````````````//
! Blanket access is provided to the following modules as they don't
! contain variables which are likely to be modified 'by accident.'
      use compar
      use indices
      use geometry

! Gas phase global variables:
!`````````````````````````````````````````````````````````````````````//
! Species mass fractions.
      use fldvar, only: X_g
! Temperature.
      use fldvar, only: T_g
! Rate of production.
      use rxns, only: R_gp
! Rate of consumption divided by species mass fraction.
      use rxns, only: RoX_gc
! Net rate of change of gas phase mass.
      use rxns, only: SUM_R_g
! Heat of reaction source term
      use energy, only: HOR_g

! Gas phase global variables:
!`````````````````````````````````````````````````````````````````````//
! Species mass fractions.
      use fldvar, only: X_s
! Temperatures.
      use fldvar, only: T_s
! Rate of production.
      use rxns, only: R_sp
! Rate of consumption divided by species mass fraction.
      use rxns, only: RoX_sc
! Net rate of change of gas phase mass.
      use rxns, only: SUM_R_s
! Heat of reaction source term
      use energy, only: HOR_s


! Reaction global variables:
!`````````````````````````````````````````````````````````````````````//
! Number of TFM reactions.
      use rxns, only: NO_OF_RXNS
! Rate of interphase mass transfer.
      use rxns, only: R_phase
! Data structure storing reaction information.
      use rxns, only: Reaction


! Global run-time parameters.
!`````````````````````````````````````````````````````````````````````//
! Number of solids phases (or one if MMAX is zero)
      use param, only: DIMENSION_M
! Number of gas species.
      use param, only: DIMENSION_N_g
! Number of solids phase species.
      use param, only: DIMENSION_N_s
      use param1, only: dimension_lm
      use param1, only: zero
! Indicates that the energy equations are solved.
      use run, only: ENERGY_EQ
! Indicates that the species equations are solved.
      use run, only: SPECIES_EQ
! Simulation units (CGS or SI)
      use run, only: UNITS
! Number of solids phases.
      use physprop, only: MMAX
! Number of species comprising each phase.
      use physprop, only: NMAX
! Indicate that DEM is used.
      use discretelement, only: DISCRETE_ELEMENT
! Indicate that TFM/DEM hybrid is used.
      use discretelement, only: DES_CONTINUUM_HYBRID
! Small value for species mass fractions
      use toleranc
      use functions

      implicit none

! Local variables:
!`````````````````````````````````````````````````````````````````````//
! Index/loop counters
      INTEGER :: IJK    ! fluid cell index
      INTEGER :: H      ! reaction loop counter
      INTEGER :: L, M   ! global phase index loop counters
      INTEGER :: NN      ! global species index
      INTEGER :: lN     ! reaction speices index
      INTEGER :: lM     ! reaction phase index
      INTEGER :: mXfr   ! global phase index for mass transfer

! This is a mask for the number of solids phases. The mask is used to
! avoid issues for DEM based reactions where MMAX and DES_MMAX may be
! in conflict.
      INTEGER :: lMMAX
      LOGICAL :: noSolids

! User-defined reaction rates returned from USR_RATES
      DOUBLE PRECISION RATES(NO_OF_RXNS)
! Reaction and species specific reaction rate.
      DOUBLE PRECISION lRate

! Single reaction formation/consumption rate of each gas species.
      DOUBLE PRECISION :: r_Rgp(DIMENSION_N_g)
      DOUBLE PRECISION :: r_ROXgc(DIMENSION_N_g)
! Single reaction formation/consumption rate of each solids species.
      DOUBLE PRECISION :: r_Rsp(DIMENSION_M, DIMENSION_N_s)
      DOUBLE PRECISION :: r_ROXsc(DIMENSION_M, DIMENSION_N_s)

! Rate of interphase enthalpy transfer
      DOUBLE PRECISION :: r_RxH(0:MMAX, 0:MMAX)
! Local heat of reactions
      DOUBLE PRECISION :: r_HoRg, r_HoRs(1:MMAX)
! Conversion factor for HoR (CGS to SI)
      DOUBLE PRECISION :: l_2SI

! Reaction limiters. If a species mass fraction is less than this
! value, then the reaction is suppressed.
      DOUBLE PRECISION :: speciesLimiter

! External functions:
!`````````````````````````````````````````````````````````````````````//
! Calculates enthalpy (cal/gram)
      DOUBLE PRECISION, EXTERNAL :: CALC_H

! Initialize global storage arrays to zero
!---------------------------------------------------------------------//
!     Gas                Solids
      R_GP    = ZERO;    R_SP    = ZERO
      ROX_GC  = ZERO;    ROX_SC  = ZERO
      SUM_R_G = ZERO;    SUM_R_S = ZERO
      HOR_G   = ZERO;    HOR_S   = ZERO

      R_PHASE = ZERO

! Set the species limiter:
      speciesLimiter = ZERO_X_gs

! Set the conversion factor for heat of reaction.
      l_2SI = merge(4.183925d3, ONE, UNITS(1:2) == 'SI')

! Set the number of TFM solids phases.
! the following 2 lines could be replaced by assing lmmax=smax
      noSolids = DISCRETE_ELEMENT .AND. (.NOT.DES_CONTINUUM_HYBRID)
      lMMAX = merge(0, MMAX, noSolids)

! Loop over each fluid cell.
      DO IJK = ijkstart3, ijkend3
      IF (FLUID_AT(IJK)) THEN

      RATES(:) = ZERO

! Calculate user defined reaction rates.
      CALL USR_RATES(IJK, RATES)

! Loop over reactions.
      RXN_LP: DO H = 1, NO_OF_RXNS

! Skip empty reactions
         IF(Reaction(H)%nSpecies == 0) CYCLE RXN_LP
         IF(COMPARE(RATES(H),ZERO)) CYCLE RXN_LP

! Initialize reaction loop arrays
         r_Rgp = ZERO; r_ROXgc = ZERO; r_HoRg = ZERO
         r_Rsp = ZERO; r_ROXsc = ZERO; r_HoRs = ZERO

         r_RxH = ZERO

! Calculate the rate of formation/consumption for each species.
!---------------------------------------------------------------------//
         DO lN = 1, Reaction(H)%nSpecies
! Global phase index.
            M = Reaction(H)%Species(lN)%pMap
! Global species index.
            NN = Reaction(H)%Species(lN)%sMap
! Index for interphase mass transfer. For a gas/solid reaction, the
! index is stored with the gas phase. For solid/solid mass transfer
! the index is stored with the source phase.
            mXfr = Reaction(H)%Species(lN)%mXfr
            lRate = RATES(H) * Reaction(H)%Species(lN)%MWxStoich
! Gas Phase:
            IF(M == 0) THEN
! Consumption of gas phase species.
               IF(lRate < ZERO) THEN
                  IF(X_g(IJK,NN) > speciesLimiter) THEN
                     r_ROXgc(NN) = r_ROXgc(NN) - lRate/X_g(IJK,NN)
! Enthalpy transfer associated with mass transfer. (gas/solid)
                     IF(M /= mXfr) r_RxH(M,mXfr) = r_RxH(M,mXfr) +     &
                        lRate * CALC_H(T_G(IJK),0,NN)
                  ELSE
! There is an insignificant amount of reactant. Skip this reaction.
                     RATES(H) = ZERO
                     CYCLE RXN_LP
                  ENDIF
               ELSE
! Formation of gas phase species.
                  r_Rgp(NN) = r_Rgp(NN) + lRate
! Enthalpy transfer associated with mass transfer. (gas/solid)
                  IF(M /= mXfr) r_RxH(M,mXfr) = r_RxH(M,mXfr) +        &
                     lRate * CALC_H(T_s(IJK,mXfr),0,NN)
               ENDIF


! Solids Phase M:
            ELSE
! Consumption of solids phase species.
               IF(lRate < ZERO) THEN
                  IF(X_s(IJK,M,NN) > speciesLimiter) THEN
                     r_ROXsc(M,NN) = r_ROXsc(M,NN) - lRate/X_s(IJK,M,NN)
! Enthalpy transfer associated with mass transfer. (solid/solid) This
! is only calculated from the source (reactant) material.
                     IF(M /= mXfr) THEN
                        IF(M < mXfr) THEN
                           r_RxH(M,mXfr) =  r_RxH(M,mXfr) + lRate *    &
                              Reaction(H)%Species(lN)%xXfr *           &
                              CALC_H(T_s(IJK,M),M,NN)
                        ELSE
                           r_RxH(mXfr,M) =  r_RxH(mXfr,M) - lRate *    &
                              Reaction(H)%Species(lN)%xXfr *           &
                              CALC_H(T_s(IJK,M),M,NN)
                        ENDIF
                     ENDIF
                  ELSE
! There is an insignificant amount of reactant. Skip this reaction.
                     RATES(H) = ZERO
                     CYCLE RXN_LP
                  ENDIF
               ELSE
! Formation of solids phase species.
                  r_Rsp(M,NN) = r_Rsp(M,NN) + lRate
               ENDIF
            ENDIF
         ENDDO ! Loop of species



! Map the local reaction production/consumption into global variables.
!---------------------------------------------------------------------//
! The outer reaction loop was cycled if there was insufficient amount
! of species available.

! Gas phase production
         R_gp(IJK,:NMAX(0)) = R_gp(IJK,:NMAX(0)) +                     &
            r_Rgp(:NMAX(0))
! Gas phase consumption divided by species mass fraction
         ROX_gc(IJK,:NMAX(0)) = ROX_gc(IJK,:NMAX(0)) +                 &
            r_ROXgc(:NMAX(0))

         DO M=1,lMMAX
! Solids phase species production
            R_sp(IJK,M,:NMAX(M)) = R_sp(IJK,M,:NMAX(M)) +              &
               r_Rsp(M,:NMAX(M))
! Solids phase species consumption divided by species mass fraction
            ROX_sc(IJK,M,:NMAX(M)) = ROX_sc(IJK,M,:NMAX(M)) +          &
               r_ROXsc(M,:NMAX(M))
         ENDDO


! Calculate and store the heat of reaction.
!---------------------------------------------------------------------//
         IF(ENERGY_EQ) THEN
! Automated heat of reaction calculations
            IF(Reaction(H)%Calc_DH) THEN
! Loop over reaction species.
               DO lN = 1, Reaction(H)%nSpecies
! Global phase index.
                  M = Reaction(H)%Species(lN)%pMap
! Global species index.
                  NN = Reaction(H)%Species(lN)%sMap
! Rate of formation/consumption for speices N
                  lRate = RATES(H) * Reaction(H)%Species(lN)%MWxStoich
! Gas phase enthalpy change from energy equation derivation.
                  IF(M == 0) THEN
                     r_HORg = r_HORg + CALC_H(T_g(IJK),0,NN) * lRate
! Solid phase enthalpy change from energy equation derivation.
                  ELSE
                     r_HORs(M) = r_HORs(M) + CALC_H(T_s(IJK,M),M,NN)*lRate
                  ENDIF
               ENDDO

! Complete the skew-symettric for enthalpy transfer with mass transfer
               DO M=1, lMMAX
                   DO L=0, M-1
                    r_RxH(M,L) = - r_RxH(L,M)
                  ENDDO
               ENDDO
! Apply enthalpy transfer associated with mass transfer to get the
! complete heat of reaction of heat phse for Reaction H.
               DO L=0, lMMAX
                  DO M = 0, lMMAX
                     IF(L == M) CYCLE
                     IF(L == 0) THEN
                        r_HORg = r_HORg - r_RxH(L,M)
                     ELSE
                        r_HORs(L) = r_HORs(L) - r_RxH(L,M)
                     ENDIF
                  ENDDO
               ENDDO
! Apply unit conversion and store heats of reaction in the global array.
               HOR_g(IJK) = HOR_g(IJK) + l_2SI*r_HORg
               DO M=1,lMMAX
                  HOR_s(IJK,M) = HOR_s(IJK,M) + l_2SI*r_HORs(M)
               ENDDO
            ELSE
! User-defined heat of reaction.
               HOR_g(IJK) = HOR_g(IJK) + Reaction(H)%HoR(0) * RATES(H)
               DO M=1, lMMAX
                  HOR_s(IJK,M) = HOR_s(IJK,M) +                        &
                     Reaction(H)%HoR(M) * RATES(H)
               ENDDO
            ENDIF
         ENDIF  ! ENERGY_EQ

! Update rate of interphase mass transfer.
!---------------------------------------------------------------------//
          DO LM=1, (DIMENSION_LM+DIMENSION_M-1)
             R_PHASE(IJK,LM) = R_PHASE(IJK,LM) +                       &
                RATES(H) * Reaction(H)%rPHASE(LM)
          ENDDO
      ENDDO RXN_LP ! Loop over reactions.


! Calculate the toal rate of formation and consumption for each species.
!---------------------------------------------------------------------//
      IF(SPECIES_EQ(0)) THEN
         SUM_R_G(IJK) = SUM( R_gp(IJK,:NMAX(0)) -                      &
            ROX_gc(IJK,:NMAX(0))*X_g(IJK,:NMAX(0)))
      ELSE
         DO H=1, NO_OF_RXNS
            IF(Reaction(H)%nPhases <= 0) CYCLE
            DO M=1, lMMAX
               LM = 1 + ((M-1)*M)/2
               SUM_R_G(IJK) = SUM_R_G(IJK) +                           &
                  RATES(H) * Reaction(H)%rPHASE(LM)
            ENDDO
         ENDDO
      ENDIF

      DO M=1, lMMAX
         IF(SPECIES_EQ(M)) THEN
            SUM_R_S(IJK,M) = SUM( R_sp(IJK,M,:NMAX(M)) -               &
               RoX_sc(IJK,M,:NMAX(M))*X_s(IJK,M,:NMAX(M)))
         ELSE
            DO H=1, NO_OF_RXNS
               IF(Reaction(H)%nPhases <= 0) CYCLE
               DO L=0, M-1
                  LM = 1 + L + ((M-1)*M)/2
                  SUM_R_S(IJK,M) = SUM_R_S(IJK,M) -                    &
                     RATES(H) * Reaction(H)%rPHASE(LM)
               ENDDO
               DO L=M+1, lMMAX
                  LM = 1 + M + ((L-1)*L)/2
                  SUM_R_S(IJK,M) = SUM_R_S(IJK,M) +                    &
                     RATES(H) * Reaction(H)%rPHASE(LM)
               ENDDO
            ENDDO
         ENDIF
      ENDDO


      ENDIF  ! Fluid_At(IJK)
      END DO ! IJK

      RETURN

      END SUBROUTINE RRATES0
