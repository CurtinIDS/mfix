!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RRATES(IER)                                            C
!  Purpose: Calculate reaction rates for various reactions in cell ijk C
!                                                                      C
!  Author:                                            Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!  O. Levenspiel, "Chemical Reaction Engineering," Wiley Eastern,
!       New Delhi (1974)
!                                                                      C
!                                                                      C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
!
      SUBROUTINE RRATES(IER)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar        !//d
      USE constant
      USE energy
      USE exit, ONLY: mfix_exit
      USE fldvar
      USE functions
      USE funits
      USE geometry
      USE indices
      USE machine
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run, only: energy_eq, units, species_eq
      USE rxns
      USE sendrecv      !// 400
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Local phase and species indices
      INTEGER          L, LM, M, N

!                      cell index
      INTEGER          IJK

      DOUBLE PRECISION R_tmp(0:MMAX, 0:MMAX), RxH_xfr(0:MMAX, 0:MMAX)
      DOUBLE PRECISION X_tmp(0:MMAX, 0:MMAX, Dimension_n_all)
      DOUBLE PRECISION RXNA, Trxn
      DOUBLE PRECISION, EXTERNAL ::calc_h
!
      DOUBLE PRECISION :: refT_M
      DOUBLE PRECISION :: refT_L

!-----------------------------------------------
      IER =0
      R_tmp = UNDEFINED
!
!  ---  Remember to include all the local variables here for parallel
!  ---- processing
!$omp  parallel do private(ijk, R_tmp, L, LM, M, N)

      DO IJK = IJKSTART3, IJKEND3

         IF (FLUID_AT(IJK)) THEN
!
!
!  User input is required in sections 1 through 4.
!
!1111111111111111111111111111111111111111111111111111111111111111111111111111111
!
! 1. Write the rates of various reactions:
!    Write the reaction rates for each of the reactions as RXNxF and RXNxB (both
!    quantities >= 0), where x identifies the reaction, F stands for forward
!    rate, and B stands for the backward rate.  The rates can be in
!    g-mole/(cm^3.s) or g/(cm^3.s).  For the sake of clarity, give the reaction
!    scheme and the units in a comment statement above the rate expression.
!    The volume (cm^3) is that of the computational cell.  Therefore, for
!    example, the rate term of a gas phase reaction will have a multiplicative
!    factor of epsilon. Note that X_g and X_s are mass fractions
!
!              No   1   2    3
!     GAS Species   A,  R, Inert

!  a)   A -> 3R
!       Levenspiel (p. 111, Example 4)
!            converting the rate expression given by Levenspiel for use here
!           -rA = 1E-2 * C_A^0.5 g-mol/L.s    (C_A is in mol/L)
!               = 1E-3 * 1E-2 * C_A^0.5 g-mol/cm^3.s
!               = EP_g * 1E-3 * 1E-2 * C_A^0.5 g-mol/cm^3.s (in reactor volume)
!               = EP_g * 1E-3 * 1E-2 * C_A^0.5 g-mol/cm^3.s (in reactor volume)
!               = EP_g * 1E-3 * 1E-2 * (1.E3*RO_g*X_A/M_A)^0.5 g-mol/cm^3.s (in reactor volume)

        Rxna = EP_g(IJK) * C(1) * (RO_g(IJK)*X_g(IJK,1)/Mw_g(1))**C(2)

!
!
!2222222222222222222222222222222222222222222222222222222222222222222222222222222
!
! 2. Write the formation and consumption rates of various species:
!    Obtain the rates of formation and consumption of various species
!    in g/(cm^3.s) from the rate expressions RXNxF and RXNxB obtained in the
!    previous section.  Pay attention to the units of RXNxF and RXNxB.
!    the formation rates for gas species n are added to get R_gp (IJK, n).
!    All the consumption rates are added and then divided by X_g(IJK, n) to
!    get RoX_gc(IJK, n).  If X_g(IJK, n) is zero and species n is likely
!    to be consumed in a reaction then it is recommended that RoX_gc (IJK, n)
!    be initialized to the derivative of the consumption rate with respect to
!    X_g at X_g=0.
!    If the slope is not known analytically a small value such as 1.0e-9 may
!    instead be used.  A similar procedure is used for all the species in the
!    solids phases also.

!  (1) A
      R_gp(IJK, 1) = ZERO
      IF(X_g(IJK, 1) .GT. ZERO) THEN
        RoX_gc(IJK, 1) = RXNA  * MW_g(1) / X_g(IJK, 1)
      ELSE
        RoX_gc(IJK, 1) = 1.0e-9
      ENDIF

!  (2) R
      R_gp(IJK, 2) = 3. * RXNA * MW_g(2)
      RoX_gc(IJK, 2) = ZERO

!  (3) Inert
      R_gp(IJK, 3) = ZERO
      RoX_gc(IJK, 3) = ZERO

!
!  GAS SPECIES
!
!
!  SOLIDS SPECIES
!
!
!3333333333333333333333333333333333333333333333333333333333333333333333333333333
!
! 3.  Determine the g/(cm^3.s) transferred from one phase to the other.
!          R_tmp(To phase #, From phase #)
!     e.g. R_tmp(0,1) -  mass generation of gas phase from solids-1,
!          R_tmp(0,2) -  mass generation of gas phase from solids-2,
!          R_tmp(1,0) -  mass generation of solid-1 from gas = -R_tmp(0,1)
!          R_tmp(1,2) -  mass generation of solid-1 from solids-2.
!     Note, for example, that if gas is generated from solids-1 then
!     R_tmp(0,1) > 0.
!     The R-phase matrix is skew-symmetric and diagonal elements are not needed.
!     Only one of the two skew-symmetric elements -- e.g., R_tmp(0,1) or
!     R_tmp(1,0) -- needs to be specified.
!
!     X_tmp(M,L, N) is the mass fraction of species N in the material transferred
!     between phase-M and phase-L. If the destination phase is M, then N is the
!     index of the species in phase-M, otherwise N is the index of the species
!     in phase-L; e.g. (1) In the reaction C+1/2O2 --> CO,
!     the destination phase is gas phase. Then N will be equal to the index of CO
!     in gas phase. (2) If H2O is trasferred between liquid and gas phases either
!     evaporation or condensation, then the index must change depending upon the
!     direction of mass transfer. For condensation, N is the species index
!     of H2O in the liquid phase and for evaporation it is the species index
!     of H2O in the gas phase. Also Sum_over_N (X_tmp(M,L, N)) should be equal to 1.
!
      if(MMAX > 0) then
        R_tmp(0,1) =  ZERO
	X_tmp(0, 1, :) = ZERO
	if(R_tmp(0,1) > 0) then
	  X_tmp(0, 1, 1) = 1.0
	else
	  X_tmp(0, 1, 1) = 1.0
	endif
      endif
!
!4444444444444444444444444444444444444444444444444444444444444444444444444444444
!
! 4.  Determine the heat of reactions in cal/(cm^3.s) at the
!     temperature T_g or T_s1.  Note that for exothermic reactions
!     HOR_g (or HOR_s) will be negative. The assignment of heat of reaction
!     is user defined as it depends upon the microphysics near the interface,
!     which is averaged out in the multiphase flow equations.  For example,
!     heat of Reaction for the C + O2 reaction is split into parts;
!     CO formation is assigned to the solid phase and CO2 formation from CO to
!     the gas phase.
!     *** This section is no longer needed as the heats of reactions are
!         calculated below.  If you need to override the automatic calculation,
!         comment out the calculations below.
!
!==============================================================================
!
!     No user input is required below this line
!-----------------------------------------------------------------------------
!   Determine g/(cm^3.s) of mass generation for each of the phases by adding
!   the reaction rates of all the individual species.

            SUM_R_G(IJK) = ZERO
            IF (SPECIES_EQ(0)) THEN
               IF (NMAX(0) > 0) THEN
                  SUM_R_G(IJK) = SUM_R_G(IJK) + SUM(R_GP(IJK,:NMAX(0))-ROX_GC(&
                     IJK,:NMAX(0))*X_G(IJK,:NMAX(0)))
               ENDIF
	    ELSE
	      DO M = 1, MMAX
	        IF(R_tmp(0,M) .NE. UNDEFINED)THEN
		  SUM_R_G(IJK) = SUM_R_G(IJK) + R_tmp(0,M)
		ELSEIF(R_tmp(M,0) .NE. UNDEFINED)THEN
		  SUM_R_G(IJK) = SUM_R_G(IJK) - R_tmp(M,0)
		ENDIF
	      ENDDO
            ENDIF
!
            DO M = 1, MMAX
               SUM_R_S(IJK,M) = ZERO
               IF (SPECIES_EQ(M)) THEN
                  IF (NMAX(M) > 0) THEN
                     SUM_R_S(IJK,M) = SUM_R_S(IJK,M) + SUM(R_SP(IJK,M,:NMAX(M))&
                        -ROX_SC(IJK,M,:NMAX(M))*X_S(IJK,M,:NMAX(M)))
                  ENDIF
	       ELSE
 	         DO L = 0, MMAX
	           IF(R_tmp(M,L) .NE. UNDEFINED)THEN
		     SUM_R_s(IJK,M) = SUM_R_s(IJK,M) + R_tmp(M,L)
		   ELSEIF(R_tmp(L,M) .NE. UNDEFINED)THEN
		     SUM_R_s(IJK,M) = SUM_R_s(IJK,M) - R_tmp(L,M)
		   ENDIF
	         ENDDO
               ENDIF
            END DO

!
!           Calculate the enthalpy of transferred material
!

            IF(ENERGY_EQ) THEN ! calculate heat of reactions only if energy eq. are solved
            DO M = 0, MMAX-1
               if(m==0) then
                  refT_M = T_G(IJK)
               else
                  refT_M = T_S(IJK,M)
               endif
            DO L = M+1, MMAX
               refT_L = T_S(IJK,L)


	        RxH_xfr(M, L) = zero
	        IF(R_tmp(M,L) .NE. UNDEFINED)THEN
		   IF(R_tmp(M,L) > ZERO) then ! phase-M is generated from phase-L
                     DO N = 1, NMAX(M)
                       RxH_xfr(M, L) =  RxH_xfr(M, L) + R_tmp(M,L) * X_tmp(M,L, N) * &
			                              CALC_H(refT_M, M, N)

                     END DO
		   else    ! phase-L is generated from phase-M
                     DO N = 1, NMAX(L)
                       RxH_xfr(M, L) =  RxH_xfr(M, L) + R_tmp(M,L) * X_tmp(M,L, N) * &
			                              CALC_H(refT_L, L, N)
                     END DO
		   endif
		ELSEIF(R_tmp(L,M) .NE. UNDEFINED)THEN
		   IF(R_tmp(L,M)> ZERO) then ! phase-L is generated from phase-M
                     DO N = 1, NMAX(L)
                       RxH_xfr(M, L) =  RxH_xfr(M, L) - R_tmp(L,M) * X_tmp(L,M, N) * &
			                              CALC_H(refT_L, L, N)
                     END DO
		   else ! phase-M is generated from phase-L
                     DO N = 1, NMAX(M)
                       RxH_xfr(M, L) =  RxH_xfr(M, L) - R_tmp(L,M) * X_tmp(L,M, N) * &
			                              CALC_H(refT_M, M, N)
                     END DO
		   endif
	        ENDIF
	      ENDDO
            END DO

            DO M = 1, MMAX
 	      DO L = 0, M-1
	        RxH_xfr(M, L) = -RxH_xfr(L, M)
	      ENDDO
            END DO

!
!           Calculate heats of reactions
!
            HOR_G(IJK) = zero
            DO N = 1, NMAX(0)
              HOR_G(IJK) = HOR_G(IJK) + &
	         (R_gp(IJK, N) - RoX_gc(IJK, N) * X_g(IJK, N)) * CALC_H(T_G(IJK), 0, N)
            END DO
            DO L = 1, MMAX
	      HOR_G(IJK) = HOR_G(IJK) - RxH_xfr(0, L)
	    ENDDO
            IF (UNITS == 'SI') HOR_G(IJK) = 4183.925D0*HOR_G(IJK)    !in J/kg K

            DO M = 1, MMAX
              HOR_s(IJK, M) = zero
              DO N = 1, NMAX(M)
                HOR_s(IJK, M) = HOR_s(IJK, M) + &
		  (R_sp(IJK, M, N) - RoX_sc(IJK, M, N) * X_s(IJK, M, N)) * &
              CALC_H(T_S(IJK,M), M, N)
              END DO
              DO L = 0, MMAX
	        if(M .NE. L) HOR_s(IJK, M) = HOR_s(IJK, M) - RxH_xfr(M, L)
	      ENDDO
              IF (UNITS == 'SI') HOR_s(IJK, M) = 4183.925D0*HOR_s(IJK, M)    !in J/kg K
            END DO
            ENDIF ! for energy_eq
!
!
!     Store R_tmp values in an array.  Only store the upper triangle without
!     the diagonal of R_tmp array.
!
            DO L = 0, MMAX
               DO M = L + 1, MMAX
                  LM = L + 1 + (M - 1)*M/2
                  IF (R_TMP(L,M) /= UNDEFINED) THEN
                     R_PHASE(IJK,LM) = R_TMP(L,M)
                  ELSE IF (R_TMP(M,L) /= UNDEFINED) THEN
                     R_PHASE(IJK,LM) = -R_TMP(M,L)
                  ELSE
                     CALL START_LOG
                     IF(.not.DMP_LOG)call open_pe_log(ier)
                     WRITE (UNIT_LOG, 1000) L, M
                     CALL END_LOG
                     call mfix_exit(myPE)
                  ENDIF
               END DO
            END DO

         ENDIF
      END DO

 1000 FORMAT(/1X,70('*')//' From: RRATES',/&
         ' Message: Mass transfer between phases ',I2,' and ',I2,&
         ' (R_tmp) not specified',/1X,70('*')/)
      RETURN
      END SUBROUTINE RRATES
