!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CN_EXTRAPOL                                            C
!                                                                      C
!  Purpose: Performs the explicit extrapolation phase when             C
!           Crank-Nicholson based second order accurate time           C
!           integration is activated (CN_ON = TRUE)                    C
!                                                                      C
!           Make sure dt is set to dt/2 when CN_ON=T                   C
!                                                                      C
!  Author: Aeolus Research, Inc. (A.Gel)             Date: APR-30-01   C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: ROP_g, EP_g, ROP_s, IJKMAX2, MMAX, U_s, V_s,  C
!                        W_s                                           C
!                                                                      C
!  Variables modified: 19 field variables are time integrated
!                      EP_g , EP_go  : Void fraction at t & t-dt time  C
!                      P_G , P_go    : Gas pressure at t & t-dt time   C
!                      P_star,P_staro: Solids pressure that maintains  C
!                                      EP_g >= EP_star                 C
!                      RO_g , RO_go  : Gas density                     C
!                      ROP_g , ROP_go: Macroscopic gas density         C
!                      U_g , U_go    : x-component of gas velocity     C
!                      V_g , V_go    : y-component of gas velocity     C
!                      W_g , W_go    : z-component of gas velocity     C
!                      T_g , T_go    : Gas phase temperature           C
!                      X_g , X_go    : Gas species mass fraction       C
!                      Scalar,Scalaro: User-defined Scalars            C
!                      ROP_s , ROP_so: Macroscopic density of solids   C
!                                      phases                          C
!                      T_s , T_so    : Solid phase temperature         C
!                      THETA_m,THETA_mo: Granular temperature of m^th  C
!                                        phase                         C
!                      trD_S_C, trD_S_Co: trace of D_s                 C
!                      U_s , U_so    : x-component of solid velocity   C
!                      V_s , V_so    : y-component of solid velocity   C
!                      W_s , W_so    : z-component of solid velocity   C
!                      X_s , X_so    : Solid species mass fraction     C
!                                                                      C
!                                                                      C
!  Local variables: NONE                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CN_EXTRAPOL
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE scalars
      USE trace
      USE run
      USE geometry
      USE indices
      USE compar
      USE mpi_utility
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                    Indices
      INTEGER :: M, IJK
!-----------------------------------------------
!
!!!$omp    parallel do private(IJK)
      DO ijk = IJKSTART3, IJKEND3
!          IF(.NOT.IS_ON_myPE_OWNS(I_OF(IJK), J_OF(IJK), K_OF(IJK))) CYCLE
          EP_G(ijk)   = 2.*EP_G(ijk) - EP_GO(ijk)
!AE TIME Pressure is not postprocessed/calculated in the current implementation
!          P_G(ijk)    = 2.*P_G(ijk)  - P_GO(ijk)
!          P_STAR(ijk) = 2.*P_STAR(ijk) - P_STARO(ijk)
          RO_G(ijk)   = 2.*RO_G(ijk) - RO_GO(ijk)
          ROP_G(ijk)  = 2.*ROP_G(ijk)- ROP_GO(ijk)
          U_G(ijk) = 2.*U_G(ijk) - U_GO(ijk)
          V_G(ijk) = 2.*V_G(ijk) - V_GO(ijk)
          W_G(ijk) = 2.*W_G(ijk) - W_GO(ijk)
          IF (ENERGY_EQ) T_G(ijk) = 2.*T_G(ijk) - T_GO(ijk)

!AE TIME Update P_g with state equation based on extrapolated ro_g, T_g
!          EOSG = UNSCALE(PG)*MW/(GAS_CONST*TG)
!          IF (.NOT.WALL_AT(IJK)) THEN
!             P_G(ijk) = ((RO_G(ijk)*(GAS_CONST*T_G(ijk))/MW_MIX_G(IJK))- P_ref)/P_scale
!             P_G(ijk) = ((RO_G(ijk)*(GAS_CONST*T_G(ijk))/MW_MIX_G(IJK)))
!         ENDIF

          IF (SPECIES_EQ(0)) THEN
            IF (NMAX(0) > 0) THEN
              X_G(ijk,:NMAX(0)) = 2.*X_G(ijk,:NMAX(0)) - X_GO(ijk,:NMAX(0))
            ENDIF
          ENDIF

          IF (NScalar > 0) THEN
            Scalar(ijk,:NScalar) = 2.*Scalar(ijk,:NScalar) - ScalarO(ijk,:NScalar)
          ENDIF

          DO M = 1, MMAX
            ROP_S(ijk,M) = 2.*ROP_S(ijk,M) - ROP_SO(ijk,M)
            IF (ENERGY_EQ) T_S(ijk,M) = 2.*T_S(ijk,M) - T_SO(ijk,M)
            IF (GRANULAR_ENERGY) THEN
              THETA_M(ijk,M) = 2.*THETA_M(ijk,M) - THETA_MO(ijk,M)
              TRD_S_C(ijk,M) = 2.*TRD_S_C(ijk,M) - TRD_S_CO(ijk,M)
            ENDIF
            U_S(ijk,M) = 2.*U_S(ijk,M) - U_SO(ijk,M)
            V_S(ijk,M) = 2.*V_S(ijk,M) - V_SO(ijk,M)
            W_S(ijk,M) = 2.*W_S(ijk,M) - W_SO(ijk,M)
            IF (SPECIES_EQ(M)) THEN
              IF (NMAX(M) > 0) THEN
                X_S(ijk,M,:NMAX(M)) = 2.*X_S(ijk,M,:NMAX(M)) - X_SO(ijk,M,:NMAX(M))
              ENDIF
            ENDIF
          END DO


      END DO


      RETURN
      END SUBROUTINE CN_EXTRAPOL

!// Comments on the modifications for DMP version implementation
!// 120 Replaced the index for initialization: (:IJKMAX2) to just (:)
