!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_DES_2FLUID                                         !
!                                                                      !
!  Purpose: This subroutine is only called from the CONTINUUM side. It !
!  is only called at the start of each time step for explicitly coupled!
!  cases. Otherwise, it called every iteration.                        !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_DES_2FLUID

      use discretelement, only: DES_CONTINUUM_COUPLED
      use discretelement, only: DES_CONTINUUM_HYBRID
      use discretelement, only: DES_EXPLICITLY_COUPLED

      use particle_filter, only: DES_INTERP_SCHEME_ENUM
      use particle_filter, only: DES_INTERP_GARG

      use discretelement, only: DES_CONTINUUM_COUPLED
      use particle_filter, only: DES_DIFFUSE_MEAN_FIELDS

! Contribution to gas momentum equation due to drag
      use discretelement, only: DRAG_BM
! Scalar cell center total drag force
      use discretelement, only: F_GDS
! Flag for 3D simulatoins.
      use geometry, only: DO_K
      use run, only: ANY_SPECIES_EQ

      use des_thermo, only: CALC_CONV_DES
      use rxns, only: RRATE

      IMPLICIT NONE

      IF(.NOT.DES_CONTINUUM_COUPLED) RETURN

      IF(DES_EXPLICITLY_COUPLED) THEN
! Bin particles to the fluid grid.
         CALL PARTICLES_IN_CELL
! Calculate interpolation weights
         CALL CALC_INTERP_WEIGHTS
! Calculate mean fields (EPg).
         CALL COMP_MEAN_FIELDS

! Calculate gas phase source terms: gas-solids heat transfer
         IF(CALC_CONV_DES) CALL CONV_GS_GAS1
! Calculate gas phase source terms: gas-solids mass transfer
         IF(RRATE) CALL RXNS_GS_GAS1
      ENDIF

! Calculate gas phase source terms: gas-solids drag force.
      SELECT CASE(DES_INTERP_SCHEME_ENUM)
      CASE(DES_INTERP_GARG) ; CALL DRAG_GS_GAS0
      CASE DEFAULT; CALL DRAG_GS_GAS1
      END SELECT

! Calculate solids phase source terms: solids-solids drag force.
      IF(DES_CONTINUUM_HYBRID) THEN
         SELECT CASE(DES_INTERP_SCHEME_ENUM)
         CASE DEFAULT; CALL DRAG_SS_TFM_NONINTERP
         END SELECT
      ENDIF

! Apply the diffusion filter.
      IF(DES_DIFFUSE_MEAN_FIELDS) THEN
         CALL DIFFUSE_MEAN_FIELD(F_GDS,'F_GDS')
         CALL DIFFUSE_MEAN_FIELD(DRAG_BM(:,1),'DRAG_BM(1)')
         CALL DIFFUSE_MEAN_FIELD(DRAG_BM(:,2),'DRAG_BM(2)')
         IF(DO_K) CALL DIFFUSE_MEAN_FIELD(DRAG_BM(:,3),'DRAG_BM(3)')
!         IF(ENERGY_EQ) THEN
!            CALL DIFFUSE_MEAN_FIELD(CONV_Sc,'CONV_Sc')
!            CALL DIFFUSE_MEAN_FIELD(CONV_Sp,'CONV_Sp')
!         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE CALC_DES_2FLUID


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: DES_2FLUID_CONV                                         !
!  Author: J.Musser                                   Date: 15-Jan-11  !
!                                                                      !
!  Purpose: This routine is called from the continuum phase and        !
!  calculates the source term from the particles to the fluid.         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE  DES_2FLUID_CONV(S_P, S_C)

      Use discretelement, only: DES_EXPLICITLY_COUPLED
      Use des_thermo, only: CONV_Sp, CONV_Sc
      USE geometry, only: FLAG
      Use param, only: DIMENSION_3

      use run, only: ODT
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED

      IMPLICIT NONE

! Passed Variables
!---------------------------------------------------------------------//
! Source term on LHS
      DOUBLE PRECISION, INTENT(INOUT) :: S_P(DIMENSION_3)
! Source term on RHS
      DOUBLE PRECISION, INTENT(INOUT) :: S_C(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
      IF(DES_ONEWAY_COUPLED) RETURN

      IF(DES_EXPLICITLY_COUPLED) THEN
         WHERE(FLAG==1)
            S_P = S_P + CONV_Sp ! GAMMA
            S_C = S_C + CONV_Sc ! GAMMA*Tp
         END WHERE

! Redistribute the energy over the fluid time step. Note that by the
! time this routine is called, S_C and S_P have already been multiplied
! by the fluid cell volume. Thus, the mapping should result in units
! of energy per time.
      ELSE
         WHERE(FLAG==1) &
            S_C = S_C + oDT*CONV_Sc ! GAMMA*(Tg-Ts)
      ENDIF

      RETURN
      END SUBROUTINE  DES_2FLUID_CONV


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: DES_2FLUID_RXNS                                         !
!  Author: J.Musser                                   Date: 15-Jan-11  !
!                                                                      !
!  Purpose: This routine is called from the continuum phase and        !
!  calculates the source term from the particles to the fluid.         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE  DES_2FLUID_RXNS

      Use discretelement, only: DES_EXPLICITLY_COUPLED
      Use des_thermo, only: CONV_Sp, CONV_Sc
      USE geometry, only: FLAG
      Use param, only: DIMENSION_3
      USE rxns, only : RRATE
      USE rxns
      USE des_rxns
      USE energy, only: HOR_g
      use fldvar, only: X_g
      use run, only: DT
      use param1, only: ZERO, SMALL_NUMBER, ONE
      use stiff_chem, only: STIFF_CHEMISTRY
      use compar, only: IJKSTART3, IJKEND3
      use geometry, only: VOL
      use functions, only: FLUID_AT
      use toleranc, only: ZERO_X_gs
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Flag to use stiff chemistry solver
      use stiff_chem, only: stiff_chemistry

      IMPLICIT NONE

! Passed Variables
!---------------------------------------------------------------------//
! Source term on LHS
!      DOUBLE PRECISION, INTENT(INOUT) :: S_P(DIMENSION_3)
! Source term on RHS
!      DOUBLE PRECISION, INTENT(INOUT) :: S_C(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: IJK
      DOUBLE PRECISION :: toTFM, lDT

      IF(DES_ONEWAY_COUPLED) RETURN
      IF(STIFF_CHEMISTRY) RETURN

! For DEM simulations that do not have a homogeneous gas phase reaction,
! the gas phase arrays need to be cleared.
      IF(.NOT.RRATE .OR. STIFF_CHEMISTRY) THEN
         SUM_R_G = ZERO
         HOR_G = ZERO
         R_GP = ZERO
         ROX_GC = ZERO
         R_PHASE = ZERO
      ENDIF

! Redistribute the energy over the fluid time step. Note that by the
! time this routine is called, S_C and S_P have already been multiplied
! by the fluid cell volume. Thus, the mapping should result in units
! of energy per time.
      lDT = merge(ONE, DT, DES_EXPLICITLY_COUPLED)

      DO IJK=IJKSTART3,IJKEND3
      IF(.NOT.FLUID_AT(IJK)) CYCLE
         toTFM = ONE/(lDT * VOL(IJK))
         R_gp(IJK,:) = R_gp(IJK,:) + DES_R_gp(IJK,:)*toTFM
         R_PHASE(IJK,:) = R_PHASE(IJK,:) + DES_R_PHASE(IJK,:)*toTFM
         SUM_R_g(IJK) = SUM_R_g(IJK) + DES_SUM_R_g(IJK)*toTFM
         HOR_g(IJK) = HOR_g(IJK) + DES_HOR_g(IJK)*toTFM
         WHERE(X_g(IJK,:) > ZERO_X_gs) RoX_gc(IJK,:) = &
            RoX_gc(IJK,:)+DES_R_gc(IJK,:)*toTFM/X_g(IJK,:)
      ENDDO

      RETURN
      END SUBROUTINE  DES_2FLUID_RXNS
