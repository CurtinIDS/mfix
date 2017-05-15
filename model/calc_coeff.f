!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      !
!  Subroutine: CALC_COEFF_ALL                                          !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, exchange rates, and reaction rates.  !
!                                                                      !
!  Author: M. Syamlal                                 Date: 25-AUG-05  !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced:                                               !
!  Variables modified:                                                 !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_COEFF_ALL(FLAG, IER)

! Global variables:
!-----------------------------------------------------------------------
! Double precision: 1.0d0
      use param1, only: ONE
! Under relaxation factor for gas-solids drag coefficient
      use ur_facs, only: UR_F_gs
! Under relaxation factor solid conductivity coefficient for IA theory
      use ur_facs, only: UR_kth_sml
! Flag for DES coupled simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
! Flag for explcit coupling between the fluid and particles.
      use discretelement, only: DES_EXPLICITLY_COUPLED

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! FLAG = 0, overwrite the coeff arrays, (e.g. start of a time step)
! FLAG = 1, do not overwrite
      INTEGER, intent(in) :: FLAG
! Error index
      INTEGER, intent(inout) :: IER

! Local variables
!-----------------------------------------------
! Under relaxation factor for gas-solids drag coefficient
      DOUBLE PRECISION :: loc_UR_F_gs ! Local copy
! Under relaxation factor solid conductivity coefficient for IA theory
      DOUBLE PRECISION :: loc_UR_kth_sml ! Local copy

!-----------------------------------------------------------------------

! 1) Backup user-defined coefficient relaxation factors.
! 2) Set user-defined coefficient relaxation factors to 1.
! Note that 'FLAG' is hard coded to 0 in time march and reset_new.
      IF(FLAG == 0) THEN
        loc_UR_F_gs = UR_F_gs;          UR_F_gs = ONE
        loc_UR_Kth_sml = UR_Kth_sml;    UR_Kth_sml = ONE
      ENDIF

! Calculate all physical properties, transport properties, and exchange
! rates.
      CALL CALC_COEFF(IER, 2)

! Calculate reaction rates and interphase mass transfer.
      CALL CALC_RRATE(IER)

! Restore all coefficient underrelaxation factors to original values.
      IF(FLAG == 0) THEN
        UR_F_gs = loc_UR_F_gs
        UR_Kth_sml = loc_UR_Kth_sml
      ENDIF

! DES interaction for explictly coupled simulations
      IF(DES_CONTINUUM_COUPLED .AND. DES_EXPLICITLY_COUPLED) &
         CALL CALC_DES_2FLUID

      RETURN
      END SUBROUTINE CALC_COEFF_ALL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_COEFF                                              !
!  Purpose: This routine directs the calculation of all physical and   !
!           transport properties, and exchange rates.                  !
!                                                                      !
!  Author: M. Syamlal                                 Date: 25-AUG-05  !
!  Reviewer:                                          Date:            !
!                                                                      !
!                                                                      !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!  Variables referenced:                                               !
!  Variables modified:                                                 !
!  Local variables:                                                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_COEFF(IER, pLevel)

! Flag for DES coupled simulation
      use discretelement, only: DES_CONTINUUM_COUPLED
! Flag for explcit coupling between the fluid and particles.
      use discretelement, only: DES_EXPLICITLY_COUPLED

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! Error index
      INTEGER, intent(inout) :: IER
! Level to calculate physical properties.
! 0) Only density
! 1) Everything but density
! 2) All physical properties
      INTEGER, intent(in) :: pLevel
!-----------------------------------------------------------------------

! Calculate physical properties: (density, specific heat, diameter)
      CALL PHYSICAL_PROP(IER, pLevel)

! Calculate transport properties: (conductivity, diffusivity, ect)
      CALL TRANSPORT_PROP(IER)

! Calculate interphase coeffs: (momentum and energy)
      CALL EXCHANGE(IER)

! Calculate DES coupled quantities.
      IF(DES_CONTINUUM_COUPLED .AND. .NOT.DES_EXPLICITLY_COUPLED) &
         CALL CALC_DES_2FLUID

      RETURN
      END SUBROUTINE CALC_COEFF

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_RRATE                                              C
!  Purpose: if rrate then calculate reaction rates and interphase      C
!           mass transfer. if present, calculate discrete reactions    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_RRATE(IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement, only : DISCRETE_ELEMENT
      USE rxns,only : RRATE, USE_RRATES
      use run, only: ANY_SPECIES_EQ

      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Error index
      INTEGER, INTENT(INOUT) :: IER

!-----------------------------------------------

! Calculate reaction rates and interphase mass transfer
      IF(RRATE) THEN
! Legacy hook: Calculate reactions from rrates.f.
         IF(USE_RRATES) THEN
            CALL RRATES (IER)
            IF(IER .EQ. 1) THEN
               CALL INIT_ERR_MSG('CALC_RRATE')
               WRITE(ERR_MSG, 1000)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ELSE
            CALL RRATES0 (IER)
         ENDIF

! DES Chemical reactions
         IF(ANY_SPECIES_EQ .AND. DISCRETE_ELEMENT) &
            CALL DES_2FLUID_RXNS
      ENDIF

      RETURN

 1000 FORMAT('Species balance equations are being solved; but chemical',/, &
         ' reactions are not specified in mfix.dat or in rrates.f.',/,  &
         ' Copy the file mfix/model/rrates.f into the run directory ',/, &
         ' and remove the initial section that returns IER=1.')

      END SUBROUTINE CALC_RRATE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_TRD_AND_TAU                                        !
!  Purpose: Calculate various terms in the gas and solids phase        !
!  stress tensor as indicated below                                    !
!                                                                      !
!  Author: M. Syamlal                                 Date: 25-AUG-05  !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_TRD_AND_TAU()

      USE run, only: jackson
! Stress tensor trace.
      USE visc_g, only : TRD_g
      USE visc_s, only : TRD_S
! Stress tensor cross terms.
      USE tau_g, only : TAU_U_G, TAU_V_G, TAU_W_G
      USE tau_g, only : cTAU_U_G, cTAU_V_G, cTAU_W_G
      USE tau_s, only : TAU_U_S, TAU_V_S, TAU_W_S
! Runtime flag for DEM model.
      USE discretelement, only: DISCRETE_ELEMENT
! Runtime flag for TFM-DEM hybrid model.
      USE discretelement, only: DES_CONTINUUM_HYBRID

      USE param1, only: zero
      implicit none

!-----------------------------------------------------------------------

! Calculate the trace of the stress tensor (gas phase; m=0)
      CALL CALC_TRD_G (TRD_G)

! Calculate the cross terms of the stress tensor (gas phase; m=0)
      CALL CALC_TAU_U_G (TAU_U_G, CTAU_U_G)
      CALL CALC_TAU_V_G (TAU_V_G, CTAU_V_G)
      CALL CALC_TAU_W_G (TAU_W_G, CTAU_W_G)

      IF (.NOT. JACKSON) THEN
         CTAU_U_G = ZERO
         CTAU_V_G = ZERO
         CTAU_W_G = ZERO
      ENDIF

! Bypass the following calculations if there are no TFM solids.
      IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN
! Calculate the cross terms of the stress tensor (solids phases; m>0)
         CALL CALC_TRD_S (TRD_S)
! Calculate the trace of the stress tensor (solids phases; m>0)
         CALL CALC_TAU_U_S (TAU_U_S)
         CALL CALC_TAU_V_S (TAU_V_S)
         CALL CALC_TAU_W_S (TAU_W_S)
      ENDIF

      RETURN
      END SUBROUTINE CALC_TRD_AND_TAU
