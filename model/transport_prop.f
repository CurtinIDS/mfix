!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: TRANSPORT_PROP                                          C
!  Purpose: Calculate the indicated transport properties that vary     C
!           with time if directed to do so by the corresponding flag   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Mods for MFIX 2.0 (old name CALC_PHYSPROP)                 C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE TRANSPORT_PROP()

! Global Variables:
!----------------------------------------------------------------------
! Number of solids phases.
      use physprop, only: MMAX
! Flags for calculating viscosity.
      use coeff, only: VISC
! Flags for calculating conductivity.
      use coeff, only: COND
! Flags for calculating diffusivity.
      use coeff, only: DIFF
! Flags for calculating particle-particle energy dissipation.
      use coeff, only: GRAN_DISS
! Kinetic theory model.
      use run, only: KT_TYPE_enum
      use run, only: ia_2005, gd_1999, gtsh_2012
      use kintheory, only: CALC_IA_ENERGY_DISSIPATION_SS
      use kintheory, only: CALC_GD_99_ENERGY_DISSIPATION_SS
      use kintheory, only: CALC_GTSH_ENERGY_DISSIPATION_SS

      implicit none

! Local variables
!-----------------------------------------------------------------------
! Loop counter
      INTEGER :: M ! Solids phase


      IF (VISC(0)) CALL CALC_MU_G()    ! Fluid viscosity
      IF (COND(0)) CALL CALC_K_G()     ! Fluid conductivity
      IF (DIFF(0)) CALL CALC_DIF_G()   ! Fluid diffusivity

      DO M = 1, MMAX
! Particle-Particle Energy Dissipation
! for gtsh theory this call needs to be done before calc_mu_s so that
! the cooling rate is available for mu_s
         IF (GRAN_DISS(M)) THEN
         SELECT CASE (KT_TYPE_ENUM)
            CASE (IA_2005)
               CALL CALC_IA_ENERGY_DISSIPATION_SS(M)
            CASE(GD_1999)
               CALL CALC_GD_99_ENERGY_DISSIPATION_SS(M)
            CASE(GTSH_2012)
               CALL CALC_GTSH_ENERGY_DISSIPATION_SS(M)
            END SELECT
         ENDIF
! these were moved after gran_diss since some quantities above are
! needed in the subsequent gtsh calculations performed via calc_mu_s
         IF (COND(M)) CALL CALC_K_S (M)   ! Solids conductivity
         IF (VISC(M)) CALL CALC_MU_S (M)  ! Solids viscosity
         IF (DIFF(M)) CALL CALC_DIF_S (M) ! Solids diffusivity
      ENDDO

      RETURN
      END SUBROUTINE TRANSPORT_PROP
