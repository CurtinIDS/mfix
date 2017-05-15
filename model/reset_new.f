!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: RESET_NEW                                              C
!  Purpose: Reset the new variables with the stored previous-time-step C
!           values of field variables.                                 C
!    *****Remember to modify update_old also
!                                                                      C
!  Author: M. Syamlal                                 Date: FEB-6-97   C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: ROP_g, EP_g, ROP_s, IJKMAX2, MMAX, U_s, V_s,  C
!                        W_s                                           C
!                                                                      C
!  Variables modified: ROP_go, ROP_so, IJK, M, U_so, V_so, W_so C
!                                                                      C
!  Local variables: NONE                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE RESET_NEW

!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------

      USE fldvar
      USE physprop, only: mmax, nmax
      USE run, only: energy_eq, call_dqmom, granular_energy, k_epsilon, species_eq
      USE scalars, only: nscalar
      USE trace, only: trd_s_c, trd_s_co

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
      INTEGER :: M
!
!                    error index
      INTEGER :: IER

!-----------------------------------------------

      EP_G(:) = EP_GO(:)
      P_G(:) = P_GO(:)
      P_STAR(:) = P_STARO(:)
      RO_G(:) = RO_GO(:)
      ROP_G(:) = ROP_GO(:)
      U_G(:) = U_GO(:)
      V_G(:) = V_GO(:)
      W_G(:) = W_GO(:)
      IF (ENERGY_EQ) T_G(:) = T_GO(:)
      IF (SPECIES_EQ(0)) THEN
        IF (NMAX(0) > 0) THEN
          X_G(:,:NMAX(0)) = X_GO(:,:NMAX(0))
        ENDIF
      ENDIF

      IF (NScalar > 0) THEN
        Scalar(:,:NScalar) = ScalarO(:,:NScalar)
      ENDIF

      IF (K_Epsilon) THEN
        K_Turb_G(:) = K_Turb_GO(:)
        E_Turb_G(:) = E_Turb_GO(:)
      ENDIF

      DO M = 1, MMAX
        ROP_S(:,M) = ROP_SO(:,M)
! add by rong
        If (Call_DQMOM) D_P(:,M)=D_Po(:,M)
!       If (NScalar>0) ome(:,M)=ome_o(:,M)
! add by rong
        IF (ENERGY_EQ) T_S(:,M) = T_SO(:,M)
        IF (GRANULAR_ENERGY) THEN
          THETA_M(:,M) = THETA_MO(:,M)
          TRD_S_C(:,M) = TRD_S_CO(:,M)
        ENDIF
        U_S(:,M) = U_SO(:,M)
        V_S(:,M) = V_SO(:,M)
        W_S(:,M) = W_SO(:,M)
        IF (SPECIES_EQ(M)) THEN
          IF (NMAX(M) > 0) THEN
            X_S(:,M,:NMAX(M)) = X_SO(:,M,:NMAX(M))
          ENDIF

          RO_S(:,M) = RO_SO(:,M)
        ENDIF
      END DO

!     Recalculate all coefficients
      CALL CALC_COEFF_ALL (0, IER)

      RETURN
      END SUBROUTINE RESET_NEW

!// Comments on the modifications for DMP version implementation
!// 120 Replaced the index for initialization: (:IJKMAX2) to just (:)
