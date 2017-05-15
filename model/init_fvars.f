!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: INIT_FVARS                                             !
!  Purpose: Initialize all field variables.                            !
!                                                                      !
!  Author: M. Syamlal                                 Date: 23-JAN-94  !
!  Reviewer: J.Musser                                 Date:  8-Oct-13  !
!                                                                      !
!  Literature/Document References:                                     !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE INIT_FVARS

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase volume farction
      USE fldvar, only: EP_G
! Pressures
      USE fldvar, only: P_G    ! Gas
      USE fldvar, only: P_s    ! Solids
      USE fldvar, only: P_STAR ! Solids at EP_star
! Densities
      USE fldvar, only: RO_G  ! Gas
      USE fldvar, only: RO_s  ! Solids
! Bulk densities: RO*EP
      USE fldvar, only: ROP_G ! Gas
      USE fldvar, only: ROP_s ! Solids
! Gas velocity components:
      USE fldvar, only: U_G  ! x-axis
      USE fldvar, only: V_G  ! y-axis
      USE fldvar, only: W_G  ! z-axis
! Solids velocity components:
      USE fldvar, only: U_S  ! x-axis
      USE fldvar, only: V_S  ! y-axis
      USE fldvar, only: W_S  ! z-axis
! Temperature
      USE fldvar, only: T_G  ! Gas
      USE fldvar, only: T_S  ! Solids
! Species mass fractions
      USE fldvar, only: X_G  ! Gas
      USE fldvar, only: X_S  ! Solids
! Solids particle diameter
      USE fldvar, only: D_P
! Granular Energy
      USE fldvar, only: THETA_M
! Granular energy conductivity
      USE physprop, only: KTH_S
! Gas turbulence
      USE fldvar, only: K_Turb_G
      USE fldvar, only: E_Turb_G
! User-scalar equation field
      USE fldvar, only: Scalar
! Reaction Rates
      USE rxns, only: ReactionRates

! Global Parameters:
!---------------------------------------------------------------------//
      USE param1, only: UNDEFINED
      USE param1, only: ZERO


      IMPLICIT NONE

! Passed Variables:
!---------------------------------------------------------------------//
! NONE

! Local Variables:
!---------------------------------------------------------------------//
! NONE

      IF(allocated(EP_G)) EP_G = UNDEFINED

      IF(allocated(P_G)) P_G  = UNDEFINED
      IF(allocated(P_S)) P_s = UNDEFINED
      IF(allocated(P_STAR)) P_STAR = UNDEFINED

      IF(allocated(RO_G)) RO_G = UNDEFINED
      IF(allocated(RO_S)) RO_s = UNDEFINED
      IF(allocated(ROP_G)) ROP_G = UNDEFINED
      IF(allocated(ROP_S)) ROP_s = UNDEFINED

      IF(allocated(U_G)) U_G = UNDEFINED
      IF(allocated(V_G)) V_G = UNDEFINED
      IF(allocated(W_G)) W_G = UNDEFINED

      IF(allocated(U_S)) U_S = UNDEFINED
      IF(allocated(V_S)) V_S = UNDEFINED
      IF(allocated(W_S)) W_S = UNDEFINED

      IF(allocated(T_G)) T_G = ZERO
      IF(allocated(T_S)) T_S = ZERO

      IF(allocated(X_G)) X_G = ZERO
      IF(allocated(X_S)) X_S = ZERO

      IF(allocated(D_P)) D_P = ZERO
      IF(allocated(THETA_M)) THETA_M = ZERO
      IF(allocated(KTH_S)) KTH_S = UNDEFINED

      IF(allocated(K_Turb_G)) K_Turb_G = ZERO
      IF(allocated(E_Turb_G)) E_Turb_G = ZERO

      IF(allocated(Scalar)) Scalar = ZERO
      IF(allocated(ReactionRates)) ReactionRates = ZERO

      RETURN
      END SUBROUTINE INIT_FVARS
