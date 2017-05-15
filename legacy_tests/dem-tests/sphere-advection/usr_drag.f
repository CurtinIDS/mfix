!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_DRAG                                               !
!                                                                      !
!  Purpose: Provide a hook for user defined drag law implementation.   !
!                                                                      !
!  This routine is called from inside fluid (TFM) and particle (DES)   !
!  loops. The fluid cell index (IJK) and phase (TFM) or particle index !
!  (DES) is passed.                                                    !
!                                                                      !
!  ***************************   WARNING   **************************  !
!  *----------------------------------------------------------------*  !
!  * The dummy arguments changed in the 2015-1 MFIX Release.        *  !
!  *                                                                *  !
!  *   1) Phase index (M) is now particle index (NP) for DES. This  *  !
!  *      is reflected in the name change M --> M_NP.               *  !
!  *                                                                *  !
!  *   2) The fluid velocity was added as a dummy argument. This    *  !
!  *      provides access to the interpolated gas velocity for      *  !
!  *      coupled DES simulations.                                  *  !
!  *                                                                *  !
!  ******************************************************************  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_USR(IJK, M_NP, lDgA, EPg, Mug, ROg, VREL, DPM, &
         ROs, lUg, lVg, lWg)

! Fluid cell I, J, K, IJK containing particle and phase
      use discretelement, only: DES_VEL_NEW

      use error_manager

      IMPLICIT NONE

! Index of fluid cell:
      INTEGER, INTENT(IN) :: IJK
! TFM SOLIDS --> Index of phase (M)
! DES SOLIDS --> Index of particle (NP); M = PIJK(NP,5)
      INTEGER, INTENT(IN) :: M_NP

! drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDgA
! gas volume fraction
      DOUBLE PRECISION, INTENT(IN) :: EPg
! gas laminar viscosity
      DOUBLE PRECISION, INTENT(IN) :: Mug
! gas density
      DOUBLE PRECISION, INTENT(IN) :: ROg
! Magnitude of gas-solids relative velocity
      DOUBLE PRECISION, INTENT(IN) :: VREL
! particle diameter of solids phase M or
! average particle diameter if PCF
      DOUBLE PRECISION, INTENT(IN) :: DPM
! particle density of solids phase M
      DOUBLE PRECISION, INTENT(IN) :: ROs
! fluid velocity components:
! o TFM: Averaged from faces to cell center
! o DES: Interpolated to the particle's position
      DOUBLE PRECISION, INTENT(IN) :: lUg, lVg, lWg


! This case doesn't need a drag force.
      lDgA = 0.0d0

! Set the particle velocity to the gas velocity
      DES_VEL_NEW(M_NP,1) = lUg
      DES_VEL_NEW(M_NP,2) = lVg
      DES_VEL_NEW(M_NP,3) = lWg


      RETURN
      END SUBROUTINE DRAG_USR
