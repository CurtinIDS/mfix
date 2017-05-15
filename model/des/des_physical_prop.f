!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine name: PHYSICAL_PROP                                      !
!                                                                      !
!  Purpose: Calculate physical properties that vary with time.         !
!                                                                      !
!  Author: J.Musser                                   Date: 09-May-11  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DES_PHYSICAL_PROP

! Modules
!---------------------------------------------------------------------//
      Use des_rxns
      Use des_thermo
      Use discretelement
      Use funits
      Use param
      Use param1
      Use physprop
      Use run

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Logical indicating to write additional data about the particle for
! debugging purposes.

! Local Variables
!---------------------------------------------------------------------//
! Index of particle
      INTEGER :: NP
!......................................................................!


! Specific heat
!-----------------------------------------------------------------------
! This only needs calculated when solving the energy equations.
      IF(ENERGY_EQ) THEN
      ENDIF

      RETURN
      END SUBROUTINE DES_PHYSICAL_PROP
