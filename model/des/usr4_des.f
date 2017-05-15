!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS4_DES                                               !
!                                                                      !
!  Purpose: This routine is called before thermodynamic calculations   !
!  and is user-definable. The user may insert code in this routine or  !
!  call appropriate user defined subroutines.                          !
!                                                                      !
!  This routien is called from an IJK loop, hence the IJK dummy        !
!  argument is the only defined index. This routine is usefule for     !
!  calculations based on fluid grid properties needed for reaction     !
!  rate calculations.                                                  !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR4_DES(IJK)

      Use des_rxns
      Use des_thermo
      Use discretelement
      Use fldvar
      USE geometry
      USE indices
      Use run
      Use usr

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJK ! Fluid cell index


      RETURN
      END SUBROUTINE USR4_DES
