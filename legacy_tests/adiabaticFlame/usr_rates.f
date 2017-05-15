!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: USR_RATES                                              !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 10-Oct-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR_RATES(IJK, RATES)

      USE param
      USE param1
      USE parallel
      USE fldvar
      USE rxns
      USE energy
      USE geometry
      USE run
      USE indices
      USE physprop
      USE constant
      USE funits
      USE compar
      USE sendrecv

      USE toleranc
      USE usr
      USE fun_avg
      USE functions

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: IJK

      DOUBLE PRECISION, DIMENSION(NO_OF_RXNS), INTENT(OUT) :: RATES

!-----------------------------------------------
      INCLUDE 'species.inc'
      INCLUDE 'usrnlst.inc'

! Reaction specific variables:
!`````````````````````````````````````````````````````````````````````//
      DOUBLE PRECISION c_O2    ! Oxygen concentration mol/cm^3
      DOUBLE PRECISION c_CH4   ! Methane concentration mol/cm^3


! CH4_Comb:  CH4 + 2O2 --> CO2 + 2H2O  (mol/cm^3.s)
!---------------------------------------------------------------------//
! Note: The CH4 combustion rate is artificial and used for the
! adiabatic flame test case.


      c_O2  = (RO_g(IJK)*X_g(IJK,O2)/MW_g(O2))
      c_CH4 = (RO_g(IJK)*X_g(IJK,CH4)/MW_g(CH4))

      RATES(CH4_Comb) = C(1) * EP_g(IJK) * c_O2 * c_CH4


      RETURN

      END SUBROUTINE USR_RATES
