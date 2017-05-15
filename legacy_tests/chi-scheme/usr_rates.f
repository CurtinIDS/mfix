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

! O2_to_O3:  3.0*O2 --> 2.0*O3         (mol/cm^3.s)
!---------------------------------------------------------------------//
! Note: The rate is artificial and used for the chi-scheme test case.

      c_O2 = (RO_g(IJK)*X_g(IJK,O2)/MW_g(O2))

      RATES(O2_to_O3) = EP_g(IJK) * c(1) * c_O2 / 3.0d0


      RETURN

      END SUBROUTINE USR_RATES
