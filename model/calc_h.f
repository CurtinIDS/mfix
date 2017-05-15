!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function name: CALC_H                                               C
!  Purpose: Calculate specific enthalpy of species N in phase M        C
!                                                                      C
!  Author: M. Syamlal                                Date: 27-DEC-2007 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION CALC_H(refT, M, NN)

! Modules
!---------------------------------------------------------------------//
      USE physprop, only: mw_g, mw_s, HfrefoR
      USE constant, only: RGAS => GAS_CONST_cal
      USE read_thermochemical, only: calc_ICpoR
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! cell, phase and species indices
      DOUBLE PRECISION, INTENT(IN) :: refT   ! Temperature
      INTEGER, INTENT(IN) :: M ! Phase index
      INTEGER, INTENT(IN) :: NN ! Species index

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: ICpoR
      DOUBLE PRECISION :: lMW
      INTEGER :: IER
!---------------------------------------------------------------------//

      IER = 0

      if(M == 0)then
         lMW = MW_g(NN)
      else
         lMW = MW_s(M,NN)
      endif

! Integrate the specific heat from zero to refT
      ICpoR = calc_ICpoR(refT, M, NN, IER)

! Evaluate the enthalpy of species N at refT
      CALC_H = (HfrefoR(M,NN)  + ICpoR) * (RGAS / lMW)

      RETURN
      END FUNCTION CALC_H
