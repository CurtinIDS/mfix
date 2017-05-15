#include "version.inc"
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!  Subroutine: DES_CALC_GAMMA                                          !
!                                                                      !
!  Purpose: Calculate the heat transfer coefficient (GAMMAxSA) for     !
!  particle-fluid heat transfer.                                       !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!  REF: Ranz, W.E. and Marshall, W.R., "Friction and transfer          !
!       coefficients for single particles and packed beds," Chemical   !
!       Engineering Science, Vol. 48, No. 5, pp 247-253, 1925.         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_GAMMA_DES(NP, pGAMMA)

      USE compar
      USE constant
      USE des_thermo
      USE discretelement
      USE exit, only: mfix_exit
      USE fldvar
      USE geometry
      USE indices
      USE param1
      USE physprop
      USE fun_avg
      USE functions

      IMPLICIT NONE

! Passed variables
!---------------------------------------------------------------------//
! Index value of particle
      INTEGER, INTENT(IN) :: NP
! Convective heat transfer coefficient
      DOUBLE PRECISION, INTENT(OUT) :: pGAMMA

! Local variables
!---------------------------------------------------------------------//
! Fluid cell indices
      INTEGER IMJK, IJMK, IJKM
! Double precision value for 1/3
      DOUBLE PRECISION, PARAMETER  :: THIRD = (1.0d0/3.0d0)

      DOUBLE PRECISION N_Pr  ! Prandtl Number
      DOUBLE PRECISION N_Re  ! Reynolds Number
      DOUBLE PRECISION N_Nu  ! Nusselt Number

! Magnitude of slip velocity
      DOUBLE PRECISION SLIP
! Fluid velocity
      DOUBLE PRECISION cUg, cVg, cWg
      DOUBLE PRECISION Us, Vs, Ws
! Index value of fluid cell
      INTEGER :: IJK
!---------------------------------------------------------------------//

      IJK = PIJK(NP,4)
! Initialization
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

      SELECT CASE(DES_CONV_CORR_ENUM)

         CASE (RANZ_1952) ! (Ranz and Mrshall, 1952)
! Initialize variables
            SLIP = ZERO
            N_Re = ZERO
            N_Nu = ZERO
! Gas velocity in fluid cell IJK
            cUg = AVG_X_E(U_g(IMJK), U_g(IJK), 1)
            cVg = AVG_Y_N(V_g(IJMK), V_g(IJK))
! Particle Velocity
            Us = DES_VEL_NEW(NP,1)
            Vs = DES_VEL_NEW(NP,2)

! Calculate the magnitude of the slip velocity
            IF(NO_K) THEN
               SLIP = SQRT((cUg-Us)**2 + (cVg-Vs)**2)
            ELSE
               cWg = AVG_Z_T(W_g(IJKM), W_g(IJK))
               Ws = DES_VEL_NEW(NP,3)
               SLIP = SQRT((cUg-Us)**2 + (cVg-Vs)**2 + (cWg-Ws)**2)
            ENDIF

! Calculate the Prandtl Number
            IF(K_G(IJK) > ZERO) THEN
               N_Pr = (C_PG(IJK)*MU_G(IJK))/K_G(IJK)
            ELSE
               N_Pr = LARGE_NUMBER
            ENDIF

! Calculate the particle Reynolds Number
            IF(MU_G(IJK) > ZERO) THEN
               N_Re = (2.0d0*DES_RADIUS(NP)*SLIP*RO_g(IJK)) / MU_g(IJK)
            ELSE
               N_Re = LARGE_NUMBER
            ENDIF

! Calculate the Nusselt Number
            N_Nu = 2.0d0 + 0.6d0 *((N_Re)**HALF * (N_Pr)**THIRD)

! Calculate the convective heat transfer coefficient
            pGAMMA = (N_Nu * K_G(IJK))/(2.0d0 * DES_RADIUS(NP))

         CASE DEFAULT
            WRITE(*,*)'INVALID DES CONVECTION MODEL'
            ERROR_STOP 'INVALID DES CONVECTION MODEL'
            CALL MFIX_EXIT(myPE)
      END SELECT

      RETURN
      END SUBROUTINE CALC_GAMMA_DES
