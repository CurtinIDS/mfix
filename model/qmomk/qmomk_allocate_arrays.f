!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_ALLOCATE_ARRAYS                                  C
!  Purpose: DES - Allocate QMOMK arrays                                C
!                                                                      C
!                                                                      C
!  Author: Alberto Passalacqua                        Date:            C
!  Reviewer:                                          Date:            C
!  Comments:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

SUBROUTINE qmomk_allocate_arrays

  USE param
  USE param1
  USE indices
  USE geometry
  USE physprop
  USE qmom_kinetic_equation

  IMPLICIT NONE

  !     Allocation of QMOMK arrays
  !     See qmomb_mod.f for a description of the variables

  Allocate( QMOMK_N0(QMOMK_NN, DIMENSION_3, MMAX) )
  Allocate( QMOMK_N1(QMOMK_NN, DIMENSION_3, MMAX) )

  Allocate( QMOMK_M0(QMOMK_NMOM, DIMENSION_3, MMAX) )
  Allocate( QMOMK_M1(QMOMK_NMOM, DIMENSION_3, MMAX) )

  Allocate( QMOMK_U0(QMOMK_NN, DIMENSION_3, MMAX) )
  Allocate( QMOMK_U1(QMOMK_NN, DIMENSION_3, MMAX) )

  Allocate( QMOMK_V0(QMOMK_NN, DIMENSION_3, MMAX) )
  Allocate( QMOMK_V1(QMOMK_NN, DIMENSION_3, MMAX) )

  Allocate( QMOMK_W0(QMOMK_NN, DIMENSION_3, MMAX) )
  Allocate( QMOMK_W1(QMOMK_NN, DIMENSION_3, MMAX) )

  Allocate( QMOMK_U_S(DIMENSION_3, MMAX) )
  Allocate( QMOMK_V_S(DIMENSION_3, MMAX) )
  Allocate( QMOMK_W_S(DIMENSION_3, MMAX) )

  Allocate( QMOMK_F_GS(QMOMK_NN, DIMENSION_3, MMAX) )
  Allocate( QMOMK_TAU_DRAG(QMOMK_NN, DIMENSION_3, MMAX) )

  Allocate( QMOMK_COLLISION_TIME(DIMENSION_3, MMAX) )

  RETURN
END SUBROUTINE qmomk_allocate_arrays
