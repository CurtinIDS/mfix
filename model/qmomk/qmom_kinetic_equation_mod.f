!                                                                       C
!   Module name: QMOMKIN                                                C
!   Purpose: QMOMK mod file                                           C
!                                                                       C
!                                                                       C
!   Author: Alberto Passalacqua                        Date:            C
!   Reviewer:                                          Date:            C
!   Comments:                                                           C
!                                                                       C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE qmom_kinetic_equation

  USE qmomk_parameters

  !     QMOMK variables

  !     Old values of quadrature weights (NN, Nx, Ny, Nz, phase)
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_N0
  !     Old values of moments
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_M0
  !     Old values of abscissas
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_U0
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_V0
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_W0
  !     Current values of weights
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_N1
  !     Current values of moments
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_M1
  !     Current values of abscissas
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_U1
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_V1
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_W1

  !    Mean velocities
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: QMOMK_U_S
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: QMOMK_V_S
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: QMOMK_W_S
  !    Drag term
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_TAU_DRAG
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: QMOMK_F_GS

  !    Collision time
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: QMOMK_COLLISION_TIME

  !     QMOMB logicals

  LOGICAL QMOMK
  LOGICAL QMOMK_COUPLED
  LOGICAL PRINT_QMOMK_DATA


  CHARACTER(64) QMOMK_TYPE
  !     Strings

  !     Collision operator
  !     Accepted values: BGK, Boltzmann
  CHARACTER(64) QMOMK_COLLISIONS

  !     Wall BC type
  !     Accepted values: SPECULAR_REFLECTIVE
  CHARACTER(64) QMOMK_WALL_BC_TYPE
  !     Integers

  INTEGER QMOMK_COLLISIONS_ORDER

  !     QMOMB time step
  DOUBLE PRECISION :: QMOMK_DT

  !     QMOMB CFL number
  DOUBLE PRECISION :: QMOMK_CFL

END MODULE qmom_kinetic_equation
