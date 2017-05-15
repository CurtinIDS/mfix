!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_SPECIES                                            !
!                                                                      !
!  Purpose: Common elements for MFIX-DEM species transfer.             !
!  condition.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 16-Jun-10  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE DES_RXNS

      ! USE param
      USE rxn_com, only: reaction_block

! Data Storage:
!---------------------------------------------------------------------//
! discrete solids species mass fractions (PARTICLES, 0:MAX_DES_NMAX))
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_X_s

! Rate of production/consumption of solids species
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: DES_R_s

! Numerical integration:
!---------------------------------------------------------------------//
! Previous time step's rate of change. Used for Adams-Bashforth
! time integration scheme.
! 1) particle mass
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dMdt_OLD
! 2) particle species mass
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dXdt_OLD

! Interphase transfer variables.
!---------------------------------------------------------------------//
! Amount produced of gas phase species
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_R_gp
! Amount consumed of gas phase species
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_R_gc
! Net production (+) or consumption (-)  of gas phase species
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_SUM_R_g
! Amount of interphase mass transfer.
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_R_PHASE
! Amount of gas phase enthalpy change due to phase change/chemical reaction
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_HOR_g

! Phase change/reactive chemistry:
!---------------------------------------------------------------------//
! Actual number of Reactions
      INTEGER NO_OF_DES_RXNS

! Array linking all of the reaction data.
      TYPE(REACTION_BLOCK), DIMENSION(:), TARGET, ALLOCATABLE :: DES_Reaction


      END MODULE DES_RXNS
