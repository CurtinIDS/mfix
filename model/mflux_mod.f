!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: mflux                                                  C
!  Purpose: Module for mass fluxes and densities at faces              C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE mflux

! x-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gE
! y-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gN
! z-component of gas mass flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_gT

! y-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sN
! x-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sE
! z-component of solids mass flux
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Flux_sT

! Added mass Flux Components to be used for scalar eq.
! Note: added mass apply only to one solids phase (M=M_AM)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Flux_gSE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Flux_sSE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Flux_gSN
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Flux_sSN
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Flux_gST
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Flux_sST

! macroscopic gas density at east face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gE
! macroscopic gas density at north face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gN
! macroscopic gas density at top face
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  ROP_gT
! macroscopic solids density at north face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sN
! macroscopic solids density at east face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sE
! macroscopic solids density at top face
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  ROP_sT


! for GHD Theory
! x-component of solids total number density flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_nE
! y-component of solids total number density flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_nN
! z-component of solids total number density flux
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  Flux_nT
! end GHD Theory modification


      END MODULE mflux
