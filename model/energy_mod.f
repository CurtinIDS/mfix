! -*- f90 -*-
MODULE energy


!   Gas-phase heat of reaction
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  HOR_g

!   Solids-phase heat of reaction
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  HOR_s

!   Gas-solids heat transfer coefficient
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  GAMA_gs

!   Gas-phase radiation coefficient
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  GAMA_Rg

!   Solids-phase radiation coefficient
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  GAMA_Rs

!   Gas-phase radiation temperature
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  T_Rg

!   Solids-phase radiation temperature
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  T_Rs

CONTAINS

  !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
  !  Source terms for the energy equations.  The source term is linearized as
  !  S = S_c - S_p * T,  where S_c and S_p must be positive.
  !
  !  By default the source terms have been coded for radiation sources.

  !
  !  S_p for gas phase at i, j, k
  DOUBLE PRECISION FUNCTION S_Rpg(IJK)
    USE fldvar, only: t_g
    IMPLICIT NONE
    INTEGER IJK
    S_Rpg = 4.d0 * GAMA_Rg(IJK) *  T_g(IJK)**3
  END FUNCTION S_Rpg

  !  S_c for gas phase at i, j, k
  DOUBLE PRECISION FUNCTION S_Rcg(IJK)
    USE fldvar, only: t_g
    IMPLICIT NONE
    INTEGER IJK
    S_Rcg = GAMA_Rg(IJK) * ( T_Rg(IJK)**4 + 3.d0 * T_g(IJK)**4 )
  END FUNCTION S_Rcg

  !  S_p for solids phase at i, j, k
  DOUBLE PRECISION FUNCTION S_Rps(IJK, M)
    USE fldvar, only: t_s
    IMPLICIT NONE
    INTEGER IJK, M
    S_Rps = 4.d0 * GAMA_Rs(IJK, M) *  T_s(IJK, M)**3
  END FUNCTION S_Rps

  !  S_c for solids phase at i, j, k
  DOUBLE PRECISION FUNCTION S_Rcs(IJK, M)
    USE fldvar, only: t_s
    IMPLICIT NONE
    INTEGER IJK, M
    S_Rcs = GAMA_Rs(IJK, M) * ( T_Rs(IJK, M)**4 &
         + 3.d0 * T_s(IJK, M)**4 )
  END FUNCTION S_Rcs

END MODULE energy
