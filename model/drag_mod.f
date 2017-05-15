! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: drag                                                   C
!  Purpose: Common block containing drag arrays                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-92  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE drag

! Gas-solids drag
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_gs

! Solids-solids drag
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  F_ss

! Off diagonal friction coefficient in HYS drag relation
  DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  beta_ij

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function(s): C_DsxRe                                                C
!  Purpose:                                                            C
!     Calculate single sphere drag correlation multiplied by           C
!     the Reynolds number or                                           C
!     Calculate the single sphere drag correlation                     C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC

! Dalla Valle (1948)
!----------------------------------------------------------------->>>
  DOUBLE PRECISION FUNCTION C_DSXRE_DV(RE)
    USE param
    USE param1
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number

    C_DSXRE_DV = (0.63D0*SQRT(RE) + 4.8D0)**2
    RETURN
  END FUNCTION C_DSXRE_DV

! Schiller and Naumann (1933)
!----------------------------------------------------------------->>>
  DOUBLE PRECISION FUNCTION C_DS_SN(RE)
    USE param
    USE param1
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: RE ! Reynolds number

    C_DS_SN = 24.D0*(1.D0 + 0.15D0*RE**0.687D0)/(RE+SMALL_NUMBER)
    RETURN
  END FUNCTION C_DS_SN
!-----------------------------------------------------------------<<<

END MODULE drag
