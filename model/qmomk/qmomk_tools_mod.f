!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_TOOLS                                            C
!  Purpose: Helper functions for QMOMK implementation                  C
!                                                                      C
!  Author: Alberto Passalacqua                        Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

#include "version.inc"

MODULE qmomk_tools

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: det3
  PUBLIC :: inv3
  PUBLIC :: diag3
  PUBLIC :: cholesky3
  PUBLIC :: transpose3
  PUBLIC :: multiplyMatrix3

CONTAINS

  SUBROUTINE DET3 (A, det)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,3) :: A
    DOUBLE PRECISION, INTENT(OUT) :: det

    det = a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) +  a(1,3)*a(2,1)*a(2,3) &
         - a(3,1)*a(2,2)*a(1,3) - a(3,2)*a(2,3)*a(1,1) - a(3,3)*a(2,1)*a(1,2)

  END SUBROUTINE DET3

  SUBROUTINE INV3 (A, B)

    IMPLICIT none

    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,3) :: A
    DOUBLE PRECISION, iNTENT(OUT), DIMENSION(3,3) :: B

    DOUBLE PRECISION :: detA;

    CALL det3(A, detA)

    IF (detA == 0.D0) THEN
       PRINT *,'QMOMK: Null determinant in matrix'
       ERROR_STOP
    END IF

    B(1,1) = (A(3,3)*A(2,2) - A(3,2)*A(2,3))/detA
    B(1,2) = -(A(3,3)*A(1,2) - A(3,2)*A(1,3))/detA
    B(1,3) = (A(2,3)*A(1,2) - A(2,2)*A(1,3))/detA

    B(2,1) = -(A(3,3)*A(2,1) - A(3,1)*A(2,3))/detA
    B(2,2) = (A(3,3)*A(1,1) - A(3,1)*A(1,3))/detA
    B(2,3) = -(A(2,3)*A(1,1) - A(2,1)*A(1,3))/detA

    B(3,1) = (A(3,2)*A(2,1) - A(3,1)*A(2,2))/detA
    B(3,2) = -(A(3,2)*A(1,1) - A(3,1)*A(1,2))/detA
    B(3,3) = (A(2,2)*A(1,1) - A(2,1)*A(1,2))/detA
  END SUBROUTINE INV3

  SUBROUTINE DIAG3 (a, b, c, mat)

    IMPLICIT NONE

    INTEGER :: i, j
    DOUBLE PRECISION, INTENT(IN) :: a, b, c
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,3) :: mat

    DO i = 1, 3
       DO j = 1, 3
          mat(i,j) = 0.D0
       END DO
    END DO

    mat(1,1) = a
    mat(2,2) = b
    mat(3,3) = c
  END SUBROUTINE DIAG3

  SUBROUTINE CHOLESKY3 (A, L)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,3) :: A
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,3) :: L

    INTEGER :: i, j

    DO i = 1, 3
       DO j = 1, 3
          L(i,j) = 0.D0
       END DO
    END DO

    L(1,1) = SQRT(A(1,1))
    IF (L(1,1) == 0.D0) THEN
       PRINT *,'Impossible to find Cholesky decomposition.'
       ERROR_STOP
    ELSE
       L(2,1) = A(2,1)/L(1,1)
       L(2,2) = SQRT(A(2,2) - (L(2,1))**2)
       L(3,1) = A(3,1)/L(1,1)

       IF (L(2,2) ==  0.D0) THEN
          PRINT *,'Impossible to find Cholesky decomposition.'
          ERROR_STOP
       ELSE
          L(3,2) = (A(3,2) - L(3,1)*L(2,1))/L(2,2)
          L(3,3) = SQRT(A(3,3) - (L(3,1))**2 - (L(3,2))**2)
       END IF
    END IF
  END SUBROUTINE CHOLESKY3

  SUBROUTINE TRANSPOSE3 (A, T)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,3) :: A
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,3) :: T

    INTEGER :: i, j

    DO i = 1, 3
       DO j = 1, 3
          T(i,j) = A(j,i)
       END DO
    END DO

  END SUBROUTINE TRANSPOSE3

  SUBROUTINE MULTIPLYMATRIX3 (A, B, P)

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,3) :: A
    DOUBLE PRECISION, INTENT(IN), DIMENSION(3,3) :: B
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,3) :: P

    INTEGER :: i, j, k
    DOUBLE PRECISION ::sum

    sum = 0.D0

    DO i = 1, 3
       DO j = 1, 3
          DO k = 1, 3
             sum = sum + A(i,k)*B(k,j)
          END DO
          P(i,j) = sum
          sum = 0.D0
       END DO
    END DO
  END SUBROUTINE MULTIPLYMATRIX3

END MODULE qmomk_tools
