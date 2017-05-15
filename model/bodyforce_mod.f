!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module:  BODYFORCE                                                  C
!  Purpose: Include file for all body force statement functions        C
!                                                                      C
!  Author: M. Syamlal                                 Date:  6-MAR-92  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE bodyforce

      USE constant, only: gravity_x, gravity_y, gravity_z
      IMPLICIT NONE

      CONTAINS

! Body force on gas at i+1/2, j, k
      DOUBLE PRECISION FUNCTION BFX_g(IJK)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ijk
      BFX_g = GRAVITY_X
      END FUNCTION BFX_g

! Body force on gas at i, j+1/2, k
      DOUBLE PRECISION FUNCTION BFY_g(IJK)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ijk
      BFY_g = GRAVITY_Y
      END FUNCTION BFY_g

! Body force on gas at i, j, k+1/2
      DOUBLE PRECISION FUNCTION BFZ_g(IJK)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ijk
      BFZ_g = GRAVITY_Z
      END FUNCTION BFZ_g

! Body force on solids m at i+1/2, j, k
      DOUBLE PRECISION FUNCTION BFX_s(IJK,M)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ijk,m
      BFX_s = GRAVITY_X
      END FUNCTION BFX_s

! Body force on solids m at i, j+1/2, k
      DOUBLE PRECISION FUNCTION BFY_s(IJK,M)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ijk,m
      BFY_s = GRAVITY_Y
      END FUNCTION BFY_s

! Body force on solids m at i, j, k+1/2
      DOUBLE PRECISION FUNCTION BFZ_s(IJK,M)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ijk,m
      BFZ_s = GRAVITY_Z
      END FUNCTION BFZ_s

      END MODULE bodyforce
