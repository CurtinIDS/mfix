!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_DISTANCE                                           C
!  Purpose: Calculate the scalar and vector distance for the passed    C
!           direction.                                                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 07-JUL-92  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DISTANCE (DIST_MIN, D, DIM_MAX2, DIST_SC, DIST_VEC)

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! starting point of axis
      DOUBLE PRECISION :: DIST_MIN
! cell size across domain axis (dx,dy,dz)
      DOUBLE PRECISION :: D(*)
! maximum number of cells along axis (imax2,jmax2,kmax2)
      INTEGER :: DIM_MAX2
! distance to scalar cell center
      REAL :: DIST_SC(*)
! distance to 'vector' cell face
      REAL :: DIST_VEC(*)

! Local variables
!---------------------------------------------------------------------//
      INTEGER :: I
!---------------------------------------------------------------------//

! the array D is allocated 0:dimension, but the values are stored
! 1 to dim_max2. so the index must be shifted here to account
! for this offset in where the values are stored
      DIST_SC(1)  = DIST_MIN - 0.5 * D(2)
      DIST_VEC(1) = DIST_MIN

      DO I = 2,DIM_MAX2
         DIST_SC(I)  = DIST_VEC(I-1) + 0.5 * D(I+1)
         DIST_VEC(I) = DIST_VEC(I-1) + D(I+1)
      ENDDO

      RETURN
      END
