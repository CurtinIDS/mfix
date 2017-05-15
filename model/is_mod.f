!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: IS                                                          C
!  Purpose: Internal surface specifications                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 10-OCT-92  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE is

! Modules
!---------------------------------------------------------------------//
      use param, only: DIM_M, DIMENSION_IS
!---------------------------------------------------------------------//


! x coordinate of the west face of a region where internal surfaces
! are specified
      DOUBLE PRECISION :: IS_X_w (DIMENSION_IS)

! x coordinate of the east face of a region where internal surfaces
! are specified
      DOUBLE PRECISION :: IS_X_e (DIMENSION_IS)

! y coordinate of the south face of a region where internal surfaces
! are specified
      DOUBLE PRECISION :: IS_Y_s (DIMENSION_IS)

! y coordinate of the north face of a region where internal surfaces
! are specified
      DOUBLE PRECISION :: IS_Y_n (DIMENSION_IS)

! z coordinate of the bottom face of a region where internal surfaces
! are specified
      DOUBLE PRECISION :: IS_Z_b (DIMENSION_IS)

! z coordinate of the top face of a region where internal surfaces
! are specified
      DOUBLE PRECISION :: IS_Z_t (DIMENSION_IS)

! i index of the west face of a region where internal surfaces
! are specified
      INTEGER :: IS_I_w (DIMENSION_IS)

! i index of the east face of a region where internal surfaces
! are specified
      INTEGER :: IS_I_e (DIMENSION_IS)

! j index of the south face of a region where internal surfaces
! are specified
      INTEGER :: IS_J_s (DIMENSION_IS)

! j index of the north face of a region where internal surfaces
! are specified
      INTEGER :: IS_J_n (DIMENSION_IS)

! k index of the bottom face of a region where internal surfaces
! are specified
      INTEGER :: IS_K_b (DIMENSION_IS)

! k index of the top face of a region where internal surfaces
! are specified
      INTEGER :: IS_K_t (DIMENSION_IS)

! Type of internal surface:
! IMPERMEABLE - no gas or solids flow through the surface
! SEMIPERMEABLE - only gas flows through the surface
      CHARACTER(LEN=16) :: IS_TYPE (DIMENSION_IS)

! Logical variable to determine whether an IS is defined
       LOGICAL :: IS_DEFINED (DIMENSION_IS)

! Are there any IS defined?
       LOGICAL :: ANY_IS_DEFINED

! Character variable with values E, N, and T to determine the
! IS plane of a flow cell
       CHARACTER :: IS_PLANE (DIMENSION_IS)

! Permeability coefficients for semipermeable internal surface:
! 1- Darcy coefficient
! 2- Inertial resistance factor
      DOUBLE PRECISION :: IS_PC (DIMENSION_IS, 2)

! Solids velocity at the semipermeable surface
      DOUBLE PRECISION :: IS_VEL_s (DIMENSION_IS, DIM_M)

      END MODULE is
