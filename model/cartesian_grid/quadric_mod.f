      MODULE quadric

!     Maximum of the number of quadrics that can be read
      INTEGER, PARAMETER          :: DIM_QUADRIC = 500
!     Nnumber of quadrics
      INTEGER                     :: N_QUADRIC
!     Current Quadric
      INTEGER :: QUADRIC_ID
!     form of quadric : 'normal' or one of the pre-defined quadrics
      CHARACTER (LEN=12), DIMENSION(DIM_QUADRIC) :: quadric_form
!     Scale factor for quadrics
      DOUBLE PRECISION :: quadric_scale
!     Characteristic values of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: Lambda_x,Lambda_y,Lambda_z
!     d - coefficient of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: dquadric
!     Translation components of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: t_x,t_y,t_z
!     Rotation angles (Deg) of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: theta_x,theta_y,theta_z
!     Radius for either Spere or Cylinder (pre-defined quadrics)
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: Radius
!     Radii for Torus
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: Torus_R1, Torus_R2
!     Radii for U-coil
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: UCOIL_R1, UCOIL_R2
      !     Y-location of bends for U-coil
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: UCOIL_Y1, UCOIL_Y2
!     Radii for Bend
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: BEND_R1, BEND_R2
!     Angles for Bend
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: BEND_THETA1, BEND_THETA2
!     Y-locations of cylinder-cone-cylinder
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: C2C_Y1, C2C_Y2
!     Radii of cylinder-cone-cylinder
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: C2C_R1, C2C_R2
!     Half-angle for cone (pre-defined quadrics)
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: Half_angle
!     Reactor 1 lower, upper cylinder radii
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: REACTOR1_R1,REACTOR1_R2
!     Reactor 1 lower, upper conical transition between cylinders
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: REACTOR1_Y1,REACTOR1_Y2
!     Reactor 1 lower, upper rounding locations
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: REACTOR1_YR1,REACTOR1_YR2
!     Reactor 1 lower, upper rounding radii
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: REACTOR1_RR1,REACTOR1_RR2
!     Reactor 1 lower, upper rounding angles
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: REACTOR1_THETA1,REACTOR1_THETA2
!     Normal vector components for plane (pre-defined quadrics)
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: n_x,n_y,n_z
!     A-matrices of the quadrics
      DOUBLE PRECISION, DIMENSION(3,3,DIM_QUADRIC) :: A_QUADRIC
!     Translation-matrices of the quadrics
      DOUBLE PRECISION, DIMENSION(1,3,DIM_QUADRIC) :: T_QUADRIC
!     Clipping range  of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: clip_xmin,clip_xmax,clip_ymin,clip_ymax,clip_zmin,clip_zmax
!     Piecewise range  of the quadrics
      DOUBLE PRECISION, DIMENSION(DIM_QUADRIC) :: piece_xmin,piece_xmax,piece_ymin,piece_ymax,piece_zmin,piece_zmax
!     Clip flag
      LOGICAL, DIMENSION(DIM_QUADRIC) :: FLUID_IN_CLIPPED_REGION
!     Boundary condition ID
      INTEGER :: BC_ID_Q(DIM_QUADRIC)
!     Maximum number of groups
      INTEGER, PARAMETER :: DIM_GROUP = 50
!     Number of groups
      INTEGER :: N_GROUP
!     Number of quadric in each group
      INTEGER,DIMENSION(DIM_GROUP) :: GROUP_SIZE
!     Quadric ID list in each group
      INTEGER,DIMENSION(DIM_GROUP,DIM_QUADRIC) :: GROUP_Q
!     Quadric relation in each group
      CHARACTER(LEN=9),DIMENSION(DIM_GROUP) :: GROUP_RELATION
!     Relation between groups
      CHARACTER(LEN=9),DIMENSION(DIM_GROUP) :: RELATION_WITH_PREVIOUS
!     Tolerance for intersection between quadrics and planes
      DOUBLE PRECISION :: TOL_F
!     Maximum number of iterations while finding intersection between geometry and grid
      INTEGER :: ITERMAX_INT

    CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Function name: CROSS_PRODUCT                                        C
!  Purpose: Performs the cross product between two vectors             C
!           C = A x B                                                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      FUNCTION CROSS_PRODUCT(A,B)

        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(3) :: CROSS_PRODUCT
        DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: A,B

        CROSS_PRODUCT(1) = A(2) * B(3) - A(3) * B(2)
        CROSS_PRODUCT(2) = A(3) * B(1) - A(1) * B(3)
        CROSS_PRODUCT(3) = A(1) * B(2) - A(2) * B(1)

      END FUNCTION CROSS_PRODUCT

    END MODULE quadric
