      MODULE polygon

!     Maximum of the number of polygons that can be read
      INTEGER, PARAMETER          :: DIM_POLYGON = 10
!     Nnumber of polygons
      INTEGER                     :: N_POLYGON
!     Maximum of the number of Vertices per polygon that can be read
      INTEGER, PARAMETER          :: DIM_VERTEX = 500
!     Vertex Coordinates X and Y
      DOUBLE PRECISION, DIMENSION(DIM_POLYGON,DIM_VERTEX) :: X_VERTEX,Y_VERTEX
!     Sign of polygon interior
      DOUBLE PRECISION, DIMENSION(DIM_POLYGON) :: POLY_SIGN
!     Number of vertices for each polygon
      INTEGER, DIMENSION(DIM_POLYGON) :: N_VERTEX
!     Tolerance for polygone edge detection
      DOUBLE PRECISION :: TOL_POLY
!     Boundary condition ID
      INTEGER, DIMENSION(DIM_POLYGON,DIM_VERTEX) :: BC_ID_P

      END MODULE polygon
