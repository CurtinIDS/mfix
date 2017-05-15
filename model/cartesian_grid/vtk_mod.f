      MODULE vtk

      Use param1
! Maximum number of solids phases.
      use param, only: DIM_M
! Maximum number of gas phase species
      use param, only: DIM_N_g
! Maximum number of solids phase species
      use param, only: DIM_N_s
! Maximum number of scalar equations
      use param, only: DIM_Scalar
! Maximum number of DEM solids phase species
use param, only: DIMENSION_N_S



      INTEGER NUMBER_OF_CELLS
      INTEGER NUMBER_OF_CUT_CELLS
      INTEGER NUMBER_OF_BLOCKED_CELLS
      INTEGER NUMBER_OF_STANDARD_CELLS
      INTEGER NUMBER_OF_VTK_CELLS
      INTEGER NUMBER_OF_VTK_PARTS

      LOGICAL             :: WRITE_VTK_FILES
      LOGICAL             :: TIME_DEPENDENT_FILENAME
      LOGICAL             :: RESET_FRAME_AT_TIME_ZERO=.TRUE.
      CHARACTER (LEN=255)  :: VTU_DIR
      CHARACTER (LEN=255)  :: VTK_FILENAME,FRAME_CHAR,VTU_FILENAME,PVD_FILENAME,PVTU_FILENAME
      CHARACTER (LEN=255)  :: VTU_FRAME_FILENAME='VTU_FRAME_INDEX.TXT'
      CHARACTER (LEN=512) :: BUFFER

      CHARACTER (LEN=1), PARAMETER  :: END_REC = CHAR(10)

      INTEGER :: BOUNDARY_UNIT=122
      INTEGER :: VTK_UNIT=123,VTU_UNIT=124,PVD_UNIT=125
      INTEGER :: PVTU_UNIT=126,VTU_FRAME_UNIT=127
      INTEGER, PARAMETER :: DIM_VTK_VAR = 20
      INTEGER, DIMENSION(DIM_VTK_VAR) :: VTK_VAR


      INTEGER :: POLY_COUNTER,NUMBER_OF_POINTS

      integer, allocatable :: CLEANED_CONNECTIVITY(:,:)
      REAL, allocatable :: COORDS_OF_POINTS(:,:)

      integer, allocatable :: GLOBAL_I_OF(:)
      integer, allocatable :: GLOBAL_J_OF(:)
      integer, allocatable :: GLOBAL_K_OF(:)
      integer, allocatable :: GLOBAL_CONNECTIVITY(:,:)
      integer, allocatable :: GLOBAL_CLEANED_CONNECTIVITY(:,:)
      integer, allocatable :: GLOBAL_NUMBER_OF_NODES(:)

      REAL, allocatable :: GLOBAL_COORDS_OF_POINTS(:,:)

      LOGICAL, allocatable :: GLOBAL_INTERIOR_CELL_AT(:)
      LOGICAL, allocatable :: GLOBAL_BLOCKED_CELL_AT(:)
      LOGICAL, allocatable :: GLOBAL_STANDARD_CELL_AT(:)
      LOGICAL, allocatable :: GLOBAL_CUT_CELL_AT(:)
      LOGICAL, allocatable :: GLOBAL_SNAP(:)
      DOUBLE PRECISION, allocatable :: GLOBAL_F_AT(:)

      double precision, allocatable :: GLOBAL_X_NEW_POINT(:)
      double precision, allocatable :: GLOBAL_Y_NEW_POINT(:)
      double precision, allocatable :: GLOBAL_Z_NEW_POINT(:)

      INTEGER :: GLOBAL_NUMBER_OF_NEW_POINTS


      LOGICAL :: GLOBAL_VAR_ALLOCATED

      LOGICAL :: GRID_INFO_PRINTED_ON_SCREEN

      LOGICAL :: WRITE_ANI_CUTCELL


      INTEGER :: VTU_offset

      LOGICAL, allocatable :: BELONGS_TO_VTK_SUBDOMAIN(:)
      ! LOGICAL, allocatable :: PART_BELONGS_TO_VTK_SUBDOMAIN(:)

      INTEGER, PARAMETER :: DIMENSION_VTK = 100
! Max size for VTK_nRR
      INTEGER, PARAMETER :: VTK_nRRmax = 100

! Current VTK region
      INTEGER :: VTK_REGION

! Time interval at which vtk files are saved
      DOUBLE PRECISION :: VTK_DT(DIMENSION_VTK)

! Flag to distiguish between regular file (VTK_DBG_FILE=.FALSE., default value)
!                        and debug file (VTK_DBG_FILE=.TRUE.)
      LOGICAL :: VTK_DBG_FILE(DIMENSION_VTK)

! Current vtk time
      DOUBLE PRECISION :: VTK_TIME(DIMENSION_VTK)

! Type of data in vtk region: 'C':cell data, 'P': particle data
      CHARACTER(LEN=1) :: VTK_DATA(DIMENSION_VTK)

! FRAME index of vtk file
      INTEGER :: FRAME(DIMENSION_VTK)

! PVD file initialization flag
      LOGICAL :: PVD_FILE_INITIALIZED(DIMENSION_VTK)=.FALSE.

! Logical variable to determine whether an vtk region is defined
      LOGICAL :: VTK_DEFINED (DIMENSION_VTK)

! VTK region West face, X-coordinate
      DOUBLE PRECISION :: VTK_X_w (DIMENSION_VTK)

! VTK region East face, X-coordinate
      DOUBLE PRECISION :: VTK_X_e (DIMENSION_VTK)

! VTK region South face, Y-coordinate
      DOUBLE PRECISION :: VTK_Y_s (DIMENSION_VTK)

! VTK region North face, Y-coordinate
      DOUBLE PRECISION :: VTK_Y_n (DIMENSION_VTK)

! VTK region Bottom face, Z-coordinate
      DOUBLE PRECISION :: VTK_Z_b (DIMENSION_VTK)

! VTK region Top face, Z-coordinate
      DOUBLE PRECISION :: VTK_Z_t (DIMENSION_VTK)

! VTK number of slices in x-direction
      INTEGER :: VTK_NXS(DIMENSION_VTK)

! VTK number of slices in y-direction
      INTEGER :: VTK_NYS(DIMENSION_VTK)

! VTK number of slices in z-direction
      INTEGER :: VTK_NZS(DIMENSION_VTK)

! VTK slice tolerance
      DOUBLE PRECISION :: VTK_SLICE_TOL(DIMENSION_VTK)

! Flag to write only cut cell data in VTK file
      LOGICAL :: VTK_CUTCELL_ONLY(DIMENSION_VTK)

! VTK filename base
      CHARACTER(LEN=255) :: VTK_FILEBASE(DIMENSION_VTK)

! Gas phase volume fraction
      LOGICAL :: VTK_EP_g (DIMENSION_VTK)

! Gas pressure
      LOGICAL :: VTK_P_g (DIMENSION_VTK)

! Solids pressure
      LOGICAL :: VTK_P_star(DIMENSION_VTK)

! X-component of gas velocity
      LOGICAL :: VTK_U_g(DIMENSION_VTK)

! X-component of solids phase velocity
      LOGICAL :: VTK_U_s(DIMENSION_VTK, DIM_M)

! Y-component of gas velocity
      LOGICAL :: VTK_V_g(DIMENSION_VTK)

! Y-component of solids phase velocity
      LOGICAL :: VTK_V_s(DIMENSION_VTK, DIM_M)

! Z-component of gas velocity
      LOGICAL :: VTK_W_g(DIMENSION_VTK)

! Z-component of solids phase velocity
      LOGICAL :: VTK_W_s(DIMENSION_VTK, DIM_M)

! Gas velocity vector
      LOGICAL :: VTK_VEL_g(DIMENSION_VTK)

! Solids velocity vector
      LOGICAL :: VTK_VEL_s(DIMENSION_VTK, DIM_M)

! Macroscopic density of solids phases
      LOGICAL :: VTK_ROP_s(DIMENSION_VTK, DIM_M)

! Solids phase volume fraction
      LOGICAL :: VTK_EP_s (DIMENSION_VTK, DIM_M)

! Gas temperature
      LOGICAL :: VTK_T_g(DIMENSION_VTK)

! Solids temperature
      LOGICAL :: VTK_T_s(DIMENSION_VTK, DIM_M)

! Gas species mass fractions
      LOGICAL :: VTK_X_g(DIMENSION_VTK, DIM_N_g)

! Solids species mass fractions
      LOGICAL :: VTK_X_s(DIMENSION_VTK, DIM_M, DIM_N_s)

! Granular temperature
      LOGICAL :: VTK_Theta_m(DIMENSION_VTK, DIM_M)

! Scalar value
      LOGICAL :: VTK_Scalar(DIMENSION_VTK, DIM_scalar)

! Reaction rates
      LOGICAL :: VTK_RRate(DIMENSION_VTK,VTK_nRRmax)

! K & Epsilon values
      LOGICAL :: VTK_K_Turb_G(DIMENSION_VTK)
      LOGICAL :: VTK_E_Turb_G(DIMENSION_VTK)

! Vorticity magnitude
      LOGICAL :: VTK_VORTICITY(DIMENSION_VTK)

! Lambda_2 (vortex core indicator)
      LOGICAL :: VTK_LAMBDA_2(DIMENSION_VTK)

! Grid partition
      LOGICAL :: VTK_PARTITION(DIMENSION_VTK)

! Boundary condition ID
      LOGICAL :: VTK_BC_ID(DIMENSION_VTK)

! Wall distance
      LOGICAL :: VTK_DWALL(DIMENSION_VTK)

! Facet count (DES)
      LOGICAL :: VTK_FACET_COUNT_DES(DIMENSION_VTK)

! Neighboring facets (DES)
      LOGICAL :: VTK_NB_FACET_DES(DIMENSION_VTK)

! Cell IJK index
      LOGICAL :: VTK_IJK(DIMENSION_VTK)

! Cut face normal vector
      LOGICAL :: VTK_NORMAL(DIMENSION_VTK)

! Debug variable
      LOGICAL :: VTK_DEBUG(DIMENSION_VTK,15)

      INTEGER, DIMENSION(DIMENSION_VTK,DIM_VTK_VAR) :: VTK_VARLIST

! VTK Particle selection: 'C': center, 'P': entire particle, 'I': particle intersect
      CHARACTER(LEN=1) :: VTK_SELECT_MODE(DIMENSION_VTK)

! Particle radius
      LOGICAL :: VTK_PART_DIAMETER(DIMENSION_VTK)

! Particle velocity
      LOGICAL :: VTK_PART_VEL(DIMENSION_VTK)

! Particle angular velocity
      LOGICAL :: VTK_PART_ANGULAR_VEL(DIMENSION_VTK)

! Particle orientation
      LOGICAL :: VTK_PART_ORIENTATION(DIMENSION_VTK)

! Particle user-defined variable
      LOGICAL :: VTK_PART_USR_VAR(DIMENSION_VTK,3)

! Particle temperature
      LOGICAL :: VTK_PART_TEMP(DIMENSION_VTK)

! Particle species mass fraction
      LOGICAL :: VTK_PART_X_s(DIMENSION_VTK,100)

! Particle cohesion
      LOGICAL :: VTK_PART_COHESION(DIMENSION_VTK)

      END MODULE vtk

