!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: geometry                                               C
!  Purpose: Common block containing geometry and discretization data   C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE geometry

      Use param, only: DIM_I, DIM_J, DIM_K


! Coordinates: CARTESIAN, CYLINDRICAL
      CHARACTER(LEN=16)     COORDINATES

! Indicates whether x or r direction is not considered
      LOGICAL :: NO_I
! Indicates whether x or r direction is considered
      LOGICAL :: DO_I
! Indicates whether y direction is not considered
      LOGICAL :: NO_J
! Indicates whether y direction is considered
      LOGICAL :: DO_J
! Indicates whether z or theta direction is not considered
      LOGICAL :: NO_K
! Indicates whether z or theta direction is considered
      LOGICAL :: DO_K

! Reactor length in the x or r direction
      DOUBLE PRECISION :: XLENGTH
! Reactor length in the y direction
      DOUBLE PRECISION :: YLENGTH
! Reactor length in the z or theta direction
      DOUBLE PRECISION :: ZLENGTH

! Starting index in the x or r direction
      INTEGER :: IMIN1
! Starting index in the y direction
      INTEGER :: JMIN1
! Starting index in the z or theta direction
      INTEGER :: KMIN1

! Number of cells in the x or r direction
      INTEGER :: IMAX
! Number of cells in the y direction
      INTEGER :: JMAX
! Number of cells in the z or theta direction
      INTEGER :: KMAX

! Number of cells in the x or r direction + 1
      INTEGER :: IMAX1
! Number of cells in the y direction + 1
      INTEGER :: JMAX1
! Number of cells in the z or theta direction + 1
      INTEGER :: KMAX1

! Number of cells in the x or r direction + 2
      INTEGER :: IMAX2
! Number of cells in the y direction + 2
      INTEGER :: JMAX2
! Number of cells in the z or theta direction + 2
      INTEGER :: KMAX2

! Cell sizes in the x or r direction
      DOUBLE PRECISION :: DX (0:DIM_I)
! Cell sizes in the y direction
      DOUBLE PRECISION :: DY (0:DIM_J)
! Cell sizes in the z or theta direction
      DOUBLE PRECISION :: DZ (0:DIM_K)
! Starting value of X.  This quantity is useful for
! simulating an annular cylindrical region.
      DOUBLE PRECISION :: XMIN

! IMAX2 * JMAX2
      INTEGER :: IJMAX2
! IMAX2 * JMAX2 * KMAX2
      INTEGER :: IJKMAX2
! IMAX2 * JMAX2 * KMAX2
      INTEGER :: IJKMAX3
! IJMAX2 + 1
      INTEGER :: IJKMIN1
! IJKMAX2 - IJMAX2
      INTEGER :: IJKMAX1

! For discretization in parallel
      INTEGER :: IMIN2, JMIN2, KMIN2
      INTEGER :: IMIN3, JMIN3, KMIN3
      INTEGER :: IMAX3, JMAX3, KMAX3

! For 4th order discretization in parallel
      INTEGER :: IMIN4, JMIN4, KMIN4
      INTEGER :: IMAX4, JMAX4, KMAX4
      INTEGER :: IJKMAX4, IJKMIN4

! Cell flags.
      INTEGER, DIMENSION(:), ALLOCATABLE :: FLAG
! Cell flags with 3rd layer.
      INTEGER, DIMENSION(:), ALLOCATABLE :: FLAG3
! Flag for the East surface
      INTEGER, DIMENSION(:), ALLOCATABLE :: FLAG_E
! Flag for North surface
      INTEGER, DIMENSION(:), ALLOCATABLE :: FLAG_N
! Flag for Top surface
      INTEGER, DIMENSION(:), ALLOCATABLE :: FLAG_T
! Cell flags (bc/ic conditions)
! Allocatable type causes PG internal error, Ed's soln: pointers
!      CHARACTER(LEN=3), DIMENSION(:), ALLOCATABLE :: ICBC_FLAG
      character(LEN=3),  dimension(:), pointer :: icbc_flag

! 1 / dx_i
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oDX
! 1 / dy_j
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oDY
! 1 / dz_k
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oDZ

! 1 / dx_i+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oDX_E
! 1 / dy_j+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oDY_N
! 1 / dz_k+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oDZ_T

! Radial location at cell center (x_i).
! X = 1 in Cartesian coordinates.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X
! Azimuthal location at cell center (z_k).
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z
! Radial location at East face (x_i+1/2).
! X_E = 1 in Cartesian coordinates.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  X_E
!  Azimuthal location at top face (z_k+1/2).
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z_T

! Reciprocal of radial location at cell center (1/x_i).
! oX = 1 in Cartesian coordinates.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oX
! Reciprocal of radial location at East face (1/x_i+1/2).
! oX_E = 1 in Cartesian coordinates.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  oX_E


!  one or more periodic boundary condition is used
      LOGICAL :: CYCLIC
! Variable to flag periodic boundary condition in X
      LOGICAL :: CYCLIC_X
! Variable to flag periodic boundary condition in Y
      LOGICAL :: CYCLIC_Y
! Variable to flag periodic boundary condition in Z
      LOGICAL :: CYCLIC_Z
! Variable to flag periodic bc with pressure drop in X
      LOGICAL :: CYCLIC_X_PD
! Variable to flag periodic bc with pressure drop in Y
      LOGICAL :: CYCLIC_Y_PD
! Variable to flag periodic bc with pressure drop in Z
      LOGICAL :: CYCLIC_Z_PD
! Variable to flag periodic bc with mass flux in X
      LOGICAL :: CYCLIC_X_MF
! Variable to flag periodic bc with mass flux in Y
      LOGICAL :: CYCLIC_Y_MF
! Variable to flag periodic bc with mass flux in Z
      LOGICAL :: CYCLIC_Z_MF

! Variable to flag cylindrical coordinates
      LOGICAL :: CYLINDRICAL

! Variables for cylindrical_2d simulation
!------------------------------------------------------------>>>
! Turn on the cylindrical_2d simulation
      logical :: CYLINDRICAL_2D
! Variables for cylindrical_2d simulation
! Half width of the plate in term of cell count
      integer :: I_CYL_NUM
! Variables for cylindrical_2d simulation
! Cell number used to smooth the transition from plate to wedge
      integer :: I_CYL_TRANSITION
! For cylindrical_2d simulation
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cyl_X
! For cylindrical_2d simulation
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cyl_X_E
!------------------------------------------------------------<<<

! Factor for x direction averaging of U velocity: FX_i
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FX
! 1 - FX_i
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FX_bar
! Factor for x direction averaging of scalars: FX_i+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FX_E
! 1 - FX_i+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FX_E_bar

! Factor for y direction averaging of scalars: FY_j+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FY_N
! 1 - FY_j+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FY_N_bar

! Factor for z direction averaging of scalars: FZ_k+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FZ_T
! 1 - FZ_k+1/2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FZ_T_bar

! East face area - scalar cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AYZ
! North face area - scalar cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AXZ
! Top face area - scalar cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AXY
! Cell volume - scalar cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VOL

! Total volume of cell's DES stencil neighbors
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VOL_SURR

! East face area - U cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AYZ_U
! North face area - U cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AXZ_U
! Top face area - U cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AXY_U
! Cell volume - U cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VOL_U

! East face area - V cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AYZ_V
! North face area - V cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AXZ_V
! Top face area - V cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AXY_V
! Cell volume - V cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VOL_V

! East face area - W cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AYZ_W
! North face area - W cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AXZ_W
! Top face area - W cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AXY_W
! Cell volume - W cell
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VOL_W

! Indicates whether to use loop over "core" cells (with same
! CELL_CLASS) in LEQSOL_MOD for vectorization/performance
      LOGICAL :: USE_CORECELL_LOOP

      INTEGER :: CORE_ISTART, CORE_IEND
      INTEGER :: CORE_JSTART, CORE_JEND
      INTEGER :: CORE_KSTART, CORE_KEND

      END MODULE geometry
