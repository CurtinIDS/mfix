!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_INLET                                              !
!                                                                      !
!  Purpose: Common elements needed for the des mass inflow boundary    !
!  condition.                                                          !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      MODULE DES_BC

      USE param, only: dimension_bc, dim_m

      INTEGER :: DEM_BCMI
      INTEGER :: DEM_BCMO

      LOGICAL DEM_MIO  ! either inlet or outlet exists

! Map between DEM MI/MO IDs and the user input BC index.
      INTEGER :: DEM_BCMI_MAP(DIMENSION_BC)
      INTEGER :: DEM_BCMO_MAP(DIMENSION_BC)

! This array contains integers representing the mass/solid phase indices
! present at a specific boundary condtion in proportion to their
! respective number fraction at the inlet (i.e., it represents the
! particle number distribution of incoming solids at the inlet).  The
! array is scaled in size according to the parameter NUMFRAC_LIMIT.
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: DEM_BC_POLY_LAYOUT

! Particle injection time scale; used when pi_factor > 1 to keep track
! of time needed for next injection
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DEM_MI_TIME

! Logical that can be flagged in the mfix.dat file to force the inlet
! to operate with an ordered boundary condition.  This may be useful
! during long simulations or if the inlet appears to be taking a long
! time to randomly place particles.
      LOGICAL :: FORCE_ORD_BC

! Particle injection factor; how many solid time steps (dtsolid) pass
! before the next injection of a particle. if pi_count is greater than
! 1, then pi_factor is set to 1 (i.e. multiple particles enter every
! solids time step).
      INTEGER, DIMENSION(:), ALLOCATABLE :: PI_FACTOR   !(DES_BCMI)

! Particle injection count (injection number); how many particles are
! injected in one solids time step. pi_count is set to one if
! less than 1 particle enters per solids time step.
      INTEGER, DIMENSION(:), ALLOCATABLE :: PI_COUNT   !(DES_BCMI)


! Limit on the total number of divisions (fineness) used to represent
! the particle number distribution at an inlet.
      INTEGER, PARAMETER :: NUMFRAC_LIMIT = 10000


! the dimension of this variable is equal to the number of grid
! cells in the inlet edge/face
      TYPE DEM_MI_
! Array position of next seed location.
         INTEGER :: VACANCY
! Number of positions in the layout grid.
         INTEGER :: OCCUPANTS
! Flag for polydisperse inlets.
         LOGICAL :: POLYDISPERSE
! Uniform grid dimension (width and height).
         DOUBLE PRECISION :: WINDOW
! Offset for placing particles in ghost cell.
         DOUBLE PRECISION :: OFFSET
! Fluid cell index associated with each grid. (I/J/K like)
         INTEGER :: L
         INTEGER, ALLOCATABLE :: W(:)
         INTEGER, ALLOCATABLE :: H(:)
! Spatial location of each grid cell's lower, bottom corder.
         DOUBLE PRECISION, ALLOCATABLE :: P(:)
         DOUBLE PRECISION, ALLOCATABLE :: Q(:)
! The rank of the owning process owning the indexed grid cell.
         INTEGER, ALLOCATABLE :: OWNER(:)
      END TYPE DEM_MI_

! Construct an array of integers in values from 1 to a calculated factor
! in a random order, which is used when placing new particles.
!      TYPE(DEM_MI_DATA), DIMENSION(:), ALLOCATABLE :: MI_ORDER

! Array linking all of the reaction data.
      TYPE(DEM_MI_), DIMENSION(:), TARGET, ALLOCATABLE :: DEM_MI

      INTEGER, ALLOCATABLE :: DEM_BCMO_IJKSTART(:)
      INTEGER, ALLOCATABLE :: DEM_BCMO_IJKEND(:)

      INTEGER, ALLOCATABLE :: DEM_BCMO_IJK(:)


      INTEGER, ALLOCATABLE :: DEM_BCMI_IJKSTART(:)
      INTEGER, ALLOCATABLE :: DEM_BCMI_IJKEND(:)

      INTEGER, ALLOCATABLE :: DEM_BCMI_IJK(:)


!----------------------------------------------------------------------!


! DES specification for solids phase velocity for WALL boundary
! conditions. The current setup is fairly limited. The specified
! boundary velocities are assigned to the indicated wall where a wall
! corresponds to one of the six planes in a cubic domain. Each wall
! corresponds to a number as follows west=1, east=2, bottom=3, top=4,
! south=5, north=6. See cfwallposvel for details. To specify a y or z
! velocity to the west wall set des_bc_vw_s(1,M) or des_bc_ww_s(1,M),
! respectively (note an x velocity is not valid for a west or east wall).
! Since these are user input, they are allocated here with a constant
! preset size, but their actual size is represented by &
! (nwalls, des_mmax)
      DOUBLE PRECISION DES_BC_Uw_s(DIMENSION_BC, DIM_M)
      DOUBLE PRECISION DES_BC_Vw_s(DIMENSION_BC, DIM_M)
      DOUBLE PRECISION DES_BC_Ww_s(DIMENSION_BC, DIM_M)

      CONTAINS
!----------------------------------------------------------------------!
!  Function to exclude cells from DEM mass inlet.                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION EXCLUDE_DEM_MI_CELL(lI, lJ, lK)

      use functions, only: FUNIJK
      use functions, only: FLUID_AT
      use functions, only: IS_ON_myPE_plus2layers

      use compar, only: DEAD_CELL_AT

! Indicies of cell to check
      INTEGER, INTENT(IN) :: lI, lJ, lK
! Local value for IJK
      INTEGER :: IJK

      EXCLUDE_DEM_MI_CELL = .TRUE.

      IF(.NOT.IS_ON_myPE_plus2layers(lI,lJ,lK)) RETURN
      IF(DEAD_CELL_AT(lI,lJ,lK)) RETURN
      IJK = FUNIJK(lI,lJ,lK)
      IF(.NOT.FLUID_AT(IJK)) RETURN

      EXCLUDE_DEM_MI_CELL = .FALSE.
      RETURN
      END FUNCTION EXCLUDE_DEM_MI_CELL

      END MODULE DES_BC

