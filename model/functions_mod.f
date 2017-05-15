MODULE functions

   use compar
   use discretelement
   use geometry
   use indices
   use param1

   implicit none

   PRIVATE


! Functions for generating IJK indices for indicated basis:
!---------------------------------------------------------------------//
! INTEGER :: FUNIJK, FUNIJK_0
! INTEGER :: FUNIJK_PROC
! INTEGER :: FUNIJK_GL
! INTEGER :: FUNIJK_IO

   PUBLIC :: funijk, funijk_0, FUNIJK_PROC, FUNIJK_GL, FUNIJK_IO
   PUBLIC :: BOUND_FUNIJK

! Logical functions to determine whether index is on my PE's domain or
! indicated subset:
!---------------------------------------------------------------------//
! LOGICAL :: IS_ON_myPE_owns
! LOGICAL :: IS_ON_myPE_wobnd
! LOGICAL :: IS_ON_myPE_plus1layer
! LOGICAL :: IS_ON_myPE_plus2layers

   PUBLIC :: IS_ON_myPE_OWNS, IS_ON_myPE_wobnd
   PUBLIC :: IS_ON_myPE_plus1layer, IS_ON_myPE_plus2layers

! Functions for calculating indicated directional shift in given IJK
! index. This will return the ijk index of the computational cell
! corresponding to the indicated shift when that computational cell
! is NOT a wall cell. If the computational cell is a wall cell then
! this will return its own ijk index. For example, east_of will return
! IPJK when IPJK is a fluid or flow at cell. However, if IPJK is a
! wall cell east_of will return IJK.
!---------------------------------------------------------------------//
! INTEGER :: EAST_OF, EAST_OF_0
! INTEGER :: WEST_OF, WEST_OF_0
! INTEGER :: SOUTH_OF, SOUTH_OF_0
! INTEGER :: NORTH_OF, NORTH_OF_0
! INTEGER :: BOTTOM_OF, BOTTOM_OF_0
! INTEGER :: TOP_OF, TOP_OF_0

   PUBLIC :: EAST_OF, NORTH_OF, TOP_OF
   PUBLIC :: WEST_OF, SOUTH_OF, BOTTOM_OF
   PUBLIC :: EAST_OF_0, NORTH_OF_0, TOP_OF_0
   PUBLIC :: WEST_OF_0, SOUTH_OF_0, BOTTOM_OF_0


! Functions for calculating indicated directional shift in given IJK
! index. This will generally return the ijk index of the computational
! cell corresponding to the indicated shift regardless of the wall
! status of that computational cell. It may not return corner cells
! unless the ijk cell itself is a corner cell.
!---------------------------------------------------------------------//
! INTEGER :: IP_OF, IP_OF_0
! INTEGER :: IM_OF, IM_OF_0
! INTEGER :: JP_OF, JP_OF_0
! INTEGER :: JM_OF, JM_OF_0
! INTEGER :: KP_OF, KP_OF_0
! INTEGER :: KM_OF, KM_OF_0

   PUBLIC :: IM_OF, IP_OF, JM_OF, JP_OF, KM_OF, KP_OF
   PUBLIC :: IM_OF_0, IP_OF_0, JM_OF_0, JP_OF_0, KM_OF_0, KP_OF_0

! Logical functions to identify indicated condition:
!---------------------------------------------------------------------//
! logical function to identify various fluid/flow cells
! LOGICAL :: FLUID_AT
! LOGICAL :: P_FLOW_AT
! LOGICAL :: P_OUTFLOW_AT
! LOGICAL :: MASS_OUTFLOW_AT
! LOGICAL :: OUTFLOW_AT
! LOGICAL :: FLOW_AT
! LOGICAL :: FLUIDorP_FLOW_AT

   PUBLIC :: FLUID_AT, P_FLOW_AT, P_OUTFLOW_AT
   PUBLIC :: MASS_OUTFLOW_AT, OUTFLOW_AT
   PUBLIC :: FLOW_AT, FLUIDorP_FLOW_AT


! logical function to identify various wall cells
! LOGICAL :: WALL_AT
! LOGICAL :: NS_WALL_AT, FS_WALL_AT, PS_WALL_AT
! LOGICAL :: DEFAULT_WALL_AT
! LOGICAL :: WALL_ICBC_FLAG

   PUBLIC :: WALL_AT, NS_WALL_AT, FS_WALL_AT, PS_WALL_AT
   PUBLIC :: DEFAULT_WALL_AT, WALL_ICBC_FLAG

! Logical function to identify a cyclic cell and different
! cyclic flow boundaries
! LOGICAL :: CYCLIC_AT
! LOGICAL :: CYCLIC_AT_E, CYCLIC_AT_N, CYCLIC_AT_T

   PUBLIC :: CYCLIC_AT, CYCLIC_AT_E, CYCLIC_AT_N, CYCLIC_AT_T


! logical function to identify different flow at boundaries
! LOGICAL :: FLOW_AT_E, FLOW_AT_N, FLOW_AT_T
! LOGICAL :: PFLOW_AT_E, PFLOW_AT_N, PFLOW_AT_T
! LOGICAL :: MFLOW_AT_E, MFLOW_AT_N, MFLOW_AT_T

   PUBLIC :: FLOW_AT_E, FLOW_AT_N, FLOW_AT_T
   PUBLIC :: PFLOW_AT_E, PFLOW_AT_N, PFLOW_AT_T
   PUBLIC :: MFLOW_AT_E, MFLOW_AT_N, MFLOW_AT_T

! Logical functions to identify different impermeable/semipermeable
! surface boundaries (specific type of internal surface)
! LOGICAL :: IP_AT_E, IP_AT_N, IP_AT_T
! LOGICAL :: SIP_AT_E, SIP_AT_N, SIP_AT_T
! LOGICAL :: SP_AT_E, SP_AT_N, SP_AT_T

   PUBLIC :: IP_AT_E, IP_AT_N, IP_AT_T
   PUBLIC :: SIP_AT_E, SIP_AT_N, SIP_AT_T
   PUBLIC :: SP_AT_E, SP_AT_N, SP_AT_T

! Logical functions concerning general internal surfaces
! LOGICAL :: IS_AT_E, IS_AT_N, IS_AT_T
! LOGICAL :: NO_IS_AT_E, NO_IS_AT_N, NO_IS_AT_T
! Integer function to return internal surface ID
! INTEGER :: IS_ID_AT_E, IS_ID_AT_N, IS_ID_AT_T


   PUBLIC :: IS_AT_E, IS_AT_N, IS_AT_T
   PUBLIC :: NO_IS_AT_E, NO_IS_AT_N, NO_IS_AT_T
   PUBLIC :: IS_ID_AT_E,  IS_ID_AT_N, IS_ID_AT_T

! Additional functions
!---------------------------------------------------------------------//
! DOUBLE PRECISION :: ZMAX
! INTEGER FUNCTION :: FUNLM

   PUBLIC :: ZMAX, FUNLM


! DEM Functions
!---------------------------------------------------------------------//
   PUBLIC :: IS_NONEXISTENT, IS_NORMAL, IS_ENTERING, IS_EXITING
   PUBLIC :: IS_GHOST, IS_ENTERING_GHOST, IS_EXITING_GHOST
   PUBLIC :: IS_ANY_GHOST
   PUBLIC :: SET_NONEXISTENT, SET_NORMAL, SET_ENTERING, SET_EXITING
   PUBLIC :: SET_GHOST, SET_ENTERING_GHOST, SET_EXITING_GHOST


CONTAINS

  INCLUDE 'functions.inc'

END MODULE functions
