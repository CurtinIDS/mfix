
      SUBROUTINE ALLOCATE_DUMMY_CUT_CELL_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: ALLOCATE_DUMMY_CUT_CELL_ARRAYS
!  Purpose: allocate dummy cut cell arrays used to maintain original mfix code
!  when CARTESIAN_GRID = .FALSE.
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      Use indices

      USE cutcell

      IMPLICIT NONE

      Allocate(  CUT_TREATMENT_AT (DIMENSION_3) )
      Allocate(  CUT_U_TREATMENT_AT (DIMENSION_3) )
      Allocate(  CUT_V_TREATMENT_AT (DIMENSION_3) )
      Allocate(  CUT_W_TREATMENT_AT (DIMENSION_3) )

      Allocate(  CUT_CELL_AT (DIMENSION_3) )
      Allocate(  CUT_U_CELL_AT (DIMENSION_3) )
      Allocate(  CUT_V_CELL_AT (DIMENSION_3) )
      Allocate(  CUT_W_CELL_AT (DIMENSION_3) )

      Allocate(  BC_ID (DIMENSION_3) )
      Allocate(  BC_U_ID (DIMENSION_3) )
      Allocate(  BC_V_ID (DIMENSION_3) )
      Allocate(  BC_W_ID (DIMENSION_3) )

      Allocate(  BLOCKED_U_CELL_AT (DIMENSION_3) )
      Allocate(  BLOCKED_V_CELL_AT (DIMENSION_3) )
      Allocate(  BLOCKED_W_CELL_AT (DIMENSION_3) )

      CUT_TREATMENT_AT = .FALSE.
      CUT_U_TREATMENT_AT = .FALSE.
      CUT_V_TREATMENT_AT = .FALSE.
      CUT_W_TREATMENT_AT = .FALSE.

      CUT_CELL_AT = .FALSE.
      CUT_U_CELL_AT = .FALSE.
      CUT_V_CELL_AT = .FALSE.
      CUT_W_CELL_AT = .FALSE.

      BLOCKED_U_CELL_AT = .FALSE.
      BLOCKED_V_CELL_AT = .FALSE.
      BLOCKED_W_CELL_AT = .FALSE.

      BC_ID = 0
      BC_U_ID = 0
      BC_V_ID = 0
      BC_W_ID = 0



      RETURN
      END SUBROUTINE ALLOCATE_DUMMY_CUT_CELL_ARRAYS


