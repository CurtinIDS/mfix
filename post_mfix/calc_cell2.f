!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_CELL2 (RMIN, REACTOR_LOC,D_DIR,N_DIR,CELL_LOC)    C
!  Purpose: calculate the cell index for a reactor location            C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:  None                                         C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC, CELL_START, CELL_END                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_CELL2 (RMIN, REACTOR_LOC,D_DIR,N_DIR,CELL_LOC)
!
        Use param
        Use param1
!
! passed arguments
!
!    RMIN       - starting value of the axis -- XMIN need not be zero,
!                   YMIN and ZMIN are assumed to be zero
!    REACTOR_LOC - location along one of the axis for which the cell
!                  index is to be found
!    D_DIR       - cell lengths (DX,DY,DZ)
!    N_DIR       - number of cells in this direction (IMAX2,JMAX2,KMAX2)
!    CELL_LOC    - cell index corresponding to REACTOR_LOC
!
! local variables
!    LC         -  loop counter
!    CELL_START -  start coordinate for cell
!    CELL_END   -  end   coordinate for cell
!
      IMPLICIT NONE
      INTEGER           N_DIR , CELL_LOC , LC
      DOUBLE PRECISION  REACTOR_LOC , D_DIR(*) , CELL_START , CELL_END, RMIN
!
      CELL_LOC = - 1
      CELL_START = RMIN
      DO 100 LC = 2,N_DIR+1
         CELL_END   = CELL_START + D_DIR(LC)
         IF (REACTOR_LOC .GE.CELL_START .AND. &
                                 REACTOR_LOC .LE. CELL_END)  THEN
            CELL_LOC = LC
            RETURN
         END IF
      CELL_START = CELL_END
100   CONTINUE
!
      RETURN
      END
