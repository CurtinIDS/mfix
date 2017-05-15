!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_TIME (VAR_INDEX,LOC_X,LOC_Y,LOC_Z,M_USE,           C
!                         L,TIME_REAL)                                 C
!  Purpose: write out a "time,var(loc_x,loc_y,loc_z)" record to a file C
!                                                                      C
!  Author: P. Nicoletti                               Date: 20-FEB-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: EP_g, P_g, P_star, U_g, V_g, W_g, U_s, V_s    C
!                        W_s, ROP_s, T_g, T_s1, T_s2                   C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: L1                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_TIME(VAR_INDEX,LOC_X,LOC_Y,LOC_Z,M_USE,L, &
                          TIME_REAL)
!
 !
      Use param
      Use param1
      Use fldvar
      Use run
      Use physprop
      Use indices
      Use geometry
      Use compar
      Use functions
      IMPLICIT NONE
!
!     Passed Arguments
!
      REAL    TIME_REAL
      INTEGER VAR_INDEX,LOC_X,LOC_Y,LOC_Z,L,M_USE
!
!     Local Variables
!
      INTEGER IJK
!
      IJK = FUNIJK(LOC_X,LOC_Y,LOC_Z)
      IF (VAR_INDEX.EQ.01) WRITE (39+L,*) TIME_REAL,EP_g(IJK)
      IF (VAR_INDEX.EQ.02) WRITE (39+L,*) TIME_REAL,P_g(IJK)
      IF (VAR_INDEX.EQ.03) WRITE (39+L,*) TIME_REAL,P_star(IJK)
      IF (VAR_INDEX.EQ.04) WRITE (39+L,*) TIME_REAL,U_g(IJK)
      IF (VAR_INDEX.EQ.05) WRITE (39+L,*) TIME_REAL,V_g(IJK)
      IF (VAR_INDEX.EQ.06) WRITE (39+L,*) TIME_REAL,W_g(IJK)
      IF (VAR_INDEX.EQ.07) WRITE (39+L,*) TIME_REAL,U_s(IJK,M_USE)
      IF (VAR_INDEX.EQ.08) WRITE (39+L,*) TIME_REAL,V_s(IJK,M_USE)
      IF (VAR_INDEX.EQ.09) WRITE (39+L,*) TIME_REAL,W_s(IJK,M_USE)
      IF (VAR_INDEX.EQ.10) WRITE (39+L,*) TIME_REAL,ROP_s(IJK,M_USE)
      IF (VAR_INDEX.EQ.11) WRITE (39+L,*) TIME_REAL,T_g(IJK)
      IF (VAR_INDEX.EQ.12) WRITE (39+L,*) TIME_REAL,T_s(IJK,M_USE)
      IF (VAR_INDEX.EQ.13) WRITE (39+L,*) TIME_REAL,T_s(IJK,M_USE)
!
      RETURN
      END
