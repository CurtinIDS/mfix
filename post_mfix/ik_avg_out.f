!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: IK_AVG_OUT(VAR_INDEXN,LOC_YN,LUNIT,TIME_REAL,M_USE)    C
!  Purpose: driver for "time, variable(@specified height - averaged    C
!                       over X & Z)" output                            C
!                                                                      C
!  Author: P. Nicoletti                               Date: 15-JUL-92  C
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
!                        W_s, ROP_s, T_g, T_s1, T_s2, VOL, VOL_U    C
!                        VOL_V, VOL_W                                  C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: AVG_OUT                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE IK_AVG_OUT(VAR_INDEXN,LOC_YN,LUNIT,TIME_REAL,M_USE)
!
      Use param
      Use param1
      Use fldvar
      Use indices
      Use run
      Use physprop
      Use geometry
      Use post3d
!
      IMPLICIT NONE
!
!     Passed arguments
!
      INTEGER LUNIT , VAR_INDEXN , LOC_YN , M_USE
      REAL    TIME_REAL(*)
!
!     local variables
!
      REAL    AVG_OUT
!
      IF (VAR_INDEXN.EQ.01) THEN
         CALL IK_AVG(EP_G,VOL,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(1) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.02) THEN
         CALL IK_AVG(P_G,VOL,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(2) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.03) THEN
         CALL IK_AVG(P_STAR,VOL,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(2) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.04) THEN
         CALL IK_AVG(U_G,VOL_U,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(3) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.05) THEN
         CALL IK_AVG(V_G,VOL_V,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(3) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.06) THEN
         CALL IK_AVG(W_G,VOL_W,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(3) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.07) THEN
         CALL IK_AVG(U_S(1,M_USE),VOL_U,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(4) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.08) THEN
         CALL IK_AVG(V_S(1,M_USE),VOL_V,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(4) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.09) THEN
         CALL IK_AVG(W_S(1,M_USE),VOL_W,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(4) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.10) THEN
         CALL IK_AVG(ROP_S(1,M_USE),VOL,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(5) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.11) THEN
         CALL IK_AVG(T_G,VOL,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(6) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.12) THEN
         CALL IK_AVG(T_S(1,M_USE),VOL,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(6) , AVG_OUT
      ELSE IF (VAR_INDEXN.EQ.13) THEN
         CALL IK_AVG(T_S(1,M_USE),VOL,AVG_OUT,LOC_YN)
         WRITE (39+LUNIT,*) TIME_REAL(6) , AVG_OUT
      END IF
!
!
      RETURN
      END
