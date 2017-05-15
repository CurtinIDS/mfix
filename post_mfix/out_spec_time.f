!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_SPEC_TIME (SPEC_TIME,L,VAR_INDEX,PLOT_TYPE,M_USE)  C
!  Purpose: Output a variable (or 3 velocities) at a specified time    C
!                                                                      C
!  Author: P. Nicoletti                               Date: 01-APR-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMIN1, IMAX1, JMIN1, JMAX2, KMIN1, KMAX1      C
!                        EP_g, P_g, P_star, U_g, V_g, W_g, U_s, V_s    C
!                        W_s, ROP_s, T_g, T_s1, T_s2
!  Variables modified: I,J,K                                           C
!                                                                      C
!  Local variables: L1, NX, NY, NZ                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_SPEC_TIME(SPEC_TIME,L,VAR_INDEX,PLOT_TYPE,M_USE)
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
      REAL      SPEC_TIME
      INTEGER   L, VAR_INDEX , PLOT_TYPE , M_USE
!
!     Local variables
!
      INTEGER   NX, NY, NZ

      INTEGER   I, J, K, IJK
!
      NX = IMAX1 - IMIN1 + 1
      NY = JMAX1 - JMIN1 + 1
      NZ = KMAX1 - KMIN1 + 1
      WRITE (39+L,*) NX , NY , NZ
      DO K = KMIN1,KMAX1
         DO J = JMIN1,JMAX1
            DO I = IMIN1,IMAX1
               IJK = FUNIJK(I,J,K)
               IF (PLOT_TYPE.EQ.5) THEN
                  WRITE (39+L,*) U_g(IJK) , V_g(IJK) , W_g(IJK)
               ELSE IF (PLOT_TYPE.EQ.6) THEN
                  WRITE (39+L,*) U_s(IJK,M_USE) , V_s(IJK,M_USE) ,&
                                 W_s(IJK,M_USE)
               ELSE
                  IF (VAR_INDEX.EQ.01) WRITE (39+L,*) EP_g(IJK)
                  IF (VAR_INDEX.EQ.02) WRITE (39+L,*) P_g(IJK)
                  IF (VAR_INDEX.EQ.03) WRITE (39+L,*) P_star(IJK)
                  IF (VAR_INDEX.EQ.04) WRITE (39+L,*) U_g(IJK)
                  IF (VAR_INDEX.EQ.05) WRITE (39+L,*) V_g(IJK)
                  IF (VAR_INDEX.EQ.06) WRITE (39+L,*) W_g(IJK)
                  IF (VAR_INDEX.EQ.07) WRITE (39+L,*) U_s(IJK,M_USE)
                  IF (VAR_INDEX.EQ.08) WRITE (39+L,*) V_s(IJK,M_USE)
                  IF (VAR_INDEX.EQ.09) WRITE (39+L,*) W_s(IJK,M_USE)
                  IF (VAR_INDEX.EQ.10) WRITE (39+L,*) ROP_s(IJK,M_USE)
                  IF (VAR_INDEX.EQ.11) WRITE (39+L,*) T_g(IJK)
                  IF (VAR_INDEX.EQ.12) WRITE (39+L,*) T_s(IJK,M_USE)
                  IF (VAR_INDEX.EQ.13) WRITE (39+L,*) T_s(IJK,M_USE)
              END IF
           END DO
         END DO
      END DO
      WRITE (39+L,*) SPEC_TIME
!
      RETURN
      END
