!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_EQ(A, B, VEL, DTxFA, IJK1)                         C
!  Purpose: Complete the filling of coefficient matrix A               C
!                                                                      C
!  Author: M. Syamlal                                 Date: 23-JUL-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_EQ(A, BB, VEL, DTXFA, IJK1)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE geometry
      USE fldvar
      USE indices
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Local index of the solids phases
      INTEGER          L, M
!
!                      IMJK etc.
      INTEGER          IJK1
!
!                      Average DTxF_gs and DTxF_ss
      DOUBLE PRECISION DTxFA(0:DIMENSION_M, 0:DIMENSION_M)
!
!                      Coefficients in the discretized momentum equations
!                        Av = B, where v(0) = V_g, v(1) = V_s1 etc.
      DOUBLE PRECISION A(0:DIMENSION_M, 0:DIMENSION_M), BB(0:DIMENSION_M)
!
!                      Solids velocity array passed
      DOUBLE PRECISION VEL (DIMENSION_3, DIMENSION_M)
!
!-----------------------------------------------
!
!
!     Compute off-diagonal coefficients
!
      DO M = 1, MMAX
         IF (A(M,M) /= ZERO) THEN
            A(0,M) = -DTXFA(0,M)
            A(M,0) = -DTXFA(0,M)
            DO L = M + 1, MMAX
               IF (A(L,L) /= ZERO) THEN
                  A(L,M) = -DTXFA(L,M)
                  A(M,L) = -DTXFA(L,M)
               ELSE
                  A(L,M) = ZERO
                  A(M,L) = ZERO
               ENDIF
            END DO
         ELSE
            A(0,M) = ZERO
            A(M,0) = ZERO
            L = M + 1
            IF (MMAX - M > 0) THEN
               A(M+1:MMAX,M) = ZERO
               A(M,M+1:MMAX) = ZERO
               L = MMAX + 1
            ENDIF
            IF (BB(M) == UNDEFINED) THEN
               A(0,0) = A(0,0) + DTXFA(0,M)
               BB(0) = BB(0) + DTXFA(0,M)*VEL(IJK1,M)
               DO L = 1, MMAX
                  IF (L/=M .AND. A(L,L)/=ZERO) THEN
                     A(L,L) = A(L,L) + DTXFA(L,M)
                     BB(L) = BB(L) + DTXFA(L,M)*VEL(IJK1,M)
                  ENDIF
               END DO
            ENDIF
         ENDIF
      END DO
      DO M = 0, MMAX
         DO L = 0, MMAX
            IF (L /= M) A(M,M) = A(M,M) - A(M,L)
         END DO
      END DO
      RETURN
      END SUBROUTINE GET_EQ
