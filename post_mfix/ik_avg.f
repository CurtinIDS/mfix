!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: IK_AVG (VAR,VOLUME,AVG,JUSE)                           C
!  Purpose: Calculate the value of a variable at a specified height -  C
!           averaged over X & Z (weighted by volume of the cell)       C
!                                                                      C
!  Author: P. Nicoletti                               Date: 20-JUL-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMIN1, IMAX1, KMIN1, KMAX1, FLAG              C
!  Variables modified: I,K,IJK                                         C
!                                                                      C
!  Local variables: TOTVOL                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE IK_AVG (VAR,VOLUME,AVG,JUSE)
!
      Use param
      Use param1
      Use fldvar
      Use physprop
      Use indices
      Use geometry
      Use compar
      Use functions
!
      IMPLICIT NONE
!
!     passed arguments
!
      DOUBLE PRECISION  VAR(*)
      REAL              VOLUME(*) , AVG
      INTEGER           JUSE
      INTEGER           I, K, IJK
!
!     local variables
!
      REAL              TOTVOL
!
      TOTVOL = 0.0
      AVG    = 0.0
      DO K = KMIN1,KMAX1
         DO I = IMIN1,IMAX1
            IJK = FUNIJK(I,JUSE,K)
            IF (FLAG(IJK).EQ.1) THEN
               AVG = AVG + VAR(IJK) * VOLUME(IJK)
               TOTVOL = TOTVOL + VOLUME(IJK)
            END IF
         END DO
      END DO
!
      AVG = AVG / TOTVOL
!
      RETURN
      END
