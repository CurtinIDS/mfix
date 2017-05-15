!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: F_INIT_DATA                                            C
!  Purpose:                                                            C
!                                                                      C
!  Author: P.Nicoletti                                Date: 05-JUN-95  C
!  Reviewer:                                          Date: dd-mmm-yy  C
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
      SUBROUTINE F_INIT_DATA
!
!
      Use param
      Use param1
      Use geometry
      Use physprop
      Use constant

      IMPLICIT NONE
      INCLUDE 'xforms.inc'
!
      INTEGER I
!
      I_MAX = IMAX2
      J_MAX = JMAX2
      K_MAX = KMAX2
      MM_MAX = MMAX
!
      DO I = 0,MMAX
         N_MAX(I) = NMAX(I)
      END DO
!
      RN_SPX = REAL(N_SPX)
!
      IF (C_E.EQ.UNDEFINED) THEN
         E_PASS = -1.0
         HAVE_E = .FALSE.
      ELSE
         E_PASS = C_E
         HAVE_E = .TRUE.
      END IF
!
      IF (RO_g0.EQ.UNDEFINED) THEN
         HAVE_RO_g0 = .FALSE.
      ELSE
         HAVE_RO_g0 = .TRUE.
      END IF
!
      IF (MW_avg.EQ.UNDEFINED) THEN
         HAVE_MW_avg = .FALSE.
      ELSE
         HAVE_MW_avg = .TRUE.
      END IF
!
      RETURN
      END
      subroutine do_pan
      write (*,*) ' hello'
      return
      end
