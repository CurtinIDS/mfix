!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_EP_g(IJK)                                         C
!  Purpose: Calculate EP_g from known solids volume fractions          C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: MMAX, IJK                                     C
!  Variables modified: M, EP_g                                         C
!                                                                      C
!  Local variables: SUM_EPS                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_EP_g(IJK)
!
!
!  Include param.inc file to specify parameter values
!
      Use param
      Use param1
!
!     Physical and Numerical Parameters Section
!
      Use physprop
!
!     Field Variables
!
      Use fldvar
      Use geometry
!
!     Indices
!
      Use indices
      Use compar
      Use functions

      IMPLICIT NONE
!
!  Local variables
!
!                      Sum of EP_s
      DOUBLE PRECISION SUM_EPS
      INTEGER          IJK, M

      SUM_EPS = ZERO
      DO 50 M = 1, MMAX
        SUM_EPS = SUM_EPS + EP_s(IJK, M)
50    CONTINUE
      EP_g(IJK) = ONE - SUM_EPS
!
      RETURN
      END
