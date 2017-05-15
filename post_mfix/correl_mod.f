!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: correlations.inc                                       C
!  Purpose: Include file containing variables used for computing       C
!           averages, variances, and correlations.                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 01-AUG-92  C
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
        MODULE correl

        Use param
        Use param1

!                      Array for accumulating the sum of EP_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  SUM_EP_g
!
!                      Average value of EP_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AVG_EP_g
!
!                      Array for accumulating the sum of EP_g*EP_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  SUM_EPxEP_g
!
!                      Standard deviation of EP_g.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  SDV_EP_g
!
!                      Array for accumulating the sum of V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  SUM_V_g
!
!                      Average value of V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  AVG_V_g
!
!                      Array for accumulating the sum of V_g*V_g
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  SUM_VxV_g
!
!                      Standard deviation of V_g.
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  SDV_V_g
!
!                      A flag to check whether the subroutine is called for
!                      the first time
      INTEGER          STARTED
!
!                      Number of time step data summed
      INTEGER          NSUM


      END  MODULE correl
