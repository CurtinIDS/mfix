!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_WRITE_RESTART                                    C
!  Purpose: Writing QMOMK data for restart                             C
!                                                                      C
!                                                                      C
!  Author: Alberto Passalacqua                        Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


SUBROUTINE QMOMK_WRITE_RESTART

  USE param1
  USE qmom_kinetic_equation
  USE run

  IMPLICIT NONE

  OPEN (UNIT=901, FILE=TRIM(RUN_NAME)//'_QMOMK.RES', FORM='Unformatted', STATUS='unknown')

  REWIND (901)

  ! Only weights and abscissas are necessary for restarting a calculation
  ! Moments are NOT necessary because they are calculated directly from
  ! weights and abscissas - Alberto Passalacqua

  WRITE (901) QMOMK_N1
  WRITE (901) QMOMK_U1
  WRITE (901) QMOMK_V1
  WRITE (901) QMOMK_W1

END SUBROUTINE QMOMK_WRITE_RESTART

