!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_READ_RESTART                                     C
!  Purpose: Reading QMOMK data for restart                             C
!                                                                      C
!                                                                      C
!  Author: Alberto Passalacqua                        Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
SUBROUTINE QMOMK_READ_RESTART

  USE param
  USE param1
  USE constant
  USE fldvar
  USE cont
  USE geometry
  USE indices
  USE run
  USE compar
  USE physprop
  USE qmom_kinetic_equation
  USE qmomk_quadrature
  USE functions

  IMPLICIT NONE

  INTEGER :: M, IJK

  OPEN (UNIT=901, FILE=TRIM(RUN_NAME)//'_QMOMK.RES', FORM='Unformatted', STATUS='unknown')

  REWIND (901)

  PRINT *,'QMOMK: Reading restart info...'

  READ (901) QMOMK_N1
  READ (901) QMOMK_U1
  READ (901) QMOMK_V1
  READ (901) QMOMK_W1

  PRINT *,'QMOMK: Updating moments after restart...'

  DO M = 1, MMAX
    DO IJK = ijkstart3, ijkend3
     IF (FLUID_AT(IJK))  THEN
       CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
            QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
       CALL EIGHT_NODE_3D (QMOMK_M1(:,IJK,M), QMOMK_N1(:,IJK,M), &
            QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M))
       CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N1(:,IJK,M), &
            QMOMK_U1(:,IJK,M), QMOMK_V1(:,IJK,M), QMOMK_W1(:,IJK,M), QMOMK_M1(:,IJK,M))
     END IF
    END DO
  END DO

  QMOMK_N0 = QMOMK_N1
  QMOMK_U0 = QMOMK_U1
  QMOMK_V0 = QMOMK_V1
  QMOMK_W0 = QMOMK_W1
  QMOMK_M0 = QMOMK_M1

  PRINT *,'QMOMK: Restart successful!'

END SUBROUTINE QMOMK_READ_RESTART
