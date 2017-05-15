!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_INIT_BC                                           C
!  Author: Alberto Passalacqua                        Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

SUBROUTINE QMOMK_INIT_BC

!-----------------------------------------------
!     M o d u l e s
!-----------------------------------------------
  USE param
  USE param1
  USE constant
  USE physprop
  USE fldvar
  USE geometry
  USE compar
  USE indices
  USE bc
  USE qmom_kinetic_equation
  USE qmomk_quadrature
  USE qmomk_bc
  USE functions

  IMPLICIT NONE

  INTEGER :: L
  INTEGER :: I1, I2, J1, J2, K1, K2

  DO L = 1, DIMENSION_BC
      IF (BC_DEFINED(L)) THEN

        I1 = BC_I_W(L)
        I2 = BC_I_E(L)
        J1 = BC_J_S(L)
        J2 = BC_J_N(L)
        K1 = BC_K_B(L)
        K2 = BC_K_T(L)

        ! Wall BC's
        IF (BC_TYPE_ENUM(L)==FREE_SLIP_WALL .OR. BC_TYPE_ENUM(L)==NO_SLIP_WALL&
            .OR. BC_TYPE_ENUM(L)==PAR_SLIP_WALL) THEN
          IF (QMOMK_WALL_BC_TYPE == 'SPECULAR_REFLECTIVE') THEN
            CALL QMOMK_REFLECTIVE_WALL_BC(L, I1, I2, J1, J2, K1, K2, .TRUE.)
          ENDIF
        ! Outlet BC's
        ELSEIF (BC_TYPE_ENUM(L) == MASS_OUTFLOW .OR. BC_TYPE_ENUM(L)==P_OUTFLOW&
            .OR. BC_TYPE_ENUM(L)==OUTFLOW) THEN
          CALL QMOMK_OUTLET_BC(L, I1, I2, J1, J2, K1, K2, .TRUE.)
        ! Inlet BC's
        ELSEIF (BC_TYPE_ENUM(L) == MASS_INFLOW) THEN
          CALL QMOMK_INLET_BC(L, .TRUE.)
        ! Cyclic BC's
        ELSEIF (CYCLIC) THEN
          CALL QMOMK_CYCLIC_BC(.TRUE.)
        END IF
      END IF
  END DO

END SUBROUTINE QMOMK_INIT_BC
