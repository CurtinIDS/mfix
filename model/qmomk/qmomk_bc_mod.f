!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_BC                                               C
!  Author: Alberto Passalacqua (A.P.)                 Date:            C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

MODULE qmomk_bc

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
  USE functions

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: QMOMK_REFLECTIVE_WALL_BC
  PUBLIC :: QMOMK_OUTLET_BC
  PUBLIC :: QMOMK_INLET_BC
  PUBLIC :: QMOMK_CYCLIC_BC

CONTAINS

 ! A.P. Purely reflective boundary condition with restitution coefficien e_w
 SUBROUTINE QMOMK_REFLECTIVE_WALL_BC(L, I1, I2, J1, J2, K1, K2, INIT)
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: L, I1, I2, J1, J2, K1, K2
   LOGICAL, INTENT(IN) :: INIT

   INTEGER :: I, J, K, M, IJK, LFLUID

   IF (INIT) THEN
     DO K = K1, K2
       DO J = J1, J2
         DO I = I1, I2
  !//SP Check if current i,j,k resides on this PE
          IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
          IJK = FUNIJK(I,J,K)
  !
  ! Fluid cell at West
  !
          IF (FLUID_AT(IM_OF(IJK))) THEN
            LFLUID = IM_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)/e_w
              QMOMK_U0(:,IJK,M) = -e_w*QMOMK_U0(:,LFLUID,M)
              QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
              QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
            END DO
          ENDIF

  !
  ! Fluid cell at East
  !
          IF (FLUID_AT(IP_OF(IJK))) THEN
            LFLUID = IP_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)/e_w
              QMOMK_U0(:,IJK,M) = -e_w*QMOMK_U0(:,LFLUID,M)
              QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
              QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
            END DO
          END IF

  !
  ! Fluid cell at South
  !
          IF (FLUID_AT(JM_OF(IJK))) THEN
            LFLUID = JM_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)/e_w
              QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
              QMOMK_V0(:,IJK,M) = -e_w*QMOMK_V0(:,LFLUID,M)
              QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
            END DO
          ENDIF
  !
  ! Fluid cell at North
  !
          IF (FLUID_AT(JP_OF(IJK))) THEN
            LFLUID = JP_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)/e_w
              QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
              QMOMK_V0(:,IJK,M) = -e_w*QMOMK_V0(:,LFLUID,M)
              QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
            END DO
          ENDIF
  !
  ! Fluid cell at Bottom
  !
          IF (FLUID_AT(KM_OF(IJK))) THEN
            LFLUID = KM_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)/e_w
              QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
              QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
              QMOMK_W0(:,IJK,M) = -e_w*QMOMK_W0(:,LFLUID,M)
            END DO
          ENDIF
  !
  ! Fluid cell at Top
  !
          IF (FLUID_AT(KP_OF(IJK))) THEN
            LFLUID = KP_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)/e_w
              QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
              QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
              QMOMK_W0(:,IJK,M) = -e_w*QMOMK_W0(:,LFLUID,M)
            END DO
          ENDIF
         END DO
       END DO
     END DO

   ! Running...
   ELSE
     DO K = K1, K2
       DO J = J1, J2
         DO I = I1, I2
  !//SP Check if current i,j,k resides on this PE
          IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
          IJK = FUNIJK(I,J,K)
  !
  ! Fluid cell at West
  !
          IF (FLUID_AT(IM_OF(IJK))) THEN
            LFLUID = IM_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)/e_w
              QMOMK_U1(:,IJK,M) = -e_w*QMOMK_U1(:,LFLUID,M)
              QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
              QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
            END DO
          ENDIF

  !
  ! Fluid cell at East
  !
          IF (FLUID_AT(IP_OF(IJK))) THEN
            LFLUID = IP_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)/e_w
              QMOMK_U1(:,IJK,M) = -e_w*QMOMK_U1(:,LFLUID,M)
              QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
              QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
            END DO
          END IF

  !
  ! Fluid cell at South
  !
          IF (FLUID_AT(JM_OF(IJK))) THEN
            LFLUID = JM_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)/e_w
              QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
              QMOMK_V1(:,IJK,M) = -e_w*QMOMK_V1(:,LFLUID,M)
              QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
            END DO
          ENDIF
  !
  ! Fluid cell at North
  !
          IF (FLUID_AT(JP_OF(IJK))) THEN
            LFLUID = JP_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)/e_w
              QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
              QMOMK_V1(:,IJK,M) = -e_w*QMOMK_V1(:,LFLUID,M)
              QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
            END DO
          ENDIF
  !
  ! Fluid cell at Bottom
  !
          IF (FLUID_AT(KM_OF(IJK))) THEN
            LFLUID = KM_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)/e_w
              QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
              QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
              QMOMK_W1(:,IJK,M) = -e_w*QMOMK_W1(:,LFLUID,M)
            END DO
          ENDIF
  !
  ! Fluid cell at Top
  !
          IF (FLUID_AT(KP_OF(IJK))) THEN
            LFLUID = KP_OF(IJK)
            DO M = 1, MMAX
              QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)/e_w
              QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
              QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
              QMOMK_W1(:,IJK,M) = -e_w*QMOMK_W1(:,LFLUID,M)
            END DO
          ENDIF
         END DO
       END DO
     END DO
    END IF
   RETURN

 END SUBROUTINE QMOMK_REFLECTIVE_WALL_BC

 ! A.P. Zero-gradient outlet
 SUBROUTINE QMOMK_OUTLET_BC(L, I1, I2, J1, J2, K1, K2, INIT)
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: I1, I2, J1, J2, K1, K2, L
   LOGICAL, INTENT(IN) :: INIT

   INTEGER :: I, J, K, M, IJK, LFLUID

   IF (INIT) THEN
     DO K = K1, K2
       DO J = J1, J2
         DO I = I1, I2
!//SP Check if current i,j,k resides on this PE
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IJK = FUNIJK(I,J,K)
!
! Fluid cell at West
!
            IF (FLUID_AT(IM_OF(IJK))) THEN
              LFLUID = IM_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)
                QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
                QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
                QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
              END DO
            ENDIF

!
! Fluid cell at East
!
            IF (FLUID_AT(IP_OF(IJK))) THEN
              LFLUID = IP_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)
                QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
                QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
                QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
              END DO
            END IF

!
! Fluid cell at South
!
            IF (FLUID_AT(JM_OF(IJK))) THEN
              LFLUID = JM_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)
                QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
                QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
                QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
              END DO
            ENDIF
!
! Fluid cell at North
!
            IF (FLUID_AT(JP_OF(IJK))) THEN
              LFLUID = JP_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)
                QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
                QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
                QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
              END DO
            ENDIF
!
! Fluid cell at Bottom
!
            IF (FLUID_AT(KM_OF(IJK))) THEN
              LFLUID = KM_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)
                QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
                QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
                QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
              END DO
            ENDIF
!
! Fluid cell at Top
!
            IF (FLUID_AT(KP_OF(IJK))) THEN
              LFLUID = KP_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N0(:,IJK,M) = QMOMK_N0(:,LFLUID,M)
                QMOMK_U0(:,IJK,M) = QMOMK_U0(:,LFLUID,M)
                QMOMK_V0(:,IJK,M) = QMOMK_V0(:,LFLUID,M)
                QMOMK_W0(:,IJK,M) = QMOMK_W0(:,LFLUID,M)
              END DO
            ENDIF
         END DO
       END DO
     END DO

   ELSE
     DO K = K1, K2
       DO J = J1, J2
         DO I = I1, I2
!//SP Check if current i,j,k resides on this PE
            IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
            IJK = FUNIJK(I,J,K)
!
! Fluid cell at West
!
            IF (FLUID_AT(IM_OF(IJK))) THEN
              LFLUID = IM_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)
                QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
                QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
                QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
              END DO
            ENDIF

!
! Fluid cell at East
!
            IF (FLUID_AT(IP_OF(IJK))) THEN
              LFLUID = IP_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)
                QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
                QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
                QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
              END DO
            END IF

!
! Fluid cell at South
!
            IF (FLUID_AT(JM_OF(IJK))) THEN
              LFLUID = JM_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)
                QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
                QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
                QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
              END DO
            ENDIF
!
! Fluid cell at North
!
            IF (FLUID_AT(JP_OF(IJK))) THEN
              LFLUID = JP_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)
                QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
                QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
                QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
              END DO
            ENDIF
!
! Fluid cell at Bottom
!
            IF (FLUID_AT(KM_OF(IJK))) THEN
              LFLUID = KM_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)
                QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
                QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
                QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
              END DO
            ENDIF
!
! Fluid cell at Top
!
            IF (FLUID_AT(KP_OF(IJK))) THEN
              LFLUID = KP_OF(IJK)
              DO M = 1, MMAX
                QMOMK_N1(:,IJK,M) = QMOMK_N1(:,LFLUID,M)
                QMOMK_U1(:,IJK,M) = QMOMK_U1(:,LFLUID,M)
                QMOMK_V1(:,IJK,M) = QMOMK_V1(:,LFLUID,M)
                QMOMK_W1(:,IJK,M) = QMOMK_W1(:,LFLUID,M)
              END DO
            ENDIF
         END DO
       END DO
     END DO
   END IF

  RETURN

 END SUBROUTINE QMOMK_OUTLET_BC


 ! A.P. Velocity inlet
 SUBROUTINE QMOMK_INLET_BC(L, INIT)
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: L
   LOGICAL, INTENT(IN) :: INIT

   INTEGER :: I, J, K, M, IJK, IJK2
   DOUBLE PRECISION :: InitVal

   IF (INIT) THEN
    DO K = BC_K_B(L), BC_K_T(L)
      DO J = BC_J_S(L), BC_J_N(L)
        DO I = BC_I_W(L), BC_I_E(L)
          IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
          IJK = FUNIJK(I,J,K)
          SELECT CASE (TRIM(BC_PLANE(L)))
          CASE ('W')
            IJK2 = IM_OF(IJK)
            QMOMK_N0 (:, IJK2, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK2,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U0(1, IJK2, M) = -InitVal + BC_U_S(L, M)
            QMOMK_U0(2, IJK2, M) = +InitVal + BC_U_S(L, M)
            QMOMK_U0(3, IJK2, M) = -InitVal + BC_U_S(L, M)
            QMOMK_U0(4, IJK2, M) = +InitVal + BC_U_S(L, M)
            QMOMK_U0(5, IJK2, M) = -InitVal + BC_U_S(L, M)
            QMOMK_U0(6, IJK2, M) = +InitVal + BC_U_S(L, M)
            QMOMK_U0(7, IJK2, M) = -InitVal + BC_U_S(L, M)
            QMOMK_U0(8, IJK2, M) = +InitVal + BC_U_S(L, M)

            QMOMK_V0(1, IJK2, M) = -InitVal
            QMOMK_V0(2, IJK2, M) = -InitVal
            QMOMK_V0(3, IJK2, M) = +InitVal
            QMOMK_V0(4, IJK2, M) = +InitVal
            QMOMK_V0(5, IJK2, M) = -InitVal
            QMOMK_V0(6, IJK2, M) = -InitVal
            QMOMK_V0(7, IJK2, M) = +InitVal
            QMOMK_V0(8, IJK2, M) = +InitVal

            QMOMK_W0(1, IJK2, M) = -InitVal
            QMOMK_W0(2, IJK2, M) = -InitVal
            QMOMK_W0(3, IJK2, M) = -InitVal
            QMOMK_W0(4, IJK2, M) = -InitVal
            QMOMK_W0(5, IJK2, M) = +InitVal
            QMOMK_W0(6, IJK2, M) = +InitVal
            QMOMK_W0(7, IJK2, M) = +InitVal
            QMOMK_W0(8, IJK2, M) = +InitVal

          CASE ('E')
            QMOMK_N0 (:, IJK, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U0(1, IJK, M) = -InitVal + BC_U_S(IJK, M)
            QMOMK_U0(2, IJK, M) = +InitVal + BC_U_S(IJK, M)
            QMOMK_U0(3, IJK, M) = -InitVal + BC_U_S(IJK, M)
            QMOMK_U0(4, IJK, M) = +InitVal + BC_U_S(IJK, M)
            QMOMK_U0(5, IJK, M) = -InitVal + BC_U_S(IJK, M)
            QMOMK_U0(6, IJK, M) = +InitVal + BC_U_S(IJK, M)
            QMOMK_U0(7, IJK, M) = -InitVal + BC_U_S(IJK, M)
            QMOMK_U0(8, IJK, M) = +InitVal + BC_U_S(IJK, M)

            QMOMK_V0(1, IJK, M) = -InitVal
            QMOMK_V0(2, IJK, M) = -InitVal
            QMOMK_V0(3, IJK, M) = +InitVal
            QMOMK_V0(4, IJK, M) = +InitVal
            QMOMK_V0(5, IJK, M) = -InitVal
            QMOMK_V0(6, IJK, M) = -InitVal
            QMOMK_V0(7, IJK, M) = +InitVal
            QMOMK_V0(8, IJK, M) = +InitVal

            QMOMK_W0(1, IJK, M) = -InitVal
            QMOMK_W0(2, IJK, M) = -InitVal
            QMOMK_W0(3, IJK, M) = -InitVal
            QMOMK_W0(4, IJK, M) = -InitVal
            QMOMK_W0(5, IJK, M) = +InitVal
            QMOMK_W0(6, IJK, M) = +InitVal
            QMOMK_W0(7, IJK, M) = +InitVal
            QMOMK_W0(8, IJK, M) = +InitVal
          CASE ('S')
            IJK2 = JM_OF(IJK)
            QMOMK_N0 (:, IJK2, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U0(1, IJK2, M) = -InitVal
            QMOMK_U0(2, IJK2, M) = +InitVal
            QMOMK_U0(3, IJK2, M) = -InitVal
            QMOMK_U0(4, IJK2, M) = +InitVal
            QMOMK_U0(5, IJK2, M) = -InitVal
            QMOMK_U0(6, IJK2, M) = +InitVal
            QMOMK_U0(7, IJK2, M) = -InitVal
            QMOMK_U0(8, IJK2, M) = +InitVal

            QMOMK_V0(1, IJK2, M) = -InitVal + BC_V_S(L, M)
            QMOMK_V0(2, IJK2, M) = -InitVal + BC_V_S(L, M)
            QMOMK_V0(3, IJK2, M) = +InitVal + BC_V_S(L, M)
            QMOMK_V0(4, IJK2, M) = +InitVal + BC_V_S(L, M)
            QMOMK_V0(5, IJK2, M) = -InitVal + BC_V_S(L, M)
            QMOMK_V0(6, IJK2, M) = -InitVal + BC_V_S(L, M)
            QMOMK_V0(7, IJK2, M) = +InitVal + BC_V_S(L, M)
            QMOMK_V0(8, IJK2, M) = +InitVal + BC_V_S(L, M)

            QMOMK_W0(1, IJK2, M) = -InitVal
            QMOMK_W0(2, IJK2, M) = -InitVal
            QMOMK_W0(3, IJK2, M) = -InitVal
            QMOMK_W0(4, IJK2, M) = -InitVal
            QMOMK_W0(5, IJK2, M) = +InitVal
            QMOMK_W0(6, IJK2, M) = +InitVal
            QMOMK_W0(7, IJK2, M) = +InitVal
            QMOMK_W0(8, IJK2, M) = +InitVal

          CASE ('N')
            QMOMK_N0 (:, IJK, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U0(1, IJK, M) = -InitVal
            QMOMK_U0(2, IJK, M) = +InitVal
            QMOMK_U0(3, IJK, M) = -InitVal
            QMOMK_U0(4, IJK, M) = +InitVal
            QMOMK_U0(5, IJK, M) = -InitVal
            QMOMK_U0(6, IJK, M) = +InitVal
            QMOMK_U0(7, IJK, M) = -InitVal
            QMOMK_U0(8, IJK, M) = +InitVal

            QMOMK_V0(1, IJK, M) = -InitVal + BC_V_S(IJK, M)
            QMOMK_V0(2, IJK, M) = -InitVal + BC_V_S(IJK, M)
            QMOMK_V0(3, IJK, M) = +InitVal + BC_V_S(IJK, M)
            QMOMK_V0(4, IJK, M) = +InitVal + BC_V_S(IJK, M)
            QMOMK_V0(5, IJK, M) = -InitVal + BC_V_S(IJK, M)
            QMOMK_V0(6, IJK, M) = -InitVal + BC_V_S(IJK, M)
            QMOMK_V0(7, IJK, M) = +InitVal + BC_V_S(IJK, M)
            QMOMK_V0(8, IJK, M) = +InitVal + BC_V_S(IJK, M)

            QMOMK_W0(1, IJK, M) = -InitVal
            QMOMK_W0(2, IJK, M) = -InitVal
            QMOMK_W0(3, IJK, M) = -InitVal
            QMOMK_W0(4, IJK, M) = -InitVal
            QMOMK_W0(5, IJK, M) = +InitVal
            QMOMK_W0(6, IJK, M) = +InitVal
            QMOMK_W0(7, IJK, M) = +InitVal
            QMOMK_W0(8, IJK, M) = +InitVal

          CASE ('B')
            IJK2 = KM_OF(IJK)
            QMOMK_N0 (:, IJK2, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK2,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U0(1, IJK2, M) = -InitVal
            QMOMK_U0(2, IJK2, M) = +InitVal
            QMOMK_U0(3, IJK2, M) = -InitVal
            QMOMK_U0(4, IJK2, M) = +InitVal
            QMOMK_U0(5, IJK2, M) = -InitVal
            QMOMK_U0(6, IJK2, M) = +InitVal
            QMOMK_U0(7, IJK2, M) = -InitVal
            QMOMK_U0(8, IJK2, M) = +InitVal

            QMOMK_V0(1, IJK2, M) = -InitVal
            QMOMK_V0(2, IJK2, M) = -InitVal
            QMOMK_V0(3, IJK2, M) = +InitVal
            QMOMK_V0(4, IJK2, M) = +InitVal
            QMOMK_V0(5, IJK2, M) = -InitVal
            QMOMK_V0(6, IJK2, M) = -InitVal
            QMOMK_V0(7, IJK2, M) = +InitVal
            QMOMK_V0(8, IJK2, M) = +InitVal

            QMOMK_W0(1, IJK2, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W0(2, IJK2, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W0(3, IJK2, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W0(4, IJK2, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W0(5, IJK2, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W0(6, IJK2, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W0(7, IJK2, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W0(8, IJK2, M) = +InitVal + BC_W_S(L, M)

          CASE ('T')
            QMOMK_N0 (:, IJK, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U0(1, IJK, M) = -InitVal
            QMOMK_U0(2, IJK, M) = +InitVal
            QMOMK_U0(3, IJK, M) = -InitVal
            QMOMK_U0(4, IJK, M) = +InitVal
            QMOMK_U0(5, IJK, M) = -InitVal
            QMOMK_U0(6, IJK, M) = +InitVal
            QMOMK_U0(7, IJK, M) = -InitVal
            QMOMK_U0(8, IJK, M) = +InitVal

            QMOMK_V0(1, IJK, M) = -InitVal
            QMOMK_V0(2, IJK, M) = -InitVal
            QMOMK_V0(3, IJK, M) = +InitVal
            QMOMK_V0(4, IJK, M) = +InitVal
            QMOMK_V0(5, IJK, M) = -InitVal
            QMOMK_V0(6, IJK, M) = -InitVal
            QMOMK_V0(7, IJK, M) = +InitVal
            QMOMK_V0(8, IJK, M) = +InitVal

            QMOMK_W0(1, IJK, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W0(2, IJK, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W0(3, IJK, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W0(4, IJK, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W0(5, IJK, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W0(6, IJK, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W0(7, IJK, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W0(8, IJK, M) = +InitVal + BC_W_S(L, M)
          END SELECT
        END DO
      END DO
    END DO

   ELSE
    DO K = BC_K_B(L), BC_K_T(L)
      DO J = BC_J_S(L), BC_J_N(L)
        DO I = BC_I_W(L), BC_I_E(L)
          IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
          IJK = FUNIJK(I,J,K)
          SELECT CASE (TRIM(BC_PLANE(L)))

          CASE ('W')
            IJK2 = IM_OF(IJK)
            QMOMK_N1 (:, IJK2, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK2,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U1(1, IJK2, M) = -InitVal + BC_U_S(L, M)
            QMOMK_U1(2, IJK2, M) = +InitVal + BC_U_S(L, M)
            QMOMK_U1(3, IJK2, M) = -InitVal + BC_U_S(L, M)
            QMOMK_U1(4, IJK2, M) = +InitVal + BC_U_S(L, M)
            QMOMK_U1(5, IJK2, M) = -InitVal + BC_U_S(L, M)
            QMOMK_U1(6, IJK2, M) = +InitVal + BC_U_S(L, M)
            QMOMK_U1(7, IJK2, M) = -InitVal + BC_U_S(L, M)
            QMOMK_U1(8, IJK2, M) = +InitVal + BC_U_S(L, M)

            QMOMK_V1(1, IJK2, M) = -InitVal
            QMOMK_V1(2, IJK2, M) = -InitVal
            QMOMK_V1(3, IJK2, M) = +InitVal
            QMOMK_V1(4, IJK2, M) = +InitVal
            QMOMK_V1(5, IJK2, M) = -InitVal
            QMOMK_V1(6, IJK2, M) = -InitVal
            QMOMK_V1(7, IJK2, M) = +InitVal
            QMOMK_V1(8, IJK2, M) = +InitVal

            QMOMK_W1(1, IJK2, M) = -InitVal
            QMOMK_W1(2, IJK2, M) = -InitVal
            QMOMK_W1(3, IJK2, M) = -InitVal
            QMOMK_W1(4, IJK2, M) = -InitVal
            QMOMK_W1(5, IJK2, M) = +InitVal
            QMOMK_W1(6, IJK2, M) = +InitVal
            QMOMK_W1(7, IJK2, M) = +InitVal
            QMOMK_W1(8, IJK2, M) = +InitVal

          CASE ('E')
            QMOMK_N1 (:, IJK, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U1(1, IJK, M) = -InitVal + BC_U_S(IJK, M)
            QMOMK_U1(2, IJK, M) = +InitVal + BC_U_S(IJK, M)
            QMOMK_U1(3, IJK, M) = -InitVal + BC_U_S(IJK, M)
            QMOMK_U1(4, IJK, M) = +InitVal + BC_U_S(IJK, M)
            QMOMK_U1(5, IJK, M) = -InitVal + BC_U_S(IJK, M)
            QMOMK_U1(6, IJK, M) = +InitVal + BC_U_S(IJK, M)
            QMOMK_U1(7, IJK, M) = -InitVal + BC_U_S(IJK, M)
            QMOMK_U1(8, IJK, M) = +InitVal + BC_U_S(IJK, M)

            QMOMK_V1(1, IJK, M) = -InitVal
            QMOMK_V1(2, IJK, M) = -InitVal
            QMOMK_V1(3, IJK, M) = +InitVal
            QMOMK_V1(4, IJK, M) = +InitVal
            QMOMK_V1(5, IJK, M) = -InitVal
            QMOMK_V1(6, IJK, M) = -InitVal
            QMOMK_V1(7, IJK, M) = +InitVal
            QMOMK_V1(8, IJK, M) = +InitVal

            QMOMK_W1(1, IJK, M) = -InitVal
            QMOMK_W1(2, IJK, M) = -InitVal
            QMOMK_W1(3, IJK, M) = -InitVal
            QMOMK_W1(4, IJK, M) = -InitVal
            QMOMK_W1(5, IJK, M) = +InitVal
            QMOMK_W1(6, IJK, M) = +InitVal
            QMOMK_W1(7, IJK, M) = +InitVal
            QMOMK_W1(8, IJK, M) = +InitVal

          CASE ('S')
            IJK2 = JM_OF(IJK)
            QMOMK_N1 (:, IJK2, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK2,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U1(1, IJK2, M) = -InitVal
            QMOMK_U1(2, IJK2, M) = +InitVal
            QMOMK_U1(3, IJK2, M) = -InitVal
            QMOMK_U1(4, IJK2, M) = +InitVal
            QMOMK_U1(5, IJK2, M) = -InitVal
            QMOMK_U1(6, IJK2, M) = +InitVal
            QMOMK_U1(7, IJK2, M) = -InitVal
            QMOMK_U1(8, IJK2, M) = +InitVal

            QMOMK_V1(1, IJK2, M) = -InitVal + BC_V_S(L, M)
            QMOMK_V1(2, IJK2, M) = -InitVal + BC_V_S(L, M)
            QMOMK_V1(3, IJK2, M) = +InitVal + BC_V_S(L, M)
            QMOMK_V1(4, IJK2, M) = +InitVal + BC_V_S(L, M)
            QMOMK_V1(5, IJK2, M) = -InitVal + BC_V_S(L, M)
            QMOMK_V1(6, IJK2, M) = -InitVal + BC_V_S(L, M)
            QMOMK_V1(7, IJK2, M) = +InitVal + BC_V_S(L, M)
            QMOMK_V1(8, IJK2, M) = +InitVal + BC_V_S(L, M)

            QMOMK_W1(1, IJK2, M) = -InitVal
            QMOMK_W1(2, IJK2, M) = -InitVal
            QMOMK_W1(3, IJK2, M) = -InitVal
            QMOMK_W1(4, IJK2, M) = -InitVal
            QMOMK_W1(5, IJK2, M) = +InitVal
            QMOMK_W1(6, IJK2, M) = +InitVal
            QMOMK_W1(7, IJK2, M) = +InitVal
            QMOMK_W1(8, IJK2, M) = +InitVal

          CASE ('N')
            QMOMK_N1 (:, IJK, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U1(1, IJK, M) = -InitVal
            QMOMK_U1(2, IJK, M) = +InitVal
            QMOMK_U1(3, IJK, M) = -InitVal
            QMOMK_U1(4, IJK, M) = +InitVal
            QMOMK_U1(5, IJK, M) = -InitVal
            QMOMK_U1(6, IJK, M) = +InitVal
            QMOMK_U1(7, IJK, M) = -InitVal
            QMOMK_U1(8, IJK, M) = +InitVal

            QMOMK_V1(1, IJK, M) = -InitVal + BC_V_S(IJK, M)
            QMOMK_V1(2, IJK, M) = -InitVal + BC_V_S(IJK, M)
            QMOMK_V1(3, IJK, M) = +InitVal + BC_V_S(IJK, M)
            QMOMK_V1(4, IJK, M) = +InitVal + BC_V_S(IJK, M)
            QMOMK_V1(5, IJK, M) = -InitVal + BC_V_S(IJK, M)
            QMOMK_V1(6, IJK, M) = -InitVal + BC_V_S(IJK, M)
            QMOMK_V1(7, IJK, M) = +InitVal + BC_V_S(IJK, M)
            QMOMK_V1(8, IJK, M) = +InitVal + BC_V_S(IJK, M)

            QMOMK_W1(1, IJK, M) = -InitVal
            QMOMK_W1(2, IJK, M) = -InitVal
            QMOMK_W1(3, IJK, M) = -InitVal
            QMOMK_W1(4, IJK, M) = -InitVal
            QMOMK_W1(5, IJK, M) = +InitVal
            QMOMK_W1(6, IJK, M) = +InitVal
            QMOMK_W1(7, IJK, M) = +InitVal
            QMOMK_W1(8, IJK, M) = +InitVal

          CASE ('B')
            IJK2 = KM_OF(IJK)
            QMOMK_N1 (:, IJK2, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK2,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U1(1, IJK2, M) = -InitVal
            QMOMK_U1(2, IJK2, M) = +InitVal
            QMOMK_U1(3, IJK2, M) = -InitVal
            QMOMK_U1(4, IJK2, M) = +InitVal
            QMOMK_U1(5, IJK2, M) = -InitVal
            QMOMK_U1(6, IJK2, M) = +InitVal
            QMOMK_U1(7, IJK2, M) = -InitVal
            QMOMK_U1(8, IJK2, M) = +InitVal

            QMOMK_V1(1, IJK2, M) = -InitVal
            QMOMK_V1(2, IJK2, M) = -InitVal
            QMOMK_V1(3, IJK2, M) = +InitVal
            QMOMK_V1(4, IJK2, M) = +InitVal
            QMOMK_V1(5, IJK2, M) = -InitVal
            QMOMK_V1(6, IJK2, M) = -InitVal
            QMOMK_V1(7, IJK2, M) = +InitVal
            QMOMK_V1(8, IJK2, M) = +InitVal

            QMOMK_W1(1, IJK2, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W1(2, IJK2, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W1(3, IJK2, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W1(4, IJK2, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W1(5, IJK2, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W1(6, IJK2, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W1(7, IJK2, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W1(8, IJK2, M) = +InitVal + BC_W_S(L, M)

          CASE ('T')
            QMOMK_N1 (:, IJK, M) = BC_ROP_s (L, M)/(QMOMK_NN * RO_s(IJK,M))

            InitVal = MAX(SQRT(BC_THETA_M(L,M)), MINIMUM_THETA)

            QMOMK_U1(1, IJK, M) = -InitVal
            QMOMK_U1(2, IJK, M) = +InitVal
            QMOMK_U1(3, IJK, M) = -InitVal
            QMOMK_U1(4, IJK, M) = +InitVal
            QMOMK_U1(5, IJK, M) = -InitVal
            QMOMK_U1(6, IJK, M) = +InitVal
            QMOMK_U1(7, IJK, M) = -InitVal
            QMOMK_U1(8, IJK, M) = +InitVal

            QMOMK_V1(1, IJK, M) = -InitVal
            QMOMK_V1(2, IJK, M) = -InitVal
            QMOMK_V1(3, IJK, M) = +InitVal
            QMOMK_V1(4, IJK, M) = +InitVal
            QMOMK_V1(5, IJK, M) = -InitVal
            QMOMK_V1(6, IJK, M) = -InitVal
            QMOMK_V1(7, IJK, M) = +InitVal
            QMOMK_V1(8, IJK, M) = +InitVal

            QMOMK_W1(1, IJK, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W1(2, IJK, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W1(3, IJK, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W1(4, IJK, M) = -InitVal + BC_W_S(L, M)
            QMOMK_W1(5, IJK, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W1(6, IJK, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W1(7, IJK, M) = +InitVal + BC_W_S(L, M)
            QMOMK_W1(8, IJK, M) = +InitVal + BC_W_S(L, M)
          END SELECT
        END DO
      END DO
    END DO
   END IF
   RETURN
 END SUBROUTINE QMOMK_INLET_BC

 ! A.P. Cyclic boundary conditions
 SUBROUTINE QMOMK_CYCLIC_BC(INIT)
  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: INIT
  INTEGER :: IJK, IJK_CYCLIC, I, J, K, IJKN, IJKS, IJKE, IJKW, IJKT, IJKB, M
  DOUBLE PRECISION, DIMENSION(QMOMK_NN) :: QMOMK_N_TMP, QMOMK_U_TMP
  DOUBLE PRECISION, DIMENSION(QMOMK_NN) :: QMOMK_V_TMP, QMOMK_W_TMP

  IF (INIT) THEN
    DO M = 1, MMAX
     DO IJK = ijkstart3, ijkend3
      IF (FLUID_AT(IJK)) THEN
        IJKN = NORTH_OF(IJK)
        IJKS = SOUTH_OF(IJK)
        IJKE = EAST_OF(IJK)
        IJKW = WEST_OF(IJK)
        IJKT = TOP_OF(IJK)
        IJKB = BOTTOM_OF(IJK)

        !  x direction cyclic
        IF (CYCLIC_X .OR. CYCLIC_X_PD) THEN
          IF (CYCLIC_AT_E(IJK)) THEN
              I = I_OF(IJKE)
              J = J_OF(IJKE)
              K = K_OF(IJKE)

              IJK_CYCLIC = FUNIJK(IP1(I), J, K)

              QMOMK_N_TMP(:) = QMOMK_N0(:,IJKE,M)
              QMOMK_U_TMP(:) = QMOMK_U0(:,IJKE,M)
              QMOMK_V_TMP(:) = QMOMK_V0(:,IJKE,M)
              QMOMK_W_TMP(:) = QMOMK_W0(:,IJKE,M)

              QMOMK_N0(:,IJKE,M) = QMOMK_N0(:,IJK_CYCLIC,M)
              QMOMK_U0(:,IJKE,M) = QMOMK_U0(:,IJK_CYCLIC,M)
              QMOMK_V0(:,IJKE,M) = QMOMK_V0(:,IJK_CYCLIC,M)
              QMOMK_W0(:,IJKE,M) = QMOMK_W0(:,IJK_CYCLIC,M)

              QMOMK_N0(:,IJK_CYCLIC,M) = QMOMK_N_TMP(:)
              QMOMK_U0(:,IJK_CYCLIC,M) = QMOMK_U_TMP(:)
              QMOMK_V0(:,IJK_CYCLIC,M) = QMOMK_V_TMP(:)
              QMOMK_W0(:,IJK_CYCLIC,M) = QMOMK_W_TMP(:)
          END IF
        END IF
        ! y direction cyclic
        IF (CYCLIC_Y .OR. CYCLIC_Y_PD) THEN
          IF (CYCLIC_AT_N(IJK)) THEN
              I = I_OF(IJKN)
              J = J_OF(IJKN)
              K = K_OF(IJKN)

              IJK_CYCLIC = FUNIJK(I, JP1(J), K)

              QMOMK_N_TMP(:) = QMOMK_N0(:,IJKN,M)
              QMOMK_U_TMP(:) = QMOMK_U0(:,IJKN,M)
              QMOMK_V_TMP(:) = QMOMK_V0(:,IJKN,M)
              QMOMK_W_TMP(:) = QMOMK_W0(:,IJKN,M)

              QMOMK_N0(:,IJKN,M) = QMOMK_N0(:,IJK_CYCLIC,M)
              QMOMK_U0(:,IJKN,M) = QMOMK_U0(:,IJK_CYCLIC,M)
              QMOMK_V0(:,IJKN,M) = QMOMK_V0(:,IJK_CYCLIC,M)
              QMOMK_W0(:,IJKN,M) = QMOMK_W0(:,IJK_CYCLIC,M)

              QMOMK_N0(:,IJK_CYCLIC,M) = QMOMK_N_TMP(:)
              QMOMK_U0(:,IJK_CYCLIC,M) = QMOMK_U_TMP(:)
              QMOMK_V0(:,IJK_CYCLIC,M) = QMOMK_V_TMP(:)
              QMOMK_W0(:,IJK_CYCLIC,M) = QMOMK_W_TMP(:)
          END IF
        END IF
        ! z direction cyclic
        IF (CYCLIC_Z .OR. CYCLIC_Z_PD) THEN
          IF (CYCLIC_AT_T(IJK)) THEN
              I = I_OF(IJKT)
              J = J_OF(IJKT)
              K = K_OF(IJKT)

              IJK_CYCLIC = FUNIJK(I, J, KP1(K))

              QMOMK_N_TMP(:) = QMOMK_N0(:,IJKT,M)
              QMOMK_U_TMP(:) = QMOMK_U0(:,IJKT,M)
              QMOMK_V_TMP(:) = QMOMK_V0(:,IJKT,M)
              QMOMK_W_TMP(:) = QMOMK_W0(:,IJKT,M)

              QMOMK_N0(:,IJKT,M) = QMOMK_N0(:,IJK_CYCLIC,M)
              QMOMK_U0(:,IJKT,M) = QMOMK_U0(:,IJK_CYCLIC,M)
              QMOMK_V0(:,IJKT,M) = QMOMK_V0(:,IJK_CYCLIC,M)
              QMOMK_W0(:,IJKT,M) = QMOMK_W0(:,IJK_CYCLIC,M)

              QMOMK_N0(:,IJK_CYCLIC,M) = QMOMK_N_TMP(:)
              QMOMK_U0(:,IJK_CYCLIC,M) = QMOMK_U_TMP(:)
              QMOMK_V0(:,IJK_CYCLIC,M) = QMOMK_V_TMP(:)
              QMOMK_W0(:,IJK_CYCLIC,M) = QMOMK_W_TMP(:)
          END IF
        END IF
      END IF
     END DO
    END DO
  ELSE
    DO M = 1, MMAX
     DO IJK = ijkstart3, ijkend3
      IF (FLUID_AT(IJK)) THEN
        IJKN = NORTH_OF(IJK)
        IJKS = SOUTH_OF(IJK)
        IJKE = EAST_OF(IJK)
        IJKW = WEST_OF(IJK)
        IJKT = TOP_OF(IJK)
        IJKB = BOTTOM_OF(IJK)

        !  x direction cyclic
        IF (CYCLIC_X .OR. CYCLIC_X_PD) THEN
          IF (CYCLIC_AT_E(IJK)) THEN
              I = I_OF(IJKE)
              J = J_OF(IJKE)
              K = K_OF(IJKE)

              IJK_CYCLIC = FUNIJK(IP1(I), J, K)

              QMOMK_N_TMP(:) = QMOMK_N1(:,IJKE,M)
              QMOMK_U_TMP(:) = QMOMK_U1(:,IJKE,M)
              QMOMK_V_TMP(:) = QMOMK_V1(:,IJKE,M)
              QMOMK_W_TMP(:) = QMOMK_W1(:,IJKE,M)

              QMOMK_N1(:,IJKE,M) = QMOMK_N1(:,IJK_CYCLIC,M)
              QMOMK_U1(:,IJKE,M) = QMOMK_U1(:,IJK_CYCLIC,M)
              QMOMK_V1(:,IJKE,M) = QMOMK_V1(:,IJK_CYCLIC,M)
              QMOMK_W1(:,IJKE,M) = QMOMK_W1(:,IJK_CYCLIC,M)

              QMOMK_N1(:,IJK_CYCLIC,M) = QMOMK_N_TMP(:)
              QMOMK_U1(:,IJK_CYCLIC,M) = QMOMK_U_TMP(:)
              QMOMK_V1(:,IJK_CYCLIC,M) = QMOMK_V_TMP(:)
              QMOMK_W1(:,IJK_CYCLIC,M) = QMOMK_W_TMP(:)
          END IF
        END IF
        ! y direction cyclic
        IF (CYCLIC_Y .OR. CYCLIC_Y_PD) THEN
          IF (CYCLIC_AT_N(IJK)) THEN
              I = I_OF(IJKN)
              J = J_OF(IJKN)
              K = K_OF(IJKN)

              IJK_CYCLIC = FUNIJK(I, JP1(J), K)

              QMOMK_N_TMP(:) = QMOMK_N1(:,IJKN,M)
              QMOMK_U_TMP(:) = QMOMK_U1(:,IJKN,M)
              QMOMK_V_TMP(:) = QMOMK_V1(:,IJKN,M)
              QMOMK_W_TMP(:) = QMOMK_W1(:,IJKN,M)

              QMOMK_N1(:,IJKN,M) = QMOMK_N1(:,IJK_CYCLIC,M)
              QMOMK_U1(:,IJKN,M) = QMOMK_U1(:,IJK_CYCLIC,M)
              QMOMK_V1(:,IJKN,M) = QMOMK_V1(:,IJK_CYCLIC,M)
              QMOMK_W1(:,IJKN,M) = QMOMK_W1(:,IJK_CYCLIC,M)

              QMOMK_N1(:,IJK_CYCLIC,M) = QMOMK_N_TMP(:)
              QMOMK_U1(:,IJK_CYCLIC,M) = QMOMK_U_TMP(:)
              QMOMK_V1(:,IJK_CYCLIC,M) = QMOMK_V_TMP(:)
              QMOMK_W1(:,IJK_CYCLIC,M) = QMOMK_W_TMP(:)
          END IF
        END IF
        ! z direction cyclic
        IF (CYCLIC_Z .OR. CYCLIC_Z_PD) THEN
          IF (CYCLIC_AT_T(IJK)) THEN
              I = I_OF(IJKT)
              J = J_OF(IJKT)
              K = K_OF(IJKT)

              IJK_CYCLIC = FUNIJK(I, J, KP1(K))

              QMOMK_N_TMP(:) = QMOMK_N1(:,IJKT,M)
              QMOMK_U_TMP(:) = QMOMK_U1(:,IJKT,M)
              QMOMK_V_TMP(:) = QMOMK_V1(:,IJKT,M)
              QMOMK_W_TMP(:) = QMOMK_W1(:,IJKT,M)

              QMOMK_N1(:,IJKT,M) = QMOMK_N1(:,IJK_CYCLIC,M)
              QMOMK_U1(:,IJKT,M) = QMOMK_U1(:,IJK_CYCLIC,M)
              QMOMK_V1(:,IJKT,M) = QMOMK_V1(:,IJK_CYCLIC,M)
              QMOMK_W1(:,IJKT,M) = QMOMK_W1(:,IJK_CYCLIC,M)

              QMOMK_N1(:,IJK_CYCLIC,M) = QMOMK_N_TMP(:)
              QMOMK_U1(:,IJK_CYCLIC,M) = QMOMK_U_TMP(:)
              QMOMK_V1(:,IJK_CYCLIC,M) = QMOMK_V_TMP(:)
              QMOMK_W1(:,IJK_CYCLIC,M) = QMOMK_W_TMP(:)
          END IF
        END IF
      END IF
     END DO
    END DO
  END IF
 END SUBROUTINE QMOMK_CYCLIC_BC

END MODULE qmomk_bc
