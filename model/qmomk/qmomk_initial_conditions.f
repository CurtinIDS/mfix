!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: QMOMK_INITIAL_CONDITIONS                               C
!  Author: Alberto Passalacqua                        Date: 30-Jul-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE QMOMK_INITIAL_CONDITIONS

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
      USE qmom_kinetic_equation
      USE qmomk_quadrature
      USE qmomk_parameters
      USE ic
      USE functions

      IMPLICIT NONE

      DOUBLE PRECISION InitVal
      INTEGER :: I, M, IJK

      DO IJK = IJKSTART3, IJKEND3

        IF(.NOT.FLUID_AT(IJK)) cycle
        DO M = 1, MMAX
          ! A.P.
          ! Initializing weights as volume_fraction/number_of_nodes
          ! Note that MFIX doesn't initialize the volume fraction in
          ! boundaries, and QMOM need a positive value there.
          ! Here the boundaries are not considered, and they will
          ! be initialized separetely in the correspoding module

          DO I = 1, QMOMK_NN
             QMOMK_N0 (I, IJK, M) = ROP_S(IJK, M)/(QMOMK_NN * RO_s(IJK,M))
          END DO

          ! A.P. Granular temperature minimum value is bounded
          !      This is required by QMOM
          InitVal = MAX(SQRT(THETA_M(IJK,M)), MINIMUM_THETA)

          QMOMK_U0(1, IJK, M) = -InitVal + U_S(IJK, M)
          QMOMK_U0(2, IJK, M) = +InitVal + U_S(IJK, M)
          QMOMK_U0(3, IJK, M) = -InitVal + U_S(IJK, M)
          QMOMK_U0(4, IJK, M) = +InitVal + U_S(IJK, M)
          QMOMK_U0(5, IJK, M) = -InitVal + U_S(IJK, M)
          QMOMK_U0(6, IJK, M) = +InitVal + U_S(IJK, M)
          QMOMK_U0(7, IJK, M) = -InitVal + U_S(IJK, M)
          QMOMK_U0(8, IJK, M) = +InitVal + U_S(IJK, M)

          QMOMK_V0(1, IJK, M) = -InitVal + V_S(IJK, M)
          QMOMK_V0(2, IJK, M) = -InitVal + V_S(IJK, M)
          QMOMK_V0(3, IJK, M) = +InitVal + V_S(IJK, M)
          QMOMK_V0(4, IJK, M) = +InitVal + V_S(IJK, M)
          QMOMK_V0(5, IJK, M) = -InitVal + V_S(IJK, M)
          QMOMK_V0(6, IJK, M) = -InitVal + V_S(IJK, M)
          QMOMK_V0(7, IJK, M) = +InitVal + V_S(IJK, M)
          QMOMK_V0(8, IJK, M) = +InitVal + V_S(IJK, M)

          QMOMK_W0(1, IJK, M) = -InitVal + W_S(IJK, M)
          QMOMK_W0(2, IJK, M) = -InitVal + W_S(IJK, M)
          QMOMK_W0(3, IJK, M) = -InitVal + W_S(IJK, M)
          QMOMK_W0(4, IJK, M) = -InitVal + W_S(IJK, M)
          QMOMK_W0(5, IJK, M) = +InitVal + W_S(IJK, M)
          QMOMK_W0(6, IJK, M) = +InitVal + W_S(IJK, M)
          QMOMK_W0(7, IJK, M) = +InitVal + W_S(IJK, M)
          QMOMK_W0(8, IJK, M) = +InitVal + W_S(IJK, M)
        END DO
      END DO

      DO IJK = IJKSTART3, IJKEND3
       DO M = 1, MMAX
         IF (FLUID_AT(IJK)) THEN
           CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N0(:,IJK,M), &
                QMOMK_U0(:,IJK,M), QMOMK_V0(:,IJK,M), QMOMK_W0(:,IJK,M), QMOMK_M0(:,IJK,M))

           CALL EIGHT_NODE_3D (QMOMK_M0(:,IJK,M), QMOMK_N0(:,IJK,M), &
                QMOMK_U0(:,IJK,M), QMOMK_V0(:,IJK,M), QMOMK_W0(:,IJK,M))

           CALL MOMENTS_TWENTY_EIGHT_NODES (QMOMK_N0(:,IJK,M), &
                QMOMK_U0(:,IJK,M), QMOMK_V0(:,IJK,M), QMOMK_W0(:,IJK,M), QMOMK_M0(:,IJK,M))
         END IF
       END DO
      END DO

      END SUBROUTINE QMOMK_INITIAL_CONDITIONS
