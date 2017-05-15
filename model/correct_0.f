!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CORRECT_0                                               C
!  Purpose: Correct the fluid pressure and gas velocities              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CORRECT_0()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE pgcor
      USE ur_facs
      IMPLICIT NONE
!-----------------------------------------------

      CALL CORRECT_0G (PP_G, UR_FAC(1), D_E, D_N, D_T, P_G, &
         U_G, V_G, W_G)
!      CALL CORRECT_0S (PP_G, D_E, D_N, D_T, U_S, V_S, W_S, IER)

      RETURN
      END SUBROUTINE CORRECT_0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CORRECT_0g                                              C
!  Purpose: Correct the fluid pressure and velocities.                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CORRECT_0G(PP_G,UR_FAC,D_E,D_N,D_T,P_G,U_G,V_G,W_G)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE cutcell
      USE quadric
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Pressure correction
      DOUBLE PRECISION, INTENT(IN) :: Pp_g(DIMENSION_3)
! Under relaxation factor for Pressure correction
      DOUBLE PRECISION, INTENT(IN) :: UR_fac
! Pressure correction coefficient -- East
      DOUBLE PRECISION, INTENT(IN) :: d_e(DIMENSION_3, 0:DIMENSION_M)
! Pressure correction coefficient -- North
      DOUBLE PRECISION, INTENT(IN) :: d_n(DIMENSION_3, 0:DIMENSION_M)
! Pressure correction coefficient -- Top
      DOUBLE PRECISION, INTENT(IN) :: d_t(DIMENSION_3, 0:DIMENSION_M)
! Pressure
      DOUBLE PRECISION, INTENT(INOUT) :: P_g(DIMENSION_3)
! Velocity components
      DOUBLE PRECISION, INTENT(INOUT) :: U_g(DIMENSION_3), &
                       V_g(DIMENSION_3),&
                       W_g(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKE, IJKN, IJKT
!-----------------------------------------------

! Underrelax pressure correction.  Velocity corrections should not be
! underrelaxed, so that the continuity eq. is satisfied.

!!$omp    parallel do private(IJK,IJKE,IJKN,IJKT)
      DO IJK = ijkstart3, ijkend3

         IF (FLUIDORP_FLOW_AT(IJK)) THEN

            P_G(IJK) = P_G(IJK) + UR_FAC*PP_G(IJK)

            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IF(.NOT.CARTESIAN_GRID) THEN
               U_G(IJK) = U_G(IJK) - D_E(IJK,0)*(PP_G(IJKE)-PP_G(IJK))
               V_G(IJK) = V_G(IJK) - D_N(IJK,0)*(PP_G(IJKN)-PP_G(IJK))
               IF (DO_K) THEN
                  IJKT = TOP_OF(IJK)
                  W_G(IJK) = W_G(IJK) - D_T(IJK,0)*(PP_G(IJKT)-PP_G(IJK))
               ENDIF
            ELSE
               U_G(IJK) = U_G(IJK) - D_E(IJK,0)*&
                  (PP_G(IJKE)*A_UPG_E(IJK) - PP_G(IJK)*A_UPG_W(IJK))
               V_G(IJK) = V_G(IJK) - D_N(IJK,0)*&
                  (PP_G(IJKN)*A_VPG_N(IJK) - PP_G(IJK)*A_VPG_S(IJK))
               IF (DO_K) THEN
                  IJKT = TOP_OF(IJK)
                  W_G(IJK) = W_G(IJK) - D_T(IJK,0)*&
                     (PP_G(IJKT)*A_WPG_T(IJK) - PP_G(IJK)*A_WPG_B(IJK))
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CORRECT_0G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CORRECT_0s                                              C
!  Purpose: Correct the solids velocities.                             C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-JUN-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CORRECT_0S(PP_G, D_E, D_N, D_T, U_S, V_S, W_S)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE cutcell
      USE quadric
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Pressure correction
      DOUBLE PRECISION, INTENT(IN) :: Pp_g(DIMENSION_3)
! Pressure correction coefficient -- East
      DOUBLE PRECISION, INTENT(IN) :: d_e(DIMENSION_3, 0:DIMENSION_M)
! Pressure correction coefficient -- North
      DOUBLE PRECISION, INTENT(IN) :: d_n(DIMENSION_3, 0:DIMENSION_M)
! Pressure correction coefficient -- Top
      DOUBLE PRECISION, INTENT(IN) :: d_t(DIMENSION_3, 0:DIMENSION_M)
! Velocity components
      DOUBLE PRECISION, INTENT(INOUT) :: U_s(DIMENSION_3, DIMENSION_M),&
                       V_s(DIMENSION_3, DIMENSION_M),&
                       W_s(DIMENSION_3, DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Solids index
      INTEGER :: M
! Indices
      INTEGER :: IJK, IJKE, IJKN, IJKT
!-----------------------------------------------
! Velocity corrections should not be underrelaxed, so that
! the continuity eq. is satisfied.

      DO M = 1, MMAX
!!$omp    parallel do private(IJK,IJKE,IJKN,IJKT)
         DO IJK = ijkstart3, ijkend3
            IF (FLUIDORP_FLOW_AT(IJK)) THEN

               IJKE = EAST_OF(IJK)
               IJKN = NORTH_OF(IJK)
               IF(.NOT.CARTESIAN_GRID) THEN
                  U_S(IJK,M) = U_S(IJK,M) - D_E(IJK,M)*&
                     (PP_G(IJKE)-PP_G(IJK))
                  V_S(IJK,M) = V_S(IJK,M) - D_N(IJK,M)*&
                     (PP_G(IJKN)-PP_G(IJK))
                  IF (DO_K) THEN
                     IJKT = TOP_OF(IJK)
                     W_S(IJK,M) = W_S(IJK,M) - D_T(IJK,M)*&
                        (PP_G(IJKT)-PP_G(IJK))
                  ENDIF
               ELSE
                  U_S(IJK,M) = U_S(IJK,M) - D_E(IJK,M)*&
                     (PP_G(IJKE)*A_UPG_E(IJK) - PP_G(IJK)*A_UPG_W(IJK))
                  V_S(IJK,M) = V_S(IJK,M) - D_N(IJK,M)*&
                     (PP_G(IJKN)*A_VPG_N(IJK) - PP_G(IJK)*A_VPG_S(IJK))
                  IF (DO_K) THEN
                     IJKT = TOP_OF(IJK)
                     W_S(IJK,M) = W_S(IJK,M) - D_T(IJK,M)*&
                        (PP_G(IJKT)*A_WPG_T(IJK) - PP_G(IJK)*A_WPG_B(IJK))
                  ENDIF
               ENDIF

            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE CORRECT_0S


