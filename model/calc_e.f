! TO DO:
! Check the formulation based on MCp
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_e_e                                                C
!  Purpose: Calculate coefficients linking velocity correction to      C
!           solids volume fraction correction -- East                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 1-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_E_E(A_M, MCP, E_E)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Index of close packed solids phase. note mcp is the lowest index of
! those solids phases that have close_packed=t and of the solids
! phase that is used for the solids correction equation (when mmax=1).
      INTEGER, INTENT(IN) :: Mcp
! Coefficient for solids correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
! Note that the E_e coefficients are only used when mcp is assigned
! (when any solids phase has close_packed=T).
      DOUBLE PRECISION, INTENT(INOUT) :: e_e(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
!-----------------------------------------------

! returning if mcp is not defined.
      IF (MCP == UNDEFINED_I) RETURN

! returning if the x-momentum equation of this phase is not to be
! solved.
      IF (.NOT.MOMENTUM_X_EQ(MCP)) RETURN

!!$omp parallel do private(ijk)
      DO IJK = ijkstart3, ijkend3
         IF (SIP_AT_E(IJK) .OR. MFLOW_AT_E(IJK)) THEN
            E_E(IJK) = ZERO
         ELSE
            IF (A_M(IJK,0,MCP) /= ZERO) THEN
! calculating the correction coefficient
               E_E(IJK) = AYZ(IJK)/(-A_M(IJK,0,MCP))
            ELSE
               E_E(IJK) = LARGE_NUMBER
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_E_E

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_e_n                                                C
!  Purpose: Calculate coefficients linking velocity correction to      C
!           solids volume fraction correction -- North                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 1-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_E_N(A_M, MCP, E_N)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Index of close packed solids phase. note mcp is the lowest index of
! those solids phases that have close_packed=t and of the solids
! phase that is used for the solids correction equation (when mmax=1).
      INTEGER, INTENT(IN) :: Mcp
! Coefficient for solids correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
! Note that the E_n coefficients are only used when mcp is assigned
! (when the solids phase has close_packed=T)
      DOUBLE PRECISION, INTENT(INOUT) :: e_n(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
!-----------------------------------------------

! returning if mcp is not defined.
      IF (MCP == UNDEFINED_I) RETURN
! returning if the y-momentum equation of this phase is not to be
! solved.
      IF (.NOT.MOMENTUM_Y_EQ(MCP)) RETURN

!!$omp parallel do private(IJK)
      DO IJK = ijkstart3, ijkend3
         IF (SIP_AT_N(IJK) .OR. MFLOW_AT_N(IJK)) THEN
            E_N(IJK) = ZERO
         ELSE
            IF ((-A_M(IJK,0,MCP)) /= ZERO) THEN
! calculating the correction coefficient
               E_N(IJK) = AXZ(IJK)/(-A_M(IJK,0,MCP))
            ELSE
               E_N(IJK) = LARGE_NUMBER
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_E_N

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_e_t                                                C
!  Purpose: Calculate coefficients linking velocity correction to      C
!           solids volume fraction correction -- Top                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 1-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CALC_E_T(A_M, MCP, E_T)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE run
      USE constant
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(IN) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Index of close packed solids phase. note mcp is the lowest index of
! those solids phases that have close_packed=t and of the solids
! phase that is used for the solids correction equation (when mmax=1).
      INTEGER, INTENT(IN) :: Mcp
! Coefficient for solids correction equation. These coefficients are
! initialized to zero in the subroutine time_march before the time loop.
! Note that the E_t coefficients are only used when mcp is assigned
! (when the solids phase has close_packed=T)
      DOUBLE PRECISION, INTENT(INOUT) :: e_t(DIMENSION_3)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK
!-----------------------------------------------

! returning if mcp is not defined.
      IF (MCP == UNDEFINED_I) RETURN
! returning if the z-momentum equation of this phase is not to be
! solved
      IF (.NOT.MOMENTUM_Z_EQ(MCP)) RETURN

!!$omp parallel do private(IJK)
      DO IJK = ijkstart3, ijkend3
         IF (SIP_AT_T(IJK) .OR. MFLOW_AT_T(IJK)) THEN
            E_T(IJK) = ZERO
         ELSE
            IF ((-A_M(IJK,0,MCP)) /= ZERO) THEN
! calculating the correction coefficient
               E_T(IJK) = AXY(IJK)/(-A_M(IJK,0,MCP))
            ELSE
               E_T(IJK) = LARGE_NUMBER
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_E_T

