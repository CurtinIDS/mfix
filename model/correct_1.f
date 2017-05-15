! TO DO:
! 1. check the formulation based on MCp.

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CORRECT_1                                               C
!  Purpose: Correct the solids volume fraction                         C
!                                                                      C
!  Notes: MCP must be defined to call this routine.                    C
!         The solids volume fraction correction routines are           C
!         explicitly setup for a single solids phases.  Therefore,     C
!         even if multiple solids phases are present, only one of      C
!         the solids phases (M=MCP) can employ the correction          C
!         routines                                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-SEP-96  C
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

      SUBROUTINE CORRECT_1()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1, only: one, zero, undefined_i
      USE fldvar, only: u_s, v_s, w_s, rop_s, ro_s, ep_s
      USE physprop, only: close_packed
      USE pscor, only: mcp, epp, k_cp, e_e, e_n, e_t
      USE ur_facs, only: ur_fac
      USE indices
      USE geometry
      USE compar
      USE sendrecv
      USE cutcell, only: cartesian_grid, cut_cell_at, cg_ur_fac
      USE visc_s, only: ep_star_array
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! corrected solids volume fraction
      DOUBLE PRECISION :: EPcor
! dPodEP_s(EP_s(IJK, M)) * EPp(IJK)
      DOUBLE PRECISION :: Pp_P
! indices
      INTEGER :: IJK, IJKE, IJKN, IJKT
! solids index
      INTEGER :: M
! solids volume fraction at maximum packing
      DOUBLE PRECISION :: EP_S_CP
!-----------------------------------------------

      IF (MCP == UNDEFINED_I) THEN
! this error should be caught earlier in the routines so that this
! branch should never be entered
         RETURN
      ELSE
! the lowest solids phase index of those solids phases that can close
! pack (i.e. close_packed=T) and the index of the solids phase that is
! used to form the solids correction equation.
         M = MCP
      ENDIF


! by definition M must be close_packed (this is a redundant check)
      IF (CLOSE_PACKED(M)) THEN

! Correct solids volume fraction
! ---------------------------------------------------------------->>>
!!$omp    parallel do &
!!$omp&   private( IJK, EPCOR )
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN

               EPCOR = EP_S(IJK,M) + EPP(IJK)
               EP_S_CP = ONE - EP_STAR_ARRAY(IJK)
               IF (EPCOR>EP_S_CP .AND. EPP(IJK)>ZERO) THEN
                  EPP(IJK) = UR_FAC(2)*EPP(IJK)

                  IF(CARTESIAN_GRID) THEN
! JFD: Using a low value of CG_UR_FAC(2) for the cut cell tends to
! stabilize code (limits occurrence of negative void fraction)
                     IF(CUT_CELL_AT(IJK)) EPP(IJK) = CG_UR_FAC(2)*EPP(IJK)
                  ENDIF

                  EPCOR = EP_S(IJK,M) + EPP(IJK)
               ENDIF
               ROP_S(IJK,M) = MAX(ZERO,RO_S(IJK,M)*EPCOR)
            ENDIF
         ENDDO
! end correct solids volume fraction
! ----------------------------------------------------------------<<<

! Correct solids velocities
! ---------------------------------------------------------------->>>
!!$omp    parallel do &
!!$omp&   private( IJK, PP_P, IJKE, IJKN, IJKT )
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN

               PP_P = K_CP(IJK)*EPP(IJK)
               IF (FLOW_AT_E(IJK)) THEN
                  IJKE = EAST_OF(IJK)
                  U_S(IJK,M)=U_S(IJK,M)-E_E(IJK)*&
                     (K_CP(IJKE)*EPP(IJKE)-PP_P)
               ENDIF
               IF (FLOW_AT_N(IJK)) THEN
                  IJKN = NORTH_OF(IJK)
                  V_S(IJK,M)=V_S(IJK,M)-E_N(IJK)*&
                     (K_CP(IJKN)*EPP(IJKN)-PP_P)
               ENDIF
               IF (DO_K) THEN
                  IF (FLOW_AT_T(IJK)) THEN
                     IJKT = TOP_OF(IJK)
                     W_S(IJK,M) = W_S(IJK,M) - E_T(IJK)*&
                        (K_CP(IJKT)*EPP(IJKT)-PP_P)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
! end correct solids velocities
! ----------------------------------------------------------------<<<

      ENDIF   ! end if (close_packed)

      RETURN
      END SUBROUTINE CORRECT_1


