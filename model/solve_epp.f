!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_Epp                                               C
!  Purpose: Solve solids volume fraction correction equation.          C
!                                                                      C
!  Notes: MCP must be defined to call this routine.                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 25-SEP-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_EPP(NORMS, RESS, IER)

! Modules
!---------------------------------------------------------------------//
      use ambm, only: a_m, b_m, lock_ambm, unlock_ambm
      use geometry, only: ijkmax2
      use leqsol, only: leq_it, leq_method, leq_sweep, leq_pc, leq_tol
      use param, only: dimension_3, dimension_m
      use param1, only: undefined_i, zero, one
      use ps, only: point_source
      use pscor, only: epp, mcp
      use residual, only: resid, max_resid, ijk_resid
      use residual, only: num_resid, den_resid
      use residual, only: resid_p
      use run, only: momentum_x_eq, momentum_y_eq
      use usr_src, only: call_usr_source, calc_usr_source
      use usr_src, only: solids_correction
      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//
! Parameter to make tolerance for residual scaled with max value
! compatible with residual scaled with first iteration residual.
! Increase it to tighten convergence.
      DOUBLE PRECISION, PARAMETER :: DEN = 1.0D1 !5.0D2

! Dummy arguments
!---------------------------------------------------------------------//
! Normalization factor for solids volume fraction correction residual.
! At start of the iterate loop norms will either be 1 (i.e. not
! normalized) or a user defined value given by norm_s.  If norm_s
! was set to zero then the normalization is based on dominate
! term in the equation
      DOUBLE PRECISION, INTENT(IN) :: NORMs
! solids volume fraction correction residual
      DOUBLE PRECISION, INTENT(OUT) :: RESs
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local variables
!-----------------------------------------------
! solids phase index locally assigned to mcp
! mcp is the lowest index of those solids phases that are close_packed
! and of the solids phase that is used for the solids correction
! equation.
      INTEGER :: M
! Normalization factor for solids volume fraction correction residual
      DOUBLE PRECISION :: NORMSloc
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

! temporary use of global arrays:
! arraym1 (locally b_mmax)
! vector B_M based on dominate term in correction equation
      DOUBLE PRECISION :: B_MMAX(DIMENSION_3, DIMENSION_M)
! Septadiagonal matrix A_m, vector B_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!---------------------------------------------------------------------//

! Note that currently this subroutine is only called when MMAX=1 AND
! MCP is defined. This combintion effectively means solids phase 1
! can close pack.
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
      call lock_ambm

! Form the sparse matrix equation.  Note that the index 0 is explicitly
! used throughout this routine for creating the matrix equation.
! However, the equation is based on the index of MCP.

! initializing
      CALL INIT_AB_M (A_M, B_M, IJKMAX2, 0)
      EPP(:) = ZERO

      CALL CONV_SOURCE_EPP (A_M, B_M, B_mmax, M)

! Add point source contributions.
      IF(POINT_SOURCE) CALL POINT_SOURCE_EPP (B_M, B_mmax, M)

! Add usr source contributions
      IF(CALL_USR_SOURCE(2)) CALL CALC_USR_SOURCE(SOLIDS_CORRECTION, &
                           A_M, B_M, lB_MMAX=B_MMAX, lM=M)

!      call check_ab_m(a_m, b_m, 0, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!      call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, &
!         DO_K, ier)


! Find average residual, maximum residual and location
      NORMSloc = NORMS
      IF(NORMS == ZERO) THEN
! calculate the residual based on dominate term in correction equation
! and use this to form normalization factor
         CALL CALC_RESID_PP (B_MMAX, ONE, NUM_RESID(RESID_P,M), &
            DEN_RESID(RESID_P,M), RESID(RESID_P,M), &
            MAX_RESID(RESID_P,M), IJK_RESID(RESID_P,M))
         NORMSloc = RESID(RESID_P,M)/DEN
      ENDIF

      CALL CALC_RESID_PP (B_M, NORMSloc, NUM_RESID(RESID_P,M), &
         DEN_RESID(RESID_P,M), RESID(RESID_P,M), MAX_RESID(RESID_P,M), &
         IJK_RESID(RESID_P,M))
      RESS = RESID(RESID_P,M)
!      write(*,*) resid(resid_p, 1), max_resid(resid_p, 1), &
!         ijk_resid(resid_p, 1)


! Solve EP_s_prime equation
      CALL ADJUST_LEQ(RESID(RESID_P,M), LEQ_IT(2), LEQ_METHOD(2),&
                      LEQI, LEQM)
! note index 0 is used here since that is the index that was used for
! creating this matrix equation
      CALL SOLVE_LIN_EQ ('EPp', 2, EPP, A_M, B_M, 0, LEQI, LEQM, &
                         LEQ_SWEEP(2), LEQ_TOL(2), LEQ_PC(2), IER)

!      call out_array(EPp, 'EPp')

      call unlock_ambm

      RETURN
      END SUBROUTINE SOLVE_EPP
