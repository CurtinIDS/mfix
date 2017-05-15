!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_Pp_g                                              C
!  Purpose: Solve fluid pressure correction equation                   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-JUN-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_PP_G(NORMG, RESG, IER)

! Modules
!---------------------------------------------------------------------//
      use ambm, only: a_m, b_m, lock_ambm, unlock_ambm
      use geometry, only: ijkmax2
      use leqsol, only: leq_it, leq_method, leq_sweep, leq_pc, leq_tol
      use param, only: dimension_3, dimension_m
      use param1, only: zero, one, undefined
      use pgcor, only: pp_g
      use physprop, only: mmax, ro_g0
      use ps, only: point_source
      use residual, only: resid, max_resid, ijk_resid
      use residual, only: num_resid, den_resid
      use residual, only: resid_p
      use run, only: momentum_x_eq, momentum_y_eq
      use usr_src, only: call_usr_source, calc_usr_source
      use usr_src, only: pressure_correction
      IMPLICIT NONE

! Local parameters
!---------------------------------------------------------------------//
! Parameter to make tolerance for residual scaled with max value
! compatible with residual scaled with first iteration residual.
! Increase it to tighten convergence.
      DOUBLE PRECISION, PARAMETER :: DEN = 1.0D1   !5.0D2

! Dummy arguments
!---------------------------------------------------------------------//
! Normalization factor for gas pressure correction residual.
! At start of the iterate loop normg will either be 1 (i.e. not
! normalized) or a user defined value given by norm_g.  If norm_g
! was set to zero then the normalization is based on dominate
! term in the equation
      DOUBLE PRECISION, INTENT(IN) :: NORMg
! gas pressure correction residual
      DOUBLE PRECISION, INTENT(OUT) :: RESg
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local variables
!---------------------------------------------------------------------//
! phase index
      INTEGER :: M
! Normalization factor for gas pressure correction residual
      DOUBLE PRECISION :: NORMGloc
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

      call lock_ambm

! initializing
      PP_G(:) = ZERO
      DO M = 0, MMAX
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)
      ENDDO

! If gas momentum equations in x and y directions are not solved return
      IF (.NOT.(MOMENTUM_X_EQ(0) .OR. MOMENTUM_Y_EQ(0)) .AND.&
          RO_G0 .NE. UNDEFINED) THEN
        call unlock_ambm
        RETURN
      ENDIF

! Forming the sparse matrix equation.
      CALL CONV_PP_G (A_M, B_M)
      CALL SOURCE_PP_G (A_M, B_M, B_MMAX)
      IF(POINT_SOURCE) CALL POINT_SOURCE_PP_G (B_M, B_MMAX)
      IF(CALL_USR_SOURCE(1)) CALL CALC_USR_SOURCE(Pressure_correction,&
                           A_M, B_M, lB_MMAX=B_MMAX, lM=0)

!      call check_ab_m(a_m, b_m, 0, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)


! Find average residual, maximum residual and location
      NORMGloc = NORMG
      IF(NORMG == ZERO) THEN
! calculating the residual based on dominate term in correction equation
! and use this to form normalization factor
        CALL CALC_RESID_PP (B_MMAX, ONE, NUM_RESID(RESID_P,0), &
         DEN_RESID(RESID_P,0), RESID(RESID_P,0), MAX_RESID(RESID_P,0), &
         IJK_RESID(RESID_P,0))
         NORMGloc = RESID(RESID_P,0)/DEN
      ENDIF
      CALL CALC_RESID_PP (B_M, NORMGloc, NUM_RESID(RESID_P,0),  &
         DEN_RESID(RESID_P,0), RESID(RESID_P,0), MAX_RESID(RESID_P,0), &
         IJK_RESID(RESID_P,0))
      RESG = RESID(RESID_P,0)
!      write(*,*) resid(resid_p, 0), max_resid(resid_p, 0), &
!         ijk_resid(resid_p, 0)


! Solve P_g_prime equation
       LEQI = LEQ_IT(1)
       LEQM = LEQ_METHOD(1)
!      CALL ADJUST_LEQ(RESID(RESID_P,0),LEQ_IT(1),LEQ_METHOD(1),LEQI,LEQM,IER)

!     call check_symmetry(A_m, 0, IER)
!     call test_lin_eq(A_M, LEQ_IT(1),LEQ_METHOD(1), LEQ_SWEEP(1), LEQ_TOL(1), LEQ_PC(1),0,IER)
      CALL SOLVE_LIN_EQ ('Pp_g', 1, PP_G, A_M, B_M, 0, LEQI, LEQM, &
                         LEQ_SWEEP(1), LEQ_TOL(1), LEQ_PC(1), IER)

!      call out_array(Pp_g, 'Pp_g')

      call unlock_ambm

      RETURN
      END SUBROUTINE SOLVE_PP_G

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_Pp_g                                       C
!  Purpose: Adds point sources to the Pressure correction equation.    C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative. See          C
!         conv_Pp_g                                                    C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_PP_G(B_M, B_mmax)

! Modules 
!-----------------------------------------------
      use compar, only: dead_cell_at
      use geometry, only: vol
      use functions, only: fluid_at, funijk
      use functions, only: is_on_myPe_plus2layers
      use param, only: dimension_3, dimension_m
      use param1, only: small_number
      use ps, only: ps_defined, dimension_ps
      use ps, only: ps_massflow_g, ps_volume
      use ps, only: ps_k_b, ps_k_t
      use ps, only: ps_j_s, ps_j_n
      use ps, only: ps_i_w, ps_i_e
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax(DIMENSION_3, 0:DIMENSION_M)

! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV

! terms of bm expression
      DOUBLE PRECISION pSource

!-----------------------------------------------
      PS_LP: do PSV = 1, DIMENSION_PS

         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(PS_MASSFLOW_G(PSV) < small_number) cycle PS_LP

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

            ijk = funijk(i,j,k)
            if(fluid_at(ijk)) then
               pSource = PS_MASSFLOW_G(PSV) * (VOL(IJK)/PS_VOLUME(PSV))

               B_M(IJK,0) = B_M(IJK,0) - pSource
               B_MMAX(IJK,0) = max(abs(B_MMAX(IJK,0)), abs(B_M(IJK,0)))
            endif

         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_PP_G
