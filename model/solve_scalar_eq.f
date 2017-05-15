!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_SCALAR_EQ                                         C
!  Purpose: Solve scalar transport equations                           C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 4-12-99    C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_Scalar_EQ(IER)

! Modules
!---------------------------------------------------------------------//
      use ambm, only: a_m, b_m, lock_ambm, unlock_ambm
      use bc, only: bc_scalar, bc_hw_scalar, bc_c_scalar, bc_scalarw
      use compar, only: ijkstart3, ijkend3
      use fldvar, only: rop_g, rop_go, ep_g, u_g, v_g, w_g
      use fldvar, only: rop_s, rop_so, ep_s, u_s, v_s, w_s
      use fldvar, only: scalar, scalaro
      use functions, only: fluid_at, zmax
      use geometry, only: vol, ijkmax2
      use leqsol, only: leq_it, leq_method, leq_sweep, leq_pc, leq_tol
      use mflux, only: flux_ge, flux_gn, flux_gt
      use mflux, only: flux_se, flux_sn, flux_st
      use mflux, only: flux_gse, flux_gsn, flux_gst
      use mflux, only: flux_sse, flux_ssn, flux_sst
      use param, only: dimension_3
      use param1, only: zero
      use residual, only: resid, max_resid, ijk_resid
      use residual, only: num_resid, den_resid
      use residual, only: resid_sc
      use run, only: discretize, added_mass, M_AM, odt
      use rxns, only: sum_r_g, sum_r_s
      use scalars, only: nscalar, phase4scalar
      use scalars, only: scalar_c, scalar_p, dif_scalar
      use toleranc, only: zero_ep_s
      use ur_facs, only: ur_fac
      use usr_src, only: call_usr_source, calc_usr_source
      use usr_src, only: usr_scalar
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local variables
!---------------------------------------------------------------------//
! phase index
      INTEGER :: M
! species index
      INTEGER :: NN
! Indices
      INTEGER :: IJK
!
      DOUBLE PRECISION :: APO
!
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI
! temporary variables in residual computation
      DOUBLE PRECISION :: res1, mres1, num_res, den_res
      INTEGER :: ires1
! source vector: coefficient of dependent variable
! becomes part of a_m matrix; must be positive
      DOUBLE PRECISION :: S_P(DIMENSION_3)
! source vector: constant part becomes part of b_m vector
      DOUBLE PRECISION :: S_C(DIMENSION_3)
!
      DOUBLE PRECISION :: EPS(DIMENSION_3)
!
      character(LEN=8) :: Vname
!---------------------------------------------------------------------//

      call lock_ambm

      RESID(RESID_sc,0) = ZERO
      NUM_RESID(RESID_sc,0) = ZERO
      DEN_RESID(RESID_sc,0) = ZERO
      MAX_RESID(RESID_sc,0) = ZERO
      IJK_RESID(RESID_sc,0) = 0

! Loop over number of scalar equations
      DO NN = 1, NScalar

! Setting M to the phase assigned to convecting the indicated scalar
! equation
         M = Phase4Scalar(NN)
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)

! Gas phase
! ---------------------------------------------------------------->>>
         IF(M == 0) THEN

            DO IJK = IJKSTART3, IJKEND3
               IF (FLUID_AT(IJK)) THEN
                  APO = ROP_GO(IJK)*VOL(IJK)*ODT
                  S_P(IJK) = APO + (ZMAX(SUM_R_G(IJK)) + &
                                         Scalar_p(IJK,NN))*VOL(IJK)
                  S_C(IJK) = APO*ScalarO(IJK,NN) + &
                             Scalar(IJK,NN)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) + &
                             Scalar_c(IJK, NN)*VOL(IJK)
               ELSE
                  S_P(IJK) = ZERO
                  S_C(IJK) = ZERO
               ENDIF
            ENDDO

            IF(.NOT.ADDED_MASS) THEN
               CALL CONV_DIF_PHI (Scalar(1,NN), DIF_Scalar(1,NN), &
                                  DISCRETIZE(9), U_G, V_G, W_G, &
                                  Flux_gE, Flux_gN, Flux_gT, M, A_M, B_M)
            ELSE
               CALL CONV_DIF_PHI (Scalar(1,NN), DIF_Scalar(1,NN), &
                                  DISCRETIZE(9), U_G, V_G, W_G, &
                                  Flux_gSE, Flux_gSN, Flux_gST, M, A_M, B_M)
            ENDIF

            CALL BC_PHI (Scalar(1,NN), BC_Scalar(1,NN), &
                         BC_ScalarW(1,NN), BC_HW_Scalar(1,NN), &
                         BC_C_Scalar(1,NN), M, A_M, B_M)

            CALL SOURCE_PHI (S_P, S_C, EP_G, Scalar(1,NN), M, A_M, B_M)

! calculate any usr source terms
            IF (CALL_USR_SOURCE(9)) CALL CALC_USR_SOURCE(USR_SCALAR,&
                                  A_M, B_M, lM=M, lN=NN)

            CALL CALC_RESID_S (Scalar(1,NN), A_M, B_M, M, num_res, &
                               den_res, res1, mres1, ires1, ZERO)
            RESID(RESID_sc,0) = RESID(RESID_sc,0)+res1
            NUM_RESID(RESID_sc,0) = NUM_RESID(RESID_sc,0)+num_res
            DEN_RESID(RESID_sc,0) = DEN_RESID(RESID_sc,0)+den_res
            if(mres1 .gt. MAX_RESID(RESID_sc,0))then
               MAX_RESID(RESID_sc,0) = mres1
               IJK_RESID(RESID_sc,0) = ires1
            endif

            CALL UNDER_RELAX_S (Scalar(1,NN), A_M, B_M, M, UR_FAC(9))

!           call check_ab_m(a_m, b_m, m, .false., ier)
!           call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!           write(*,*) res1, mres1, ires1
!           call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, &
!                            DO_K, ier)

            CALL ADJUST_LEQ (res1, LEQ_IT(9), LEQ_METHOD(9), &
               LEQI, LEQM)

            write(Vname, '(A,I2)')'Scalar',NN
            CALL SOLVE_LIN_EQ (Vname, 9, Scalar(1,NN), A_M, B_M, M,&
                               LEQI, LEQM, LEQ_SWEEP(9), LEQ_TOL(9),&
                               LEQ_PC(9), IER)
!           call out_array(Scalar(1, N), Vname)

! Solids phase
! ---------------------------------------------------------------->>>
         ELSE   ! (m/=0) , i.e., solids phase

            DO IJK = IJKSTART3, IJKEND3
               IF (FLUID_AT(IJK)) THEN
                  APO = ROP_sO(IJK, M)*VOL(IJK)*ODT
                  S_P(IJK) = APO + (ZMAX(SUM_R_s(IJK, M)) + &
                                         Scalar_p(IJK, NN))*VOL(IJK)
                  S_C(IJK) = APO*ScalarO(IJK,NN) + &
                             Scalar(IJK,NN)*ZMAX((-SUM_R_s(IJK, M)))*VOL(IJK)+&
                             Scalar_c(IJK, NN)*VOL(IJK)
                  EPs(IJK) = EP_s(IJK, M)
               ELSE
                  S_P(IJK) = ZERO
                  S_C(IJK) = ZERO
                  EPS(IJK) = ZERO
               ENDIF
            ENDDO

            IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
               CALL CONV_DIF_PHI (Scalar(1,NN), DIF_Scalar(1,NN), &
                                  DISCRETIZE(9), U_s(1,m), V_s(1,m), &
                                  W_s(1,m), Flux_sE(1,M), Flux_sN(1,M),&
                                  Flux_sT(1,M), M, A_M, B_M)
            ELSE ! virtual mass term added for M = M_AM ONLY!!!!
               CALL CONV_DIF_PHI (Scalar(1,NN), DIF_Scalar(1,NN), &
                                  DISCRETIZE(9), U_s(1,m), V_s(1,m), &
                                  W_s(1,m), Flux_sSE, Flux_sSN, Flux_sST, &
                                  M, A_M, B_M)
            ENDIF

            CALL BC_PHI (Scalar(1,NN), BC_Scalar(1,NN), BC_ScalarW(1,NN), &
                         BC_HW_Scalar(1,NN), BC_C_Scalar(1,NN), M, A_M, B_M)

            CALL SOURCE_PHI (S_P, S_C, EPs, Scalar(1,NN), M, A_M, B_M)

! calculate any usr source terms
            IF (CALL_USR_SOURCE(9)) CALL CALC_USR_SOURCE(USR_SCALAR,&
                                  A_M, B_M, lM=M, lN=NN)


            CALL CALC_RESID_S (Scalar(1,NN), A_M, B_M, M, num_res, &
                              den_res, res1, mres1, ires1, ZERO)
            RESID(RESID_sc,0) = RESID(RESID_sc,0)+res1
            NUM_RESID(RESID_sc,0) = NUM_RESID(RESID_sc,0)+num_res
            DEN_RESID(RESID_sc,0) = DEN_RESID(RESID_sc,0)+den_res
            if(mres1 .gt. MAX_RESID(RESID_sc,0))then
               MAX_RESID(RESID_sc,0) = mres1
               IJK_RESID(RESID_sc,0) = ires1
            endif

            CALL UNDER_RELAX_S (Scalar(1,NN), A_M, B_M, M, UR_FAC(9))

!            call check_ab_m(a_m, b_m, m, .false., ier)
!            call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!            write(*,*) res1, mres1, ires1
!            call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, &
!                             DO_K, ier)

            CALL ADJUST_LEQ (res1, LEQ_IT(9), LEQ_METHOD(9), &
                             LEQI, LEQM)

            write(Vname, '(A,I2)')'Scalar',NN
            CALL SOLVE_LIN_EQ (Vname, 9, Scalar(1,NN), A_M, B_M, M, &
                               LEQI, LEQM, LEQ_SWEEP(9), LEQ_TOL(9), &
                               LEQ_PC(9), IER)
!            call out_array(Scalar(1, N), Vname)

         ENDIF
      ENDDO

      call unlock_ambm

      RETURN
      END SUBROUTINE SOLVE_Scalar_EQ
