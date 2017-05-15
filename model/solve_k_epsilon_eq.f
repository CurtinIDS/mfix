!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_K_Epsilon_EQ                                      C
!  Purpose: Solve K & Epsilon equations for a turbulent flow           C
!                                                                      C
!                                                                      C
!  Author: S. Benyahia                                Date: MAY-13-04  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_K_Epsilon_EQ(IER)

! Modules
!---------------------------------------------------------------------//
      use ambm, only: a_m, b_m, lock_ambm, unlock_ambm
      use bc, only: bc_k_turb_g, bc_e_turb_g
      use compar, only: ijkstart3, ijkend3
      use constant, only: to_si
      use fldvar, only: rop_g, ep_g, u_g, v_g, w_g
      use fldvar, only: k_turb_g, k_turb_go
      use fldvar, only: e_turb_g, e_turb_go
      use functions, only: fluid_at, zmax
      use indices, only: i_of, j_of, k_of
      use geometry, only: vol, ijkmax2
      use leqsol, only: leq_it, leq_method, leq_sweep, leq_pc, leq_tol
      use mflux, only: flux_ge, flux_gn, flux_gt
      use mflux, only: flux_gse, flux_gsn, flux_gst
      use param, only: dimension_bc, dimension_3
      use param1, only: zero
      use residual, only: resid, max_resid, ijk_resid
      use residual, only: num_resid, den_resid
      use residual, only: resid_ke
      use run, only: discretize, k_epsilon, added_mass, odt
      use rxns, only: sum_r_g
      use toleranc, only: zero_ep_s
      use turb, only: dif_k_turb_g, k_turb_g_c, k_turb_g_p
      use turb, only: dif_e_turb_g, e_turb_g_c, e_turb_g_p
      use ur_facs, only: ur_fac
      use usr_src, only: call_usr_source, calc_usr_source
      use usr_src, only: k_epsilon_e, k_epsilon_k
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local variables
!---------------------------------------------------------------------//
! phase index
      INTEGER :: M
! indices
      INTEGER :: IJK, I, J, K
! loop counter
      INTEGER :: LC
!
      DOUBLE PRECISION :: apo
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI
! temporary variables in residual computation
      DOUBLE PRECISION :: res1, mres1, num_res, den_res
      INTEGER :: ires1
! A default zero flux will be defined for both K & Epsilon at walls
      DOUBLE PRECISION :: BC_hw_K_Turb_G (DIMENSION_BC)
      DOUBLE PRECISION :: BC_hw_E_Turb_G (DIMENSION_BC)
      DOUBLE PRECISION :: BC_K_Turb_GW (DIMENSION_BC)
      DOUBLE PRECISION :: BC_E_Turb_GW (DIMENSION_BC)
      DOUBLE PRECISION :: BC_C_K_Turb_G (DIMENSION_BC)
      DOUBLE PRECISION :: BC_C_E_Turb_G (DIMENSION_BC)

! small value of K or E, 1 cm2/s2 = 1e-4 m2/s2 = 1e-4 m2/s3
      DOUBLE PRECISION smallTheta

! source vector: coefficient of dependent variable
! becomes part of a_m matrix; must be positive
      DOUBLE PRECISION :: S_P(DIMENSION_3)
! source vector: constant part becomes part of b_m vector
      DOUBLE PRECISION :: S_C(DIMENSION_3)
!
      character(LEN=8) :: Vname
!---------------------------------------------------------------------//

      IF( .NOT. K_Epsilon) RETURN

      call lock_ambm

      smallTheta = (to_SI)**4 * ZERO_EP_S
      RESID(RESID_ke,0) = ZERO
      NUM_RESID(RESID_ke,0) = ZERO
      DEN_RESID(RESID_ke,0) = ZERO
      MAX_RESID(RESID_ke,0) = ZERO
      IJK_RESID(RESID_ke,0) = 0

! Setting default zero flux for K & Epsilon since we use wall functions.
! If an expert user want to use Low Re K-Epilon model and needs to set
! the turbulence quantities to zero at walls, then set the hw's to UNDEFINE will
! do it. All the variables below can be changed in the same way as in the
! MFIX data file in the boundary conditions section.
      DO LC = 1, DIMENSION_BC
        BC_hw_K_Turb_G (LC) = ZERO
        BC_K_Turb_GW (LC) = ZERO
        BC_C_K_Turb_G (LC) = ZERO
        BC_hw_E_Turb_G (LC) = ZERO
        BC_E_Turb_GW (LC) = ZERO
        BC_C_E_Turb_G (LC) = ZERO
      ENDDO
! End of setting default zero flux for K & Epsilon wall boundary conditions

! Equations solved for gas phase, thus M = 0
      M = 0
      CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)

! Solve the K_Turb_G equation first
! ---------------------------------------------------------------->>>
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN
            APO = ROP_G(IJK)*VOL(IJK)*ODT
            S_P(IJK) = APO + (ZMAX(SUM_R_G(IJK)) + K_Turb_G_p(IJK))*&
                       VOL(IJK)
            S_C(IJK) = APO*K_Turb_GO(IJK) + &
                       K_Turb_G(IJK)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) + &
                       K_Turb_G_c(IJK)*VOL(IJK)
         ELSE
            S_P(IJK) = ZERO
            S_C(IJK) = ZERO
         ENDIF
      ENDDO

      IF(.NOT.ADDED_MASS) THEN
         CALL CONV_DIF_PHI (K_Turb_G, DIF_K_Turb_G, DISCRETIZE(9), &
                            U_G, V_G, W_G, Flux_gE, Flux_gN, Flux_gT, &
                            M, A_M, B_M)
      ELSE
         CALL CONV_DIF_PHI (K_Turb_G, DIF_K_Turb_G, DISCRETIZE(9), &
                            U_G, V_G, W_G, Flux_gSE, Flux_gSN, Flux_gST, &
                            M, A_M, B_M)
      ENDIF

      CALL BC_PHI (K_Turb_G, BC_K_Turb_G, BC_K_Turb_GW, BC_HW_K_Turb_G, &
                   BC_C_K_Turb_G, M, A_M, B_M)
      CALL SOURCE_PHI (S_P, S_C, EP_G, K_Turb_G, M, A_M, B_M)

! calculate any usr source terms
      IF (CALL_USR_SOURCE(9)) CALL CALC_USR_SOURCE(K_EPSILON_K,&
                            A_M, B_M, lM=0)

      CALL CALC_RESID_S (K_Turb_G, A_M, B_M, M, num_res, den_res, res1, &
                         mres1, ires1, ZERO)

      RESID(RESID_ke,0) = RESID(RESID_ke,0)+res1
      NUM_RESID(RESID_ke,0) = NUM_RESID(RESID_ke,0)+num_res
      DEN_RESID(RESID_ke,0) = DEN_RESID(RESID_ke,0)+den_res
      if(mres1 .gt. MAX_RESID(RESID_ke,0)) then
         MAX_RESID(RESID_ke,0) = mres1
         IJK_RESID(RESID_ke,0) = ires1
      endif

      CALL UNDER_RELAX_S (K_Turb_G, A_M, B_M, M, UR_FAC(9))

!      call check_ab_m(a_m, b_m, m, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!      write(*,*) res1, mres1, ires1
!      call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, &
!                       DO_K, ier)

      CALL ADJUST_LEQ (RESID(RESID_ke,0), LEQ_IT(9), LEQ_METHOD(9), &
               LEQI, LEQM)

      write(Vname, '(A,I2)')'K_Turb_G'
      CALL SOLVE_LIN_EQ (Vname, 9, K_Turb_G, A_M, B_M, M, LEQI, LEQM, &
                         LEQ_SWEEP(9), LEQ_TOL(9),  LEQ_PC(9), IER)

!      call out_array(K_Turb_G, Vname)

! remove small negative K values generated by linear solver
! same as adjust_theta.f
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN
            IF(K_Turb_G(IJK) < smallTheta) K_Turb_G(IJK) = smallTheta
         ENDIF
      ENDDO


! Now solve the E_Turb_G (dissipation) equation.
! ---------------------------------------------------------------->>>
! Initiate (again) the Am Bm matrix.
! This has to be done for every scalar equation.
      CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)

      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IF (FLUID_AT(IJK)) THEN
            APO = ROP_G(IJK)*VOL(IJK)*ODT
            S_P(IJK) = APO + (ZMAX(SUM_R_G(IJK)) + E_Turb_G_p(IJK))*&
                       VOL(IJK)
            S_C(IJK) = APO*E_Turb_GO(IJK) + &
                       E_Turb_G(IJK)*ZMAX((-SUM_R_G(IJK)))*VOL(IJK) + &
                       E_Turb_G_c(IJK)*VOL(IJK)
         ELSE
            S_P(IJK) = ZERO
            S_C(IJK) = ZERO
         ENDIF
      ENDDO

      IF(.NOT.ADDED_MASS) THEN
         CALL CONV_DIF_PHI (E_Turb_G, DIF_E_Turb_G, DISCRETIZE(9), &
                            U_G, V_G, W_G, Flux_gE, Flux_gN, Flux_gT, &
                            M, A_M, B_M)
      ELSE
         CALL CONV_DIF_PHI (E_Turb_G, DIF_E_Turb_G, DISCRETIZE(9), &
                            U_G, V_G, W_G, Flux_gSE, Flux_gSN, Flux_gST, &
                            M, A_M, B_M)
      ENDIF

      CALL BC_PHI (E_Turb_G, BC_E_Turb_G, BC_E_Turb_GW, BC_HW_E_Turb_G, &
                   BC_C_E_Turb_G, M, A_M, B_M)

      CALL SOURCE_PHI (S_P, S_C, EP_G, E_Turb_G, M, A_M, B_M)

! calculate any usr source terms
      IF (CALL_USR_SOURCE(9)) CALL CALC_USR_SOURCE(K_EPSILON_E,&
                            A_M, B_M, lM=0)

! set epsilon in fluid cells adjacent to wall cells
      CALL SOURCE_K_EPSILON_BC(A_M, B_M, M)

      CALL CALC_RESID_S (E_Turb_G, A_M, B_M, M, num_res, den_res, res1, &
                         mres1, ires1, ZERO)

      RESID(RESID_ke,0) = RESID(RESID_ke,0)+res1
      NUM_RESID(RESID_ke,0) = NUM_RESID(RESID_ke,0)+num_res
      DEN_RESID(RESID_ke,0) = DEN_RESID(RESID_ke,0)+den_res
      if(mres1 .gt. MAX_RESID(RESID_ke,0))then
        MAX_RESID(RESID_ke,0) = mres1
        IJK_RESID(RESID_ke,0) = ires1
      endif

      CALL UNDER_RELAX_S (E_Turb_G, A_M, B_M, M, UR_FAC(9))

!      call check_ab_m(a_m, b_m, m, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, m, ier)
!      write(*,*) res1, mres1, ires1
!      call test_lin_eq(ijkmax2, ijmax2, imax2, a_m(1, -3, 0), 1, &
!                       DO_K, ier)

       CALL ADJUST_LEQ (RESID(RESID_ke,0), LEQ_IT(9), LEQ_METHOD(9), &
                        LEQI, LEQM)

      write(Vname, '(A,I2)')'E_Turb_G'
      CALL SOLVE_LIN_EQ (Vname, 9, E_Turb_G, A_M, B_M, M, LEQI, LEQM, &
                         LEQ_SWEEP(9), LEQ_TOL(9), LEQ_PC(9), IER)

!      call out_array(E_Turb_G, Vname)

! remove small negative Epsilon values generated by linear solver
! same as adjust_theta.f
      DO IJK = IJKSTART3, IJKEND3
         IF (FLUID_AT(IJK)) THEN
            IF(E_Turb_G(IJK) < smallTheta) E_Turb_G(IJK) = smallTheta
         ENDIF
      ENDDO

      call unlock_ambm

      RETURN
      END SUBROUTINE SOLVE_K_Epsilon_EQ

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_K_EPSILON_BC                                     C
!  Purpose: When implementing the wall functions, the epsilon          C
!  (dissipation) value at the fluid cell near the walls needs to be    C
!  set.                                                                C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOURCE_K_EPSILON_BC(A_M, B_M, M)

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      use cutcell, only: cut_cell_at, delh_scalar
      use fldvar, only: k_turb_G, e_turb_g
      use functions, only: fluid_at, wall_at
      use functions, only: im_of, jm_of, km_of
      use functions, only: ip_of, jp_of, kp_of
      use geometry, only: dy, odz, ox, dx, cylindrical
      use indices, only: i_of, j_of, k_of
      use param, only: dimension_3, dimension_m
      use param1, only: zero, one
      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! phase index
      INTEGER, INTENT(IN) :: M
! Local variables
!---------------------------------------------------------------------//
      INTEGER :: IJK, I, J, K
!---------------------------------------------------------------------//

      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IF (FLUID_AT(IJK)) THEN
            IF(WALL_AT(JP_OF(IJK)).OR.WALL_AT(JM_OF(IJK))) THEN
               A_M(IJK,1,M) = ZERO
               A_M(IJK,-1,M) = ZERO
               A_M(IJK,2,M) = ZERO
               A_M(IJK,-2,M) = ZERO
               A_M(IJK,3,M) = ZERO
               A_M(IJK,-3,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/DY(J) &
                           *2.0D+0/0.42D+0
            ELSEIF(WALL_AT(KP_OF(IJK)).OR.WALL_AT(KM_OF(IJK))) THEN
               A_M(IJK,1,M) = ZERO
               A_M(IJK,-1,M) = ZERO
               A_M(IJK,2,M) = ZERO
               A_M(IJK,-2,M) = ZERO
               A_M(IJK,3,M) = ZERO
               A_M(IJK,-3,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)* &
                            (ODZ(K)*OX(I)*2.0D+0)/0.42D+0
            ENDIF  !for identifying wall cells in J or K direction

            IF(CYLINDRICAL) THEN
               IF (WALL_AT(IP_OF(IJK)))  THEN
                  A_M(IJK,1,M) = ZERO
                  A_M(IJK,-1,M) = ZERO
                  A_M(IJK,2,M) = ZERO
                  A_M(IJK,-2,M) = ZERO
                  A_M(IJK,3,M) = ZERO
                  A_M(IJK,-3,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/&
                               DX(I)*2.0D+0/0.42D+0

               ENDIF! for wall cells in I direction
            ELSEIF (WALL_AT(IP_OF(IJK)).OR.WALL_AT(IM_OF(IJK))) THEN
               A_M(IJK,1,M) = ZERO
               A_M(IJK,-1,M) = ZERO
               A_M(IJK,2,M) = ZERO
               A_M(IJK,-2,M) = ZERO
               A_M(IJK,3,M) = ZERO
               A_M(IJK,-3,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/DX(I) &
                             *2.0D+0/0.42D+0
            ENDIF ! for cylindrical

            IF(CUT_CELL_AT(IJK)) THEN
               A_M(IJK,1,M) = ZERO
               A_M(IJK,-1,M) = ZERO
               A_M(IJK,2,M) = ZERO
               A_M(IJK,-2,M) = ZERO
               A_M(IJK,3,M) = ZERO
               A_M(IJK,-3,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) =-((0.09D+0)**0.75*K_Turb_G(IJK)**1.5)/&
                            (0.42D+0*DELH_Scalar(IJK))
            ENDIF
         ENDIF  !for fluid at ijk
      ENDDO   ! enddo ijk
      RETURN
      END SUBROUTINE SOURCE_K_EPSILON_BC
