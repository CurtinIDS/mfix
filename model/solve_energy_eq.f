!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLVE_ENERGY_EQ                                         C
!  Purpose: Solve energy equations                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-97  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLVE_ENERGY_EQ(IER)

! Modules
!---------------------------------------------------------------------//
      use ambm, only: a_m, b_m, lock_ambm, unlock_ambm
      use bc, only: bc_t_g, bc_tw_g, bc_hw_t_g, bc_c_t_g
      use bc, only: bc_t_s, bc_tw_s, bc_hw_t_s, bc_c_t_s
      use compar, only: ijkstart3, ijkend3, mype, numpes
      use discretelement, only: des_continuum_coupled
      use energy, only: hor_s, s_rps, s_rcs, gama_gs
      use energy, only: hor_g, s_rpg, s_rcg
      use fldvar, only: ep_g, u_g, v_g, w_g, t_g, t_go, rop_go
      use fldvar, only: ep_s, u_s, v_s, w_s, t_s, t_so, rop_so
      use functions, only: fluid_at, wall_at
      use functions, only: is_on_mype_plus2layers
      use functions, only: ip_of, jp_of, kp_of
      use geometry, only: ijkmax2, vol
      use indices, only: i_of, j_of, k_of
      use indices, only: ip1, jp1, kp1
      use leqsol, only: leq_it, leq_method, leq_sweep, leq_pc, leq_tol
      use mflux, only: flux_ge, flux_gse, flux_gn, flux_gsn
      use mflux, only: flux_se, flux_sse, flux_sn, flux_ssn
      use mflux, only: flux_gt, flux_gst, flux_st, flux_sst
      use mms, only: use_mms, mms_t_g_src, mms_t_s_src
      use mpi_utility, only: global_all_sum
      use param, only: dimension_3, dimension_m
      use param1, only: zero, half
      use physprop, only: k_g, c_pg
      use physprop, only: k_s, c_ps, smax
      use ps, only: point_source
      use ps, only: ps_t_g, ps_cpxmflow_g
      use ps, only: ps_t_s, ps_cpxmflow_s
      use residual, only: resid, max_resid, ijk_resid
      use residual, only: num_resid, den_resid
      use residual, only: resid_t
      use run, only: discretize, odt, added_mass, m_am
      use toleranc, only: tmin, tmax
      use ur_facs, only: ur_fac
      use usr_src, only: call_usr_source, calc_usr_source
      use usr_src, only: gas_energy, solids_energy
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local variables
!---------------------------------------------------------------------//
! phase index
      INTEGER :: M
!  Cp * Flux
      DOUBLE PRECISION :: CpxFlux_E(DIMENSION_3), &
                          CpxFlux_N(DIMENSION_3), &
                          CpxFlux_T(DIMENSION_3)
! previous time step term
      DOUBLE PRECISION :: apo
! Indices
      INTEGER :: IJK, I, J, K
! linear equation solver method and iterations
      INTEGER :: LEQM, LEQI

! Arrays for storing errors:
! 120 - Gas phase energy equation diverged
! 121 - Solids energy equation diverged
! 12x - Unclassified
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global

! temporary use of global arrays:
! arraym1 (locally vxgama)
! the volume x average gas-solids heat transfer at cell centers
!      DOUBLE PRECISION :: VXGAMA(DIMENSION_3, DIMENSION_M)
! array1 (locally s_p)
! source vector: coefficient of dependent variable
! becomes part of a_m matrix; must be positive
      DOUBLE PRECISION :: S_P(DIMENSION_3)
! array2 (locally s_c)
! source vector: constant part becomes part of b_m vector
      DOUBLE PRECISION :: S_C(DIMENSION_3)
! array3 (locally eps)
      DOUBLE PRECISION :: EPS(DIMENSION_3)

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: VXGAMA

! Septadiagonal matrix A_m, vector b_m
!      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!      DOUBLE PRECISION B_m(DIMENSION_3, 0:DIMENSION_M)
!---------------------------------------------------------------------//

      ALLOCATE(VXGAMA(DIMENSION_3,DIMENSION_M))

      call lock_ambm         ! locks arrys a_m and b_m
                             ! (locally s_p, s_c, eps)
! Initialize error flags.
      Err_l = 0

! initializing
      DO M = 0, SMAX
         CALL INIT_AB_M (A_M, B_M, IJKMAX2, M)
      ENDDO

! Gas phase computations
! ---------------------------------------------------------------->>>
      DO IJK = IJKSTART3, IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF (IS_ON_myPE_plus2layers(IP1(I),J,K)) THEN
            IF(.NOT.ADDED_MASS) THEN
               CpxFlux_E(IJK) = HALF * (C_pg(IJK) + C_pg(IP_OF(IJK))) * Flux_gE(IJK)
            ELSE
               CpxFlux_E(IJK) = HALF * (C_pg(IJK) + C_pg(IP_OF(IJK))) * Flux_gSE(IJK)
            ENDIF
         ENDIF

         IF (IS_ON_myPE_plus2layers(I,JP1(J),K)) THEN
            IF(.NOT.ADDED_MASS) THEN
               CpxFlux_N(IJK) = HALF * (C_pg(IJK) + C_pg(JP_OF(IJK))) * Flux_gN(IJK)
            ELSE
               CpxFlux_N(IJK) = HALF * (C_pg(IJK) + C_pg(JP_OF(IJK))) * Flux_gSN(IJK)
            ENDIF
         ENDIF

         IF (IS_ON_myPE_plus2layers(I,J,KP1(K))) THEN
            IF(.NOT.ADDED_MASS) THEN
               CpxFlux_T(IJK) = HALF * (C_pg(IJK) + C_pg(KP_OF(IJK))) * Flux_gT(IJK)
            ELSE
               CpxFlux_T(IJK) = HALF * (C_pg(IJK) + C_pg(KP_OF(IJK))) * Flux_gST(IJK)
            ENDIF
         ENDIF

         IF (FLUID_AT(IJK)) THEN
            APO = ROP_GO(IJK)*C_PG(IJK)*VOL(IJK)*ODT
            S_P(IJK) = APO + S_RPG(IJK)*VOL(IJK)
            S_C(IJK) = APO*T_GO(IJK)-HOR_G(IJK)*VOL(IJK)+S_RCG(IJK)*VOL(IJK)
            IF(USE_MMS) S_C(IJK) = S_C(IJK) + MMS_T_G_SRC(IJK)*VOL(IJK)
         ELSE
            S_P(IJK) = ZERO
            S_C(IJK) = ZERO
         ENDIF
      ENDDO

! Account for convective heat transfer between gas and DES particles.
      IF(DES_CONTINUUM_COUPLED) CALL DES_2FLUID_CONV(S_P, S_C)

! calculate the convection-diffusion terms
      CALL CONV_DIF_PHI (T_g, K_G, DISCRETIZE(6), U_G, V_G, W_G, &
         CpxFlux_E, CpxFlux_N, CpxFlux_T, 0, A_M, B_M)

! calculate standard bc
      CALL BC_PHI (T_g, BC_T_G, BC_TW_G, BC_HW_T_G, BC_C_T_G, 0, A_M, B_M)

! set the source terms in a and b matrix equation form
      CALL SOURCE_PHI (S_P, S_C, EP_G, T_G, 0, A_M, B_M)

! add point sources
      IF(POINT_SOURCE) CALL POINT_SOURCE_PHI (T_g, PS_T_g, &
         PS_CpxMFLOW_g, 0, A_M, B_M)
! usr source
      IF (CALL_USR_SOURCE(6)) CALL CALC_USR_SOURCE(GAS_ENERGY, &
                            A_M, B_M, lM=0)

! Solids phases computations
! ---------------------------------------------------------------->>>
      DO M = 1, SMAX
         DO IJK = IJKSTART3, IJKEND3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)

            IF (IS_ON_myPE_plus2layers(IP1(I),J,K)) THEN
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_E(IJK) = HALF * (C_ps(IJK,M) + C_ps(IP_OF(IJK),M)) * Flux_sE(IJK,M)
               ELSE   ! M=M_AM is the only phase for which virtual mass is added
                  CpxFlux_E(IJK) = HALF * (C_ps(IJK,M) + C_ps(IP_OF(IJK),M)) * Flux_sSE(IJK)
               ENDIF
            ENDIF

            IF (IS_ON_myPE_plus2layers(I,JP1(J),K)) THEN
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_N(IJK) = HALF * (C_ps(IJK,M) + C_ps(JP_OF(IJK),M)) * Flux_sN(IJK,M)
               ELSE
                  CpxFlux_N(IJK) = HALF * (C_ps(IJK,M) + C_ps(JP_OF(IJK),M)) * Flux_sSN(IJK)
               ENDIF
            ENDIF

            IF (IS_ON_myPE_plus2layers(I,J,KP1(K))) THEN
               IF(.NOT.ADDED_MASS .OR. M /= M_AM) THEN
                  CpxFlux_T(IJK) = HALF * (C_ps(IJK,M) + C_ps(KP_OF(IJK),M)) * Flux_sT(IJK,M)
               ELSE
                  CpxFlux_T(IJK) = HALF * (C_ps(IJK,M) + C_ps(KP_OF(IJK),M)) * Flux_sST(IJK)
               ENDIF
            ENDIF

            IF (FLUID_AT(IJK)) THEN
               APO = ROP_SO(IJK,M)*C_PS(IJK,M)*VOL(IJK)*ODT
               S_P(IJK) = APO + S_RPS(IJK,M)*VOL(IJK)
               S_C(IJK) = APO*T_SO(IJK,M) - HOR_S(IJK,M)*VOL(IJK) + &
                  S_RCS(IJK,M)*VOL(IJK)
               VXGAMA(IJK,M) = GAMA_GS(IJK,M)*VOL(IJK)
               EPS(IJK) = EP_S(IJK,M)
               IF(USE_MMS) S_C(IJK) = S_C(IJK) + MMS_T_S_SRC(IJK)*VOL(IJK)
            ELSE
               S_P(IJK) = ZERO
               S_C(IJK) = ZERO
               VXGAMA(IJK,M) = ZERO
               EPS(IJK) = ZERO
               IF(USE_MMS) EPS(IJK) = EP_S(IJK,M)
            ENDIF
         ENDDO   ! end do (ijk=ijkstart3,ijkend3)

! calculate the convection-diffusion terms
         CALL CONV_DIF_PHI (T_s(1,M), K_S(1,M), DISCRETIZE(6), &
            U_S(1,M), V_S(1,M), W_S(1,M), CpxFlux_E, CpxFlux_N, &
            CpxFlux_T, M, A_M, B_M)

! calculate standard bc
         CALL BC_PHI (T_s(1,M), BC_T_S(1,M), BC_TW_S(1,M), &
            BC_HW_T_S(1,M), BC_C_T_S(1,M), M, A_M, B_M)

! set the source terms in a and b matrix equation form
         CALL SOURCE_PHI (S_P, S_C, EPS, T_S(1,M), M, A_M, B_M)

! Add point sources.
         IF(POINT_SOURCE) CALL POINT_SOURCE_PHI (T_s(:,M), PS_T_s(:,M),&
            PS_CpxMFLOW_s(:,M), M, A_M, B_M)

! usr source
         IF (CALL_USR_SOURCE(6)) CALL CALC_USR_SOURCE(SOLIDS_ENERGY, &
                               A_M, B_M, lM=M)

      ENDDO   ! end do (m=1,smax)

! Solve gas and solids phase equations
! ---------------------------------------------------------------->>>
! use partial elimination on interphase heat transfer term
      IF (SMAX > 0 .AND. .NOT.USE_MMS) &
        CALL PARTIAL_ELIM_S (T_G, T_S, VXGAMA, A_M, B_M)

      CALL CALC_RESID_S (T_G, A_M, B_M, 0, NUM_RESID(RESID_T,0),&
         DEN_RESID(RESID_T,0), RESID(RESID_T,0), MAX_RESID(RESID_T,&
         0), IJK_RESID(RESID_T,0), ZERO)

      CALL UNDER_RELAX_S (T_G, A_M, B_M, 0, UR_FAC(6))

!      call check_ab_m(a_m, b_m, 0, .false., ier)
!      call write_ab_m(a_m, b_m, ijkmax2, 0, ier)
!      write(*,*) &
!         resid(resid_t, 0), max_resid(resid_t, 0), &
!         ijk_resid(resid_t, 0)


      DO M = 1, SMAX
         CALL CALC_RESID_S (T_S(1,M), A_M, B_M, M, NUM_RESID(RESID_T,M), &
            DEN_RESID(RESID_T,M), RESID(RESID_T,M), MAX_RESID(&
            RESID_T,M), IJK_RESID(RESID_T,M), ZERO)

         CALL UNDER_RELAX_S (T_S(1,M), A_M, B_M, M, UR_FAC(6))
      ENDDO

! set/adjust linear equation solver method and iterations
      CALL ADJUST_LEQ(RESID(RESID_T,0), LEQ_IT(6), LEQ_METHOD(6), &
         LEQI, LEQM)
!      call test_lin_eq(a_m(1, -3, 0), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6),  0, ier)

      CALL SOLVE_LIN_EQ ('T_g', 6, T_G, A_M, B_M, 0, LEQI, LEQM, &
         LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6), IER)
! Check for linear solver divergence.
      IF(ier == -2) Err_l(myPE) = 120

! bound temperature in any fluid or flow boundary cells
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.WALL_AT(IJK))&
            T_g(IJK) = MIN(TMAX, MAX(TMIN, T_g(IJK)))
      ENDDO
!      call out_array(T_g, 'T_g')

      DO M = 1, SMAX
         CALL ADJUST_LEQ (RESID(RESID_T,M), LEQ_IT(6), LEQ_METHOD(6), &
            LEQI, LEQM)
!         call test_lin_eq(a_m(1, -3, M), LEQI, LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6),  0, ier)
         CALL SOLVE_LIN_EQ ('T_s', 6, T_S(1,M), A_M, B_M, M, LEQI, &
            LEQM, LEQ_SWEEP(6), LEQ_TOL(6), LEQ_PC(6), IER)
! Check for linear solver divergence.
         IF(ier == -2) Err_l(myPE) = 121

! bound temperature in any fluid or flow boundary cells
         DO IJK = IJKSTART3, IJKEND3
            IF(.NOT.WALL_AT(IJK))&
               T_s(IJK, M) = MIN(TMAX, MAX(TMIN, T_s(IJK, M)))
         ENDDO
      ENDDO   ! end do (m=1, smax)

      call unlock_ambm

! If the linear solver diverged, temperatures may take on unphysical
! values. To prevent them from propogating through the domain or
! causing failure in other routines, force an exit from iterate and
! reduce the time step.
      CALL global_all_sum(Err_l, Err_g)
      IER = maxval(Err_g)

      RETURN
      END SUBROUTINE SOLVE_ENERGY_EQ
