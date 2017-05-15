!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: set_constprop                                           C
!  Purpose: This routine serves two purposes:                          C
!    1) initializes various variables everywhere in the domain with    C
!       a zero value. the real need for this is unclear. undefined     C
!       may be a better approach...                                    C
!    2) if defined, sets physical properties to their specified        C
!       constant value in the fluid domain. cannot set in flow         C
!       boundaries or later checks will also complain (may be          C
!       overly strict check)                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_CONSTPROP

! Modules
!-----------------------------------------------
      USE param1, only: zero, half, one, undefined

      use discretelement, only: des_mmax

      USE fldvar, only: ro_g
      USE fldvar, only: ro_s, d_p
      use fldvar, only: p_s

      USE visc_g, only: mu_gt, epmu_gt, lambda_gt, eplambda_gt
      USE visc_s, only: mu_s, epmu_s, lambda_s, eplambda_s, lambda_s_c
      USE visc_s, only: ep_star_array
      USE visc_s, only: ep_g_blend_start, ep_g_blend_end

      USE physprop, only: nmax, mmax, smax
      USE physprop, only: ro_s0, ro_g0, d_p0
      USE physprop, only: mu_g0, mu_g, mu_s0
      USE physprop, only: mw_avg, mw_mix_g
      USE physprop, only: c_pg0, k_g0, c_pg, k_g, dif_g0, dif_g
      USE physprop, only: c_ps0, k_s0, c_ps, k_s, dif_s0, dif_s
      USE physprop, only: cv

      USE constant, only: ep_s_max_ratio, d_p_ratio, ep_s_max, m_max
      use constant, only: ep_star

      USE run, only: call_dqmom
      USE run, only: yu_standish, fedors_landel
      USE run, only: kt_type_enum, ia_2005, gd_1999, gtsh_2012
      USE run, only: blending_stress, sigm_blend, tanh_blend

      USE drag, only: f_gs, f_ss

      use kintheory, only: mu_sm_ip, mu_sl_ip, xi_sm_ip, xi_sl_ip
      use kintheory, only: fnu_s_ip, ft_sm_ip, ft_sl_ip
      use kintheory, only: kth_sl_ip, knu_sm_ip, knu_sl_ip, kvel_s_ip
      use kintheory, only: ed_ss_ip, edvel_sl_ip, edt_s_ip, edvel_sm_ip
      use kintheory, only: a2_gtsh, xsi_gtsh

      use mms, only: use_mms

      USE compar, only: ijkstart3, ijkend3
      use functions, only: wall_at, fluid_at

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! indices
      INTEGER :: IJK, M, I, J
      DOUBLE PRECISION :: old_value, DP_TMP(SMAX)
!-----------------------------------------------

! First, initialize certain transport coefficients, physical
! properties, and other variable types everywhere in the
! domain to zero. Then, set these to their specified value
! (if defined) in fluid cells. Some are also set in flow cells.
! NOTE: DO NOT simply zero existing field variables.
      RO_g = ZERO
      MU_g = ZERO
      MU_gt = ZERO
      EPMU_GT = ZERO
      LAMBDA_GT = ZERO
      EPLAMBDA_GT = ZERO
      K_g = ZERO
      C_PG = ZERO
      MW_MIX_G = zero
      DIF_g = ZERO
      F_GS = ZERO

      RO_S = ZERO
      D_P = zero     ! this is done in init_fvars
      MU_s = ZERO
      EPMU_s = ZERO
      LAMBDA_s_c = ZERO
      LAMBDA_s = ZERO
      EPLAMBDA_s = ZERO
! not defining p_s in flow cells can become a problem for certain
! cases as gradients in P_s are used (e.g., when a PO is west of a
! MI...?). Simply assigning it a zero value may not be the best
! approach but it is currently used....
      P_S = ZERO
      K_s = ZERO
      C_PS = ZERO
      DIF_S = ZERO
      F_SS = ZERO

! Set the flag for recalculating gas viscosity.
!      RECALC_VISC_G = (MU_g0==UNDEFINED .OR. L_SCALE0/=ZERO .OR.&
!                       K_EPSILON .OR. ISHII)

! Set default value for virtual mass coefficient
      Cv = HALF

! Variables specific to various kinetic theory models
      IF (KT_TYPE_ENUM == IA_2005) THEN
         MU_sM_ip = ZERO
         MU_sL_ip = ZERO
         XI_sM_ip = ZERO
         XI_sL_ip = ZERO
         Fnu_s_ip = ZERO
         FT_sM_ip = ZERO
         FT_sL_ip = ZERO
         Kth_sL_ip = ZERO
         Knu_sM_ip = ZERO
         Knu_sL_ip = ZERO
         Kvel_s_ip = ZERO
         ED_ss_ip = ZERO
         EDvel_sL_ip = ZERO
      ENDIF

      IF (KT_TYPE_ENUM == IA_2005 .OR. &
          KT_TYPE_ENUM == GD_1999 .OR.  &
          KT_TYPE_ENUM == GTSH_2012) THEN
         EDT_s_ip = ZERO
         EDvel_sM_ip = ZERO
      ENDIF

      IF(KT_TYPE_ENUM == GTSH_2012) THEN
         A2_gtsh = ZERO
         xsi_gtsh = zero
      ENDIF

! Set specified constant physical properties values
      DO IJK = ijkstart3, ijkend3

         IF (.NOT.WALL_AT(IJK)) THEN
! Fluid and inflow/outflow cells: FLAG < 100
            IF (RO_G0 /= UNDEFINED) RO_G(IJK) = RO_G0
            IF (C_PG0 /= UNDEFINED) C_PG(IJK) = C_PG0
            IF (MW_AVG /= UNDEFINED) MW_MIX_G(IJK) = MW_AVG
         ENDIF

         IF (FLUID_AT(IJK)) THEN
! Strictly Fluid cells: FLAG = 1
            IF (MU_G0 /= UNDEFINED) THEN
               MU_G(IJK) = MU_G0
               MU_GT(IJK) = MU_G0
               EPMU_GT(IJK) = MU_G0
               LAMBDA_GT(IJK) = -(2.0d0/3.0d0)*MU_G0
               EPLAMBDA_GT(IJK) = -(2.0d0/3.0d0)*MU_G0
            ENDIF
            IF (K_G0 /= UNDEFINED) K_G(IJK) = K_G0
            IF (DIF_G0 /= UNDEFINED) DIF_G(IJK,:NMAX(0)) = DIF_G0
         ENDIF

         IF (USE_MMS) THEN
            IF (MU_G0 /= UNDEFINED) THEN
               MU_G(IJK) = MU_G0
               MU_GT(IJK) = MU_G0
               EPMU_GT(IJK) = MU_G0
               LAMBDA_GT(IJK) = -(2.0d0/3.0d0)*MU_G0
               EPLAMBDA_GT(IJK) = -(2.0d0/3.0d0)*MU_G0
            ENDIF
            IF (K_G0 /= UNDEFINED) K_G(IJK) = K_G0
            IF (DIF_G0 /= UNDEFINED) DIF_G(IJK,:NMAX(0)) = DIF_G0
         ENDIF

      ENDDO

! ensure ro_s(ijk,m) is assigned to ro_s0(m) or ro_s(np) so that the
! function ep_s works for discrete phases. might be able to move this
! to set_ic_dem and set_ic_mppic. also required d_p(ijk,m) for hybrid
! use.  note check_solids_common_all ensures d_p0 is set for all
! solids defined also either ro_s0 must be set or ro_s0
      DO M = 1, MMAX+DES_MMAX
         DO IJK = ijkstart3, ijkend3
            IF(.NOT.WALL_AT(IJK)) THEN
! Fluid and inflow/outflow cells: FLAG < 100
               IF (RO_S0(M) /= UNDEFINED) RO_S(IJK,M) = RO_S0(M)
               IF (C_PS0(M) /= UNDEFINED) C_PS(IJK,M) = C_PS0(M)
               IF (D_P0(M) /= UNDEFINED) D_P(IJK,M) = D_P0(M)
            ENDIF

            IF (FLUID_AT(IJK)) THEN
! Strictly fluid cells: FLAG = 1
               IF (MU_S0(M) /= UNDEFINED) THEN
                  MU_S(IJK,M) = MU_S0(M)
                  EPMU_S(IJK,M) = MU_S0(M)
                  LAMBDA_S(IJK,M) = (-2./3.)*MU_S(IJK,M)
                  EPLAMBDA_S(IJK,M) = (-2./3.)*MU_S(IJK,M)
               ENDIF
               IF (K_S0(M) /= UNDEFINED) K_S(IJK,M) = K_S0(M)
               IF (DIF_S0(M) /= UNDEFINED) DIF_S(IJK,M,:NMAX(M)) = DIF_S0(M)
            ENDIF

            IF (USE_MMS) THEN
               IF (MU_S0(M) /= UNDEFINED) THEN
                  MU_S(IJK,M) = MU_S0(M)
                  EPMU_S(IJK,M) = MU_S0(M)
                  LAMBDA_S(IJK,M) = (-2./3.)*MU_S(IJK,M)
                  EPLAMBDA_S(IJK,M) = (-2./3.)*MU_S(IJK,M)
               ENDIF
               IF (K_S0(M) /= UNDEFINED) K_S(IJK,M) = K_S0(M)
               IF (DIF_S0(M) /= UNDEFINED) DIF_S(IJK,M,:NMAX(M)) = DIF_S0(M)
            ENDIF

! set ep_star_array to user input ep_star in all cells.
            EP_star_array(ijk) = ep_star
! initializing blending stress parameters
            IF(BLENDING_STRESS.AND.TANH_BLEND) THEN
               ep_g_blend_start(ijk) = ep_star_array(ijk) * 0.99d0
               ep_g_blend_end(ijk)   = ep_star_array(ijk) * 1.01d0
            ELSE IF(BLENDING_STRESS.AND.SIGM_BLEND) THEN
               ep_g_blend_start(ijk) = ep_star * 0.97d0
               ep_g_blend_end(ijk) = ep_star * 1.01d0
            ELSE
               ep_g_blend_start(ijk) = ep_star_array(ijk)
               ep_g_blend_end(ijk)   = ep_star_array(ijk)
            ENDIF

         ENDDO   ! end loop over ijk
      ENDDO   ! end loop over MMAX


! Initializing parameters needed if a correlation is used to compute
! ep_star: initializing the indexing system.
      IF(YU_STANDISH .OR. FEDORS_LANDEL) THEN
         DO M = 1, SMAX
            IF(EP_S_MAX(M) == UNDEFINED) EP_S_MAX(M) = ONE-EP_STAR
         ENDDO

         IF (.NOT.CALL_DQMOM) THEN

! refer to Syam's dissertation
            IF (SMAX == 2) THEN
               ep_s_max_ratio(1,2) = ep_s_max(1)/ &
                  (ep_s_max(1)+(1.-ep_s_max(1))*ep_s_max(2))
            ENDIF

! initialize local variables
            DO I = 1, SMAX
               DP_TMP(I) = D_P0(I)
               M_MAX(I) = I
            ENDDO

! Rearrange the indices from coarsest particles to finest to be
! used in CALC_ep_star. Done here because it may need to be done
! for auto_restart
            DO I = 1, SMAX
               DO J = I, SMAX
                  IF(DP_TMP(I) < DP_TMP(J)) THEN
                     old_value = DP_TMP(I)
                     DP_TMP(I) = DP_TMP(J)
                     DP_TMP(J) = old_value
                  ENDIF
               ENDDO
            ENDDO

            DO I = 1, SMAX
               DO J = 1, SMAX
                  IF(DP_TMP(I) == D_P0(J) .AND. D_P0(I) .NE. D_P0(J)) THEN
                     M_MAX(I) = J
                  ENDIF
               ENDDO
            ENDDO
         ENDIF    ! if .not. call_dqmom
      ELSE   ! if .not. Yu-standish or Fedors-Landel
         EP_S_MAX(:) = ZERO
         EP_S_MAX_RATIO(:,:) = ZERO
         D_P_RATIO(:,:) = ZERO
         M_MAX(:) = 0
      ENDIF

      RETURN
      END SUBROUTINE SET_CONSTPROP
