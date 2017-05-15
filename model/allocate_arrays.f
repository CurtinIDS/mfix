! -*- f90 -*-
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutinee: ALLOCATE_ARRAYS                                        C
!  Purpose: allocate arrays                                            C
!                                                                      C
!  Author: M. Syamlal                                Date: 17-DEC-98   C
!  Reviewer:                                                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ALLOCATE_ARRAYS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use ambm
      use cdist
      use cont
      use des_rxns
      use drag
      use energy
      use fldvar
      use generate_particles, only: particle_count
      use geometry
      use ghdtheory
      use indices
      use kintheory
      use mflux
      use param
      use param1
      use pgcor
      use physprop
      use pscor
      use residual
      use run
      use rxns
      use scalars
      use iterate, only: errorpercent
      use tau_g
      use tau_s
      use trace
      use turb
      use visc_g
      use visc_s
      use vshear

      IMPLICIT NONE

!-----------------------------------------------
! Variables
!-----------------------------------------------

!ambm
      Allocate( A_m(DIMENSION_3, -3:3, 0:DIMENSION_M) )
      Allocate( B_m(DIMENSION_3, 0:DIMENSION_M) )

!cont
      Allocate( DO_CONT(0:DIMENSION_M) )

!drag
      Allocate(  F_gs(DIMENSION_3, DIMENSION_M) )
      Allocate(  F_ss(DIMENSION_3, 0:DIMENSION_LM) )

!Off diagonal friction coefficient in HYS drag relation
      IF(DRAG_TYPE_ENUM.EQ.HYS) &
         Allocate(  beta_ij(DIMENSION_3, 0:DIMENSION_M, 0:DIMENSION_M) )

!energy
      Allocate(  HOR_g (DIMENSION_3) )
      Allocate(  HOR_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  GAMA_gs (DIMENSION_3, DIMENSION_M) )
      Allocate(  GAMA_Rg (DIMENSION_3) )
      Allocate(  GAMA_Rs (DIMENSION_3, DIMENSION_M) )
      Allocate(  T_Rg (DIMENSION_3) )
      Allocate(  T_Rs (DIMENSION_3, DIMENSION_M) )

!fldvar
      Allocate(  EP_g (DIMENSION_3) )
      Allocate(  epg_jfac (DIMENSION_3p) )
      Allocate(  epg_ifac (DIMENSION_3p) )
      Allocate(  eps_ifac (DIMENSION_3p, DIMENSION_M) )
      Allocate(  EP_go (DIMENSION_3p) )
      Allocate(  P_g (DIMENSION_3) )
      Allocate(  P_go (DIMENSION_3p) )
      Allocate(  RO_g (DIMENSION_3) )
      Allocate(  RO_go (DIMENSION_3p) )
      Allocate(  ROP_g (DIMENSION_3) )
      Allocate(  ROP_go (DIMENSION_3p) )
      Allocate(  RO_S (DIMENSION_3, DIMENSION_M) )
      Allocate(  RO_So (DIMENSION_3p, DIMENSION_M) )
      Allocate(  ROP_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  ROP_so (DIMENSION_3p, DIMENSION_M) )

      Allocate(  EP_SS(DIMENSION_3,DIMENSION_M,DIMENSION_N_S) )
      Allocate(  ERR_ARRAY(DIMENSION_3,DIMENSION_M) )

      Allocate(  T_g (DIMENSION_3) )
      Allocate(  T_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  T_go (DIMENSION_3p) )
      Allocate(  T_so (DIMENSION_3p, DIMENSION_M) )
      Allocate(  X_g (DIMENSION_3, DIMENSION_N_g) )
      Allocate(  X_s (DIMENSION_3, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  X_go (DIMENSION_3p, DIMENSION_N_g) )
      Allocate(  X_so (DIMENSION_3p, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  U_g (DIMENSION_3) )
      Allocate(  U_go (DIMENSION_3p) )
      Allocate(  U_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  U_so (DIMENSION_3p, DIMENSION_M) )
      Allocate(  V_g (DIMENSION_3) )
      Allocate(  V_go (DIMENSION_3p) )
      Allocate(  V_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  V_so (DIMENSION_3p, DIMENSION_M) )
      Allocate(  W_g (DIMENSION_3) )
      Allocate(  W_go (DIMENSION_3p) )
      Allocate(  W_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  W_so (DIMENSION_3p, DIMENSION_M) )
      Allocate(  P_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  P_s_c (DIMENSION_3, DIMENSION_M) )
      Allocate(  P_s_v (DIMENSION_3) )
      Allocate(  P_s_f (DIMENSION_3) )
      Allocate(  P_s_p (DIMENSION_3) )
      Allocate(  P_star (DIMENSION_3) )
      Allocate(  P_staro (DIMENSION_3p) )
      Allocate(  THETA_m (DIMENSION_3, DIMENSION_M) )
      Allocate(  THETA_mo (DIMENSION_3p, DIMENSION_M) )

      IF(K_Epsilon)THEN
        Allocate(  K_Turb_G (DIMENSION_3) )
        Allocate(  K_Turb_Go (DIMENSION_3p) )
        Allocate(  E_Turb_G (DIMENSION_3) )
        Allocate(  E_Turb_Go (DIMENSION_3p) )
      ENDIF

      IF(DIMENSION_Scalar /= 0) THEN
        Allocate(  Scalar (DIMENSION_3,  DIMENSION_Scalar) )
        Allocate(  Scalaro (DIMENSION_3p, DIMENSION_Scalar) )
      ENDIF


!pgcor
      Allocate(  d_e(DIMENSION_3p, 0:DIMENSION_M) )
      Allocate(  d_n(DIMENSION_3p, 0:DIMENSION_M) )
      Allocate(  d_t(DIMENSION_3p, 0:DIMENSION_M) )
      Allocate(  Pp_g(DIMENSION_3p) )
      Allocate(  PHASE_4_P_g(DIMENSION_3p) )

!physprop
      Allocate(  MU_g (DIMENSION_3) )
      Allocate(  C_pg (DIMENSION_3) )
      Allocate(  C_ps (DIMENSION_3, DIMENSION_M) )
      Allocate(  K_g (DIMENSION_3) )
      Allocate(  K_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  Kth_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  Kphi_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DIF_g (DIMENSION_3p, DIMENSION_N_g) )
      Allocate(  DIF_s (DIMENSION_3p, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  MW_MIX_g (DIMENSION_3) )

!pscor
      Allocate(  e_e(DIMENSION_3p) )
      Allocate(  e_n(DIMENSION_3p) )
      Allocate(  e_t(DIMENSION_3p) )
      Allocate(  K_cp(DIMENSION_3p) )
      Allocate(  EPp(DIMENSION_3p) )
      Allocate(  PHASE_4_P_s(DIMENSION_3p) )

!residual
      Allocate( RESID(NRESID, 0:DIMENSION_M) )
      Allocate( MAX_RESID(NRESID, 0:DIMENSION_M) )
      Allocate( IJK_RESID(NRESID, 0:DIMENSION_M) )
      Allocate( NUM_RESID(NRESID, 0:DIMENSION_M) )
      Allocate( DEN_RESID(NRESID, 0:DIMENSION_M) )
      Allocate( RESID_PACK(NRESID*2*(DIMENSION_M+1)))

!rxns
      if (nRR .gt. 0) Allocate( ReactionRates(DIMENSION_3,nRR) )
      Allocate(  R_gp (DIMENSION_3p, DIMENSION_N_g) )
      Allocate(  R_sp (DIMENSION_3p, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  RoX_gc (DIMENSION_3p, DIMENSION_N_g) )
      Allocate(  RoX_sc (DIMENSION_3p, DIMENSION_M, DIMENSION_N_s) )
      Allocate(  SUM_R_g (DIMENSION_3p) )
      Allocate(  SUM_R_s (DIMENSION_3p, DIMENSION_M) )
      Allocate(  R_phase (DIMENSION_3, DIMENSION_LM+DIMENSION_M-1) )

!scalars
      IF(DIMENSION_Scalar /= 0) then
        Allocate(  Scalar_c (DIMENSION_3p,  DIMENSION_Scalar) )
        Allocate(  Scalar_p (DIMENSION_3p,  DIMENSION_Scalar) )
        Allocate(  Dif_Scalar (DIMENSION_3p, DIMENSION_Scalar) )
      ENDIF

! add by rong for dqmom
      Allocate(  D_p  (DIMENSION_3, DIMENSION_M) )
      Allocate(  D_po (DIMENSION_3, DIMENSION_M) )
!      Allocate(  ome  (DIMENSION_3, DIMENSION_M) )
!      Allocate(  ome_o (DIMENSION_3, DIMENSION_M) )
      Allocate(  Source_a(DIMENSION_3, DIMENSION_M) )
      Allocate(  S_bar( 0:DIM_Scalar2-1 ))
      Allocate(  Matrix_a(DIM_Scalar2,DIM_scalar2))
      Allocate(  Matrix_b(DIM_Scalar2,DIM_scalar2))
      Allocate(  Matrix_c(DIM_Scalar2,DIM_scalar2))
      Allocate(  Inv_a(DIM_Scalar2,DIM_scalar2))
      Allocate(  A( 1:DIMENSION_Scalar))
      Allocate(  omega(1:DIMENSION_m))
      ALLocate(  beta_a( DIM_Scalar,DIM_Scalar))
      ALLocate(  ystart( 1:DIM_Scalar2))
!     ALLocate(  g_a( 1:DIMENSION_Scalar))

! K-Epsilon Turbulence model
      IF(K_Epsilon) THEN
        Allocate(  K_Turb_G_c   (DIMENSION_3p) )
        Allocate(  K_Turb_G_p   (DIMENSION_3p) )
        Allocate(  Dif_K_Turb_G (DIMENSION_3p) )
        Allocate(  E_Turb_G_c   (DIMENSION_3p) )
        Allocate(  E_Turb_G_p   (DIMENSION_3p) )
        Allocate(  Dif_E_Turb_G (DIMENSION_3p) )
      ENDIF

! Simonin or Ahmadi model
      IF(KT_TYPE_ENUM==SIMONIN_1996 .OR.&
         KT_TYPE_ENUM==AHMADI_1995) THEN
        Allocate(  K_12 (DIMENSION_3) )
        Allocate(  Tau_12 (DIMENSION_3) )
        Allocate(  Tau_1 (DIMENSION_3) )
      ENDIF

!tau_g
      Allocate(  TAU_U_g(DIMENSION_3p) )
      Allocate(  TAU_V_g(DIMENSION_3p) )
      Allocate(  TAU_W_g(DIMENSION_3p) )
      Allocate(  DF_gu(DIMENSION_3p, -3:3) )
      Allocate(  DF_gv(DIMENSION_3p, -3:3) )
      Allocate(  DF_gw(DIMENSION_3p, -3:3) )
      Allocate(  CTAU_U_G(DIMENSION_3P))
      Allocate(  CTAU_V_G(DIMENSION_3P))
      Allocate(  CTAU_W_G(DIMENSION_3P))

!tau_s
      Allocate(  TAU_U_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  TAU_V_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  TAU_W_s(DIMENSION_3p, DIMENSION_M) )

! generate_particles / particle_count
      Allocate(  PARTICLE_COUNT(DIMENSION_3) )

!trace
      Allocate(  trD_s_C (DIMENSION_3, DIMENSION_M) )
      Allocate(  trD_s2 (DIMENSION_3, DIMENSION_M) )
      Allocate(  trD_s_Co (DIMENSION_3, DIMENSION_M) )
      Allocate(  trD_s_Co2 (DIMENSION_3, DIMENSION_M) )
!visc_g
      Allocate(  trD_g(DIMENSION_3) )
      Allocate(  MU_gt (DIMENSION_3) )
      Allocate(  EPMU_gt (DIMENSION_3p) )
      Allocate(  LAMBDA_gt (DIMENSION_3p) )
      Allocate(  EPLAMBDA_gt (DIMENSION_3) )
      Allocate(  L_scale (DIMENSION_3) )

!visc_s
      Allocate(  MU_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  EPMU_s (DIMENSION_3p, DIMENSION_M) )
      Allocate(  LAMBDA_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  EPLAMBDA_s (DIMENSION_3p, DIMENSION_M) )
      Allocate(  ALPHA_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  MU_s_c (DIMENSION_3, DIMENSION_M) )
      Allocate(  LAMBDA_s_c (DIMENSION_3, DIMENSION_M) )
      Allocate(  LAMBDA_s_v (DIMENSION_3) )
      Allocate(  LAMBDA_s_f (DIMENSION_3) )
      Allocate(  LAMBDA_s_p (DIMENSION_3) )
      Allocate(  MU_s_v (DIMENSION_3) )
      Allocate(  MU_s_f (DIMENSION_3) )
      Allocate(  MU_s_p (DIMENSION_3) )
      Allocate(  MU_b_v (DIMENSION_3) )
      Allocate(  EP_star_array (DIMENSION_3) )
      Allocate(  EP_g_blend_start (DIMENSION_3) )
      Allocate(  EP_g_blend_end (DIMENSION_3) )
      Allocate(  trD_s(DIMENSION_3, DIMENSION_M) )
      Allocate(  I2_devD_s (DIMENSION_3) )
      Allocate(  TrM_s (DIMENSION_3) )
      Allocate(  TrDM_s (DIMENSION_3) )

!shear quantities
      Allocate(  VSH(DIMENSION_3) )
      Allocate(  VSHE(DIMENSION_3) )

!mflux
      Allocate( Flux_gE(DIMENSION_3p) )
      Allocate( Flux_sE(DIMENSION_3p, DIMENSION_M) )
      Allocate( Flux_gN(DIMENSION_3p) )
      Allocate( Flux_sN(DIMENSION_3p, DIMENSION_M) )
      Allocate( Flux_gT(DIMENSION_3p) )
      Allocate( Flux_sT(DIMENSION_3p, DIMENSION_M) )
      IF(ADDED_MASS) THEN ! Fluxes calculated for just one 'bubble' species (M=M_AM)
         Allocate( Flux_gSE(DIMENSION_3p) )
         Allocate( Flux_sSE(DIMENSION_3p) )
         Allocate( Flux_gSN(DIMENSION_3p) )
         Allocate( Flux_sSN(DIMENSION_3p) )
         Allocate( Flux_gST(DIMENSION_3p) )
         Allocate( Flux_sST(DIMENSION_3p) )
      ENDIF
      Allocate( ROP_gE(DIMENSION_3p) )
      Allocate( ROP_sE(DIMENSION_3p, DIMENSION_M) )
      Allocate( ROP_gN(DIMENSION_3p) )
      Allocate( ROP_sN(DIMENSION_3p, DIMENSION_M) )
      Allocate( ROP_gT(DIMENSION_3p) )
      Allocate( ROP_sT(DIMENSION_3p, DIMENSION_M) )

! allocate variables for GHD Theory
      IF (KT_TYPE_ENUM == GHD_2007) THEN
        Allocate(  Flux_nE(DIMENSION_3p) )
        Allocate(  Flux_nN(DIMENSION_3p) )
        Allocate(  Flux_nT(DIMENSION_3p) )
        Allocate(  Zeta0(DIMENSION_3p) )   ! zeroth rate of cooling
        Allocate(  ZetaU(DIMENSION_3p) )   ! 1st order cooling rate transport coefficient
        Allocate(  DiT(DIMENSION_3p, DIMENSION_M) )   ! thermal diffusivity
        Allocate(  DijF(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )   ! mass mobility
        Allocate(  Lij(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )   ! thermal mobility
        Allocate(  Dij(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )   ! ordinary diffusion
        Allocate(  DijQ(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )   ! Dufour coeff.
        Allocate(  JoiX(DIMENSION_3p, DIMENSION_M) )   ! X- species mass flux
        Allocate(  JoiY(DIMENSION_3p, DIMENSION_M) )   ! Y- species mass flux
        Allocate(  JoiZ(DIMENSION_3p, DIMENSION_M) )   ! Z- species mass flux
        Allocate(  FiX(DIMENSION_3p, DIMENSION_M) )   ! X- external force
        Allocate(  FiY(DIMENSION_3p, DIMENSION_M) )   ! Y- external force
        Allocate(  FiZ(DIMENSION_3p, DIMENSION_M) )   ! Z- external force
        Allocate(  FiXvel(DIMENSION_3p, DIMENSION_M) )   ! X- external force
        Allocate(  FiYvel(DIMENSION_3p, DIMENSION_M) )   ! Y- external force
        Allocate(  FiZvel(DIMENSION_3p, DIMENSION_M) )   ! Z- external force
        Allocate(  DELTAU(DIMENSION_3p, DIMENSION_M) )
        Allocate(  DELTAV(DIMENSION_3p, DIMENSION_M) )
        Allocate(  DELTAW(DIMENSION_3p, DIMENSION_M) )
        Allocate(  dragFx(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  dragFy(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  dragFz(DIMENSION_3p, DIMENSION_M) )   ! Z- drag force
        Allocate(  dragFxflux(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  dragFyflux(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  dragFzflux(DIMENSION_3p, DIMENSION_M) )   ! Z- drag force
        Allocate(  FiMinusDragX(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  JoiMinusDragX(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  FiMinusDragY(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  JoiMinusDragY(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  FiMinusDragZ(DIMENSION_3p, DIMENSION_M) )   ! Z- drag force
        Allocate(  JoiMinusDragZ(DIMENSION_3p, DIMENSION_M) )   ! Z- drag force
        Allocate(  beta_cell_X(DIMENSION_3p, DIMENSION_M) )   ! X- drag force
        Allocate(  beta_cell_Y(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  beta_cell_Z(DIMENSION_3p, DIMENSION_M) )   ! Y- drag force
        Allocate(  beta_ij_cell_X(DIMENSION_3p, DIMENSION_M,DIMENSION_M) )   ! X- drag force
        Allocate(  beta_ij_cell_Y(DIMENSION_3p, DIMENSION_M,DIMENSION_M) )   ! Y- drag force
        Allocate(  beta_ij_cell_Z(DIMENSION_3p, DIMENSION_M,DIMENSION_M) )   ! Y- drag force
        Allocate(  DEL_DOT_J(DIMENSION_3p, DIMENSION_M) )
        Allocate(  DiT_HarmE(DIMENSION_3p) )
        Allocate(  DiT_HarmN(DIMENSION_3p) )
        Allocate(  DiT_HarmT(DIMENSION_3p) )
        Allocate(  Dij_HarmE(DIMENSION_3p, DIMENSION_M) )
        Allocate(  Dij_HarmN(DIMENSION_3p, DIMENSION_M) )
        Allocate(  Dij_HarmT(DIMENSION_3p, DIMENSION_M) )
        Allocate(  DijF_HarmE(DIMENSION_3p, DIMENSION_M) )
        Allocate(  DijF_HarmN(DIMENSION_3p, DIMENSION_M) )
        Allocate(  DijF_HarmT(DIMENSION_3p, DIMENSION_M) )
      ENDIF


! We need to set this even when KT_TYPE is not set to IA_NONEP - at
! least in the current version of the code and needs to be revisited
      Allocate(  KTMOM_U_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  KTMOM_V_s(DIMENSION_3p, DIMENSION_M) )
      Allocate(  KTMOM_W_s(DIMENSION_3p, DIMENSION_M) )

! allocate variables for Iddir & Arastoopour (2005) kinetic theory
! EDvel_sM_ip & EDT_s_ip are also used for Garzy & Dufty (1999) kinetic theory
      IF (KT_TYPE_ENUM == IA_2005) THEN
         Allocate(  trD_s2_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  MU_sM_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  MU_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  XI_sM_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  XI_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  Fnu_s_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  FT_sM_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  FT_sL_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  Kth_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  Knu_sM_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  Knu_sL_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  Kvel_s_ip(DIMENSION_3, DIMENSION_M, DIMENSION_M) )
         Allocate(  EDvel_sL_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  ED_ss_ip(DIMENSION_3p, 0:DIMENSION_LM) )
      ENDIF
      IF (KT_TYPE_ENUM == GTSH_2012) THEN
         Allocate(  A2_gtsh(DIMENSION_3) )
         Allocate(  xsi_gtsh(DIMENSION_3) )
      ENDIF
      IF (KT_TYPE_ENUM == IA_2005 .OR. &
          KT_TYPE_ENUM == GD_1999 .OR. &
          KT_TYPE_ENUM == GTSH_2012) THEN
         Allocate(  EDT_s_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
         Allocate(  EDvel_sM_ip(DIMENSION_3p, DIMENSION_M, DIMENSION_M) )
      ENDIF

      Allocate(errorpercent(0:MMAX))

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_ARRAYS_GEOMETRY                               !
!  Author: M. Syamlal                                 Date: 21-JAN-92  !
!                                                                      !
!  Purpose: Calculate X, X_E,  oX, oX_E                                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_ARRAYS_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Domain decomposition and dimensions
      use geometry, only: oDX, oDX_E
      use geometry, only: oDZ, oDZ_T
      use geometry, only: oDY, oDY_N
      use geometry, only: X, X_E, oX, oX_E, cyl_X, cyl_X_E
      use geometry, only: Z, Z_T
! Averaging factors.
      use geometry, only: FX_E, FX_E_bar, FX, FX_bar
      use geometry, only: FY_N, FY_N_bar
      use geometry, only: FZ_T, FZ_T_bar
! Domain flags.
      use geometry, only: ICBC_FLAG
      use geometry, only: FLAG, FLAG3
      use geometry, only: FLAG_E, FLAG_N, FLAG_T
! Domain volumes and areas.
      use geometry, only: VOL, VOL_SURR, AYZ, AXZ, AXY! Scalar grid
      use geometry, only: VOL_U, AYZ_U, AXZ_U, AXY_U  ! X-Momentum
      use geometry, only: VOL_V, AYZ_V, AXZ_V, AXY_V  ! Y-Momentum
      use geometry, only: VOL_W, AYZ_W, AXZ_W, AXY_W  ! Z-Momentum
! Axis decomposition
      USE param, only: DIMENSION_I, DIMENSION_J, DIMENSION_K
      USE param, only: DIMENSION_3, DIMENSION_4
      USE param, only: DIMENSION_3L, DIMENSION_3P
! Flag for POST_MFIX
      use cdist, only: bDoing_postmfix

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager

      IMPLICIT NONE

! Local Variables:
!---------------------------------------------------------------------//
! Error Flag
      INTEGER :: IER
! Flag indicating that the arrays were previously allocated.
      INTEGER, SAVE :: CALLED = -1
!......................................................................!

      CALLED = CALLED + 1

      IF(CALLED > 0) THEN
         IF(.NOT.bDoing_postmfix) THEN
            RETURN
         ELSEIF(mod(CALLED,2) /= 0) THEN
            RETURN
         ENDIF
      ENDIF

! Initialize the error manager.
      CALL INIT_ERR_MSG("ALLOCATE_ARRAYS_GEOMETRY")

! Allocate geometry components related to the mesh. Check the
! allocation error status and abort if any failure is detected.
      ALLOCATE( X     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( cyl_X     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( X_E   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( cyl_X_E   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oX    (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oX_E  (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oDX   (0:DIMENSION_I), STAT=IER)
      ALLOCATE( oDX_E (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( oDY   (0:DIMENSION_J), STAT=IER )
      ALLOCATE( oDY_N (0:DIMENSION_J), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( Z     (0:DIMENSION_K), STAT=IER )
      ALLOCATE( Z_T   (0:DIMENSION_K), STAT=IER )
      ALLOCATE( oDZ   (0:DIMENSION_K), STAT=IER )
      ALLOCATE( oDZ_T (0:DIMENSION_K), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( FX     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( FX_bar (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( FX_E     (0:DIMENSION_I), STAT=IER)
      ALLOCATE( FX_E_bar (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      ALLOCATE( FY_N     (0:DIMENSION_J), STAT=IER )
      ALLOCATE( FY_N_bar (0:DIMENSION_J), STAT=IER )
      IF(IER /= 0) goto 500

      ALLOCATE( FZ_T     (0:DIMENSION_K), STAT=IER )
      ALLOCATE( FZ_T_bar (0:DIMENSION_K), STAT=IER )
      IF(IER /= 0) goto 500

! Flags for the scalar grid.
      Allocate( FLAG  (DIMENSION_3), STAT=IER )
      Allocate( FLAG3 (DIMENSION_4), STAT=IER )
      IF(IER /= 0) goto 500

! Flags for the momentum grids.
      Allocate( FLAG_E (DIMENSION_3), STAT=IER )
      Allocate( FLAG_N (DIMENSION_3), STAT=IER )
      Allocate( FLAG_T (DIMENSION_3), STAT=IER )
      IF(IER /= 0) goto 500

! Text flags for scalar grid.
      Allocate( ICBC_FLAG (DIMENSION_3L), STAT=IER )
      IF(IER /= 0) goto 500

! Volume and face-areas of scalar grid.
      Allocate( VOL (DIMENSION_3),  STAT=IER )
      Allocate( AYZ (DIMENSION_3P), STAT=IER )
      Allocate( AXZ (DIMENSION_3P), STAT=IER )
      Allocate( AXY (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

      ! total volume of each cell's surrounding stencil cells
      Allocate( VOL_SURR (DIMENSION_3), STAT=IER )

! Volume and face-areas of X-Momentumn grid.
      Allocate( VOL_U (DIMENSION_3),  STAT=IER )
      Allocate( AYZ_U (DIMENSION_3P), STAT=IER )
      Allocate( AXZ_U (DIMENSION_3P), STAT=IER )
      Allocate( AXY_U (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

! Volume and face-areas of Y-Momentum grid.
      Allocate( VOL_V (DIMENSION_3),  STAT=IER )
      Allocate( AYZ_V (DIMENSION_3P), STAT=IER )
      Allocate( AXZ_V (DIMENSION_3P), STAT=IER )
      Allocate( AXY_V (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

! Volume and face-areas of Z-Momentum grid.
      Allocate( VOL_W (DIMENSION_3),  STAT=IER )
      Allocate( AYZ_W (DIMENSION_3P), STAT=IER )
      Allocate( AXZ_W (DIMENSION_3P), STAT=IER )
      Allocate( AXY_W (DIMENSION_3P), STAT=IER )
      IF(IER /= 0) goto 500

! Collect the error flags from all ranks. If all allocaitons were
! successfull, do nothing. Otherwise, flag the error and abort.
! Note that the allocation status is checked in groups. This can
! be increase if tracking the source of an allocation failure.
  500 CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Failure during array allocation.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_ARRAYS_INCREMENTS                             !
!  Author: M. Syamlal, W. Rogers                      Date: 10-DEC-91  !
!                                                                      !
!  Purpose: The purpose of this module is to create increments to be   !
!           stored in the array STORE_INCREMENT which will be added    !
!           to cell index ijk to find the effective indices of its     !
!           neighbors. These increments are found using the 'class'    !
!           of cell ijk. The class is determined based on the          !
!           neighboring cell type, i.e. wall or fluid.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE ALLOCATE_ARRAYS_INCREMENTS

      USE param
      USE param1
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE fldvar
      USE funits

! Module procedures
!---------------------------------------------------------------------//
      use mpi_utility, only: GLOBAL_ALL_SUM
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Error flag.
      INTEGER :: IER
! Flag indicating that the arrays were previously allocated.
      LOGICAL, SAVE :: ALREADY_ALLOCATED = .FALSE.
!......................................................................!

      IF(ALREADY_ALLOCATED) RETURN

! Initialize the error manager.
      CALL INIT_ERR_MSG("ALLOCATE_ARRAYS_INCREMENTS")

! Allocate increment arrays and report an allocation errors.
      Allocate( I_OF (DIMENSION_3), STAT=IER)
      Allocate( J_OF (DIMENSION_3), STAT=IER)
      Allocate( K_OF (DIMENSION_3), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Im1 (0:DIMENSION_I), STAT=IER)
      Allocate( Ip1 (0:DIMENSION_I), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Jm1 (0:DIMENSION_J), STAT=IER)
      Allocate( Jp1 (0:DIMENSION_J), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Km1 (0:DIMENSION_K), STAT=IER)
      Allocate( Kp1 (0:DIMENSION_K), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( STORE_LM (DIMENSION_M, DIMENSION_M), STAT=IER)
      Allocate( CELL_CLASS (DIMENSION_3), STAT=IER)
      IF(IER /= 0) goto 500


! Allocate increment arrays and report an allocation errors.
      Allocate( I3_OF (DIMENSION_4), STAT=IER)
      Allocate( J3_OF (DIMENSION_4), STAT=IER)
      Allocate( K3_OF (DIMENSION_4), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Im1_3 (-1:DIMENSION_I+1), STAT=IER)
      Allocate( Ip1_3 (-1:DIMENSION_I+1), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Jm1_3 (-1:DIMENSION_J+1), STAT=IER)
      Allocate( Jp1_3 (-1:DIMENSION_J+1), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( Km1_3 (-1:DIMENSION_K+1), STAT=IER)
      Allocate( Kp1_3 (-1:DIMENSION_K+1), STAT=IER)
      IF(IER /= 0) goto 500

      Allocate( CELL_CLASS3 (DIMENSION_4), STAT=IER)
      IF(IER /= 0) goto 500

! Collect the error flags from all ranks. If all allocaitons were
! successfull, do nothing. Otherwise, flag the error and abort.
! Note that the allocation status is checked in groups. This can
! be increase if tracking the source of an allocation failure.
  500 CALL GLOBAL_ALL_SUM(IER)

      IF(IER /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Failure during array allocation.')

      ALREADY_ALLOCATED = .TRUE.

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE ALLOCATE_ARRAYS_INCREMENTS
