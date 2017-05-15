!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module: ic                                                          !
!  Author: M. Syamlal                                 Date: dd-mmm-yy  !
!                                                                      !
!  Purpose: Global initial conditions variables.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      MODULE ic

! Maximum number of IC regions.
      use param, only: DIMENSION_IC
! Maximum number of solids phases.
      use param, only: DIM_M
! Maximum number of gas phase species
      use param, only: DIM_N_g
! Maximum number of solids phase species
      use param, only: DIM_N_s
! Maximum number of scalar equations
      use param, only: DIM_Scalar

! IC region West face, X-coordinate
      DOUBLE PRECISION :: IC_X_w (DIMENSION_IC)

! IC region East face, X-coordinate
      DOUBLE PRECISION :: IC_X_e (DIMENSION_IC)

! IC region South face, Y-coordinate
      DOUBLE PRECISION :: IC_Y_s (DIMENSION_IC)

! IC region North face, Y-coordinate
      DOUBLE PRECISION :: IC_Y_n (DIMENSION_IC)

! IC region Bottom face, Z-coordinate
      DOUBLE PRECISION :: IC_Z_b (DIMENSION_IC)

! IC region Top face, Z-coordinate
      DOUBLE PRECISION :: IC_Z_t (DIMENSION_IC)

! IC region, West face, I Index
      INTEGER :: IC_I_w (DIMENSION_IC)

! IC region, East face, I Index
      INTEGER :: IC_I_e (DIMENSION_IC)

! IC region, South face, J Index
      INTEGER :: IC_J_s (DIMENSION_IC)

! IC region, North face, J Index
      INTEGER :: IC_J_n (DIMENSION_IC)

! IC region, Bottom face, K Index
      INTEGER :: IC_K_b (DIMENSION_IC)

! IC region, Top face, K Index
      INTEGER :: IC_K_t (DIMENSION_IC)

! Type of initial condition: PATCH
      CHARACTER(LEN=16) :: IC_TYPE(DIMENSION_IC)

! Initial gas phase volume fraction
      DOUBLE PRECISION :: IC_EP_g (DIMENSION_IC)

! Initial gas pressure
      DOUBLE PRECISION :: IC_P_g (DIMENSION_IC)

! Initial gas pressure
      DOUBLE PRECISION :: IC_P_star(DIMENSION_IC)

! Initial turbulence length scale
      DOUBLE PRECISION :: IC_L_scale(DIMENSION_IC)

! Initial macroscopic density of solids phases
      DOUBLE PRECISION :: IC_ROP_s(DIMENSION_IC, DIM_M)

! Initial solids phase volume fraction
      DOUBLE PRECISION :: IC_EP_s (DIMENSION_IC, DIM_M)

! Initial gas phase temperature
      DOUBLE PRECISION :: IC_T_g(DIMENSION_IC)

! Initial solids phase temperature
      DOUBLE PRECISION :: IC_T_s(DIMENSION_IC, DIM_M)

! Initial granular temperature
      DOUBLE PRECISION :: IC_Theta_m(DIMENSION_IC, DIM_M)

! Initial x-component of gas velocity
      DOUBLE PRECISION :: IC_U_g(DIMENSION_IC)

! Initial x-component of solids phase velocity
      DOUBLE PRECISION :: IC_U_s(DIMENSION_IC, DIM_M)

! Initial y-component of gas velocity
      DOUBLE PRECISION :: IC_V_g(DIMENSION_IC)

! Initial y-component of solids phase velocity
      DOUBLE PRECISION :: IC_V_s(DIMENSION_IC, DIM_M)

! Initial z-component of gas velocity
      DOUBLE PRECISION :: IC_W_g(DIMENSION_IC)

! Initial z-component of solids phase velocity
      DOUBLE PRECISION :: IC_W_s(DIMENSION_IC, DIM_M)

! Logical variable to determine whether an ic is defined
      LOGICAL :: IC_DEFINED (DIMENSION_IC)

! Initial gas species mass fractions
      DOUBLE PRECISION :: IC_X_g(DIMENSION_IC, DIM_N_g)

! Initial solids species mass fractions
      DOUBLE PRECISION :: IC_X_s(DIMENSION_IC, DIM_M, DIM_N_s)

! Gas phase radiation coefficient
      DOUBLE PRECISION :: IC_GAMA_Rg (DIMENSION_IC)

! Gas phase radiation temperature
      DOUBLE PRECISION :: IC_T_Rg(DIMENSION_IC)

! Solids phase-1 radiation coefficient
      DOUBLE PRECISION :: IC_GAMA_Rs(DIMENSION_IC, DIM_M)

! Solids phase-1 radiation temperature
      DOUBLE PRECISION :: IC_T_Rs(DIMENSION_IC, DIM_M)

! Initial scalar value in a region
      DOUBLE PRECISION :: IC_Scalar(DIMENSION_IC, DIM_scalar)

! Initial K & Epsilon values in a region
      DOUBLE PRECISION :: IC_K_Turb_G(DIMENSION_IC)
      DOUBLE PRECISION :: IC_E_Turb_G(DIMENSION_IC)

! Initial conditions for DES cases (such as, DEM, MPPIC, hybrid)

! Flag to extend the lattice distribution in a given IC to available area
      LOGICAL :: IC_DES_FIT_TO_REGION (DIMENSION_IC)

! Flag to specify the initial constant number of particles per cell
! for the MPPIC method initialization.
! Statistical weight of parcels will be calculated by the code
      INTEGER :: IC_PIC_CONST_NPC(DIMENSION_IC, DIM_M)

! Flag to specify the initial constant statistical weight.
! for the MPPIC method initialization.
! Number of computational particles/parcels will be calculated by the code
      DOUBLE PRECISION :: IC_PIC_CONST_STATWT(DIMENSION_IC, DIM_M)


      END MODULE ic
