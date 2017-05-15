!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: physprop                                                    C
!  Purpose: Common block containing physical property data             C
!                                                                      C
!  Author: M. Syamlal                                 Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      MODULE physprop

! Modules
!---------------------------------------------------------------------//
      Use param, only: dim_m, dim_n, dim_n_g, dim_n_s
!---------------------------------------------------------------------//


! Number of solids phases
      INTEGER :: MMAX

! Real number of solids phases for GHD theory
      INTEGER :: SMAX

! Particle diameters
      DOUBLE PRECISION :: D_p0(DIM_M)

! Constant or baseline solids phase densities.
      DOUBLE PRECISION :: RO_s0(DIM_M)

! Constant solids phase species mass fractions. These values delinate
! the baseline/initial solids phase composition for ariable density
      DOUBLE PRECISION :: X_S0(DIM_M, DIM_N_s)

! Density of solid species (constant)
      DOUBLE PRECISION :: RO_Xs0(DIM_M, DIM_N_s)

! The index of an inert solids phase species. This is needed for
! calculating the variable solids phase density.
      INTEGER :: INERT_SPECIES(DIM_M)

! Inert solids phase species mass fraction in dilute region.  This is needed for
! calculating the variable solids phase density.
      DOUBLE PRECISION :: DIL_INERT_X_VSD(DIM_M)

! Factor to define dilute region for special treatment of
! the variable solids phase density.
      DOUBLE PRECISION :: DIL_FACTOR_VSD

! Particle shape factor
      DOUBLE PRECISION :: SHAPE_FACTOR(DIM_M)

! Specified constant solids viscosity
      DOUBLE PRECISION MU_s0(DIM_M)

! Flag indicates whether the phase becomes close-packed at ep_star
      LOGICAL :: CLOSE_PACKED (DIM_M)

! Specified constant gas density
      DOUBLE PRECISION RO_g0

! Specified constant gas viscosity
      DOUBLE PRECISION MU_g0

! Virtual (added) mass coefficient Cv
      DOUBLE PRECISION Cv

! Gas viscosity
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  MU_g

! Average molecular weight of gas
      DOUBLE PRECISION MW_AVG

! Constant constant-pressure specific heat of gas
      DOUBLE PRECISION C_pg0

! Reference temperature for enthalpy calculations (K)
      DOUBLE PRECISION, PARAMETER :: T_ref = 298

! Constant pressure specific heat of gas
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  C_pg

! Constant constant-pressure specific heat of solids
      DOUBLE PRECISION C_ps0(DIM_M)

! Constant pressure specific heat of solids
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  C_ps

! Specified constant gas conductivity
      DOUBLE PRECISION :: K_g0

! Conductivity of gas
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  K_g

! Specified constant solids conductivity
      DOUBLE PRECISION K_s0(DIM_M)

! Conductivity of solids
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  K_s

! Granular Temperature Conductivity (associated with temperature grad)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Kth_s

! Granular Temperature Conductivity (associated with volume fraction grad)
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  Kphi_s

! Specified constant gas diffusivity
      DOUBLE PRECISION DIF_g0

! Diffusivity of gas species N
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  DIF_g

! Specified constant solids diffusivity
      DOUBLE PRECISION DIF_s0(DIM_M)

! Diffusivity of solids species N
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE ::  DIF_s

! Total number of gas or solids species
      INTEGER :: NMAX(0:DIM_M) ! Runtime (all phases)
      INTEGER :: NMAX_g        ! Number of gas phase species
      INTEGER :: NMAX_s(DIM_M) ! Number of solids phase species

! Molecular weight of gas species
      DOUBLE PRECISION :: MW_g (DIM_N_g)

! Molecular weight of solids species
      DOUBLE PRECISION :: MW_s (DIM_M, DIM_N_s)

! Molecular weight of gas mixture
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MW_MIX_g

! Logical for reading thermochemical database.
      LOGICAL :: DATABASE_READ = .TRUE.

! Polynomical coefficients for calculating specific heat.
      DOUBLE PRECISION Alow (7,0:DIM_M, DIM_N) ! Tlow --> Tcom
      DOUBLE PRECISION Ahigh(7,0:DIM_M, DIM_N) ! Tcom --> Thigh

! Range where the polynomials are valid.
      DOUBLE PRECISION Thigh(0:DIM_M, DIM_N) ! Upper bound
      DOUBLE PRECISION Tlow (0:DIM_M, DIM_N) ! Lower bound
      DOUBLE PRECISION Tcom (0:DIM_M, DIM_N) ! Switch from low to high

! Heat of formation at Tref divided by the gas constant.
      DOUBLE PRECISION HfrefoR(0:DIM_M, DIM_N)

! Reference values.
      DOUBLE PRECISION ICpoR_l(0:DIM_M, DIM_N)
      DOUBLE PRECISION ICpoR_h(0:DIM_M, DIM_N)

CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_MW (X_g, DIM, L, NMAX, MW_g)                      C
!  Purpose: Calculate average molecular weight of gas                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-OCT-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:None                                           C
!  Variables modified:None                                             C
!                                                                      C
!  Local variables: SUM, N                                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION CALC_MW (X_G, DIM, L, NMAX, MW_G)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE toleranc
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Mass fraction array's Ist dimension
      INTEGER          DIM
!
!
!                      Max of X_g array 2nd index and MW_g array index
      INTEGER          NMAX
!
!                      Mass fraction array
!
      DOUBLE PRECISION X_g(DIM, NMAX)
!
!                      Moleculare weight array
!
      DOUBLE PRECISION MW_g(NMAX)
!
!                      Mass fraction array Ist index
      INTEGER          L
!
!  Local variable
!
!                      local sum
      DOUBLE PRECISION SUM
!
!                      local index
      INTEGER           NN
!-----------------------------------------------
!
      SUM = ZERO
      DO NN = 1, NMAX
         SUM = SUM + X_G(L,NN)/MW_G(NN)
      END DO
      CALC_MW = ONE/MAX(SUM,OMW_MAX)
!
      RETURN
      END FUNCTION CALC_MW

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  FUNCTION: blend_function                                            C
!  Purpose: To calculate blending function                             C
!                                                                      C
!  Author: S. Pannala                                 Date: 28-FEB-06  C
!  Reviewer:                                          Date:            C
!  Modified:                                          Date: 24-OCT-06  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      DOUBLE PRECISION FUNCTION blend_function(IJK)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE fldvar
      USE fun_avg
      USE functions
      USE geometry
      USE indices
      USE param
      USE param1
      USE run
      USE toleranc
      USE visc_s
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! IJK index
      INTEGER, INTENT(IN) :: IJK
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
! Blend Factor
      Double Precision:: blend, blend_right
! Scale Factor
      Double Precision, Save:: scale
! Midpoint
      Double Precision, Save:: ep_mid_point
!-----------------------------------------------

! Tan hyperbolic blending of stresses
      IF(TANH_BLEND) THEN
         IF(EP_g(IJK) .LT. ep_g_blend_end(ijk).AND. EP_g(IJK) .GT. ep_g_blend_start(ijk)) THEN
            ep_mid_point = (ep_g_blend_end(IJK)+ep_g_blend_start(IJK))/2.0d0
            blend = tanh(2.0d0*pi*(ep_g(IJK)-ep_mid_point)/ &
            (ep_g_blend_end(IJK)-ep_g_blend_start(IJK)))
            blend = (blend+1.0d0)/2.0d0
         ELSEIF(EP_g(IJK) .GE. ep_g_blend_end(ijk)) THEN
            blend = 1.0d0
         ELSEIF(EP_g(IJK) .LE. ep_g_blend_start(ijk)) THEN
            blend = 0.0d0
         ENDIF

! Truncated and Scaled Sigmoidal blending of stresses
      ELSEIF(SIGM_BLEND) THEN
         IF(FIRST_PASS) THEN
            blend_right =  1.0d0/(1+0.01d0**((ep_g_blend_end(IJK)-ep_star_array(IJK))&
            /(ep_g_blend_end(IJK)-ep_g_blend_start(IJK))))
            blend_right = (blend_right+1.0d0)/2.0d0
            scale = 1.0d0/blend_right
            write(*,*) 'Blending value at end and scaling factor', blend_right, scale
            FIRST_PASS = .FALSE.
         ENDIF
         IF(EP_g(IJK) .LT. ep_g_blend_end(ijk)) THEN
            blend =  scale/(1+0.01d0**((ep_g(IJK)-ep_star_array(IJK))&
            /(ep_g_blend_end(IJK)-ep_g_blend_start(IJK))))
         ELSEIF(EP_g(IJK) .GE. ep_g_blend_end(ijk)) THEN
            blend = 1.0d0
         ENDIF

      ENDIF

      blend_function = blend

      RETURN
      END FUNCTION blend_function

      END MODULE physprop
