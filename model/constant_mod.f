!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: constant                                                    C
!  Purpose: Common block containing physical constants and constants   C
!           used in the numerical technique                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 5-FEB-92   C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!    Gera, D., Syamlal, M., and O'Brien, T. J., "Hydrodynamics of      C
!      particle segregation in fluidized beds", Int. J. of Multiphase  C
!      Flow, Vol 30, 2004, pp. 419-428.                                C
!    Johnson, P. C., and Jackson, R., "Frictional-collisional          C
!      constitutive relations for granluar materials, with application C
!      to plane shearing", JFM, Vol. 176, 1987, pp. 67-93.             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE constant


! Modules
!---------------------------------------------------------------------//
      Use param, only: dim_m, dimension_c
!---------------------------------------------------------------------//

! Packed bed (close packed) void fraction
      DOUBLE PRECISION :: EP_star

! parameters used in the correlation to calculate the local maximum
! solids volume fraction for a polydisperse powder: ep_s_max_ratio,
! d_p_ratio and ep_s_max, m_max
      DOUBLE PRECISION :: ep_s_max_ratio(DIM_M, DIM_M), &
                          d_p_ratio(DIM_M, DIM_M)
! maximum packing volume fraction for indicate particulate phase
! its value will default to 1-ep_star
      DOUBLE PRECISION :: ep_s_max(DIM_M)
! Index to rearrange particles from coarsest to finest for use in
! function CALC_ep_star(IJK,IER)
      INTEGER :: M_MAX(DIM_M)

! SWITCH enables us to turn on/off modifications to certain kinetic
! theory models for granular solids (i.e. no gas) that have been
! adjusted to account for the presence of a fluid phase. If one wants
! to simulate gas-particle flow then set SWITCH=1. As a result, the
! effects of drag on particle viscosity/conductivity will be
! incorporated. Additional gas-solids terms may also have been
! introduced into the granular energy balance depending on the KT
! model (see source_granular_energy for details). If we want to
! simulate pure granular flow without the effects of an interstitial
! gas, set SWITCH=0.
      DOUBLE PRECISION, PARAMETER :: SWITCH=1.d0

! ALPHA is a parameter introduced into the theory of Lun_1984 for
! calculating solids viscosity. It also appears when invoking the
! solids frictional model FRICTION, which uses the Lun et al.
! theory. The factor (2+alpha)/3 was eliminated in the complete
! analysis of Lun et al. but was introduced as an adjustable
! parameter. To recover the original theory alpha should be set to
! 1. For details see Johnson and Jackson, 1987.
      DOUBLE PRECISION, PARAMETER :: ALPHA = 1.6d0


! parameter used in the solids-solids drag model invoked in the
! default KT (Lun_1984). For details see Gera et al., 2004
      DOUBLE PRECISION :: SEGREGATION_SLOPE_COEFFICIENT

! SWITCH_IA enforces consistency in the solids viscosity and
! conductivity so that the results using 2 or more identical
! solids phases are the same as an equivalent single solids
! phase. Set to false to use original (published) theory of
! Iddir-Arastoopour.
      LOGICAL, PARAMETER :: SWITCH_IA = .TRUE.


! PHIP = Specularity coefficient associated with particle wall
! collisions
      DOUBLE PRECISION :: PHIP
! PHIP0 specularity coefficient for r->0
      double precision :: phip0
! k4phi k=7/2*mu*(1+e_w)
      double precision :: k4phi
! e_w = particle-wall coefficient of restitution
      DOUBLE PRECISION :: e_w

! Parameters used in the solids frictional model FRICTION:
! - Fr, N_Pc, D_Pc, and EPS_F_min are all used in the equation for
!   Pc, the critical solids pressure:
!     Fr = Constant with dyne/cm2 units of pressure. It will be
!          automatically converted to Pa in calc_mu_s.f
!     N_Pc = exponent in numerator
!     D_Pc = exponent in denominator
!     EPS_f_min = minimum solids fraction above which friction
!                 kicks in
! - N_Pf appears as an exponent in the equation of state for Pf, the
!   frictional pressure:
! - delta is a small deviation in void fraction near packing where
!   Pc and dPc/deps are calculated.
      DOUBLE PRECISION :: EPS_f_min
      DOUBLE PRECISION :: Fr, N_Pc, D_Pc, N_Pf, delta
      PARAMETER(Fr = 0.5d0, N_Pc=2d0, D_Pc=5d0, N_Pf=1.03d0, delta=1d-2)

! Coefficient of restitution
      DOUBLE PRECISION :: C_e

! (1+C_e)/2.
      DOUBLE PRECISION :: eta

! particle-type dependent rest. coef. for use in GHD theory
      DOUBLE PRECISION :: r_p(DIM_M, DIM_M)

! Coeficient of friction
      DOUBLE PRECISION :: C_f

! Angle of internal friction (degrees)
      DOUBLE PRECISION :: Phi

! Angle of wall-particle friction (degrees)
      DOUBLE PRECISION :: Phi_w

! (k=) Sin(PHI) in frictional-flow stress formulation
      DOUBLE PRECISION :: Sin_Phi

! Sin^2(PHI) in plastic-flow stress formulation
      DOUBLE PRECISION :: Sin2_Phi

! (3-2k^2)/6k^2 in Plastic-flow stress formulation
      DOUBLE PRECISION :: F_Phi

! tan(PHI_w)
      DOUBLE PRECISION :: tan_Phi_w

! Excluded volume (Boyle-Massoudi stress tensor)
      DOUBLE PRECISION :: V_ex

! Coefficients for calibrating Syamlal-O'Brien drag correlation with
! Umf data
      DOUBLE PRECISION :: drag_c1, drag_d1

! success-factor for aggregation and breakage
      DOUBLE PRECISION :: AGGREGATION_EFF
      DOUBLE PRECISION :: BREAKAGE_EFF

! UNIT conversion factor for pressure (Barye to Pa if SI)
      DOUBLE PRECISION :: to_SI

! Gravitational acceleration
      DOUBLE PRECISION :: GRAVITY, GRAVITY_X, GRAVITY_Y, GRAVITY_Z

! Universal gas constant
      DOUBLE PRECISION :: GAS_CONST

! Universal gas constant in cal/mol.K
      DOUBLE PRECISION, PARAMETER :: GAS_CONST_cal = 1.987207D0

! Pi, the ubiquitous irrational number
      DOUBLE PRECISION, PARAMETER :: Pi = 4.D0*ATAN(1.D0)

! Square root of Pi
      DOUBLE PRECISION, PARAMETER :: SQRT_Pi = 2.D0*SQRT(ATAN(1.D0))

! Maximum pressure correction allowed in one iteration
      DOUBLE PRECISION :: MAX_DELP

! User defined constants
      DOUBLE PRECISION :: C (DIMENSION_C)

! Names of user defined constants (for output file only)
      CHARACTER(LEN=20) :: C_NAME (DIMENSION_C)

! Move these to turb at some point:
! Scale factor for gas turbulence length scale
      DOUBLE PRECISION :: K_scale

! Default value for characteristic length for turbulence
      DOUBLE PRECISION :: L_scale0

! Maximum value of turbulent viscosity
      DOUBLE PRECISION :: MU_gmax


      END MODULE constant
