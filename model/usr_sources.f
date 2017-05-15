!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: USR_SOURCES                                             C
!  Purpose: Hook for user defined source terms                         C
!                                                                      C
!  Comments:                                                           C
!  Discretized equations take form as a matrix equation Ax=b.          C
!  Terms that may be represented as a coefficient of the dependent     C
!  variable go to the center coefficient (ap) through sourcelhs        C
!  (i.e., included on the left-hand-side (lhs)).                       C
!  Terms that are a constant go to the source vector (b) through       C
!  sourcerhs (i.e., included on the right-hand-side (rhs)).            C
!                                                                      C
!  Source terms are often a function of the dependent variable. To     C
!  aid in convergence, this dependency should be acknowledged in the   C
!  equation. Incorporate the dependency of the source term on the      C
!  dependent variable through both ap and b.                           C
!                                                                      C
!  See Patankar, S. V., Numerical heat transfer and fluid flow,        C
!  Taylor and Francis, 1980, for rules and suggestions for             C
!  appropriate discretization of the source term.                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR_SOURCES(lEQ_NO, IJK, sourcelhs, sourcerhs, M, N)

! Modules 
!-----------------------------------------------
      use constant, only: pi, gas_const, gravity
      use fldvar, only: u_g, v_g, w_g
      use fldvar, only: u_s, v_s, w_s
      use fldvar, only: ep_g, rop_g, ro_g, T_g, X_g, P_g
      use fldvar, only: ep_s, rop_s, ro_s, T_s, X_s, d_p, theta_m
      use fldvar, only: k_turb_g, e_turb_g, scalar
      use functions
      use geometry
      use indices, only: i_of, j_of, k_of
      use indices, only: im1, ip1, jm1, jp1, km1, kp1
      use param1, only: zero, one, half, undefined, undefined_i
      use physprop
      use scalars, only: phase4scalar
      use usr_src
      IMPLICIT NONE

! Dummy arguments
!-----------------------------------------------
! reference equation 
      INTEGER, INTENT(IN) :: lEQ_NO
! index
      INTEGER, INTENT(IN) :: IJK
! source terms which appear appear in the 
! center coefficient (lhs) - part of a_m matrix
! source vector (rhs) - part of b_m vector
      DOUBLE PRECISION, INTENT(OUT) :: sourcelhs, sourcerhs
! Phase index 
      INTEGER, INTENT(IN) :: M
! Species index OR scalar equation number 
! (if applicable otherwise undefined_i)
      INTEGER, INTENT(IN) :: N

! Local variables
!-----------------------------------------------

!-----------------------------------------------
! initialize
      sourcelhs = zero
      sourcerhs = zero

      SELECT CASE(lEQ_NO)

! source for pressure correction equation: pp_g
      CASE (PRESSURE_CORRECTION)

! source for solids correctione quation:: epp
      CASE (SOLIDS_CORRECTION)

! source for gas continuity equation: rop_g
      CASE (GAS_CONTINUITY)

! source for solids continuity equation: rop_s
      CASE (SOLIDS_CONTINUITY)

! source for gas u momentum equation: U_g
      CASE (GAS_U_MOM)

! source for solids u momentum equation: U_s
      CASE (SOLIDS_U_MOM)

! source for gas v momentum equation: V_g
      CASE (GAS_V_MOM)

! source for solids v momentum equation: V_s
      CASE (SOLIDS_V_MOM)

! source for gas w momentum equation: W_g
      CASE (GAS_W_MOM)

! source for solids w momentum equation: W_s
      CASE (SOLIDS_W_MOM)

! source for gas temperature equation: T_g
      CASE (GAS_ENERGY)

! source for solids temperature equation: T_s
      CASE (SOLIDS_ENERGY)

! source for gas species equation: X_g
      CASE (GAS_SPECIES)

! source for solids species equation: X_s
      CASE (SOLIDS_SPECIES)

! source for granular energy equation: Theta_m
      CASE (GRAN_ENERGY)

! source for user scalar equations: scalar
      CASE (USR_SCALAR)

! source for k_epsilon turbulence equations: k_turb_g
      CASE (K_EPSILON_K)

! source for k_epsilon turbulence equations: e_turb_g
      CASE (K_EPSILON_E)

      END SELECT
      RETURN
      END SUBROUTINE USR_SOURCES



