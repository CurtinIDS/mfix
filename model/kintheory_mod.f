!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module: kintheory                                                   C
!  Purpose: Common block containing constants, variables, functions    C
!  used by various kinetic theory models                               C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Garzo, V., Tenneti, S., Subramaniam, S., and Hrenya, C. M.,         C
!      "Enskog kinetic theory for monodisperse gas-solid flows", JFM,  C
!      Vol. 712, 2012, pp. 129-168                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      MODULE kintheory


! coefficient terms needed for stress (in addition to
! mu_s and lambda_s which are defined in visc_s_mod, allocated
! in allocate_arrays and initialized in set_constprop; and
! P_s which is defined in fldvar_mod, allocated in
! allocate_arrays, and initialized in init_fvars)
!     stress term with gradient in particle M velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: MU_sM_ip
!     stress term with gradient in particle L velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: MU_sL_ip
!     stress term with trace in particle M velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: XI_sM_ip
!     stress term with trace in particle L velocity
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: XI_sL_ip


! coefficient terms needed for momentum source (in addition to
! F_SS which is defined in drag_mod, allocated in allocate_arrays
! and initialized in set_constprop)
!     momentum source term with gradient in number density
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: Fnu_s_ip
!     momentum source term with gradient in mixture temperature or
!     with the gradient in temperature of species M
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: FT_sM_ip
!     momentum source term with gradient in temperature of species L
      DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: FT_sL_ip


! coefficient terms needed for heat flux (in addition to
! kth_s which is defined in set_constprop, allocated
! allocate_arrays and initialized in init_fvars)
!     heat flux term with gradient in granular temperature of species L
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Kth_sL_ip
!     heat flux term with gradient in number density of species M
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Knu_sM_ip
!     heat flux term with gradient in number density of species L
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Knu_sL_ip
!     heat flux term with velocity difference
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Kvel_s_ip


! coefficient terms needed for energy dissipation
!     energy dissipation with difference in species granular
!     temperature: transfer between solid solid phases
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ED_ss_ip
!     energy dissipation term
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: EDT_s_ip
!     energy dissipation with divergence of velocity of species M
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: EDvel_sM_ip
!     energy dissipation with divergence of velocity of species L
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: EDvel_sL_ip
!     coefficient A2, xsi used in multiple places in GTSH theory
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: A2_gtsh
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xsi_gtsh

! Solids source terms needed for Iddir & Arastoopour (2005)
! kinetic theory model
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  KTMOM_U_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  KTMOM_V_s
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE ::  KTMOM_W_s

! parameter in the theory of GTSH that is related to length scale
! of lubrication effects. For details see GTSH, 2012.
      DOUBLE PRECISION, PARAMETER :: EpM = 0.01d0



      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the magnitude of the gas-solids relative         C
!  velocity at i, j, k                                                 C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION KT_RVEL(IJK, m)

! Modules
!---------------------------------------------------------------------//
      use fldvar, only: u_s, v_s, w_s
      use fldvar, only: u_g, v_g, w_g
      use run, only: shear
      use vshear, only: vsh
      use indices, only: i_of
      use functions, only: im_of, jm_of, km_of
      use functions, only: fluid_at
      use fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! ijk index
      INTEGER, INTENT(IN) :: ijk
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: I, IMJK, IJMK, IJKM
! Cell center value of solids and gas velocities
      DOUBLE PRECISION :: USCM, VSCM, WSCM, &
                          UGC, VGC, WGC
! y-component of velocity with 'shear' applied
      DOUBLE PRECISION :: vs_j, vs_jm
!---------------------------------------------------------------------//
      I = I_OF(IJK)
      IMJK  = IM_OF(IJK)
      IJMK  = JM_OF(IJK)
      IJKM  = KM_OF(IJK)

! Awkward here but captures what was done in early version of
! calc_mu_s. Otherwise any call to this routine must be done with
! 'shear' already applied to v_s.
      vs_j = V_S(IJK,M)
      vs_jm = V_S(IJMK,M)
      IF (SHEAR) THEN
         vs_j = V_S(IJK,M)+VSH(IJK)
         IF(FLUID_AT(IJMK)) THEN
            vs_jm = V_S(IJMK,M)+VSH(IJMK)
         ELSE
            vs_jm = v_s(IJMK,M)
         ENDIF
      ENDIF

! Start calculation for relative velocity
! Calculate velocity components at i, j, k
      UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
      VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
      WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

      USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
      VSCM = AVG_Y_N(vs_jm, vs_j)
      WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

! Magnitude of gas-solids relative velocity
      KT_RVEL = SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + &
                     (WGC - WSCM)**2)

      RETURN
      END FUNCTION KT_RVEL

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION KT_COS_THETA(ijk, M)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, small_number
      USE fldvar, only: u_g, v_g, w_g
      USE fldvar, only: u_s, v_s, w_s
      USE fldvar, only: ep_s
      USE run, only: shear
      USE vshear, only: vsh
      USE toleranc, only: zero_ep_s
      USE indices, only: i_of
      USE functions, only: im_of, jm_of, km_of
      USE functions, only: fluid_at
      USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! ijk index
      INTEGER, INTENT(IN) :: ijk
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: I, IMJK, IJMK, IJKM
! Cell center value of solids and gas velocities
      DOUBLE PRECISION :: USCM, VSCM, WSCM, &
                          UGC, VGC, WGC
      DOUBLE PRECISION :: vs_j, vs_jm, speed, rvel
!---------------------------------------------------------------------//

      I = I_OF(IJK)
      IMJK  = IM_OF(IJK)
      IJMK  = JM_OF(IJK)
      IJKM  = KM_OF(IJK)

! Awkward here but captures what was done in early version of
! calc_mu_s. Otherwise any call to this routine must be done with
! 'shear' already applied to v_s.
      vs_j = V_S(IJK,M)
      vs_jm = V_S(IJMK,M)
      IF (SHEAR) THEN
         vs_j = V_S(IJK,M)+VSH(IJK)
         IF(FLUID_AT(IJMK)) THEN
            vs_jm = V_S(IJMK,M)+VSH(IJMK)
         ELSE
            vs_jm = v_s(IJMK,M)
         ENDIF
      ENDIF

      UGC = AVG_X_E(U_G(IMJK),U_G(IJK),I)
      VGC = AVG_Y_N(V_G(IJMK),V_G(IJK))
      WGC = AVG_Z_T(W_G(IJKM),W_G(IJK))

      USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
      VSCM = AVG_Y_N(vs_jm, vs_j)
      WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

      RVEL =  SQRT((UGC - USCM)**2 + (VGC - VSCM)**2 + &
                   (WGC - WSCM)**2)

      SPEED = SQRT(USCM**2+VSCM**2+WSCM**2)

! motion viewed by the particles (crossing trajectory effect)
      IF(SPEED > Small_Number .AND. RVEL > Small_Number .AND. &
         EP_S(IJK,M) > ZERO_EP_S) THEN
         KT_Cos_Theta = ( (UGC-USCM)*USCM + (VGC-VSCM)*VSCM + &
                       (WGC-WSCM)*WSCM )/(RVEL*SPEED)
      ELSE
         KT_Cos_Theta = ZERO
      ENDIF

      RETURN
      END FUNCTION  KT_COS_THETA


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate a single particle drag coefficient/ep_s, which   C
!  is used in evaluting limiting values of certain granular kinetic    C
!  theory terms                                                        C
!                                                                      C
!  Comments: This is currently based on wen-yu but in the the future   C
!  we may want to make this consistent with drag_gs (that is, replace  C
!  this function call and assign appropriate variable within drag_gs   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION KT_DGA(IJK, m)

! Modules
!---------------------------------------------------------------------//
      use param1, only: zero, one, small_number, large_number
      use fldvar, only: rop_g
      use fldvar, only: d_p
      use physprop, only: mu_g
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! ijk index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! single particle drag coefficient, reynolds number
      DOUBLE PRECISION :: C_d, Re
! local value for relative velocity
      DOUBLE PRECISION :: rvel
!---------------------------------------------------------------------//

! initialization
      kt_dga = zero
      rvel = zero
      RVEL = KT_RVEL(IJK, M)

! Defining single particle drag coefficient dgA based on Wen-Yu
! correlation
      RE = D_p(IJK,M)*RVEL*ROP_G(IJK)/(MU_G(IJK) + SMALL_NUMBER)
      IF(RE .LE. 1000.d0)THEN
         C_d = (24.d0/(Re+SMALL_NUMBER)) * &
            (ONE + 0.15d0 * Re**0.687D0)
      ELSE
         C_d = 0.44d0
      ENDIF

! dga_s is local to this routine as is lrvel
      kt_dga = 0.75d0 * C_d * RVEL * ROP_g(IJK) / D_p(IJK,M)

! set value for 1st iteration and 1st time step
      IF(RVEL == ZERO) kt_dga = LARGE_NUMBER

      RETURN
      END FUNCTION KT_DGA

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_IA_ENERGY_DISSIPATION_SS                           C
!                                                                      C
!  Purpose: Implement kinetic theory of Iddir & Arastoopour (2005)     C
!     for calculation of source terms in granular energy equation      C
!                                                                      C
!  Author: Janine E. Galvin, Univeristy of Colorado                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_IA_ENERGY_DISSIPATION_SS(M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: pi
      USE constant, only: C_e
      USE fldvar, only: ro_s, rop_s, d_p
      USE fldvar, only: theta_m
      USE physprop, only: mmax
      USE rdf, only: g_0
      USE param1, only: zero

      USE functions, only: fluid_at, funlm
      USE compar, only: ijkstart3, ijkend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! Index
      INTEGER :: IJK
! Solids phase index
      INTEGER :: L
! Index for storing solids-solids drag coefficients
! in the upper triangle of the matrix
      INTEGER :: LM
! variables for IA theory
      DOUBLE PRECISION :: ED_common_term
      DOUBLE PRECISION :: EDvel_sL, EDvel_sM
      DOUBLE PRECISION :: M_PM, M_PL, MPSUM, NU_PL, NU_PM, D_PM, &
                          D_PL, DPSUMo2
      DOUBLE PRECISION :: Ap_lm, Dp_lm, R1p_lm, R10p_lm, R3p_lm, &
                          R4p_lm, R5p_lm, Bp_lm
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
          IF ( FLUID_AT(IJK) ) THEN

             D_PM = D_P(IJK,M)
             M_PM = (PI/6.d0)*(D_PM**3)*RO_S(IJK,M)
             NU_PM = ROP_S(IJK,M)/M_PM

             DO L = 1, MMAX
                LM = FUNLM(L,M)
                D_PL = D_P(IJK,L)
                M_PL = (PI/6.d0)*(D_PL**3)*RO_S(IJK,L)

                MPSUM = M_PM + M_PL
                DPSUMo2 = (D_PM+D_PL)/2.d0
                NU_PL = ROP_S(IJK,L)/M_PL

                ED_common_term = (3.d0/4.d0)*(DPSUMo2*DPSUMo2)*&
                   (1.d0+C_E)*G_0(IJK,M,L)*NU_PM*NU_PL*(M_PM*&
                   M_PL/MPSUM)*((M_PM*M_PL)**1.5)

                IF (M .eq. L) THEN
                   Ap_lm = MPSUM/(2.d0)
                   Dp_lm = M_PL*M_PM/(2.d0*MPSUM)
                   R1p_lm = 1.d0/( (Ap_lm**1.5)*(Dp_lm**3) )
                   R3p_lm = 1.d0/( (Ap_lm**1.5)*(Dp_lm**3.5) )

! Dissipation associated with the difference in temperature
! (e.g. interphase transfer term).  For the case case (L=M)
! the term ED_s * (TL-TM) cancels. Therefore, explicity set
! the term to zero for this case.
                   ED_ss_ip(IJK,LM) = ZERO

! Dissipation associated with temperature
                   EDT_s_ip(IJK,M,L) = -ED_common_term* (1.d0-C_E)*&
                      (M_PL/MPSUM)*(DSQRT(PI)/6.d0)*R1p_lm*&
                      (Theta_m(IJK,M)**1.5)

! Dissipation associated with divergence of
! velocity of solid phase L: do not explicity include
! terms that cancel when EDvel_sL is summed with EDvel_sM
                   EDvel_sL = ED_common_term*((1.d0-C_E)*(M_PL/MPSUM)*&
                      (DPSUMo2*PI/48.d0)*(M_PM*M_PL/MPSUM)*R3p_lm)
                      EDvel_sL_ip(IJK,M,L) = EDvel_sL*(Theta_m(IJK,L))

! Dissipation associated with divergence of
! velocity of solid phase M: do not explicity include
! terms that cancel when EDvel_sL is summed with EDvel_sM
! commented by sof, no need to re-do computation
                   !EDvel_sM = ED_common_term*((1.d0-C_E)*(M_PL/MPSUM)*&
                   !   (DPSUMo2*PI/48.d0)*(M_PM*M_PL/MPSUM)*R3p_lm)
                   !EDvel_sM_ip(IJK,M,L) = EDvel_sM*(Theta_m(IJK,M)) ! M=L

                   EDvel_sM_ip(IJK,M,L) =  EDvel_sL_ip(IJK,M,L)

                ELSE

                   Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*Theta_m(IJK,M))/&
                      2.d0
                   Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-&
                      Theta_m(IJK,M) ))/(2.d0*MPSUM)
                   Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+M_PL*&
                      Theta_m(IJK,L) ))/(2.d0*MPSUM*MPSUM)

                   R1p_lm = (1.d0/((Ap_lm**1.5)*(Dp_lm**3)))+ &
                      ((9.d0*Bp_lm*Bp_lm)/(Ap_lm**2.5 * Dp_lm**4))+&
                      ((30.d0*Bp_lm**4)/(2.d0*Ap_lm**3.5 * Dp_lm**5))

                   R3p_lm = (1.d0/((Ap_lm**1.5)*(Dp_lm**3.5)))+&
                      ((21.d0*Bp_lm*Bp_lm)/(2.d0 * Ap_lm**2.5 * Dp_lm**4.5))+&
                      ((315.d0*Bp_lm**4)/(8.d0 * Ap_lm**3.5 *Dp_lm**5.5))

                   R4p_lm = (3.d0/( Ap_lm**2.5 * Dp_lm**3.5))+&
                      ((35.d0*Bp_lm*Bp_lm)/(2.d0 * Ap_lm**3.5 * Dp_lm**4.5))+&
                      ((441.d0*Bp_lm**4)/(8.d0 * Ap_lm**4.5 * Dp_lm**5.5))

                   R5p_lm = (1.d0/(Ap_lm**2.5 * Dp_lm**3))+ &
                      ((5.d0*Bp_lm*Bp_lm)/(Ap_lm**3.5 * Dp_lm**4))+&
                      ((14.d0*Bp_lm**4)/( Ap_lm**4.5 * Dp_lm**5))

                   R10p_lm = (1.d0/(Ap_lm**2.5 * Dp_lm**2.5))+&
                      ((25.d0*Bp_lm*Bp_lm)/(2.d0* Ap_lm**3.5 * Dp_lm**3.5))+&
                      ((1225.d0*Bp_lm**4)/(24.d0* Ap_lm**4.5 * Dp_lm**4.5))

! Dissipation associated with the difference in temperature
! (e.g. interphase transfer term). to solved using PEA
                   ED_ss_ip(IJK,LM) = ED_common_term*DSQRT(PI)*(M_PM*M_PL/&
                      (2.d0*MPSUM))*R5p_lm*( (Theta_M(IJK,M)*Theta_M(IJK,L))**3 )

! Dissipation associated with temperature
                   EDT_s_ip(IJK,M,L) = -ED_common_term* (1.d0-C_E)*(M_PL/MPSUM)*&
                      (DSQRT(PI)/6.d0)*R1p_lm*( (Theta_m(IJK,M)*Theta_m(IJK,L))**3 )

! Dissipation associated with divergence of
! velocity of solid phase L
                   EDvel_sL = ED_common_term*( ((3.d0*DPSUMo2*PI/40.d0)*M_PL*&
                      R10p_lm)+( (DPSUMo2*PI/4.d0)*(M_PM*M_PL/MPSUM)*Bp_lm*&
                      R4p_lm)-( (1.d0-C_E)*(M_PL/MPSUM)*(DPSUMo2*PI/16.d0)*&
                      M_PL*Bp_lm*R4p_lm)+( (1.d0-C_E)*(M_PL/MPSUM)*&
                      (DPSUMo2*PI/48.d0)*(M_PM*M_PL/MPSUM)*R3p_lm))

                   EDvel_sL_ip(IJK,M,L) = EDvel_sL*( Theta_m(IJK,M)**3.5 *&
                      Theta_m(IJK,L)**2.5 )

! Dissipation associated with divergence of
! velocity of solid phase M
                   EDvel_sM = ED_common_term*( (-(3.d0*DPSUMo2*PI/40.d0)*M_PM*&
                      R10p_lm)+( (DPSUMo2*PI/4.d0)*(M_PM*M_PL/MPSUM)*Bp_lm*&
                      R4p_lm)+( (1.d0-C_E)*(M_PL/MPSUM)*(DPSUMo2*PI/16.d0)*&
                       M_PM*Bp_lm*R4p_lm)+( (1.d0-C_E)*(M_PL/MPSUM)*&
                      (DPSUMo2*PI/48.d0)*(M_PM*M_PL/MPSUM)*R3p_lm))

                   EDvel_sM_ip(IJK,M,L) = EDvel_sM*( Theta_m(IJK,M)**2.5 *&
                      Theta_m(IJK,L)**3.5 )

                ENDIF   ! end if/else (m.eq.l)

             ENDDO   ! end do (l=1,mmax)

         ENDIF   ! end if(fluid_at)
      ENDDO   ! end do ijk

      RETURN
      END SUBROUTINE CALC_IA_ENERGY_DISSIPATION_SS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_GD_99_ENERGY_DISSIPATION_SS                        C
!                                                                      C
!  Purpose: Implement kinetic theory of Garzo & Dufty (1999)           C
!  for calculation of source terms in granular energy equation         C
!                                                                      C
!  Author: Janine E. Galvin                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_GD_99_ENERGY_DISSIPATION_SS(M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: pi
      USE constant, only: C_e
      USE fldvar, only: rop_s, d_p
      USE fldvar, only: theta_m, ep_s
      USE rdf, only: g_0

      USE functions, only: fluid_at
      USE compar, only: ijkstart3, ijkend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
! variables for GD model
      DOUBLE PRECISION :: press_star, c_star, zeta0_star, &
                          nu_gamma_star, &
                          lambda_num, cd_num, zeta1
      DOUBLE PRECISION :: D_PM, EP_SM
      DOUBLE PRECISION :: nu0, Chi
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
          IF ( FLUID_AT(IJK) ) THEN

! Note: k_boltz = M_PM

! local aliases
             Chi = G_0(IJK,M,M)
             EP_SM = EP_s(IJK,M)
             D_PM = D_P(IJK,M)

!            nu0=p_k/eta0
             nu0 = (96.d0/5.d0)*(EP_SM/D_PM)*DSQRT(Theta_m(IJK,M)/PI)

             press_star = 1.d0 + 2.d0*(1.d0+C_E)*EP_SM*Chi

             c_star = 32.0d0*(1.0d0 - C_E)*(1.d0 - 2.0d0*C_E*C_E) &
                / (81.d0 - 17.d0*C_E + 30.d0*C_E*C_E*(1.0d0-C_E))

             zeta0_star = (5.d0/12.d0)*Chi*(1.d0 - C_E*C_E) &
                * (1.d0 + (3.d0/32.d0)*c_star)

             nu_gamma_star = ((1.d0+C_E)/48.d0)*Chi*(128.d0-96.d0*C_E + &
                15.d0*C_E*C_E-15.d0*C_E*C_E*C_E+ (c_star/64.d0) * (15.d0* &
                C_E*C_E*C_E-15.d0*C_E*C_E+498.d0*C_E-434.d0))

             lambda_num = (3.d0/8.d0)*( (1.d0-C_E)*(5.d0*C_E*C_E+4.d0*C_E-1.d0)+ &
                (c_star/12.d0)*(159.d0*C_E+3.d0*C_E*C_E-19.d0*C_E- &
                15.d0*C_E*C_E*C_E) ) * (1.d0+C_E)

! does not include factor of 1/nu0.  the factor 1/nu0 in cD will cancel
! with the factor of nu0 used to obtain zeta1 from zeta1_star
             cd_num = ( (4.d0/15.d0)*lambda_num*EP_SM*Chi + &
                (press_star-1.d0)*(2.d0/3.d0-C_E)*c_star ) / &
                ( 0.5d0*zeta0_star+nu_gamma_star + (5.d0*c_star/64.d0) * &
                (1.d0+(3.d0*c_star/64.d0))*Chi*(1.d0-C_E*C_E))

! does not include factor of 1/nu0 in the first term.  the factor 1/nu0
! will cancel with the factor of nu0 used to obtain zeta1 from zeta1_star
! zeta1 = nu0*zeta1_star
             zeta1 = -(1.d0-C_E)*(press_star-1.d0) + (5.d0/32.d0) * &
                (1.d0-C_E*C_E)*(1.d0+(3.d0*c_star/64.d0))*Chi*cd_num

! in the energy equation this term is multiplied by 3/2*n*kboltz*T
! leave multiplication of theta for source routine
             EDvel_sM_ip(IJK,M,M) = (3.d0/2.d0)*ROP_s(IJK,M)*zeta1

! in the energy equation this term is multiplied by 3/2*n*kboltz*T
! leave multiplication of theta for source routine
             EDT_s_ip(IJK,M,M) = (3.d0/2.d0)*ROP_s(IJK,M)*nu0*zeta0_star

         ENDIF   ! end if (fluid_at)
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      RETURN
      END SUBROUTINE CALC_GD_99_ENERGY_DISSIPATION_SS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_GTSH_ENERGY_DISSIPATION_SS                         C
!                                                                      C
!  Purpose: Implement kinetic theory of Garzo, Tenneti, Subramaniam    C
!  Hrenya (2012) for calculation of source terms in granular           C
!  energy equation                                                     C
!                                                                      C
!  Author: Sofiane Benyahia                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_GTSH_ENERGY_DISSIPATION_SS(M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: pi
      USE constant, only: C_e
      USE fldvar, only: ep_s
      USE fldvar, only: ro_g, rop_g
      USE fldvar, only: ro_s, d_p
      USE fldvar, only: theta_m
      USE physprop, only: mu_g
      USE rdf, only: g_0
      USE param1, only: zero, one, small_number

      USE functions, only: fluid_at
      USE compar, only: ijkstart3, ijkend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
! variables
      DOUBLE PRECISION :: D_PM, EP_SM, V_p, N_p, M_p
      DOUBLE PRECISION :: nu0, Chi
      DOUBLE PRECISION :: VREL
      DOUBLE PRECISION :: Re_m, Re_T
      DOUBLE PRECISION :: zeta_star, mu2_0, mu4_0, mu4_1
      DOUBLE PRECISION :: omega, nu_j, rho_10, rho_11

!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

! Local aliases
            Chi = G_0(IJK,M,M)
            EP_SM = EP_s(IJK,M)
            D_PM = D_P(IJK,M)
            V_p = pi*D_PM**3/6d0
            n_p = EP_SM/V_p
            M_p = V_p * ro_s(ijk,m)

            nu0 = (96.d0/5.d0)*(EP_SM/D_PM)*DSQRT(Theta_m(IJK,M)/PI)

! First calculate the Re number (Re_m, Re_T) in eq. 3.1 in GTSH theory
            VREL = KT_RVEL(IJK, M)
! Note: rop_g = ro_g * (1-phi)
            Re_m = D_PM*VREL*ROP_g(ijk)/Mu_g(ijk)
            Re_T = ro_g(ijk)*D_PM*dsqrt(theta_m(ijk,m)) / mu_g(ijk)

! Now calculate and store xsi (used in many spots) in eq. 8.2 in
! GTSH theory.
            xsi_gtsh(ijk) = one/6d0*D_PM*VREL**2*(3d0*pi*mu_g(ijk)*&
               D_PM/M_p)**2 / dsqrt(pi*theta_m(ijk,m)) * &
               S_star(EP_SM, Chi)

! eq. (6.22) GTSH theory
            mu2_0 = dsqrt(2d0*pi) * Chi * (one-C_E**2)

! eq. (6.23) GTSH theory
            mu4_0 = (4.5d0+C_E**2) * mu2_0

! eq. (6.24) GTSH theory
!            mu4_1 = (6.46875d0+0.3125d0*C_E**2 + 2d0/(one-C_E)) * mu2_0
! this is done to avoid /0 in case c_e = 1.0
            mu4_1 = (6.46875d0+0.9375d0*C_E**2)*mu2_0 + &
               2d0*dsqrt(2d0*pi)* Chi*(one+C_E)

            A2_gtsh(ijk) = zero ! for EP_SM = zero
            if(EP_SM> small_number) then ! avoid singularity
! Now calculate zeta_star in eq. 8.10 in GTSH theory
               zeta_star = 4.5d0*dsqrt(2d0*Pi)*&
                  (ro_g(ijk)/ro_s(ijk,m))**2*Re_m**2 * &
                  S_star(EP_SM,Chi) / (EP_SM*(one-EP_SM)**2 * Re_T**4)

! Now calculate important parameter A2. This is used in many transport
! coefficients so that it's best to store it in an array instead of
! having many function calls.
               A2_gtsh(ijk) = (5d0*mu2_0 - mu4_0) / (mu4_1 - 5d0* &
                  (19d0/16d0*mu2_0 - 1.5d0*zeta_star))
            endif

! Now calculate first term rho_0 in cooling rate rho, eq (6.26) in GTSH
! theory.
! note that theta_(ijk,m) does't contain mass m, so that T/m = theta_m.
! note that edt_s_ip will need to be multiplied by 3/2 rop_s(ijk,m) in
! source_granular_energy to have the same meaning as in GD_99 theory.
            EDT_s_ip(ijk,M,M) = 4d0/3d0*dsqrt(pi)*(one-C_E**2)*Chi* &
               (one+0.1875d0*A2_gtsh(ijk))*n_p*D_PM**2*dsqrt(theta_m(ijk,m))

! Calculate the second term in cooling rate, eq (7.25) in GTSH theory.
! First calculate eq (7.28, 7.29) of GTSH theory.
            omega = (one+C_E)*nu0*((one-C_E**2)*(5d0*C_E-one) - &
               A2_gtsh(ijk)/ 6d0 * (15d0*C_E**3-3d0*C_E**2+81d0*C_E-61d0))
            nu_j = (one+C_E)/192d0*Chi*nu0* &
               (241d0-177d0*C_E+30d0*C_E**2-30d0*C_E**3)

! Now calculate eq (7.26, 7.27) of GTSH theory. ! corrected by W. Fullmer
            rho_10 = 2d0*Chi*EP_SM*(C_E**2-one)
            rho_11 = 25d0/1024d0*EP_SM*Chi**2*(one-C_E**2)* &
               (one+3d0/128d0*A2_gtsh(ijk)) * (omega/10d0 - &
               (one+C_E)*nu0*(one/3d0-C_E)*A2_gtsh(ijk)/2d0) / &
               (nu_j + G_gtsh(EP_SM, Chi, IJK, M)/m_p + 1.5d0* &
               xsi_gtsh(ijk)/theta_m(ijk,m) -1.5d0*EDT_s_ip(ijk,M,M))

! note that EDvel_sM_ip will later need to be multiplied by 3/2 rop_s(ijk,m).
            EDvel_sM_ip(IJK,M,M) = rho_10 + rho_11

         ENDIF   ! end if (fluid_at)
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      RETURN
      END SUBROUTINE CALC_GTSH_ENERGY_DISSIPATION_SS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Functions: G_tsh, K_phi, R_d, S_star                                C
!                                                                      C
!  Purpose: Implement kinetic theory of Garzo, Tenneti, Subramaniam    C
!  Hrenya (2012) for calculation of source terms in granular           C
!  energy equation                                                     C
!                                                                      C
!  Author: Sofiane Benyahia                                            C
!                                                                      C
!  Comments:                                                           C
!     these functions are needed only for gtsh theory                  C
!     Function gamma, eq. (8.1) in GTSH theory                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION G_gtsh (EP_SM, Chi, IJK, M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: pi
      USE param1, only: one
      USE physprop, only: mu_g
      USE fldvar, only: ro_g
      USE fldvar, only: d_p, theta_m
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids volume fraction of phase m at ijk
      DOUBLE PRECISION, INTENT(IN) :: EP_SM
! radial distribution function of phase M at ijk
      DOUBLE PRECISION, INTENT(IN) :: Chi
! Index
      INTEGER, INTENT(IN) :: IJK
! Solids phase index. Note M should be equal to 1 since theory
! valid for only mmax = 1.
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: Re_T
      DOUBLE PRECISION :: Rdiss, RdissP

!---------------------------------------------------------------------//

      if(EP_SM <= 0.1d0) then
         RdissP = one+3d0*dsqrt(EP_SM/2d0)
      else
         RdissP = &
         one + 3d0*dsqrt(EP_SM/2d0) + 135d0/64d0*EP_SM*dlog(EP_SM) + &
         11.26d0*EP_SM*(one-5.1d0*EP_SM+16.57d0*EP_SM**2-21.77d0*    &
         EP_SM**3) - EP_SM*Chi*dlog(epM)
      endif

      Re_T = ro_g(ijk)*d_p(ijk,m)*dsqrt(theta_m(ijk,m)) / mu_g(ijk)

      Rdiss = RdissP + Re_T * K_phi(EP_SM)

      G_gtsh = 3d0*Pi*mu_g(ijk)*d_p(ijk,m)*Rdiss  ! eq. (8.1) in GTSH theory

      RETURN
      END FUNCTION G_gtsh

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Comments:                                                           C
!     Function K(phi), eq. (2.7) in Yi/Zenk/Mitrano/Hrenya JFM (2013)  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION K_phi (phi)
      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! solids volume fraction
      DOUBLE PRECISION, INTENT(IN) :: phi
!---------------------------------------------------------------------//

      K_phi = (0.096d0 + 0.142d0*phi**0.212d0) / (1d0-phi)**4.454d0
      K_phi = 0.0d0 ! set to zero for compatibility with GTSH JFM (2012)
      RETURN
      END FUNCTION K_phi

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Comments:                                                           C
!     Function R_d, eq. (8.6) in GTSH theory                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION R_d (phi)
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids volume fraction
      DOUBLE PRECISION, INTENT(IN) :: phi
!---------------------------------------------------------------------//

      R_d = 1.0d0  ! this avoids singularity at phi = 0.0
      if((phi > 1d-15) .and. (phi <= 0.4d0)) then
        R_d = (1d0+3d0*dsqrt(phi/2d0)+135d0/64d0*phi*dlog(phi)+17.14d0*phi) / &
              (1d0+0.681d0*phi-8.48*phi**2+8.16d0*phi**3)
      elseif(phi > 0.4d0) then
        R_d = 10d0*phi/(1d0-phi)**3 + 0.7d0
      endif
      RETURN
      END FUNCTION R_d

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Comments:                                                           C
!     Function S_Star(phi), eq. (8.3, 8.5) in GTSH theory              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      DOUBLE PRECISION FUNCTION S_star (phi, Chi)
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids volume fraction
      DOUBLE PRECISION, INTENT(IN) :: phi
! radial distribution function
      DOUBLE PRECISION, INTENT(IN) :: Chi

!---------------------------------------------------------------------//

      S_star = 1.0d0
      if(phi >= 0.1d0) &
        S_star = R_d(phi)**2/(Chi*(1d0+3.5d0*dsqrt(phi)+5.9*phi))
      RETURN
      END FUNCTION S_Star

      END MODULE kintheory
