!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MU_s                                               C
!  Purpose: Calculate the vicosity, second viscosity, and pressure     C
!  associated with the Mth 'solids' phase. In the 'default' or         C
!  undefined case closure for 'granular conductivity' is also          C
!  calculated.                                                         C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MU_s(M)

! Modules
!---------------------------------------------------------------------//
! some constants
      USE param1, only: undefined
! specified constant solids phase viscosity
      USE physprop, only: mu_s0
! invoke user defined quantity
      USE usr_prop, only: usr_mus, calc_usr_prop
      USE usr_prop, only: solids_viscosity
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

!---------------------------------------------------------------------//

      IF (USR_MUS(M)) THEN
         CALL CALC_USR_PROP(Solids_Viscosity,lm=M)
      ELSEIF (MU_S0(M) == UNDEFINED) THEN
! this is a necessary check as one may have mu_s0 defined but still 
! need to call this routine (for set_epmus)
         CALL CALC_DEFAULT_MUs(M)
      ENDIF

      CALL SET_EPMUS_VALUES(M)
      RETURN
      END SUBROUTINE CALC_MU_s


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_EPMUS_VALUES                                        C
!  Purpose: This routine sets the internal variables epmu_s and        C
!  eplambda_s that are used in the stress calculations. If the         C
!  keyword Ishii is invoked then these quantities represent the        C
!  viscosity and second viscosity multiplied by the volume fraction    C
!  otherwise they are simply viscosity/second viscosity (i.e. are      C
!  multipled by one).                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SET_EPMUS_VALUES(M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE fldvar, only: eps_ifac
      USE functions, only: fluid_at
      USE param1, only: zero
      USE visc_s, only: mu_s, epmu_s, lambda_s, eplambda_s
      USE mms, only: use_mms

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell index
      INTEGER :: IJK
!---------------------------------------------------------------------//

!      EPMU_S(:,M) = EPS_IFAC(:,M)*MU_s(:,M)
!      EPLAMBDA_S(:,M) = EPS_IFAC(:,M)*LAMBDA_S(:,M)

      DO IJK=ijkstart3,ijkend3
         IF(FLUID_AT(IJK) .OR. USE_MMS) THEN
! if ishii then multiply by volume fraction otherwise multiply by 1
            EPMU_S(IJK,M) = EPS_IFAC(IJK,M)*Mu_s(IJK,M)
            EPLAMBDA_S(IJK,M) = EPS_IFAC(IJK,M)*Lambda_s(IJK,M)
         ELSE
            EPMU_S(IJK,M) = ZERO
            EPLAMBDA_S(IJK,M) = ZERO
         ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SET_EPMUS_VALUES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine CALC_DEFAULT_MUS                                         C
!  Purpose: By default MFIX assumes each Mth phase is a granular       C
!  solids material. Accordingly, various theories are invoked to       C
!  provide closures for the solids stress term. These include          C
!  different kinetic theory and soil mechanic models, which may be     C
!  specified by the user through various keyword selections.           C
!  The granular stress terms (granular viscosity, second/bulk          C
!  viscosity, solids pressure), and if required, a granular            C
!  conductivity term are closed here.                                  C
!                                                                      C
!  Comments:                                                           C
!     GRANULAR_ENERGY = .FALSE.                                        C
!        EP_g < EP_star   -->    friction_schaeffer                    C
!        EP_g >= EP_star  -->    viscous (algebraic)                   C
!                                                                      C
!     GRANULAR_ENERGY = .TRUE.                                         C
!        FRICTION = .TRUE.                                             C
!           EP_s(IJK,M) > EPS_f_min  -->  friction + viscous(pde)      C
!           EP_s(IJK,M) < EP_f_min   -->  viscous (pde)                C
!        FRICTION = .FALSE.                                            C
!           EP_g < EP_star  -->  friction_schaeffer + viscous(pde)     C
!           EP_g >= EP_star -->  viscous (pde)                         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_DEFAULT_MUS(M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
!
      USE error_manager, only: err_msg, init_err_msg, finl_err_msg
      USE error_manager, only: flush_err_msg
! solids pressure
      USE fldvar, only: p_s, p_s_c, p_s_v
      USE fldvar, only: p_s_f
! some constants
      USE param1, only: one
! runtime flag indicates whether the solids phase becomes close-packed
! at ep_star
      USE physprop, only: close_packed
! number of solids phases
      USE physprop, only: mmax
      USE physprop, only: blend_function

! runtime flag to use qmomk
      USE qmom_kinetic_equation, only: qmomk

! runtime flag to solve partial differential granular energy eqn(s)
      USE run, only: granular_energy
! kinetic theories
      USE run, only: kt_type_enum
      USE run, only: lun_1984
      USE run, only: simonin_1996
      USE run, only: ahmadi_1995
      USE run, only: gd_1999
      USE run, only: gtsh_2012
      USE run, only: ia_2005
      USE run, only: ghd_2007
      USE run, only: kt_type
! filtered subgrid model
      USE run, only: subgrid_type_enum
      USE run, only: milioli, igci, undefined_subgrid_type
! frictional theories
      USE run, only: friction, schaeffer
! runtime flag for blending stress (only with schaeffer)
      USE run, only: blending_stress
! solids transport coefficients
      USE visc_s, only: mu_s, epmu_s, mu_s_c, mu_s_v
      USE visc_s, only: mu_s_p, mu_s_f
      USE visc_s, only: lambda_s, eplambda_s, lambda_s_c, lambda_s_v
      USE visc_s, only: lambda_s_p, lambda_s_f
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell index
      INTEGER :: IJK
! blend factor
      DOUBLE PRECISION :: BLEND

! General initializations
!---------------------------------------------------------------------
! Zero out various quantities
      CALL INIT0_MU_S(M)

! Calculate quantities that are functions of velocity only and may be
! needed for subsequent calculations by various model options
      CALL INIT1_MU_S(M)


! Viscous-flow stress tensor
!---------------------------------------------------------------------
      IF (.NOT. QMOMK) THEN
! if QMOMK then do not solve algebraic or PDE form of granular
! temperature governing equation
         IF(.NOT.GRANULAR_ENERGY) THEN
            IF(SUBGRID_TYPE_ENUM .ne. UNDEFINED_SUBGRID_TYPE) THEN
               IF (SUBGRID_TYPE_ENUM .EQ. IGCI) THEN
                  CALL SUBGRID_STRESS_IGCI(M)
               ELSEIF (SUBGRID_TYPE_ENUM .EQ. MILIOLI) THEN
                  CALL SUBGRID_STRESS_MILIOLI(M)
               ENDIF
            ELSE
              call gt_algebraic(M)   ! algebraic granular energy equation
            ENDIF
         ELSE   ! granular energy transport equation
            SELECT CASE(KT_TYPE_ENUM)
               CASE (IA_2005)   ! polydisperse theory
                  CALL gt_pde_ia(M)
               CASE (GD_1999)   ! strictly monodisperse theory
                  CALL gt_pde_gd(M)
               CASE (GTSH_2012)   ! strictly monodisperse theory
                  CALL gt_pde_gtsh(M)
               CASE (GHD_2007)   ! polydisperse GHD theory for mixture temperature
                  CALL TRANSPORT_COEFF_GHD(M)
               CASE(LUN_1984)   ! monodisperse/ad-hoc polydisperse theory
                  CALL gt_pde_lun(M)
               CASE(SIMONIN_1996)   ! monodisperse/ad-hoc polydisperse theory
                  CALL gt_pde_simonin(M)
               CASE (AHMADI_1995)   ! monodisperse/ad-hoc polydisperse theory
                  CALL gt_pde_ahmadi(M)
               CASE DEFAULT
! should never hit this
                  CALL INIT_ERR_MSG("CALC_MU_S")
                  WRITE(ERR_MSG, 1100) KT_TYPE
 1100 FORMAT('ERROR 1100: Invalid KT_TYPE: ', A,' The check_data ',&
         'routines should',/, 'have already caught this error and ',&
         'prevented the simulation from ',/,'running. Please notify ',&
         'the MFIX developers.')
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
                  CALL FINL_ERR_MSG
            END SELECT  ! end selection of kt_type_enum
         ENDIF
      ENDIF


! Frictional stress tensors
! Note that only one of these can be used at this time
!---------------------------------------------------------------------
! Schaeffer's frictional formulation
      IF (SCHAEFFER .AND. CLOSE_PACKED(M)) call friction_schaeffer(M)
! Princeton's frictional implementation
      IF (FRICTION .AND. CLOSE_PACKED(M)) call friction_princeton(M)


! Assign viscosity, second viscosity and pressure in mth solids phase.
! Note that a plastic pressure component is calculated in a separate
! routine (see calc_p_star) which is then directly incorporated into
! the solids momentum equations (see source_u_s, source_v_s and
! source_w_s).
!---------------------------------------------------------------------
      Mu_s_c(:,M) = Mu_s_v(:)
      LAMBDA_s_c(:,M)= Lambda_s_v(:)
      P_s_c(:,M) = P_s_v(:)

      IF(BLENDING_STRESS) THEN
! blend plastic & viscous stresses (not available with friction)
         DO IJK = ijkstart3, ijkend3
            blend =  blend_function(IJK)
            Mu_s(IJK,M) = (ONE-blend)*Mu_s_p(IJK) &
                + blend*Mu_s_v(IJK)
! is there any point in blending lambda_s_p since it is never
! assigned? the only thing this would effectively do is reduce
! lambda_s_v by the blend value
            LAMBDA_s(IJK,M) = (ONE-blend)*Lambda_s_p(IJK) + &
               blend*Lambda_s_v(IJK)
! note that p_star is the plastic pressure term. p_star is calculated
! and blended in its own routine to accomodate slightly different
! timing of the calculation to coincide with updated volume fractions.
! however, we must factor P_s_v by the blend value here.
           P_s(IJK,M) = blend*P_s_v(IJK)
         ENDDO
      ELSE
         Mu_s(:,M) = Mu_s_p(:) + Mu_s_v(:) + Mu_s_f(:)
         LAMBDA_s(:,M) = Lambda_s_p(:) + Lambda_s_v(:) + Lambda_s_f(:)
         P_s(:,M) = P_s_v(:) + P_s_f(:)
      ENDIF  ! end if/else (blending_stress)

      RETURN
      END SUBROUTINE CALC_DEFAULT_MUs



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: gt_algebraic                                            C
!  Purpose: solve an algebraic form of the granular energy equation    C
!  based on Lun et al. (1984)                                          C
!  Obtained from the granular energy equation of Lun et al. (1984) by  C
!  assuming that the granular energy is dissipated locally so that     C
!  diffusion and convection contributions are neglected and only       C
!  generation and dissipation terms are retained.                      C
!                                                                      C
!  References:                                                         C
!  Syamlal, M., 1987c, "A Review of Granular Stress Constitutive       C
!     Relations," Topical Report, DOE/MC/21353-2372, NTIS/DE87006499,  C
!     National Technical Information Service, Springfield, VA.         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine Gt_algebraic (M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE constant, only: C_e, sqrt_pi
      USE constant, only: v_ex
      USE fldvar, only: p_s_v
      USE fldvar, only: ep_g, ep_s
      USE fldvar, only: d_p, ro_s
      USE fldvar, only: theta_m
      USE functions, only: fluid_at
      USE param1, only: zero, half, one, small_number
      USE rdf, only: g_0
      USE trace, only: trd_s_c, trd_s2
      USE visc_s, only: mu_s_v, lambda_s_v
      use visc_s, only: ep_g_blend_start
      USE visc_s, only: alpha_s
      USE visc_s, only: trm_s, trdM_s
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell index
      INTEGER :: IJK
! Coefficients of quadratic equation
      DOUBLE PRECISION :: aq, bq, cq
! Constant in equation for mth solids phase pressure
      DOUBLE PRECISION :: K_1m
! Constant in equation for mth solids phase bulk viscosity
      DOUBLE PRECISION :: K_2m
! Constant in equation for mth solids phase viscosity
      DOUBLE PRECISION :: K_3m
! Constant in equation for mth solids phase dissipation
      DOUBLE PRECISION :: K_4m
! Constant in Boyle-Massoudi stress term
      DOUBLE PRECISION :: K_5m
! Value of EP_s * SQRT( THETA )for Mth solids phase continuum
      DOUBLE PRECISION :: EP_sxSQRTHETA
! Value of EP_s * EP_s * THETA for Mth solids phase continuum
      DOUBLE PRECISION :: EP_s2xTHETA, temp_local
!---------------------------------------------------------------------//

!!$omp parallel do default(shared)                                    &
!!$omp private(IJK, K_1m, K_2m, K_3m, K_4m, K_5m, temp_local,         &
!!$omp         aq, bq, cq, EP_sxSQRTHETA, EP_s2xTHETA)
       DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN
            IF(EP_g(IJK) .GE. EP_g_blend_start(IJK)) THEN

! Calculate K_1m, K_2m, K_3m, K_4m
               K_1m = 2.D0 * (ONE + C_e) * RO_S(IJK,M) * G_0(IJK, M, M)
               K_3m = HALF * D_p(IJK,M) * RO_S(IJK,M) * ( ( (SQRT_PI / &
                  (3.D0 * (3.D0 - C_e))) * ( HALF*(3d0*C_e+ONE) + &
                  0.4D0*(ONE + C_e)*(3.D0*C_e - ONE) * EP_s(IJK,M)* &
                  G_0(IJK, M,M)) ) + 8.D0*EP_s(IJK,M)*G_0(IJK, M,M)*&
                  (ONE + C_e)/ (5.D0*SQRT_PI) )
               K_2m = 4.D0 * D_p(IJK,M) * RO_S(IJK,M) * (ONE + C_e) *&
                  EP_s(IJK,M) * G_0(IJK, M,M) / (3.D0 * SQRT_PI) - &
                  2.D0/3.D0 * K_3m
               K_4m = 12.D0 * (ONE - C_e*C_e) *&
                  RO_S(IJK,M) * G_0(IJK, M,M) / (D_p(IJK,M) * SQRT_PI)
               aq = K_4m*EP_s(IJK,M)
               bq = K_1m*EP_s(IJK,M)*trD_s_C(IJK,M)
               cq = -(K_2m*trD_s_C(IJK,M)*trD_s_C(IJK,M)&
                      + 2.D0*K_3m*trD_s2(IJK,M))

! Boyle-Massoudi Stress term
               IF(V_ex .NE. ZERO) THEN
                  K_5m = 0.4 * (ONE + C_e) * G_0(IJK, M,M) * &
                     RO_S(IJK,M) * ( (V_ex * D_p(IJK,M)) / &
                     (ONE - EP_s(IJK,M) * V_ex) )**2
                  bq = bq + EP_s(IJK,M) * K_5m * &
                     (trM_s(IJK) + 2.D0 * trDM_s(IJK))
               ELSE
                  K_5m = ZERO
               ENDIF

! Calculate EP_sxSQRTHETA and EP_s2xTHETA
               temp_local = bq**2 - 4.D0 * aq * cq
               EP_sxSQRTHETA = (-bq + SQRT(temp_local))&
                  / ( 2.D0 * K_4m )
               EP_s2xTHETA = EP_sxSQRTHETA * EP_sxSQRTHETA

               IF(EP_s(IJK,M) > SMALL_NUMBER)THEN
! Find pseudo-thermal temperature in the Mth solids phase
                  THETA_m(IJK,M) = EP_s2xTHETA/(EP_s(IJK,M)*EP_s(IJK,M))
               ELSE
                  THETA_m(IJK,M) = ZERO
               ENDIF

! Find pressure in the Mth solids phase
               P_s_v(IJK) = K_1m * EP_s2xTHETA

! second viscosity in Mth solids phase
               LAMBDA_s_v(IJK) = K_2m * EP_sxSQRTHETA

! shear viscosity in Mth solids phase
               MU_s_v(IJK) = K_3m * EP_sxSQRTHETA

! Boyle-Massoudi stress coefficient
               ALPHA_s(IJK, M) = -K_5m * EP_s2xTHETA
            ENDIF

         ENDIF   ! Fluid_at
      ENDDO
!!$omp end parallel do

      RETURN
      END SUBROUTINE GT_ALGEBRAIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE_LUN                                              C
!  Purpose: Calculate granular stress terms (viscosity, bulk viscosity C
!  solids pressure) & granular conductivity                            C
!                                                                      C
!  Author: Kapil Agrawal, Princeton University      Date: 6-FEB-98     C
!                                                                      C
!  Literature/Document References:                                     C
!  Lun, C.K.K., S.B. Savage, D.J. Jeffrey, and N. Chepurniy,           C
!     Kinetic theories for granular flow - inelastic particles in      C
!     Couette-flow and slightly inelastic particles in a general       C
!     flow field. Journal of Fluid Mechanics, 1984. 140(MAR):          C
!     p. 223-256                                                       C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine gt_pde_lun (M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: switch
      USE constant, only: alpha
      USE constant, only: pi
      USE constant, only: eta

      USE drag, only: f_gs

      USE fldvar, only: ro_g
      USE fldvar, only: ep_s
      USE fldvar, only: d_p, rop_s, ro_s
      USE fldvar, only: theta_m
      USE fldvar, only: p_s_v

      USE kintheory, only: kt_dga

      USE param1, only: zero, one, small_number

      USE physprop, only: smax
      USE physprop, only: Kth_s, kphi_s

      USE rdf, only: g_0, DG_0DNU

      USE toleranc, only: dil_ep_s

      USE visc_s, only: mu_s_v, mu_b_v, lambda_s_v

      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: IJK
! solids phase index
      INTEGER :: MM
! use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION :: Mu_star, Kth, Kth_star
! sum of ep_s * g_0
      DOUBLE PRECISION :: SUM_EpsGo
! single particle drag coefficient/ep_s
      DOUBLE PRECISION :: dga_sm
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

! The following is purely an ad-hoc modification so that the underlying
! monodisperse theory can be used for polydisperse systems in a
! consistent manner.  That is, solids pressure, viscosity and
! conductivity must be additive. The non-linear terms (eps^2) are
! corrected so the stresses of two or more identical solids phases are
! equal to those of a equivalent single solids phase. sof June 15 2005.
            SUM_EpsGo = ZERO
            DO MM = 1, SMAX
               SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,MM)*G_0(IJK,M,MM)
            ENDDO

            IF(EP_S(IJK,M) < DIL_EP_S) &
               dga_sM = KT_DGA(ijk, M)

! Pressure
            P_s_v(IJK) = ROP_s(IJK,M)*(1.d0 + 4.d0 * Eta *&
                SUM_EpsGo)*Theta_m(IJK,M)

! Bulk and shear viscosity
            Mu_s_v(IJK) = (5.d0*DSQRT(Pi*Theta_m(IJK,M))*D_p(IJK,M)*&
               RO_S(IJK,M))/96.d0
            Mu_b_v(IJK) = (256.d0*Mu_s_v(IJK)*EP_s(IJK,M)*SUM_EpsGo)&
                /(5.d0*Pi)

! added Ro_g = 0 for granular flows (no gas). sof Aug-02-2005
            IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN
               Mu_star = Mu_s_v(IJK)
            ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
               Mu_star = ZERO
            ELSEIF(EP_S(IJK,M) < DIL_EP_S) THEN
               Mu_star = RO_S(IJK,M)*EP_s(IJK,M)* G_0(IJK,M,M)*&
                   Theta_m(IJK,M)* Mu_s_v(IJK)/ &
                   (RO_S(IJK,M)*SUM_EpsGo*Theta_m(IJK,M) &
                   + 2.d0*DgA_SM/RO_S(IJK,M)* Mu_s_v(IJK))
            ELSE
               Mu_star = RO_S(IJK,M)*EP_s(IJK,M)* G_0(IJK,M,M)*&
                   Theta_m(IJK,M)*Mu_s_v(IJK)/ &
                   (RO_S(IJK,M)*SUM_EpsGo*Theta_m(IJK,M)+ &
                   (2.d0*F_gs(IJK,M)*Mu_s_v(IJK)/&
                   (RO_S(IJK,M)*EP_s(IJK,M) )) )
            ENDIF

! Shear viscosity
            Mu_s_v(IJK) =&
                ((2.d0+ALPHA)/3.d0)*((Mu_star/(Eta*(2.d0-Eta)*&
                G_0(IJK,M,M)))*(ONE+1.6d0*Eta*SUM_EpsGo)&
                *(ONE+1.6d0*Eta*(3d0*Eta-2d0)*&
                SUM_EpsGo)+(0.6d0*Mu_b_v(IJK)*Eta))

! Second viscosity as defined in MFIX
            LAMBDA_S_V(IJK) = Eta*Mu_b_v(IJK) - (2d0*Mu_s_v(IJK))/3d0

            Kth=75d0*RO_S(IJK,M)*D_p(IJK,M)*DSQRT(Pi*Theta_m(IJK,M))/&
                (48d0*Eta*(41d0-33d0*Eta))

            IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN
               Kth_star=Kth
            ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
               Kth_star = ZERO
            ELSEIF(EP_S(IJK,M) < DIL_EP_S) THEN
               Kth_star = RO_S(IJK,M)*EP_s(IJK,M)* &
                  G_0(IJK,M,M)*Theta_m(IJK,M)* Kth/ &
                  (RO_S(IJK,M)*SUM_EpsGo*Theta_m(IJK,M) &
                  + 1.2d0*DgA_sM/RO_S(IJK,M)* Kth)
            ELSE
               Kth_star = RO_S(IJK,M)*EP_s(IJK,M)* &
                  G_0(IJK,M,M)*Theta_m(IJK,M)*Kth/ &
                  (RO_S(IJK,M)*SUM_EpsGo*Theta_m(IJK,M)+ &
                  (1.2d0*F_gs(IJK,M)*Kth/(RO_S(IJK,M)*EP_s(IJK,M))) )
            ENDIF

! Granular conductivity
            Kth_s(IJK,M) = Kth_star/G_0(IJK,M,M)*(&
               ( ONE+(12d0/5.d0)*Eta*SUM_EpsGo ) * &
               ( ONE+(12d0/5.d0)*Eta*Eta*(4d0*Eta-3d0)* SUM_EpsGo ) + &
               (64d0/(25d0*Pi))*(41d0-33d0*Eta)*(Eta*SUM_EpsGo)**2 )



! Dufour coefficient: granular 'conductivity' in the Mth solids phase
! associated with gradient in volume fraction
!--------------------------------------------------------------------
!     Kphi_s has been set to zero.  To activate the feature uncomment the
!     following lines and also the lines in source_granular_energy.
            Kphi_s(IJK,M) = ZERO
!     &            (Kth_star/(G_0(IJK,M,M)))*(12d0/5.)*Eta*(Eta-1.)*
!     &            (2.*Eta-1.)*(1.+(12d0/5.)*Eta*EP_s(IJK,M)*
!     &            G_0(IJK,M,M))*(EP_s(IJK,M)*
!     &            DG_0DNU(EP_s(IJK,M))
!     &            + 2*G_0(IJK,M,M))*Theta_m(IJK,M)
!--------------------------------------------------------------------

         ENDIF   ! Fluid_at
      ENDDO   ! IJK loop
      RETURN
      END SUBROUTINE GT_PDE_LUN


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE_ahmadi                                           C
!  Purpose: Calculate granular stress terms (viscosity, bulk viscosity C
!  solids pressure) & granular conductivity using ahmadi model         C
!                                                                      C
!  Author: Sofiane Benyahia, Fluent Inc.            Date: 02-01-05     C
!                                                                      C
!  Literature/Document References:                                     C
!  Cao, J. and Ahmadi, G., 1995, Gas-particle two-phase turbulent      C
!     flow in a vertical duct. Int. J. Multiphase Flow, vol. 21,       C
!     No. 6, pp. 1203-1228.                                            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine gt_pde_ahmadi (M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: C_e, eta

      USE drag, only: f_gs

      USE fldvar, only: p_s_v
      USE fldvar, only: ROP_s, ro_s, d_p, theta_m
      USE fldvar, only: ep_s
      USE fldvar, only: k_turb_g, e_turb_g
      USE kintheory, only: kt_dga

      USE param1, only: zero, half, one, small_number

      USE physprop, only: smax
      USE physprop, only: kth_s

      USE rdf, only: g_0

      USE toleranc, only: dil_ep_s

      USE turb, only: tau_1

      USE visc_s, only: lambda_s_v, mu_s_v, mu_b_v
      USE visc_s, only: ep_star_array

      USE functions, only: fluid_at
      USE compar, only: ijkstart3, ijkend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local parameters
!---------------------------------------------------------------------//
! constant in simonin model
      DOUBLE PRECISION, PARAMETER :: C_mu = 9.0D-02

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: IJK
! solids phase index
      INTEGER :: MM, L
! defining parameters for Ahmadi
      DOUBLE PRECISION :: Tau_12_st
      DOUBLE PRECISION :: Tmp_Ahmadi_Const
! sum of ep_s * g_0
      DOUBLE PRECISION :: SUM_EpsGo
! single particle drag coefficient/ep_s
      DOUBLE PRECISION :: dga_sl
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

! Define time scales and constants

! time scale of turbulent eddies
            Tau_1(ijk) = 3.d0/2.d0*C_MU*K_Turb_G(IJK)/&
               (E_Turb_G(IJK)+small_number)

! NOTE: Some calculations are based explicitly on solids phase 1!
! Parameters based on L=1: tau_12_st....
            L = 1

! Particle relaxation time. For very dilute flows avoid singularity
! by redefining the drag as single particle drag
            IF(Ep_s(IJK,M) > DIL_EP_S .AND. &
               F_GS(IJK,L) > small_number) THEN
               Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,L)
            ELSE
               dgA_sl = kt_dga(IJK, L)
               Tau_12_st = RO_S(IJK,M)/DgA_SL
            ENDIF


! The following is purely an ad-hoc modification so that the underlying
! monodisperse theory can be used for polydisperse systems in a
! consistent manner.  That is, solids pressure, viscosity and
! conductivity must be additive. THe non-linear terms (eps^2) are
! corrected so the stresses of two or more identical solids phases are
! equal to those of a equivalent single solids phase. sof June 15 2005.
            SUM_EpsGo = ZERO
            DO MM = 1, SMAX
               SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,MM)*G_0(IJK,M,MM)
            ENDDO

            P_s_v(IJK) = ROP_s(IJK,M)*Theta_m(IJK,M) * &
               ( (ONE + 4.d0*SUM_EpsGo ) + HALF*(ONE - C_e*C_e) )

            IF(EP_s(IJK,M) < (ONE-EP_star_array(ijk))) THEN
               Tmp_Ahmadi_Const = &
                  ONE/(ONE+ Tau_1(ijk)/Tau_12_st * &
                  (ONE-EP_s(IJK,M)/(ONE-EP_star_array(ijk)))**3)
            ELSE
               Tmp_Ahmadi_Const = ONE
            ENDIF

! Shear viscosity
! Note that Ahmadi coefficient 0.0853 in C_mu was replaced by 0.1567
! to include 3/2*sqrt(3/2) because K = 3/2 Theta_m
            Mu_s_v(IJK) = Tmp_Ahmadi_Const &
                *0.1045d0*(ONE/G_0(IJK,M,M)+3.2d0*EP_s(IJK,M)+12.1824d0*   &
                G_0(IJK,M,M)*EP_s(IJK,M)*EP_s(IJK,M))*D_p(IJK,M)*RO_S(IJK,M)*  &
                DSQRT(Theta_m(IJK,M))

! Bulk viscosity
! The following formulation is a guess for Mu_b be by taking 5/3 of the
! collisional viscosity contribution. In this case col. visc. is the
! eps^2 contribution to Mu_s_v(IJK). This might be changed later if
! communications with Ahmadi reveals a different appoach
            Mu_b_v(IJK) = 5.d0/3.d0* Tmp_Ahmadi_Const                  &
                *0.1045d0*(12.1824d0*G_0(IJK,M,M)*EP_s(IJK,M)*EP_s(IJK,M)) &
                *D_p(IJK,M)*RO_S(IJK,M)* DSQRT(Theta_m(IJK,M))

! Second viscosity as defined in MFIX
            LAMBDA_S_V(IJK) = Eta*Mu_b_v(IJK) - (2d0*Mu_s_v(IJK))/3d0

! Defining Ahmadi conductivity from his equation 42 in Cao and Ahmadi 1995 paper
! note the constant 0.0711 is now 0.1306 because K = 3/2 theta_m
            Kth_s(IJK,M) = 0.1306D0*RO_S(IJK,M)*D_p(IJK,M)*(ONE+C_e**2)* &
                (ONE/G_0(IJK,M,M)+4.8D0*EP_s(IJK,M)+12.1184D0 &
                *EP_s(IJK,M)*EP_s(IJK,M)*G_0(IJK,M,M) )  &
                *DSQRT(Theta_m(IJK,M))


         ENDIF   ! Fluid_at
      ENDDO   ! IJK loop
      RETURN
      END SUBROUTINE GT_PDE_AHMADI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE_simonin                                          C
!  Purpose: Calculate granular stress terms (viscosity, bulk viscosity C
!  solids pressure) & granular conductivity using simonin model        C
!                                                                      C
!  Author: Sofiane Benyahia, Fluent Inc.            Date: 02-01-05     C
!                                                                      C
!  Literature/Document References:                                     C
!  Simonin, O., 1996. Combustion and turbulence in two-phase flows,    C
!     Von Karman institute for fluid dynamics, lecture series,         C
!     1996-02                                                          C
!  Balzer, G., Simonin, O., Boelle, A., and Lavieville, J., 1996,      C
!     A unifying modelling approach for the numerical prediction       C
!     of dilute and dense gas-solid two phase flow. CFB5, 5th int.     C
!     conf. on circulating fluidized beds, Beijing, China.             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine gt_pde_simonin (M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: C_e, eta
      USE constant, only: pi

      USE drag, only: f_gs

      USE fldvar, only: p_s_v
      USE fldvar, only: ep_g, ro_g
      USE fldvar, only: ROP_s, ro_s, d_p, theta_m
      USE fldvar, only: ep_s
      USE fldvar, only: k_turb_g, e_turb_g

      USE kintheory, only: kt_rvel, kt_dga, kt_cos_theta
      USE param1, only: zero, one, large_number, small_number

      USE physprop, only: smax
      USE physprop, only: kth_s

      USE rdf, only: g_0

      USE toleranc, only: dil_ep_s

      USE turb, only: tau_1, tau_12, k_12

      USE visc_s, only: lambda_s_v, mu_s_v, mu_b_v

      USE functions, only: fluid_at
      USE compar, only: ijkstart3, ijkend3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local parameters
!---------------------------------------------------------------------//
! constant in simonin model
      DOUBLE PRECISION, PARAMETER :: C_mu = 9.0D-02

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: IJK
! solids phase index
      INTEGER :: MM, L
! defining parameters for Simonin model
      DOUBLE PRECISION :: Tau_12_st, Tau_2_c, Tau_2, Zeta_r, C_Beta
      DOUBLE PRECISION :: Sigma_c, Zeta_c, Omega_c, Zeta_c_2, X_21, Nu_t
      DOUBLE PRECISION :: MU_2_T_Kin, Mu_2_Col, Kappa_kin, Kappa_Col
      DOUBLE PRECISION :: cos_theta
! sum of ep_s * g_0
      DOUBLE PRECISION :: SUM_EpsGo
! single particle drag coefficient/ep_s
      DOUBLE PRECISION :: dga_sl
! relative velocity between gas and solids
      DOUBLE PRECISION :: rvel_l
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

! Define time scales and constants

! time scale of turbulent eddies
            Tau_1(ijk) = 3.d0/2.d0*C_MU*K_Turb_G(IJK)/&
               (E_Turb_G(IJK)+small_number)

! NOTE: Some calculations are based explicitly on solids phase 1!
! Parameters based on L=1: tau_12_st, zeta_r, cos_theta, c_beta,
! tau_12, tau_2, nu_t, k_12...
            L = 1

! Particle relaxation time. For very dilute flows avoid singularity
! by redefining the drag as single particle drag
            IF(Ep_s(IJK,M) > DIL_EP_S .AND. &
               F_GS(IJK,L) > small_number) THEN
               Tau_12_st = Ep_s(IJK,M)*RO_S(IJK,M)/F_GS(IJK,L)
            ELSE
               DgA_sL = KT_DGA(ijk, l)
               Tau_12_st = RO_S(IJK,M)/DgA_sL
            ENDIF

! This is Zeta_r**2 as defined by Simonin
            rvel_l = KT_RVEL(ijk, l)
            Zeta_r = 3.0d0 * RVEL_L**2 / &
               (2.0d0*K_Turb_G(IJK)+small_number)

            cos_theta = KT_COS_Theta(IJK, L)
            C_Beta = 1.8d0 - 1.35d0*Cos_Theta**2

! Lagrangian Integral time scale: Tau_12
! the time-scale of the fluid turbulent motion viewed by the
! particle (crossing trajectory effect):
! no solids tau_12 = tau_1
            Tau_12(ijk) = Tau_1(ijk)/sqrt(ONE+C_Beta*Zeta_r)

! Defining the inter-particle collision time
            IF(Ep_s(IJK,M) > DIL_EP_S) THEN
               Tau_2_c = D_p(IJK,M)/(6.d0*Ep_s(IJK,M)*G_0(IJK,M,M) &
               *DSQRT(16.d0*(Theta_m(ijk,m)+Small_number)/PI))
            ELSE             ! assign it a large number
               Tau_2_c = LARGE_NUMBER
            ENDIF

! Define some constants
            Sigma_c = (ONE+ C_e)*(3.d0-C_e)/5.d0
! Zeta_c: const. to be used in the K_2 Diffusion coefficient.
            Zeta_c  = (ONE+ C_e)*(49.d0-33.d0*C_e)/100.d0
            Omega_c = 3.d0*(ONE+ C_e)**2 *(2.d0*C_e-ONE)/5.d0
            Zeta_c_2= 2./5.*(ONE+ C_e)*(3.d0*C_e-ONE)

! mixed time scale in the generalized Simonin theory (switch between dilute
! and kinetic theory formulation of the stresses)
            Tau_2 = ONE/(2./Tau_12_st+Sigma_c/Tau_2_c)
! The ratio of these two time scales.
            Nu_t =  Tau_12(ijk)/Tau_12_st


! The ratio of densities
            X_21 = Ep_s(IJK,M)*RO_S(IJK,M)/(EP_g(IJK)*RO_g(IJK))

! Definition of an "algebraic" form of of Simonin K_12 PDE. This is obtained
! by equating the dissipation term to the exchange terms in the PDE and
! neglecting all other terms, i.e. production, convection and diffusion.
! This works because Tau_12 is very small for heavy particles
            K_12(ijk) = Nu_t / (ONE+Nu_t*(ONE+X_21)) * &
                (2.d+0 *K_Turb_G(IJK) + 3.d+0 *X_21*theta_m(ijk,m))

! Realizability Criteria
            IF(K_12(ijk) > DSQRT(6.0D0*K_Turb_G(IJK)*theta_m(ijk,m))) THEN
               K_12(ijk) = DSQRT(6.0D0*K_Turb_G(IJK)*theta_m(ijk,m))
            ENDIF

! The following is purely an ad-hoc modification so that the underlying
! monodisperse theory can be used for polydisperse systems in a
! consistent manner.  That is, solids pressure, viscosity and
! conductivity must be additive. The non-linear terms (eps^2) are
! corrected so the stresses of two or more identical solids phases are
! equal to those of a equivalent single solids phase. sof June 15 2005.
            SUM_EpsGo = ZERO
            DO MM = 1, SMAX
               SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,MM)*G_0(IJK,M,MM)
            ENDDO


! Solids pressure
! Note this formulation is the same as standard granular pressure
            P_s_v(IJK) = ROP_s(IJK,M) * &
                  (ONE + 4.d0*Eta*SUM_EpsGo)*Theta_m(IJK,M)

! Solids viscosity: shear and bulk
! Turbulent Kinetic (MU_2_T_Kin) and collisional (Mu_2_Col) viscosities
            MU_2_T_Kin = (2.0d0/3.0d0*K_12(ijk)*Nu_t + Theta_m(IJK,M) * &
                (ONE+ zeta_c_2*EP_s(IJK,M)*G_0(IJK,M,M)))*Tau_2
            Mu_2_Col = 8.d0/5.d0*EP_s(IJK,M)*G_0(IJK,M,M)*Eta* (MU_2_T_Kin+ &
                D_p(IJK,M)*DSQRT(Theta_m(IJK,M)/PI))
            Mu_b_v(IJK) = 5.d0/3.d0*EP_s(IJK,M)*RO_S(IJK,M)*Mu_2_Col
            Mu_s_v(IJK) = EP_s(IJK,M)*RO_S(IJK,M)*(MU_2_T_Kin + Mu_2_Col)


! Second viscosity as defined in MFIX
            LAMBDA_S_V(IJK) = Eta*Mu_b_v(IJK) - (2d0*Mu_s_v(IJK))/3d0

! Solids conductivity
! Defining Turbulent Kinetic diffusivity: Kappa
            Kappa_kin = (9.d0/10.d0*K_12(ijk)*Nu_t + 3.0D0/2.0D0 * &
                Theta_m(IJK,M)*(ONE+ Omega_c*EP_s(IJK,M)*G_0(IJK,M,M)))/&
                (9.d0/(5.d0*Tau_12_st) + zeta_c/Tau_2_c)
            Kappa_Col = 18.d0/5.d0*EP_s(IJK,M)*G_0(IJK,M,M)*Eta* &
                (Kappa_kin+ 5.d0/9.d0*D_p(IJK,M)*DSQRT(Theta_m(IJK,M)/PI))
            Kth_s(IJK,M) =  EP_s(IJK,M)*RO_S(IJK,M)*(Kappa_kin + Kappa_Col)


         ENDIF   ! Fluid_at
      ENDDO   ! IJK loop
      RETURN
      END SUBROUTINE GT_PDE_SIMONIN



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE_GD                                               C
!  Purpose: Implement kinetic theory of Garzo and Dufty (1999) for     C
!  calculation of granular stress terms and granular conductivity      C
!                                                                      C
!  Author: Janine E. Galvin                                            C
!                                                                      C
!  Literature/Document References:                                     C
!  Garzo, V., and Dufty, J., Homogeneous cooling state for a           C
!     granular mixture, Physical Review E, 1999, Vol 60 (5), 5706-     C
!     5713                                                             C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
      Subroutine gt_pde_gd (M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: C_e, pi, switch
      USE drag, only: f_gs
      USE fldvar, only: p_s_v
      USE fldvar, only: rop_s, ro_s, d_p, theta_m
      USE fldvar, only: ep_s
      USE fldvar, only: ro_g
      USE kintheory, only: kt_dga
      USE param1, only: zero, small_number
      USE physprop, only: kth_s, kphi_s
      USE rdf, only: g_0, dg_0dnu
      USE toleranc, only: dil_ep_s
      USE visc_s, only: lambda_s_v, mu_s_v, mu_b_v
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell Indices
      INTEGER :: IJK
! Use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION :: Mu_star, Kth_star
!
      DOUBLE PRECISION :: D_PM, M_PM, NU_PM, EP_SM, RO_SM, ROP_SM
!
      DOUBLE PRECISION :: chi, dChiOdphi
!
      DOUBLE PRECISION :: c_star, zeta0_star, nu_eta_star, &
                          gamma_star, eta_k_star, eta_star, eta0, &
                          kappa0, nu_kappa_star, kappa_k_star, &
                          qmu_k_star, qmu_star, kappa_star, press_star
! single particle drag coefficient
      DOUBLE PRECISION :: dga_sm
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
          IF ( FLUID_AT(IJK) ) THEN

! local aliases
             D_PM = D_P(IJK,M)
             RO_SM = RO_S(IJK,M)
             ROP_SM = ROP_S(IJK,M)
             EP_SM = EP_S(IJK,M)
             M_PM = (PI/6.d0)*D_PM**3 * RO_SM
             NU_PM = ROP_SM/M_PM
             Chi = G_0(IJK,M,M)
             dChiOdphi = DG_0DNU(EP_SM)
             IF(EP_SM <= DIL_EP_S) &
                dga_sm = KT_DGA(IJK, M)


! Pressure/Viscosity/Bulk Viscosity
! Note: k_boltz = M_PM
!-----------------------------------
! Find pressure in the Mth solids phase
             press_star = 1.d0 + 2.d0*(1.d0+C_E)*EP_SM*Chi

! n*k_boltz = n*m = ep_s*ro_s
             P_s_v(IJK) = ROP_sM*Theta_m(IJK,M)*press_star

! find bulk and shear viscosity
             c_star = 32.0d0*(1.0d0 - C_E)*(1.d0 - 2.0d0*C_E*C_E) &
                / (81.d0 - 17.d0*C_E + 30.d0*C_E*C_E*(1.0d0-C_E))

             zeta0_star = (5.d0/12.d0)*Chi*(1.d0 - C_E*C_E) &
                * (1.d0 + (3.d0/32.d0)*c_star)

             nu_eta_star = Chi*(1.d0 - 0.25d0*(1.d0-C_E)*(1.d0-C_E)) &
                * (1.d0-(c_star/64.d0))

             gamma_star = (4.d0/5.d0)*(32.d0/PI)*EP_SM*EP_SM &
                * Chi*(1.d0+C_E) * (1.d0 - (c_star/32.d0))

             eta_k_star = (1.d0 - (2.d0/5.d0)*(1.d0+C_E)*(1.d0-3.d0*C_E) &
                * EP_SM*Chi ) / (nu_eta_star - 0.5d0*zeta0_star)

             eta_star = eta_k_star*(1.d0 + (4.d0/5.d0)*EP_SM*Chi &
                * (1.d0+C_E) ) + (3.d0/5.d0)*gamma_star

             eta0 = 5.0d0*M_PM*DSQRT(Theta_m(IJK,M)/PI) / (16.d0*D_PM*D_PM)

             IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN
                Mu_star = eta0
             ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
                Mu_star = ZERO
             ELSEIF(EP_SM <= DIL_EP_S) THEN
                Mu_star = RO_SM*EP_SM*Chi*Theta_m(IJK,M)*eta0 / &
                   ( RO_S(IJK,M)*EP_SM*Chi*Theta_m(IJK,M) + &
                   2.d0*DgA_sM*eta0/RO_S(IJK,M) )
             ELSE
                Mu_star = RO_SM*EP_SM*Chi*Theta_m(IJK,M)*eta0 / &
                   ( RO_SM*EP_SM*Chi*Theta_m(IJK,M) + &
                   (2.d0*F_gs(IJK,M)*eta0/(RO_SM*EP_SM)) )
             ENDIF

! Shear and bulk viscosity
             Mu_s_v(IJK) = Mu_star * eta_star
             Mu_b_v(IJK) = Mu_star * gamma_star

! Second viscosity
             LAMBDA_S_V(IJK) = Mu_b_v(IJK) - (2.d0/3.d0)*Mu_s_v(IJK)


! Granular Conductivity/Dufour Coefficient
!-----------------------------------
             kappa0 = (15.d0/4.d0)*eta0

             nu_kappa_star = (Chi/3.d0)*(1.d0+C_E) * ( 1.d0 + &
                (33.d0/16.d0)*(1.d0-C_E) + ((19.d0-3.d0*C_E)/1024.d0)*&
                c_star)
!             nu_mu_star = nu_kappa_star

             kappa_k_star = (2.d0/3.d0)*(1.d0 +0.5d0*(1.d0+press_star)*&
                c_star + (3.d0/5.d0)*EP_SM*Chi*(1.d0+C_E)*(1.d0+C_E) * &
                (2.d0*C_E - 1.d0 + ( 0.5d0*(1.d0+C_E) - 5.d0/&
                (3*(1.d0+C_E))) * c_star ) ) / (nu_kappa_star - &
                2.d0*zeta0_star)

             kappa_star = kappa_k_star * (1.d0 + (6.d0/5.d0)*EP_SM* &
                Chi*(1.d0+C_E) ) + (256.d0/25.d0)*(EP_SM* &
                EP_SM/PI)*Chi*(1.d0+C_E)*(1.d0+(7.d0/32.d0)* &
                c_star)

             IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN
                Kth_star= kappa0
             ELSEIF(Theta_m(IJK,M) .LT. SMALL_NUMBER)THEN
                Kth_star = ZERO
             ELSEIF(EP_SM <= DIL_EP_S) THEN
                Kth_star = RO_SM*EP_SM*Chi*Theta_m(IJK,M)*kappa0/ &
                   (RO_SM*EP_SM*Chi*Theta_m(IJK,M) + 1.2d0*DgA_sM* &
                   kappa0/RO_SM)
             ELSE
                Kth_star = RO_SM*EP_SM*Chi*Theta_m(IJK,M)*kappa0/ &
                   (RO_SM*EP_SM*Chi*Theta_m(IJK,M)+ &
                   (1.2d0*F_gs(IJK,M)*kappa0/(RO_SM*EP_SM)) )
             ENDIF

! Granular conductivity
             Kth_s(IJK,M) = Kth_star * kappa_star

! transport coefficient of the Mth solids phase associated
! with gradient in volume fraction in heat flux
             qmu_k_star = 2.d0*( (1.d0+EP_SM*dChiOdphi)* &
                zeta0_star*kappa_k_star + ( (press_star/3.d0) + &
                (2.d0/3.d0)* EP_SM*(1.d0+C_E)*(Chi+EP_SM* dChiOdphi) )*&
                c_star - (4.d0/5.d0)*EP_SM*Chi* (1.d0+(EP_SM/2.d0)*&
                dChiOdphi)* (1.d0+C_E) * ( C_E*(1.d0-C_E)+0.25d0*&
                ((4.d0/3.d0)+C_E* (1.d0-C_E))*c_star ) ) / &
                (2.d0*nu_kappa_star-3.d0*zeta0_star)

             qmu_star = qmu_k_star*(1.d0+(6.d0/5.d0)*EP_SM*Chi*&
                (1.d0+C_E) )

             IF (EP_SM .LT. SMALL_NUMBER) THEN
                Kphi_s(IJK,M) = ZERO
             ELSE
                Kphi_s(IJK,M) = (Theta_m(IJK,M)*Kth_star/NU_PM)*qmu_star
             ENDIF

          ENDIF   ! Fluid_at
      ENDDO   ! IJK loop
      RETURN
      END SUBROUTINE GT_PDE_GD



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE_GTSH                                             C
!  Purpose: Implement kinetic theory of Garzo, Tenneti, Subramaniam    C
!  Hrenya (2012) for calculation of granular stress terms and          C
!                                                                      C
!  Author: Sofiane Benyahia                                            C
!                                                                      C
!  Literature/Document References:                                     C
!  Garzo, V., Tenneti, S., Subramaniam, S., and Hrenya, C. M.,         C
!      "Enskog kinetic theory for monodisperse gas-solid flows", JFM,  C
!      Vol. 712, 2012, pp. 129-168                                     C
!                                                                      C
!  Comments:                                                           C
!  And also based on C.M. Hrenya hand-notes dated Sep 2013             C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
      Subroutine gt_pde_gtsh (M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: C_e, pi

      USE fldvar, only: p_s_v
      USE fldvar, only: ro_g
      USE fldvar, only: rop_s, ro_s, d_p, theta_m
      USE fldvar, only: ep_s

      USE kintheory, only: EDT_s_ip, xsi_gtsh, a2_gtsh
      USE kintheory, only: kt_rvel
      USE kintheory, only: epm, G_gtsh, S_star, K_phi, R_d

      USE param1, only: zero, one, small_number

      USE physprop, only: mu_g
      USE physprop, only: kth_s, kphi_s

      USE rdf, only: g_0, dg_0dnu

      USE visc_s, only: lambda_s_v, mu_s_v, mu_b_v

      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell Indices
      INTEGER :: IJK
!
      DOUBLE PRECISION :: D_PM, M_PM, NU_PM, EP_SM
!
      DOUBLE PRECISION :: chi, dChiOdphi, eta0
!
      DOUBLE PRECISION :: nu0, nuN, etaK
      DOUBLE PRECISION :: dZeta_dT, dGama_dT, NuK, Kth0, KthK
      DOUBLE PRECISION :: Rdissdphi, Kphidphi, Re_T, dGamadn, dRdphi
      DOUBLE PRECISION :: denom
      DOUBLE PRECISION :: dSdphi, R_dphi, Tau_st, dPsidn, MuK
! relative velocity
      DOUBLE PRECISION :: RVEL
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
          IF ( FLUID_AT(IJK) ) THEN

! local aliases
             D_PM = D_P(IJK,M)
             EP_SM = EP_S(IJK,M)
             M_PM = (PI/6.d0)*D_PM**3 * RO_S(IJK,M)
             NU_PM = ROP_S(IJK,M)/M_PM
             Chi = G_0(IJK,M,M)
             dChiOdphi = DG_0DNU(EP_SM)
             RVEL = kt_rvel(ijk, m)


! note that T = (m_pm*theta_m)
! Pressure/Viscosity/Bulk Viscosity
!-----------------------------------
! solids pressure, eq (6.14) of GTSH theory
             P_s_v(IJK) = ROP_s(IJK,M)*Theta_m(IJK,M)* &
                (one+2d0*(one+C_E)*Chi*EP_SM)

! evaluating shear viscosity, eq (7.3) of GTSH theory
! starting with nu_0 equ (7.6-7.8)
             eta0 = 0.3125d0/(dsqrt(pi)*D_PM**2)*M_pm*&
                dsqrt(theta_m(ijk,m))

             nu0 = (96.d0/5.d0)*(EP_SM/D_PM)*DSQRT(Theta_m(IJK,M)/PI)
!             nu0 = NU_PM*M_pm*theta_m(ijk,m)/eta0

             nuN = 0.25d0*nu0*Chi*(3d0-C_E)*(one+C_E) * &
                    (one+0.4375d0*A2_gtsh(ijk))
! defining kinetic part of shear viscosity nuK  equ (7.7)
             etaK = rop_s(ijk,m)*theta_m(ijk,m) / (nuN-0.5d0*( &
                EDT_s_ip(ijk,M,M)-xsi_gtsh(ijk)/theta_m(ijk,m) - &
                2d0*G_gtsh(EP_SM, Chi, IJK, M)/M_PM)) * (one -0.4d0 * &
                (one+C_E)*(one-3d0*C_E)*EP_SM*Chi)

! bulk viscosity lambda eq. (7.5)
             Mu_b_v(IJK) = 25.6d0/pi * EP_SM**2 * Chi *(one+C_E) * &
                (one - A2_gtsh(ijk)/16d0)*eta0

! Finally shear viscosity, eq (7.9) of GTSH theory
             Mu_s_v(IJK) = etaK*(one+0.8d0*EP_SM*Chi*(one+C_E)) + &
                0.6d0*Mu_b_v(IJK)

! Second viscosity as defined in MFIX
             LAMBDA_S_V(IJK) = Mu_b_v(IJK) - (2.d0/3.d0)*Mu_s_v(IJK)


! Conductivity/Dufour coefficient
!-----------------------------------
! Calculate conductivity Kth_s(IJK,M), eq. 7.12 GTSH theory.
! Start with calculating dZeta/dT and dGama/dT
! note that 1/Tau**2 = (3d0*pi*mu_g(ijk)*D_PM/M_p)**2 defined
! under eq. 8.2 GTSH
             dZeta_dT = -0.5d0*xsi_gtsh(ijk)/(M_pm*theta_m(ijk,m))

             dGama_dT = 3d0*pi*D_PM**2*RO_g(ijk)*K_phi(EP_SM)/ &
                (2d0*M_pm*dsqrt(theta_m(ijk,m)))

! evaluating eq (7.16) in GTSH
             NuK = nu0*(one+C_E)/3d0*Chi*( one+2.0625d0*(one-C_E)+ &
                ((947d0-579*C_E)/256d0*A2_gtsh(ijk)) )

! evaluating eq. (7.13)
             Kth0 = 3.75d0*eta0/M_pm

! evaluating kinetic conductivity Kk eq. (7.14)
! note that 1/2m/T Psi and m dZeta_dT cancel out.
             KthK = zero
             IF(EP_SM > SMALL_NUMBER) KthK = 2d0/3d0*Kth0*nu0 / (NuK - &
                2d0*EDT_s_ip(ijk,M,M) - 2d0*theta_m(ijk,m)*dGama_dT) * &
                (one+2d0*A2_gtsh(ijk)+0.6d0*EP_SM*Chi* &
                (one+C_E)**2*(2*C_E-one+A2_gtsh(ijk)*(one+C_E)))

! the conductivity Kth from eq (7.17) in GTSH theory:
             Kth_s(IJK,M) = KthK*(one+1.2d0*EP_SM*Chi*(one+C_E)) + &
                (10.24d0/pi* EP_SM**2*Chi*(one+C_E)*(one+0.4375d0* &
                A2_gtsh(ijk))*Kth0)

! Finaly notice that conductivity K in eq (7.10) must be
! multiplied by m because of grad(T)
             Kth_s(IJK,M) = M_pm * Kth_s(IJK,M)

! Calculate the Dufour coefficient Kphi_s(IJK,M) in equation (7.18) of
! GTSH theory.
! First, calculate terms in 2 n/m x Gama_n, dRdiss/dphi and dK_phi/dphi.
! Notice that 2 n/m Gama_n = 2 phi/m Gama_phi, so multiply the
! deriviatives of Rdiss and K_phi by phi to avoid possible division by
! phi.
             Rdissdphi = ZERO
             IF(EP_SM > SMALL_NUMBER) Rdissdphi = &
                1.5d0*dsqrt(EP_SM/2d0)+135d0/64d0*EP_SM*(dlog(EP_SM)+one) +&
                11.26d0*EP_SM*(one-10.2*EP_SM+49.71d0*EP_SM**2-87.08d0* &
                EP_SM**3) - EP_SM*dlog(epM)*(Chi+EP_SM*dChiOdphi)
! corrections due to W. Fullmer
             Kphidphi = EP_SM*(0.212d0*0.142d0/(EP_SM**0.788d0*&
                (one-EP_SM)**4.454d0) + 4.454d0*K_phi(EP_SM)/(one-EP_SM))
             Kphidphi = zero  ! this is compatible with K_phi = zero

             Re_T = ro_g(ijk)*D_PM*dsqrt(theta_m(ijk,m)) / mu_g(ijk)

! The term phi x Gama_phi becomes
             dGamadn = 3d0*pi*D_pm*Mu_g(ijk)*(Rdissdphi+Re_T*Kphidphi)

! Second, calculate terms in 2 rho x Psi_n, which is same as ro_s x
! EP_SM x Psi_n
! Take EP_SM inside the derivative dS_star/dphi to avoid singularities
! as EP_SM -> 0

! calculating the term phi*dRd/dphi
             dRdphi = zero
             IF((EP_SM > SMALL_NUMBER) .AND. (EP_SM <= 0.4d0)) THEN
                denom = one+0.681d0*EP_SM-8.48d0*EP_SM**2+8.16d0*EP_SM**3
                dRdphi = (1.5d0*dsqrt(EP_SM/2d0)+135d0/64d0*EP_SM*&
                   (dlog(EP_SM)+one)+ 17.14d0*EP_SM)/denom - EP_SM*&
                   (one+3d0*dsqrt(EP_SM/2d0) + 135d0/64d0*EP_SM*&
                   dlog(EP_SM)+17.14*EP_SM)/denom**2 * &
                   (0.681d0-16.96d0*EP_SM+24.48d0*EP_SM**2)
             ELSEIF(EP_SM > 0.4d0) THEN
                dRdphi = 10d0*(one+2d0*EP_SM)/(one-EP_SM)**4
             ENDIF

! calculating the term phi*dS_star/dphi
             dSdphi = zero
             IF(EP_SM >= 0.1d0) THEN
                R_dphi = R_d(EP_SM)
                denom = one+3.5d0*dsqrt(EP_SM)+5.9d0*EP_SM
                dSdphi = 2d0*R_dphi*dRdphi/(Chi*denom) - &
                   EP_SM*R_dphi**2 * (dChiOdphi/(Chi**2*denom) + &
                   (1.75d0/dsqrt(EP_SM)+5.9d0)/(Chi*denom**2))
             ENDIF

! defining the relaxation time Tau_st
             Tau_st = M_pm/(3d0*pi*mu_g(ijk)*D_pm)

! The term phi x Psi_n becomes
             dPsidn = dsqrt(pi)*D_pm**4*RVEL**2 / &
                (36d0*Tau_st**2*dsqrt(theta_m(ijk,m))) * dSdphi

! Now compute the kinetic contribution to Dufour coef. Muk eq (7.20) GSTH
             Muk = ZERO  ! This is assumed to avoid /0 for EP_SM = 0
             IF(EP_SM > SMALL_NUMBER) Muk = Kth0*Nu0*M_pm*&
                theta_m(ijk,m)/NU_PM / (NuK-1.5d0*(EDT_s_ip(ijk,M,M)-&
                xsi_gtsh(ijk)/theta_m(ijk,m))) * ( KthK/(Kth0*Nu0)*&
                (2d0/M_pm*dGamadn-ro_s(IJK,m)/(M_pm* theta_m(ijk,m))*&
                dPsidn + EDT_s_ip(ijk,M,M)*(one+EP_SM/Chi*dChiOdphi)) +&
                2d0/3d0*A2_gtsh(ijk) + 0.8d0*EP_SM*Chi* (one+C_E)*&
                (one+0.5d0*EP_SM/Chi*dChiOdphi)*(C_E*(C_E-one)+ &
                A2_gtsh(ijk)/6d0*(16d0-3d0*C_E+3d0*C_E**2)))

! Finaly compute the Dufour coefficient Mu (Kphi_s(IJK,M)) from eq (7.22) GTSH
             Kphi_s(IJK,M) = Muk*(one+1.2d0*EP_SM*Chi*(one+C_E))

          ENDIF   ! Fluid_at
      ENDDO   ! outer IJK loop

      RETURN
      END SUBROUTINE GT_PDE_GTSH



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GT_PDE_IA                                               C
!  Purpose: Implement kinetic theory of Iddir & Arastoopour (2005) for C
!  calculation of granular stress terms and granular conductivity      C
!                                                                      C
!  Author: Janine E. Galvin, Univeristy of Colorado                    C
!                                                                      C
!  Literature/Document References:                                     C
!  Iddir, Y.H., PhD Modeling of the multiphase mixture of particles    C
!     using the kinetic theory approach, PhD Dissertation in           C
!     Chemical Engineering, Illinois Institute of Technology, 2004     C
!  Iddir, Y.H., H. Arastoopour, and C.M. Hrenya, Analysis of binary    C
!     and ternary granular mixtures behavior using the kinetic         C
!     theory approach. Powder Technology, 2005, p. 117-125.            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine gt_pde_ia (M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: switch, switch_ia
      USE constant, only: C_e
      USE constant, only: pi

      USE drag, only: f_gs
      USE fldvar, only: RO_g
      USE fldvar, only: ROP_s, RO_s, D_p, theta_m
      USE fldvar, only: p_s_v
      USE fldvar, only: ep_s

      USE kintheory, only: mu_sm_ip, mu_sl_ip
      USE kintheory, only: xi_sm_ip, xi_sl_ip
      USE kintheory, only: kth_sl_ip, knu_sm_ip, knu_sl_ip, kvel_s_ip
      USE kintheory, only: kt_dga

      USE param1, only: one, zero, small_number

      USE physprop, only: smax
      USE physprop, only: Kth_s

      USE rdf, only: g_0

      USE toleranc, only: dil_ep_s
      USE ur_facs, only: UR_Kth_sml
      USE visc_s, only: mu_s_v, lambda_s_v

      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: IJK
! solids phase index
      INTEGER :: L
! use to compute MU_s(IJK,M) & Kth_S(IJK,M)
      DOUBLE PRECISION :: Mu_star, Mu_s_dil, Kth_star, K_s_dil
! variables for Iddir equipartition model
      DOUBLE PRECISION :: P_s_sum, P_s_MM, P_s_LM
      DOUBLE PRECISION :: MU_common_term, K_common_term
      DOUBLE PRECISION :: Mu_sM_sum, MU_s_MM, MU_s_LM, MU_sM_LM, MU_sL_LM
      DOUBLE PRECISION :: XI_sM_sum, XI_s_v
      DOUBLE PRECISION :: M_PM, M_PL, MPSUM, NU_PL, NU_PM, D_PM, D_PL, DPSUMo2
      DOUBLE PRECISION :: Ap_lm, Dp_lm, R0p_lm, R1p_lm, R8p_lm, R9p_lm, Bp_lm,&
                          R5p_lm, R6p_lm, R7p_lm
      DOUBLE PRECISION :: K_s_sum, K_s_MM, K_s_LM
! Sum of ep_s * g_0
      DOUBLE PRECISION :: SUM_EpsGo
! Current value of Kth_sl_ip (i.e., without underrelaxation)
      DOUBLE PRECISION :: Kth_sL_iptmp
! drag coefficient/ep_s
      DOUBLE PRECISION :: dga_sm
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

! Added for consistency of IA KT (see constant_mod for details)
            IF(SWITCH_IA) THEN
               SUM_EpsGo = ZERO
               DO L = 1, SMAX
                  SUM_EpsGo =  SUM_EpsGo+EP_s(IJK,L)*G_0(IJK,M,L)
               ENDDO
            ELSE
               SUM_EpsGo =  EP_s(IJK,M)*G_0(IJK,M,M)
            ENDIF

            P_s_sum = ZERO
            Mu_sM_sum = ZERO
            XI_sM_sum = ZERO
            IF(EP_S(IJK,M) <= DIL_EP_s) &
               dga_sm = KT_DGA(ijk, M)

            D_PM = D_P(IJK,M)
            M_PM = (PI/6.d0)*D_PM**3 * RO_S(IJK,M)
            NU_PM = ROP_S(IJK,M)/M_PM

            P_s_MM = NU_PM*Theta_m(IJK,M)

            MU_s_dil = (5.d0/96.d0)*D_PM* RO_S(IJK,M)*&
                DSQRT(PI*Theta_m(IJK,M)/M_PM)

            IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN
               Mu_star = MU_s_dil
            ELSEIF(Theta_m(IJK,M)/M_PM < SMALL_NUMBER)THEN
               Mu_star = ZERO
            ELSEIF(EP_S(IJK,M) <= DIL_EP_s) THEN
               Mu_star = MU_s_dil*EP_s(IJK,M)*G_0(IJK,M,M)/ &
                  ( SUM_EpsGo + 2.0d0*DgA_sM*MU_s_dil /&
                    (RO_S(IJK,M)**2 *(Theta_m(IJK,M)/M_PM)) )
            ELSE
               Mu_star = MU_s_dil*EP_S(IJK,M)*G_0(IJK,M,M)/ &
                  ( SUM_EpsGo + 2.0d0*F_gs(IJK,M)*MU_s_dil /&
                    (RO_S(IJK,M)**2 *EP_s(IJK,M)*&
                     (Theta_m(IJK,M)/M_PM)) )
            ENDIF

            MU_s_MM = (Mu_star/G_0(IJK,M,M))*&
                (1.d0+(4.d0/5.d0)*(1.d0+C_E)*SUM_EpsGo)**2

            DO L = 1, SMAX
               D_PL = D_P(IJK,L)
               M_PL = (PI/6.d0)*D_PL**3 * RO_S(IJK,L)
               MPSUM = M_PM + M_PL
               DPSUMo2 = (D_PM+D_PL)/2.d0
               NU_PL = ROP_S(IJK,L)/M_PL

               IF ( L .eq. M) THEN
                  Ap_lm = MPSUM/(2.d0)
                  Dp_lm = M_PL*M_PM/(2.d0*MPSUM)
                  R0p_lm = ONE/( Ap_lm**1.5 * Dp_lm**2.5 )
                  R1p_lm = ONE/( Ap_lm**1.5 * Dp_lm**3 )

                  P_s_LM = PI*(DPSUMo2**3 / 48.d0)*G_0(IJK,M,L)*&
                      (M_PM*M_PL/MPSUM)* (M_PM*M_PL)**1.5 *&
                      NU_PM*NU_PL*(1.d0+C_E)*R0p_lm*Theta_m(IJK,M)

                  MU_s_LM = DSQRT(PI)*( DPSUMo2**4 / 240d0 )*&
                      G_0(IJK,M,L)*(M_PL*M_PM/MPSUM)**2 *&
                      (M_PL*M_PM)**1.5 * NU_PM*NU_PL*&
                      (1.d0+C_E) * R1p_lm * DSQRT(Theta_m(IJK,M))

! This is Mu_i_1 as defined in eq (16) of Galvin document
                  MU_sM_ip(IJK,M,L) = (MU_s_MM + MU_s_LM)

! This is Mu_ii_2 as defined in eq (17) of Galvin document
                  MU_sL_ip(IJK,M,L) = MU_s_LM

! solids phase viscosity associated with the trace of
! solids phase M (eq. 18 from Galvin theory document)
                  XI_sM_ip(IJK,M,L) = (5.d0/3.d0)*MU_s_LM

! solids phase viscosity associated with the trace
! of (sum of) all solids phases (eq. 19)
                  XI_sL_ip(IJK,M,L) = (5.d0/3.d0)*MU_s_LM

               ELSE
                  Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*&
                      Theta_m(IJK,M))/2.d0
                  Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-&
                      Theta_m(IJK,M) ))/(2.d0*MPSUM)
                  Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+&
                      M_PL*Theta_m(IJK,L) ))/&
                      (2.d0*MPSUM*MPSUM)
                  R0p_lm = (1.d0/(Ap_lm**1.5 * Dp_lm**2.5))+ &
                      ((15.d0*Bp_lm*Bp_lm)/(2.d0* Ap_lm**2.5 *&
                      Dp_lm**3.5))+&
                      ((175.d0*(Bp_lm**4))/(8.d0*Ap_lm**3.5 * &
                      Dp_lm**4.5))
                  R1p_lm = (1.d0/((Ap_lm**1.5)*(Dp_lm**3)))+ &
                      ((9.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**4))+&
                      ((30.d0*Bp_lm**4) /( 2.d0*Ap_lm**3.5 * &
                      Dp_lm**5))

                  P_s_LM = PI*(DPSUMo2**3 / 48.d0)*G_0(IJK,M,L)*&
                      (M_PM*M_PL/MPSUM)* (M_PM*M_PL)**1.5 *&
                      NU_PM*NU_PL*(1.d0+C_E)*R0p_lm* &
                      (Theta_m(IJK,M)*Theta_m(IJK,L))**2.5

                  MU_common_term = DSQRT(PI)*( DPSUMo2**4 / 240d0 )*&
                      G_0(IJK,M,L)*(M_PL*M_PM/MPSUM)**2 *&
                      (M_PL*M_PM)**1.5 * NU_PM*NU_PL*&
                       (1.d0+C_E) * R1p_lm
                  MU_sM_LM = MU_common_term * Theta_m(IJK,M)**2 *&
                      Theta_m(IJK,L)**3
                  MU_sL_LM = MU_common_term * Theta_m(IJK,L)**2 *&
                      Theta_m(IJK,M)**3

! solids phase 'viscosity' associated with the divergence
! of solids phase M. defined in eq (16) of Galvin document
                  MU_sM_ip(IJK,M,L) = MU_sM_LM

! solids phase 'viscosity' associated with the divergence
! of all solids phases. defined in eq (17) of Galvin document
                  MU_sL_ip(IJK,M,L) = MU_sL_LM

! solids phase viscosity associated with the trace of
! solids phase M
                  XI_sM_ip(IJK,M,L) = (5.d0/3.d0)*MU_sM_LM

! solids phase viscosity associated with the trace
! of all solids phases
                  XI_sL_ip(IJK,M,L) = (5.d0/3.d0)*MU_sL_LM
               ENDIF

               P_s_sum = P_s_sum + P_s_LM
               MU_sM_sum = MU_sM_sum + MU_sM_ip(IJK,M,L)
               XI_sM_sum = XI_sM_sum + XI_sM_ip(IJK,M,L)
            ENDDO

! Find the term proportional to the identity matrix
! (pressure in the Mth solids phase)
            P_s_v(IJK) = P_s_sum + P_S_MM

! Find the term proportional to the gradient in velocity
! of phase M  (shear viscosity in the Mth solids phase)
            MU_s_v(IJK) = MU_sM_sum + MU_sL_ip(IJK,M,M)
            XI_s_v = XI_sM_sum + XI_sL_ip(IJK,M,M)

! Bulk viscosity in the Mth solids phase
            LAMBDA_s_v(IJK) = -(2.d0/3.d0)*Mu_s_v(IJK) + XI_s_v


! Find the granular conductivity
            K_s_sum = ZERO

            K_s_dil = (75.d0/384.d0)*D_PM* RO_S(IJK,M)*&
                 DSQRT(PI*Theta_m(IJK,M)/M_PM)

            IF(SWITCH == ZERO .OR. RO_G(IJK) == ZERO) THEN
               Kth_star = K_s_dil
            ELSEIF(Theta_m(IJK,M)/M_PM < SMALL_NUMBER)THEN
               Kth_star = ZERO

            ELSEIF(EP_S(IJK,M) <= DIL_EP_s) THEN
               Kth_star = K_s_dil*EP_s(IJK,M)*G_0(IJK,M,M)/ &
                   (SUM_EpsGo+ 1.2d0*DgA_sM*K_s_dil &
                   / (RO_S(IJK,M)**2 *(Theta_m(IJK,M)/M_PM)))
            ELSE
               Kth_star = K_s_dil*EP_S(IJK,M)*G_0(IJK,M,M)/ &
                   (SUM_EpsGo+ 1.2d0*F_gs(IJK,M)*K_s_dil &
                   / (RO_S(IJK,M)**2 *EP_s(IJK,M)*(Theta_m(IJK,M)/M_PM)))
            ENDIF

! Kth doesn't include the mass.
            K_s_MM = (Kth_star/(M_PM*G_0(IJK,M,M)))*&
                 (1.d0+(3.d0/5.d0)*(1.d0+C_E)*(1.d0+C_E)*SUM_EpsGo)**2

            DO L = 1, SMAX
               D_PL = D_P(IJK,L)
               M_PL = (PI/6.d0)*D_PL**3 *RO_S(IJK,L)
               MPSUM = M_PM + M_PL
               DPSUMo2 = (D_PM+D_PL)/2.d0
               NU_PL = ROP_S(IJK,L)/M_PL

               IF ( L .eq. M) THEN

! solids phase 'conductivity' associated with the
! difference in velocity. again these terms cancel when
! added together so do not explicity include them
! in calculations
                  Kvel_s_ip(IJK,M,L) = ZERO
                  !     K_common_term*NU_PM*NU_PL*&
                  !     (3.d0*PI/10.d0)*R0p_lm*Theta_m(IJK,M)

! solids phase 'conductivity' associated with the
! difference in the gradient in number densities.
! again these terms cancel so do not explicity include
! them in calculations
                  Knu_sL_ip(IJK,M,L) = ZERO
                  !     K_common_term*NU_PM*&
                  !     (PI*DPSUMo2/6.d0)*R1p_lm*(Theta_m(IJK,M)**(3./2.))

                  Knu_sM_ip(IJK,M,L) = ZERO
                  !     K_common_term*NU_PL*&
                  !     (PI*DPSUMo2/6.d0)*R1p_lm*(Theta_m(IJK,M)**(3./2.))

                  K_s_sum = K_s_sum + K_s_MM

               ELSE
                  Ap_lm = (M_PM*Theta_m(IJK,L)+M_PL*Theta_m(IJK,M))/2.d0
                  Bp_lm = (M_PM*M_PL*(Theta_m(IJK,L)-&
                      Theta_m(IJK,M) ))/(2.d0*MPSUM)
                  Dp_lm = (M_PL*M_PM*(M_PM*Theta_m(IJK,M)+&
                      M_PL*Theta_m(IJK,L) ))/(2.d0*MPSUM*MPSUM)
                  R0p_lm = (1.d0/(Ap_lm**1.5 * Dp_lm**2.5))+&
                      ((15.d0*Bp_lm*Bp_lm)/(2.d0* Ap_lm**2.5 * Dp_lm**3.5))+&
                      ((175.d0*(Bp_lm**4))/(8.d0*Ap_lm**3.5 * Dp_lm**4.5))

                  R1p_lm = (1.d0/((Ap_lm**1.5)*(Dp_lm**3)))+ &
                      ((9.d0*Bp_lm*Bp_lm)/(Ap_lm**2.5 * Dp_lm**4))+&
                      ((30.d0*Bp_lm**4)/(2.d0*Ap_lm**3.5 * Dp_lm**5))

                  R5p_lm = (1.d0/(Ap_lm**2.5 * Dp_lm**3 ) )+ &
                      ((5.d0*Bp_lm*Bp_lm)/(Ap_lm**3.5 * Dp_lm**4))+&
                      ((14.d0*Bp_lm**4)/(Ap_lm**4.5 * Dp_lm**5))

                  R6p_lm = (1.d0/(Ap_lm**3.5 * Dp_lm**3))+ &
                      ((7.d0*Bp_lm*Bp_lm)/(Ap_lm**4.5 * Dp_lm**4))+&
                      ((126.d0*Bp_lm**4)/(5.d0*Ap_lm**5.5 * Dp_lm**5))

                  R7p_lm = (3.d0/(2.d0*Ap_lm**2.5 * Dp_lm**4))+ &
                      ((10.d0*Bp_lm*Bp_lm)/(Ap_lm**3.5 * Dp_lm**5))+&
                      ((35.d0*Bp_lm**4)/(Ap_lm**4.5 * Dp_lm**6))

                  R8p_lm = (1.d0/(2.d0*Ap_lm**1.5 * Dp_lm**4))+ &
                      ((6.d0*Bp_lm*Bp_lm)/(Ap_lm**2.5 * Dp_lm**5))+&
                      ((25.d0*Bp_lm**4)/(Ap_lm**3.5 * Dp_lm**6))

                  R9p_lm = (1.d0/(Ap_lm**2.5 * Dp_lm**3))+ &
                      ((15.d0*Bp_lm*Bp_lm)/(Ap_lm**3.5 * Dp_lm**4))+&
                      ((70.d0*Bp_lm**4)/(Ap_lm**4.5 * Dp_lm**5))

                  K_common_term = DPSUMo2**3 * M_PL*M_PM/(2.d0*MPSUM)*&
                       (1.d0+C_E)*G_0(IJK,M,L) * (M_PM*M_PL)**1.5

! solids phase 'conductivity' associated with the
! gradient in granular temperature of species M
                  K_s_LM = - K_common_term*NU_PM*NU_PL*(&
                   ((DPSUMo2*DSQRT(PI)/16.d0)*(3.d0/2.d0)*Bp_lm*R5p_lm)+&
                   ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*PI/6.d0)*&
                   (3.d0/2.d0)*R1p_lm)-(&
                   ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PM/8.d0)*Bp_lm*R6p_lm)+&
                   ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                   8.d0)*M_PM*R9p_lm)+&
                   ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PL*M_PM/(MPSUM*MPSUM))*&
                   M_PL*Bp_lm*R7p_lm)+&
                   ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                   2.d0)*(M_PM/MPSUM)**2 * M_PL*R8p_lm)+&
                   ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PM*M_PL/(2.d0*MPSUM))*&
                   R9p_lm)-&
                   ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*DPSUMo2*DSQRT(PI)*&
                   (M_PM*M_PL/MPSUM)*Bp_lm*R7p_lm) )*Theta_M(IJK,L) )*&
                   (Theta_M(IJK,M)**2 * Theta_M(IJK,L)**3)

! solids phase 'conductivity' associated with the
! gradient in granular temperature of species L
                    !Kth_sL_ip(IJK,M,L) = K_common_term*NU_PM*NU_PL*(&
                  Kth_sl_iptmp = K_common_term*NU_PM*NU_PL*(&
                    (-(DPSUMo2*DSQRT(PI)/16.d0)*(3.d0/2.d0)*Bp_lm*R5p_lm)-&
                    ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*PI/6.d0)*&
                    (3.d0/2.d0)*R1p_lm)+(&
                    ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PL/8.d0)*Bp_lm*R6p_lm)+&
                    ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                    8.d0)*M_PL*R9p_lm)+&
                    ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PL*M_PM/(MPSUM*MPSUM))*&
                    M_PM*Bp_lm*R7p_lm)+&
                    ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*DSQRT(PI)/&
                    2.d0)*(M_PM/MPSUM)**2 * M_PM*R8p_lm)+&
                    ((DPSUMo2*DSQRT(PI)/16.d0)*(M_PM*M_PL/(2.d0*MPSUM))*&
                    R9p_lm)+&
                    ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*DPSUMo2*DSQRT(PI)*&
                    (M_PM*M_PL/MPSUM)*Bp_lm*R7p_lm) )*Theta_M(IJK,M) )*&
                    (Theta_M(IJK,L)**2 * Theta_M(IJK,M)**3)

                  Kth_sL_ip(IJK,M,L) = (ONE-UR_Kth_sml)*Kth_sL_ip(IJK,M,L) +&
                    UR_Kth_sml * Kth_sl_iptmp

! solids phase 'conductivity' associated with the
! difference in velocity
                  Kvel_s_ip(IJK,M,L) = K_common_term*NU_PM*NU_PL*&
                    (M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(3.d0*PI/10.d0)*&
                    R0p_lm * (Theta_m(IJK,M)*Theta_m(IJK,L))**2.5

! solids phase 'conductivity' associated with the
! difference in the gradient in number densities
                  Knu_sL_ip(IJK,M,L) = K_common_term*(&
                    ((DPSUMo2*DSQRT(PI)/16.d0)*Bp_lm*R5p_lm)+&
                    ((M_PL/(8.d0*MPSUM))*(1.d0-C_E)*(DPSUMo2*PI/6.d0)*&
                    R1p_lm) )* (Theta_m(IJK,M)*Theta_m(IJK,L))**3

! to avoid recomputing Knu_sL_ip, sof.
                  Knu_sM_ip(IJK,M,L) = NU_PL * Knu_sL_ip(IJK,M,L)
                  Knu_sL_ip(IJK,M,L) = NU_PM * Knu_sL_ip(IJK,M,L)
                  K_s_sum = K_s_sum + K_s_LM
               ENDIF
            ENDDO

! granular conductivity in Mth solids phase
            Kth_s(IJK,M) = K_s_sum

         ENDIF   ! Fluid_at
      ENDDO   ! IJK loop
      RETURN
      END SUBROUTINE GT_PDE_IA



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: FRICTION_SCHAEFFER                                      C
!  Purpose: Calculate frictional-flow stress terms                     C
!                                                                      C
!  Author: M. Syamlal                               Date: 10-FEB-93    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine friction_schaeffer(M)

! Modules
!---------------------------------------------------------------------//

      USE param1, only: zero, small_number

      USE fldvar, only: ep_g, ep_s
      USE fldvar, only: p_star, theta_m

      USE physprop, only: smax, close_packed

      USE constant, only: to_SI
      USE constant, only: sin2_phi

      use run, only: granular_energy

      USE visc_s, only: ep_g_blend_end
      USE visc_s, only: mu_s_p, lambda_s_p
      USE visc_s, only: i2_devD_s

      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local parameters
!---------------------------------------------------------------------//
! maximum value of solids viscosity in poise
      DOUBLE PRECISION :: MAX_MU_s
      PARAMETER (MAX_MU_s = 1000.D0)

! Local variables
!---------------------------------------------------------------------//
! cell index
      INTEGER :: IJK
! solids phase index
      INTEGER :: MM
! sum of all solids volume fractions
      DOUBLE PRECISION :: SUM_EPS_CP
! factor in frictional-flow stress terms
      DOUBLE PRECISION :: qxP_s
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

            IF(EP_g(IJK) .LT. EP_g_blend_end(IJK)) THEN
! Tardos, PT, (1997), pp. 61-74 explains in his equation (3) that
! solids normal and shear frictional stresses have to be treated
! consistently. Therefore, add closed_packed check for consistency
! with the treatment of the normal frictional force (see for example
! source_v_s.f). sof May 24 2005.
               SUM_EPS_CP=0.0
               DO MM=1,SMAX
                  IF (CLOSE_PACKED(MM)) SUM_EPS_CP = SUM_EPS_CP + &
                     EP_S(IJK,MM)
               ENDDO

! plastic pressure (p_star) is calculated seperately to ensure
! its value is up-to-date with latest value of ep_g
!               P_star(IJK) = Neg_H(EP_g(IJK),EP_star_array(IJK))

!-----------------------------------------------------------------------
! Gray and Stiles (1988) frictional-flow stress tensor
!              IF(Sin2_Phi .GT. SMALL_NUMBER) THEN
!                 qxP_s = SQRT( (4. * Sin2_Phi) * I2_devD_s(IJK) + &
!                    trD_s_C(IJK,M) * trD_s_C(IJK,M))
!                 MU_s_p(IJK) = P_star(IJK) * Sin2_Phi/ &
!                    (qxP_s + SMALL_NUMBER)*(EP_S(IJK,M)/SUM_EPS_CP)
!                 MU_s_p(IJK) = MIN(MU_s_p(IJK), to_SI*MAX_MU_s)
!                 LAMBDA_s_p(IJK) = P_star(IJK) * F_Phi /&
!                    (qxP_s + SMALL_NUMBER)*(EP_S(IJK,M)/SUM_EPS_CP)
!                 LAMBDA_s_p(IJK) = MIN(LAMBDA_s(IJK, M), to_SI*MAX_MU_s)
!              ENDIF
!-----------------------------------------------------------------------

! Schaeffer (1987)
               qxP_s = SQRT( (4.D0 * Sin2_Phi) * I2_devD_s(IJK))
! multiply by factor ep_s/sum_eps_cp for consistency with solids
! pressure treatment
               MU_s_p(IJK) = P_star(IJK) * Sin2_Phi/&
                  (qxP_s + SMALL_NUMBER) * (EP_S(IJK,M)/SUM_EPS_CP)
               MU_s_p(IJK) = MIN(MU_s_p(IJK), to_SI*MAX_MU_s)

               LAMBDA_s_p(IJK) = ZERO

! when solving for the granular energy equation (PDE) setting
! theta = 0 is done in solve_granular_energy.f to avoid
! convergence problems. (sof)
               IF(.NOT.GRANULAR_ENERGY) THETA_m(IJK, M) = ZERO
            ENDIF
         ENDIF   ! end if (fluid_at(ijk)

      ENDDO   ! IJK loop
      RETURN
      END SUBROUTINE FRICTION_SCHAEFFER



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: FRICTION_PRINCETON                                      C
!  Purpose: Calculate frictional contributions to granular stress      C
!  terms based on model by Srivastava & Sundaresan (2003)              C
!                                                                      C
!  Author: Anuj Srivastava, Princeton University    Date: 20-APR-98    C
!                                                                      C
!  Literature/Document References:                                     C
!  Srivastava, A., and Sundaresan, S., Analysis of a frictional-       C
!     kinetic model for gas-particle flow, Powder Technology,          C
!     2003, 129, 72-85.                                                C
!  Benyahia, S., Validation study of two continuum granular            C
!     frictional flow theories, Industrial & Engineering Chemistry     C
!     Research, 2008, 47, 8926-8932.                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine Friction_princeton(M)

! Modules
!---------------------------------------------------------------------//
      USE constant, only: alpha
      USE constant, only: eta
      USE constant, only: sqrt_pi
      USE constant, only: to_si
      USE constant, only: sin_phi
      USE constant, only: eps_f_min
      USE constant, only: fr, n_pc, d_pc, n_pf, delta

      USE fldvar, only: ep_g, ep_s
      USE fldvar, only: d_p, ro_s
      USE fldvar, only: theta_m
      USE fldvar, only: p_s_f

      USE param1, only: zero, one, small_number

      USE physprop, only: smax, close_packed

      USE rdf, only: g_0

      USE run, only: kt_type_enum, ghd_2007
      USE run, only: savage
      USE trace, only: trd_s2, trd_s_c
      USE visc_s, only: mu_s_f, lambda_s_f
      USE visc_s, only: mu_s_v, mu_b_v
      USE visc_s, only: ep_star_array

      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell index
      INTEGER :: IJK
! solids phase index
      INTEGER :: MM
! used to compute frictional terms
      DOUBLE PRECISION :: Chi, Pc, Mu_zeta,PfoPc, N_Pff
      DOUBLE PRECISION :: ZETA
! sum of all solids volume fractions
      DOUBLE PRECISION :: SUM_EPS_CP
!     parameters in pressure linearization; simple averaged Dp
      DOUBLE PRECISION :: dpc_dphi, dp_avg
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

            IF (EP_g(IJK) .LT. (ONE-eps_f_min)) THEN

! This close_packed check was added for consistency with the Schaeffer
! model implementation.
               SUM_EPS_CP = ZERO
               dp_avg = ZERO
               DO MM=1,SMAX
                  dp_avg = dp_avg + D_p(IJK,MM)
                  IF (CLOSE_PACKED(MM)) SUM_EPS_CP = SUM_EPS_CP + &
                     EP_S(IJK,MM)
               ENDDO
               dp_avg = dp_avg/DBLE(SMAX)

               IF (SAVAGE.EQ.1) THEN !form of Savage (not to be used with GHD theory)
                  Mu_zeta = &
                     ((2d0+ALPHA)/3d0)*((Mu_s_v(IJK)/(Eta*(2d0-Eta)*&
                     G_0(IJK,M,M)))*(1d0+1.6d0*Eta*EP_s(IJK,M)*&
                     G_0(IJK,M,M))*(1d0+1.6d0*Eta*(3d0*Eta-2d0)*&
                     EP_s(IJK,M)*G_0(IJK,M,M))+(0.6d0*Mu_b_v(IJK)*Eta))
                  ZETA = &
                     ((48d0*Eta*(1d0-Eta)*RO_S(IJK,M)*EP_s(IJK,M)*&
                     EP_s(IJK,M)*G_0(IJK,M,M)*&
                     (Theta_m(IJK,M)**1.5d0))/&
                     (SQRT_Pi*D_p(IJK,M)*2d0*Mu_zeta))**0.5d0

               ELSEIF (SAVAGE.EQ.0) THEN !S:S form
                  ZETA = (SMALL_NUMBER +&
                     trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                     trD_s_C(IJK,M))/3.d0))**0.5d0

               ELSE          !combined form
                  IF(KT_TYPE_ENUM == GHD_2007) THEN
                    ZETA = ((Theta_m(IJK,M)/dp_avg**2) +&
                       (trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                       trD_s_C(IJK,M))/3.d0)))**0.5d0
                  ELSE
                    ZETA = ((Theta_m(IJK,M)/D_p(IJK,M)**2) +&
                       (trD_s2(IJK,M) - ((trD_s_C(IJK,M)*&
                       trD_s_C(IJK,M))/3.d0)))**0.5d0
                  ENDIF
               ENDIF


               IF ((ONE-EP_G(IJK)) .GT. ((ONE-ep_star_array(ijk))-delta)) THEN
! Linearized form of Pc; this is more stable and provides continuous function.
                  dpc_dphi = (to_SI*Fr)*((delta**5)*(2d0*(ONE-&
                      ep_star_array(IJK)-delta) - 2d0*eps_f_min)+&
                      ((ONE-ep_star_array(ijk)-delta)-eps_f_min)*&
                      (5*delta**4))/(delta**10)

                  Pc = (to_SI*Fr)*(((ONE-ep_star_array(IJK)-delta) -&
                     EPS_f_min)**N_Pc)/(delta**D_Pc)
                  Pc = Pc + dpc_dphi*((ONE-EP_G(IJK))+delta-(ONE-&
                     ep_star_array(IJK)))

                ! Pc = 1d25*(((ONE-EP_G(IJK))- (ONE-ep_star_array(ijk)))**10d0) ! old commented Pc
               ELSE
                  Pc = (to_SI*Fr)*(((ONE-EP_G(IJK)) - EPS_f_min)**N_Pc) / &
                     (((ONE-ep_star_array(ijk)) - (ONE-EP_G(IJK)) +&
                     SMALL_NUMBER)**D_Pc)
               ENDIF

               IF (trD_s_C(IJK,M) .GE. ZERO) THEN
                  N_Pff = DSQRT(3d0)/(2d0*Sin_Phi) !dilatation
               ELSE
                  N_Pff = N_Pf !compaction
               ENDIF

               IF ((trD_s_C(IJK,M)/(ZETA*N_Pff*DSQRT(2d0)&
                    *Sin_Phi)) .GT. 1d0) THEN
                  P_s_f(IJK) =ZERO
                  PfoPc = ZERO
               ELSEIF(trD_s_C(IJK,M) == ZERO) THEN
                  P_s_f(IJK) = Pc
                  PfoPc = ONE
               ELSE
                  P_s_f(IJK) = Pc*(1d0 - (trD_s_C(IJK,M)/(ZETA&
                     *N_Pff*DSQRT(2d0)*Sin_Phi)))**(N_Pff-1d0)
                  PfoPc = (1d0 - (trD_s_C(IJK,M)/(ZETA&
                     *N_Pff*DSQRT(2d0)*Sin_Phi)))**(N_Pff-1d0)
               ENDIF

               Chi = DSQRT(2d0)*P_s_f(IJK)*Sin_Phi*(N_Pff - (N_Pff-1d0)*&
                  (PfoPc)**(1d0/(N_Pff-1d0)))

               IF (Chi < ZERO) THEN
                  P_s_f(IJK) = Pc*((N_Pff/(N_Pff-1d0))**(N_Pff-1d0))
                  Chi = ZERO
               ENDIF

               Mu_s_f(IJK) = Chi/(2d0*ZETA)
               Lambda_s_f(IJK) = -2d0*Mu_s_f(IJK)/3d0

! Modification of the stresses when more than one solids phase is used:
! This is NOT done when mixture mom. eq. is solved (i.e. for GHD theory)
               IF(KT_TYPE_ENUM /= GHD_2007) THEN
                  P_s_f(IJK) = P_s_f(IJK) * (EP_S(IJK,M)/SUM_EPS_CP)
                  Mu_s_f(IJK) = Mu_s_f(IJK) * (EP_S(IJK,M)/SUM_EPS_CP)
                  Lambda_s_f(IJK) = Lambda_s_f(IJK) * (EP_S(IJK,M)/SUM_EPS_CP)
               ENDIF

            ENDIF   ! end if ep_g < 1-eps_f_min
         ENDIF   ! Fluid_at
      ENDDO   ! IJK loop
      RETURN
      END SUBROUTINE FRICTION_PRINCETON



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_STRESS_IGCI                                     C
!  Purpose: Calculate solids viscosity and pressure using subgrid      C
!  model                                                               C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Literature/Document References:                                     C
!  Igci, Y., Pannala, S., Benyahia, S., & Sundaresan S.,               C
!     Validation studies on filtered model equations for gas-          C
!     particle flows in risers, Industrial & Engineering Chemistry     C
!     Research, 2012, 51(4), 2094-2103                                 C
!                                                                      C
!  Comments:                                                           C
!  Still needs to be reviewed for accuracy with source material        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine subgrid_stress_igci(M)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, one, small_number

      USE fldvar, only: ep_s
      USE fldvar, only: d_p, ro_s
      USE fldvar, only: ro_g
      USE fldvar, only: p_s_v
! this shouldn't be necessary
      use fldvar, only: theta_m

      USE visc_s, only: mu_s_v, lambda_s_v

      USE physprop, only: mu_g

      USE run, only: filter_size_ratio
      USE run, only: subgrid_wall

      USE constant, only: gravity

      use geometry, only: vol, axy, do_k
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: IJK
! Igci models
      DOUBLE PRECISION :: pressurefac, ps_kinetic, ps_total, ps_sub
      DOUBLE PRECISION :: viscosityfac, extra_visfactor, mu_sub
      DOUBLE PRECISION :: mu_kinetic, mu_total
! factors to correct subgrid effects from the wall
      DOUBLE PRECISION :: wfactor_Ps, wfactor_mus
! the inverse Froude number, or dimensionless Filtersize
      DOUBLE PRECISION :: Inv_Froude
! one particle terminal settling velocity
      DOUBLE PRECISION :: vt
! the filter size which is a function of each grid cell volume
      DOUBLE PRECISION :: filtersize
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

! initialize
            wfactor_Ps = ONE   ! for P_S
            wfactor_mus = ONE   ! for MU_s

! particle terminal settling velocity: vt = g*d^2*(Rho_s - Rho_g) / 18 * Mu_g
            vt = GRAVITY*D_p(IJK,M)*D_p(IJK,M)*&
                 (RO_S(IJK,M) - RO_g(IJK)) / (18.0d0*MU_G(IJK))

! FilterSIZE calculation for each specific gridcell volume
            IF(DO_K) THEN
               filtersize = filter_size_ratio * (VOL(IJK)**(ONE/3.0d0))
            ELSE
               filtersize = filter_size_ratio * DSQRT(AXY(IJK))
            ENDIF


! dimensionless inverse of Froude number
            Inv_Froude =  filtersize * GRAVITY / vt**2

! various factor needed:
            pressurefac = 0.48d0*(Inv_Froude**0.86)*&
               (ONE-EXP(-Inv_Froude/1.4))
            viscosityfac = 0.37d0*(Inv_Froude**1.22)
            Extra_visfactor = ONE/(0.28d0*(Inv_Froude**0.43)+ONE)

            IF (EP_s(IJK,M) .LE. 0.0131) THEN
               Ps_kinetic = -10.4d0*(EP_s(IJK,m)**2)+0.31d0*EP_s(IJK,m)
            ELSEIF (EP_s(IJK,M) .LE. 0.290) THEN
               Ps_kinetic = -0.185d0*(EP_s(IJK,m)**3)+&
                  0.066d0*(EP_s(IJK,m)**2)-0.000183d0*EP_s(IJK,m)+&
                  0.00232d0
            ELSEIF (EP_s(IJK,M) .LE. 0.595) THEN
               Ps_kinetic = -0.00978d0*EP_s(IJK,m)+0.00615d0
            ELSE
               Ps_kinetic = -6.62d0*(EP_s(IJK,m)**3)+&
                  49.5d0*(EP_s(IJK,m)**2)-50.3d0*EP_s(IJK,m)+13.8d0
            ENDIF

            Ps_sub = pressurefac*(EP_s(IJK,M)-0.59d0)*&
               (-1.69d0*EP_s(IJK,M)-4.61d0*(EP_s(IJK,M)**2)+&
               11.d0*(EP_s(IJK,M)**3))

            IF (Ps_sub .GE. ZERO) THEN
               Ps_total=Ps_kinetic+Ps_sub
            ELSE
               Ps_total=Ps_kinetic
            ENDIF

            IF (EP_s(IJK,M) .LE. 0.02) THEN
               Mu_kinetic = 1720.d0*(EP_s(IJK,m)**4)-&
                  215.d0*(EP_s(IJK,m)**3) + 9.81d0*(EP_s(IJK,m)**2)-&
                  0.207d0*EP_s(IJK,m)+0.00254d0
            ELSEIF (EP_s(IJK,M) .LE. 0.2) THEN
               Mu_kinetic = 2.72d0*(EP_s(IJK,m)**4)-&
                  1.55d0*(EP_s(IJK,m)**3)+0.329d0*(EP_s(IJK,m)**2)-&
                  0.0296d0*EP_s(IJK,m)+0.00136d0
            ELSEIF (EP_s(IJK,M) .LE. 0.6095) THEN
               Mu_kinetic = -0.0128d0*(EP_s(IJK,m)**3)+&
                  0.0107d0*(EP_s(IJK,m)**2)-0.0005d0*EP_s(IJK,m)+&
                  0.000335d0
            ELSE
               Mu_kinetic = 23.6d0*(EP_s(IJK,m)**2)-&
                  28.0d0*EP_s(IJK,m)+8.30d0
            ENDIF

            Mu_sub = Extra_visfactor*viscosityfac*&
               (EP_s(IJK,M)-0.59d0)*(-1.22d0*EP_s(IJK,M)-&
               0.7d0*(EP_s(IJK,M)**2)-2.d0*(EP_s(IJK,M)**3))

            IF (Mu_sub .GE. ZERO) THEN
               Mu_total = Mu_kinetic+Mu_sub
            ELSE
               Mu_total = Mu_kinetic
            ENDIF

            IF (SUBGRID_WALL) THEN
               CALL SUBGRID_STRESS_WALL(wfactor_Ps,wfactor_Mus,vt,IJK)
            ENDIF

! pressure
            P_s_v(IJK) = Ps_total * wfactor_Ps * (vt**2) * &
               RO_S(IJK,M)

! shear viscosity
            Mu_s_v(IJK) = Mu_total * wfactor_mus * (vt**3) * &
               RO_S(IJK,M)/GRAVITY

! set an arbitrary value in case value gets negative (this should
! not happen unless filtersize becomes unrelastic w.r.t. gridsize)
            IF (P_s_v(IJK) .LE. SMALL_NUMBER) P_s_v(IJK) = SMALL_NUMBER
            IF (Mu_s_v(IJK) .LE. SMALL_NUMBER) Mu_s_v(IJK)= SMALL_NUMBER

! Solid second viscosity, assuming the bulk viscosity is ZERO
            lambda_s_v(IJK) = (-2.0d0/3.0d0)*Mu_s_v(IJK)

! granular temperature is zeroed in all LES/Subgrid model
! why is this necessary? theta_M shouldn't be used anyway!
            THETA_m(IJK, M) = ZERO

         ENDIF   ! endif (fluid_at(IJK))
      ENDDO   ! outer IJK loop
      RETURN
      END SUBROUTINE SUBGRID_STRESS_IGCI



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_STRESS_MILIOLI                                  C
!  Purpose: Calculate solids viscosity and pressure using subgrid      C
!  model                                                               C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Literature/Document References:                                     C
!  Milioli, C. C., et al., Filtered two-fluid models of fluidized      C
!     gas-particle flows: new constitutive relations, AICHE J,         C
!     doi: 10.1002/aic.14130                                           C
!                                                                      C
!  Comments:                                                           C
!  Still needs to be reviewed for accuracy with source material        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine subgrid_stress_MILIOLI(M)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero, one, small_number

      USE fldvar, only: ep_g, ep_s
      USE fldvar, only: d_p, ro_s
      USE fldvar, only: ro_g
      USE fldvar, only: p_s_v
! this shouldn't be necessary
      use fldvar, only: theta_m

      USE visc_s, only: mu_s_v, lambda_s_v
      USE visc_s, only: i2_devD_s

      USE physprop, only: mu_g

      USE run, only: filter_size_ratio
      USE run, only: subgrid_wall

      USE constant, only: gravity

      use geometry, only: vol, axy, do_k
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell indices
      INTEGER :: IJK
! Milioli model
      DOUBLE PRECISION :: cvisc_pot, cvisc_num, cvisc_den, Cvisc
      DOUBLE PRECISION :: cpress_pot, cpress_num, cpress_den, Cpress
! factor functions to correct subgrid effects from the wall
      DOUBLE PRECISION :: wfactor_Ps, wfactor_mus
! the inverse Froude number, or dimensionless Filtersize
      DOUBLE PRECISION :: Inv_Froude
! one particle terminal settling velocity
      DOUBLE PRECISION :: vt
! the filter size which is a function of each grid cell volume
      DOUBLE PRECISION :: filtersize
!---------------------------------------------------------------------//

      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

! initialize
            wfactor_Ps = ONE   ! for P_S
            wfactor_mus = ONE   ! for MU_s

! particle terminal settling velocity: vt = g*d^2*(Rho_s - Rho_g) / 18 * Mu_g
            vt = GRAVITY*D_p(ijk,M)*D_p(ijk,M)*&
                 (RO_S(IJK,M)-RO_g(IJK)) / (18.0d0*MU_G(IJK))

! FilterSIZE calculation for each specific gridcell volume
            IF(DO_K) THEN
               filtersize = filter_size_ratio * (VOL(IJK)**(ONE/3.0d0))
            ELSE
               filtersize = filter_size_ratio * DSQRT(AXY(IJK))
            ENDIF

! dimensionless inverse of Froude number
            Inv_Froude =  filtersize * GRAVITY / vt**2

! Cvisc:
            cvisc_pot = (0.59d0-(ONE-EP_g(IJK)))
            cvisc_num = (0.7d0*(ONE-EP_g(IJK))*cvisc_pot)
            cvisc_den = (0.8d0+(17.d0*cvisc_pot*cvisc_pot*cvisc_pot))
            IF ((ONE-EP_g(IJK)) .GE. ZERO .AND. &
                (ONE-EP_g(IJK)) .LE. 0.59) THEN
               cvisc=(cvisc_num/cvisc_den)
            ELSE
               cvisc=ZERO
            ENDIF

! aCvisc:
!            IF ((ONE-EP_g(IJK)) .GE. ZERO .AND. &
!                (ONE-EP_g(IJK)) .LE. 0.59) THEN
!               acvisc(IJK) = (0.7d0*(ONE-EP_g(IJK))*(0.59d0-&
!                  (ONE-EP_g(IJK))))/(0.8d0+17.d0*(0.59d0-&
!                  (ONE-EP_g(IJK)))*(0.59d0-(ONE-EP_g(IJK)))*&
!                  (0.59d0-(ONE-EP_g(IJK))))
!            ELSE
!               acvisc(IJK)=ZERO
!            ENDIF

! Cpress:
            cpress_pot = (0.59d0-(ONE-EP_g(IJK)))
            cpress_num = (0.4d0*(ONE-EP_g(IJK))*cpress_pot)
            cpress_den = (0.5d0+(13.d0*cpress_pot*cpress_pot*cpress_pot))
            IF ((ONE-EP_g(IJK)) .GE. ZERO .AND. &
                (ONE-EP_g(IJK)) .LE. 0.59) THEN
               cpress = (cpress_num/cpress_den)
            ELSE
               cpress = ZERO
            ENDIF

! aCpress
!            IF ((ONE-EP_g(IJK)) .GE. ZERO .AND. &
!                (ONE-EP_g(IJK)) .LE. 0.59) THEN
!                acpress(IJK) = (0.4d0*(ONE-EP_g(IJK))*(0.59d0-&
!                   (ONE-EP_g(IJK))))/(0.5d0+13.d0*(0.59d0-&
!                   (ONE-EP_g(IJK)))*(0.59d0-(ONE-EP_g(IJK)))*&
!                   (0.59d0-(ONE-EP_g(IJK))))
!            ELSE
!            acpress(IJK)=ZERO
!            ENDIF

            IF (SUBGRID_WALL) THEN
               CALL SUBGRID_STRESS_WALL(wfactor_Ps,wfactor_Mus,vt,IJK)
            ENDIF

! solid filtered pressure
            P_s_v(IJK) = RO_S(IJK,M) * Inv_froude**(2/7) * &
               filtersize**2 * DSQRT( I2_devD_s(IJK) )**2 * &
               cpress * wfactor_Ps    !16/7-2=2/7 in [Pa or kg/m.s2]

! solids filtered shear viscosity
            Mu_s_v(IJK) = RO_S(IJK,M) * filtersize**2 * &
               DSQRT( I2_devD_s(IJK) ) * cvisc * wfactor_mus  ! [kg/m.s]

! set an arbitrary value in case value gets negative (this should
! not happen unless filtersize becomes unrealistic w.r.t. gridsize)
            IF (P_s_v(IJK) .LE. SMALL_NUMBER) P_s_v(IJK) = SMALL_NUMBER
            IF (Mu_s_v(IJK) .LE. SMALL_NUMBER) Mu_s_v(IJK)= SMALL_NUMBER

! Solids second viscosity, assuming the bulk viscosity is ZERO
            lambda_s_v(IJK) = (-2.0d0/3.0d0)*Mu_s_v(IJK)

! granular temperature is zeroed in all LES/Subgrid model
! why is this necessary? theta_M shouldn't be used anyway!
            THETA_m(IJK, M) = ZERO

         ENDIF   ! endif (fluid_at(IJK))
      ENDDO   ! outer IJK loop
      RETURN
      END SUBROUTINE SUBGRID_STRESS_MILIOLI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SUBGRID_STRESS_WALL                                     C
!  Purpose: Calculate subgrid corrections arising from wall to solids  C
!  viscosity and pressure.                                             C
!                                                                      C
!  Author: Sebastien Dartevelle, LANL, May 2013                        C
!                                                                      C
!  Literature/Document References:                                     C
!  Igci, Y., and Sundaresan, S., Verification of filtered two-         C
!     fluid models for gas-particle flows in risers, AICHE J.,         C
!     2011, 57 (10), 2691-2707.                                        C
!                                                                      C
!  Comments: Currently only valid for free-slip wall but no checks     C
!  are made to ensure user has selected free-slip wall when this       C
!  option is invoked                                                   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine subgrid_stress_wall(lfactor_ps, lfactor_mus, vt, &
                 IJK)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: ONE
      USE constant, only : GRAVITY
      USE cutcell, only : DWALL
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! factor to correct the solids pressure
      DOUBLE PRECISION, INTENT(OUT) :: lfactor_ps
! factor to correct the solids viscosity
      DOUBLE PRECISION, INTENT(OUT) :: lfactor_mus
! one particle terminal settling velocity
      DOUBLE PRECISION, INTENT(IN) :: vt
! current ijk index
      INTEGER, INTENT(IN) :: IJK

! Local parameters
!---------------------------------------------------------------------//
! values are only correct for free-slip walls
      DOUBLE PRECISION, PARAMETER :: aps=9.14d0, bps=0.345d0,&
                                     amus=5.69d0, bmus=0.228d0

! Local variables
!---------------------------------------------------------------------//
! dimensionless distance to the wall
      DOUBLE PRECISION :: x_d
!---------------------------------------------------------------------//

! initialize
      lfactor_ps = ONE
      lfactor_mus = ONE

! dimensionless distance to the Wall
      x_d = DWALL(IJK) * GRAVITY / vt**2
! wall function for pressure
      lfactor_Ps = ONE / ( ONE + aps * (EXP(-bps*x_d)) )
! wall function for viscosity
      lfactor_mus = ONE / ( ONE + amus * (EXP(-bmus*x_d)) )

      RETURN
      END SUBROUTINE SUBGRID_STRESS_WALL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine add_shear(M)

! Modules
!---------------------------------------------------------------------//
! y-component of solids velocity
      USE fldvar, only: v_s
! shear velocity
      USE vshear, only: vsh
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell index
      INTEGER :: IJK
!---------------------------------------------------------------------//

!!     $omp parallel do private(IJK)
      DO IJK= ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            V_s(ijk,m)=V_s(IJK,m)+VSH(IJK)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE ADD_SHEAR



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine remove_shear(M)

! Modules
!---------------------------------------------------------------------//
! y-component of solids velocity
      USE fldvar, only: v_s
! shear velocity
      USE vshear, only: vsh
! needed for function.inc
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! cell index
      INTEGER :: IJK
!---------------------------------------------------------------------//

!!     $omp parallel do private(IJK)
      DO IJK= ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN
            V_s(IJK,m)=V_s(IJK,m)-VSH(IJK)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE REMOVE_SHEAR


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
! Subroutine: INIT0_MU_S                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine init0_mu_s (M)

! Modules
!---------------------------------------------------------------------//
      USE param1, only: zero
      USE constant, only: v_ex
      USE fldvar, only: p_s_v, p_s_p, p_s_f
      USE visc_s, only: mu_s_v, mu_s_p, mu_s_f, mu_b_v
      USE visc_s, only: lambda_s_v, lambda_s_p, lambda_s_f
      USE visc_s, only: alpha_s
!---------------------------------------------------------------------//
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M
!---------------------------------------------------------------------//

      IF (V_EX /= ZERO)  ALPHA_s(:,M) = ZERO

! viscous contributions
      Mu_s_v(:)     = ZERO
      Mu_b_v(:)     = ZERO
      LAMBDA_s_v(:) = ZERO
      P_s_v(:) = ZERO
! plastic contributions (local to this routine)
      Mu_s_p(:)     = ZERO
      LAMBDA_s_p(:) = ZERO   ! never used
      P_s_p(:) = ZERO        ! never used
! frictional contributions
      Mu_s_f(:)     = ZERO
      LAMBDA_s_f(:) = ZERO
      P_s_f(:) = ZERO

      RETURN
      END SUBROUTINE init0_mu_s


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
! Subroutine: INIT1_MU_S                                               C
! Purpose: Calculate quantities that are functions of velocity only!   C
!  - strain rate tensor: D_s (local)                                   C
!  - trace of strain rate tensor: trd_s_c (global)                     C
!    this quantity seems equivalent to the variable trd_s that is      C
!    calculated by the subroutine calc_trd_s (trd_s). formulation is   C
!    equivalent despite apparent differences:                          C
!       trd (D) = du/dx + dv/dy + 1/x dw/dz + u/x (here)               C
!               = 1/x d(xu)/dx + dv/dy + 1/x dw_dz (calc_trd_s)        C
!  - trace of square of strain rate tensor: trd_s2 (global)            C
!  - second invariant of the deviatoric stess tensor: i2_devD_s        C
!    (global)                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      Subroutine init1_mu_s (M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE constant, only: v_ex
      USE cutcell, only: cut_cell_at
      USE functions, only: fluid_at
      USE geometry, only: cylindrical
      USE physprop, only: smax
      USE param1, only: zero, one, half
      USE run, only: kt_type_enum
      USE run, only: ia_2005
      USE run, only: shear
      USE trace, only: trD_s_c
      USE trace, only: trD_s2
      USE trace, only: trD_s2_ip
      USE visc_s, only: I2_devD_s
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase index
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! solids phase index
      INTEGER :: L
! Indices
      INTEGER :: IJK
! Local DO-LOOP counters and phase index
      INTEGER :: I1, I2
! Strain rate tensor components for mth and lth solids phases
      DOUBLE PRECISION :: D_sM(3,3), D_sl(3,3)
! velocity gradient for mth and lth solids phases
      DOUBLE PRECISION :: DelV_sM(3,3), DelV_sL(3,3)
! shear related reciprocal time scale
      DOUBLE PRECISION :: SRT
!---------------------------------------------------------------------//

      IF (SHEAR) THEN
         CALL add_shear(M)
         IF (KT_TYPE_ENUM == IA_2005) THEN
            DO L = 1,SMAX
               IF (L .NE. M) THEN
                  CALL add_shear(L)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

!!$omp  parallel do default(shared) &
!!$omp  private( IJK, L, I1, I2, D_sM, D_sL, DelV_sM, DelV_sL)
      DO IJK = ijkstart3, ijkend3
         IF ( FLUID_AT(IJK) ) THEN

! Velocity derivatives (gradient and rate of strain tensor) for Mth
! solids phase at i, j, k
            CALL CALC_DERIV_VEL_SOLIDS(IJK, M, DelV_sM, D_sM)

! Calculate the trace of D_s
! friction/algebraic various kt
            trD_s_C(IJK,M) = D_sM(1,1) + D_sM(2,2) + D_sM(3,3)

! Calculate trace of the square of D_s
! friction/algebraic various kt
            trD_s2(IJK,M) = ZERO ! initialize the totalizer
            DO I1 = 1,3
               DO I2 = 1,3
                  trD_s2(IJK,M) = trD_s2(IJK,M) + D_sM(I1,I2)*D_sM(I1,I2)
               ENDDO
            ENDDO

! use this fact to prevent underflow during theta calculation
            IF (trD_s2(IJK,M) == zero) trD_s_C(IJK,M) = zero

! Frictional-flow stress tensor
! Needed by schaeffer or friction
! Calculate the second invariant of the deviator of D_s
            I2_devD_s(IJK) = ( (D_sM(1,1)-D_sM(2,2))**2 + &
                               (D_sM(2,2)-D_sM(3,3))**2 + &
                               (D_sM(3,3)-D_sM(1,1))**2 )/6.&
                 + D_sM(1,2)**2 + D_sM(2,3)**2 + D_sM(3,1)**2

            IF (V_EX /= ZERO) CALL CALC_BOYLE_MASSOUDI_STRESS(IJK,M,D_sM)

! Quantities for iddir-arastoopour theory: the trace of D_sm dot D_sl
! is required
! -------------------------------------------------------------------
            IF (KT_TYPE_ENUM == IA_2005) THEN
               DO L = 1,SMAX
                  IF (L .NE. M) THEN
                     IF (L > M) THEN  !done because trD_s2_ip(IJK,M,L) is symmetric, sof.

! Velocity derivatives (gradient and rate of strain tensor) for Mth
! solids phase at i, j, k
                        CALL CALC_DERIV_VEL_SOLIDS(IJK, L, DelV_sL, D_sL)

! Calculate trace of the D_sl dot D_sm
! (normal matrix multiplication)
                        trD_s2_ip(IJK,M,L) = ZERO ! initialize
                        DO I1 = 1,3
                           DO I2 = 1,3
                              trD_s2_ip(IJK,M,L) = trD_s2_ip(IJK,M,L)+&
                                 D_sL(I1,I2)*D_sM(I1,I2)
                           ENDDO
                        ENDDO

                     ELSE   ! elseif (m<=m)
                        trD_s2_ip(IJK,M,L) = trD_s2_ip(IJK,L,M)
                     ENDIF ! end if/else (l>m)
                  ELSE   ! elseif (L=M)
                     trD_s2_ip(IJK,M,M) = trD_s2(IJK,M)
                  ENDIF   ! end if/else (m.NE.l)
               ENDDO   ! end do (l=1,smax)
            ENDIF   ! endif (kt_type = IA theory)


         ENDIF   ! end if (fluid_at)
      ENDDO   ! end outer IJK loop
!!$omp end parallel do

      IF (SHEAR) THEN
         call remove_shear(M)
         IF (KT_TYPE_ENUM == IA_2005) THEN
            DO L = 1,SMAX
               IF (L .NE. M) THEN
                  CALL remove_shear(L)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE INIT1_MU_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: calc_boyle_massoudi_stress                              C
!  Purpose: Define quantities needed for viscosity calculations        C
!  pertaining to the boyle-massoudi stress term.                       C
!  Notes: this is applicable within the algebraic granular energy      C
!  equation only.                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_BOYLE_MASSOUDI_STRESS(IJK, M, D_S)

! Modules
!---------------------------------------------------------------------//
      USE fldvar, only: ep_g, ep_s
      USE param1, only: zero, one
      USE visc_s, only: ep_star_array
      USE visc_s, only: trm_s, trdM_s

      USE geometry, only: odx_e, ody_n, odz_t, ox, x
      USE indices, only: i_of, j_of, k_of
      USE indices, only: im1, jm1, km1
      USE functions, only: west_of, east_of, south_of, north_of
      USE functions, only: top_of, bottom_of, fluid_at
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! ijk index
      INTEGER, INTENT(IN) :: IJK
! solids phase index
      INTEGER, INTENT(IN) ::M
! rate of strain tensor at i, j, k
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: D_S

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IM, JM, KM
      INTEGER :: IJKW, IJKE, IJKS, IJKN, IJKB, IJKT
! Solids volume fraction gradient tensor
      DOUBLE PRECISION :: M_s(3,3)
! d(EP_sm)/dX
      DOUBLE PRECISION :: DEP_soDX
! d(EP_sm)/dY
      DOUBLE PRECISION :: DEP_soDY
! d(EP_sm)/XdZ
      DOUBLE PRECISION :: DEP_soXDZ
! Local DO-LOOP counters
      INTEGER :: I1, I2
!---------------------------------------------------------------------//

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IM = Im1(I)
      JM = Jm1(J)
      KM = Km1(K)
      IJKW  = WEST_OF(IJK)
      IJKE  = EAST_OF(IJK)
      IJKS  = SOUTH_OF(IJK)
      IJKN  = NORTH_OF(IJK)
      IJKB  = BOTTOM_OF(IJK)
      IJKT  = TOP_OF(IJK)

      IF(EP_g(IJK) .GE. EP_star_array(IJK)) THEN
         DEP_soDX  = ( EP_s(IJKE, M) - EP_s(IJK, M) ) * oDX_E(I)&
             * ( ONE / ( (oDX_E(IM)/oDX_E(I)) + ONE ) ) +&
             ( EP_s(IJK, M) - EP_s(IJKW, M) ) * oDX_E(IM)&
             * ( ONE / ( (oDX_E(I)/oDX_E(IM)) + ONE ) )
            DEP_soDY  = ( EP_s(IJKN, M) - EP_s(IJK, M) ) * oDY_N(J)&
             * ( ONE / ( (oDY_N(JM)/oDY_N(J)) + ONE ) ) +&
             ( EP_s(IJK, M) - EP_s(IJKS, M) ) * oDY_N(JM)&
             * ( ONE / ( (oDY_N(J)/oDY_N(JM)) + ONE ) )
         DEP_soXDZ  = (( EP_s(IJKT, M) - EP_s(IJK, M) )&
             * oX(I)*oDZ_T(K)&
             * ( ONE / ( (oDZ_T(KM)/oDZ_T(K)) + ONE ) ) +&
             ( EP_s(IJK, M) - EP_s(IJKB, M) ) * oX(I)*oDZ_T(KM)&
             * ( ONE / ( (oDZ_T(K)/oDZ_T(KM)) + ONE ) ) ) /&
             X(I)

          M_s(1,1) = DEP_soDX * DEP_soDX
          M_s(1,2) = DEP_soDX * DEP_soDY
          M_s(1,3) = DEP_soDX * DEP_soXDZ
          M_s(2,1) = DEP_soDX * DEP_soDY
          M_s(2,2) = DEP_soDY * DEP_soDY
          M_s(2,3) = DEP_soDY * DEP_soXDZ
          M_s(3,1) = DEP_soDX * DEP_soXDZ
          M_s(3,2) = DEP_soDY * DEP_soXDZ
          M_s(3,3) = DEP_soXDZ * DEP_soXDZ

          trM_s(IJK) = M_s(1,1) + M_s(2,2) + M_s(3,3)
          trDM_s(IJK) = ZERO
          DO I1 = 1,3
             DO I2 = 1,3
                trDM_s(IJK) = trDM_s(IJK) + D_s(I1,I2)*M_s(I1,I2)
             ENDDO
          ENDDO
      ENDIF   ! end if (ep_g >=ep_star_array)

      RETURN
      END SUBROUTINE CALC_BOYLE_MASSOUDI_STRESS
