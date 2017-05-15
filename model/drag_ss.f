!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_ss                                                 C
!  Purpose: This module computes the coefficient of drag between       C
!  two solids phases (M and L).                                        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAY-96  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DRAG_SS(L, M)

! Modules
!---------------------------------------------------------------------//
      use compar, only: mype
      use exit, only: mfix_exit
      use functions, only: funlm
      use run, only: granular_energy, kt_type, kt_type_enum
      use run, only: lun_1984, simonin_1996, ahmadi_1995
      use run, only: ia_2005, gd_1999, ghd_2007, gtsh_2012
      use usr_prop, only: usr_fss, calc_usr_prop
      use usr_prop, only: solidssolids_drag
      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Index of solids phases
      INTEGER, INTENT(IN) :: L, M

! Local Variables
!---------------------------------------------------------------------//
! index for storing solids-solids drag coefficients in the upper
! triangle of the matrix
      INTEGER :: LM

!---------------------------------------------------------------------//
      LM = FUNLM(L,M)

      IF (USR_FSS(LM)) THEN
         CALL CALC_USR_PROP(SolidsSolids_Drag,lL=L,lM=M)
      ELSE
         IF (GRANULAR_ENERGY) THEN
            SELECT CASE (KT_TYPE_ENUM)
               CASE(LUN_1984, SIMONIN_1996, AHMADI_1995)
                  CALL DRAG_SS_SYAM(L, M)
               CASE (IA_2005)
                  CALL DRAG_SS_IA(L, M)
               CASE (GD_1999, GTSH_2012)
! strictly speaking gd and gtsh are monodisperse theories and so
! do not have solids-solids drag
                  RETURN
! ghd theory is a polydisperse theory but is self contained and does
! not invoke this routine
               CASE (GHD_2007)
                  RETURN
            CASE DEFAULT
! should never hit this
               WRITE (*, '(A)') 'DRAG_SS'
               WRITE (*, '(A,A)') 'Unknown KT_TYPE: ', KT_TYPE
            call mfix_exit(myPE)
            END SELECT
         ELSE
            CALL DRAG_SS_SYAM(L, M)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE DRAG_SS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_SS_SYAM                                            C
!  Purpose: Calculate the solids-solids drag coefficient between a     C
!           continuous solids phase and discrete solids                C
!                                                                      C
!  Literature/Document References:                                     C
!     M. Syamlal. 1987. The particle-particle drag term in a           C
!        multiparticle model of fluidization. Technical Report.        C
!        DOE/MC/21353-2373. Office of Fossil Energy, Morgantown        C
!        Energy Technology Center, Morgantown, West Virginia.          C
!     Gera, D., Syamlal, M., O'Brien T.J. 2004. International Journal  C
!        of Multiphase Flow, 30, p419-428.                             C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
      SUBROUTINE DRAG_SS_SYAM(L,M)

! Modules
!---------------------------------------------------------------------//
      use compar, only: mype, ijkstart3, ijkend3
      USE constant, only: segregation_slope_coefficient
      USE drag, only: f_ss
      USE fldvar, only: d_p, ro_s, rop_s, theta_m, p_star
      USE fldvar, only: u_s, v_s, w_s
      USE fun_avg, only: avg_x_e, avg_y_n, avg_z_t
      use functions, only: wall_at, funlm
      use functions, only: im_of, jm_of, km_of
      USE indices, only: i_of
      USE physprop, only: close_packed
      USE rdf, only: g_0
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase indices
      INTEGER, INTENT(IN) :: M, L

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: I, IJK, IMJK, IJMK, IJKM
! index for storing solids-solids drag coefficients in the upper
! triangle of the matrix
      INTEGER :: LM
! cell center value of U_sm, U_sl, V_sm, V_sl, W_sm, W_sl
      DOUBLE PRECISION :: USCM, USCL, VSCM, VSCL, WSCM, WSCL
! relative velocity between solids phase m and l
      DOUBLE PRECISION :: VREL
! particle diameters of phase M and phase L
      DOUBLE PRECISION :: D_pm, D_pl
! particle densities of phase M and phase L
      DOUBLE PRECISION :: RO_M, RO_L
! radial distribution function between phases M and L
      DOUBLE PRECISION :: G0_ML
! solid-solid drag coefficient
      DOUBLE PRECISION :: lDss

!---------------------------------------------------------------------//
      LM = FUNLM(L,M)

      DO IJK = ijkstart3, ijkend3
         IF (WALL_AT(IJK)) CYCLE
! Evaluate at all flow boundaries and fluid cellls. This is unlike
! calculation of the fluid-solid drag coefficient, which is only
! evaluated in fluid cells and pressure inflow cells

         I = I_OF(IJK)
         IMJK = IM_OF(IJK)
         IJMK = JM_OF(IJK)
         IJKM = KM_OF(IJK)

! calculating velocity components at i, j, k (cell center)
         USCL = AVG_X_E(U_S(IMJK,L),U_S(IJK,L),I)
         VSCL = AVG_Y_N(V_S(IJMK,L),V_S(IJK,L))
         WSCL = AVG_Z_T(W_S(IJKM,L),W_S(IJK,L))

         USCM = AVG_X_E(U_S(IMJK,M),U_S(IJK,M),I)
         VSCM = AVG_Y_N(V_S(IJMK,M),V_S(IJK,M))
         WSCM = AVG_Z_T(W_S(IJKM,M),W_S(IJK,M))

! magnitude of solids-solids relative velocity
         VREL = SQRT((USCL - USCM)**2 + (VSCL - VSCM)**2 + &
                     (WSCL - WSCM)**2)

! setting aliases for easy reference
         D_PM = D_P(IJK,M)
         D_PL = D_P(IJK,L)
         RO_M = RO_S(IJK,M)
         RO_L = RO_S(IJK,L)
         G0_ML = G_0(IJK,L,M)

! evaluating the solids-solids drag coefficient
         CALL DRAG_SS_SYAM0(ldss, d_pm, d_pl, ro_m, ro_l, g0_ml, vrel)
         F_SS(IJK,LM) = lDss*ROP_S(IJK,M)*ROP_S(IJK,L)

! Gera: accounting for particle-particle drag due to enduring contact
! in a close-packed system.
         IF(CLOSE_PACKED(M) .AND. CLOSE_PACKED(L)) &
            F_SS(IJK,LM) = F_SS(IJK,LM) + &
               SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK)

      ENDDO    ! end do (ijk=ijkstart3,ijkend3)

      RETURN
      END SUBROUTINE DRAG_SS_SYAM


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DRAG_SS_SYAM0                                           C
!  Purpose: Created as a work-around for des/hybrid cases that need    C
!  to use this model while passing particle type information.          C
!                                                                      C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
      SUBROUTINE DRAG_SS_SYAM0(ldss, d_pm, d_pl, ro_m, ro_l, g0_ml, vrel)

! Modules 
!---------------------------------------------------------------------//
      use param1, only: one
      use constant, only: c_e, c_f, pi
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solid-solid drag coefficient
      DOUBLE PRECISION, INTENT(OUT) :: lDss
! particle diameters of phase M and phase L
      DOUBLE PRECISION, INTENT(IN) :: D_pm, D_pl
! particle densities of phase M and phase L
      DOUBLE PRECISION, INTENT(IN) :: RO_M, RO_L
! radial distribution function between phases M and L
      DOUBLE PRECISION, INTENT(IN) :: G0_ML
! relative velocity between solids phase m and l
      DOUBLE PRECISION, INTENT(IN) :: VREL

! Local variables
!---------------------------------------------------------------------//
! Sum of particle diameters
      DOUBLE PRECISION :: DPSUM
! Intermediate calculation
      DOUBLE PRECISION :: const
!---------------------------------------------------------------------//

      DPSUM = D_PL + D_PM
      const = 3.d0*(ONE + C_E)*(PI/2.d0 + C_F*PI*PI/8.d0)*&
              DPSUM**2/(2.d0*PI*(RO_L*D_PL**3+RO_M*D_PM**3))
      ldss = const * G0_ML * VREL

      RETURN
      END SUBROUTINE DRAG_SS_SYAM0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DRAG_SS_IA                                              !
!  Purpose: Compute collisional momentum source terms betweem solids   !
!  phase M and solids phase L using Iddir Arastoopour (2005) kinetic   !
!  theory mdoel that is proportional to the relative velocity between  !
!  between the two phases (i.e. solids-solids drag).                   !
!                                                                      !
!  Literature/Document References:                                     !
!    Iddir, Y.H., "Modeling of the multiphase mixture of particles     !
!       using the kinetic theory approach," PhD Thesis, Illinois       !
!       Institute of Technology, Chicago, Illinois, 2004               !
!    Iddir, Y.H., & H. Arastoopour, "Modeling of Multitype particle    !
!      flow using the kinetic theory approach," AIChE J., Vol 51,      !
!      no. 6, June 2005                                                !
!     Gera, D., Syamlal, M., O'Brien T.J. 2004. International Journal  !
!        of Multiphase Flow, 30, p419-428.                             !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DRAG_SS_IA(L, M)

! Modules
!---------------------------------------------------------------------//
      use compar, only: ijkstart3, ijkend3
      USE constant, only: segregation_slope_coefficient
      USE constant, only: c_e, pi
      USE drag, only: f_ss
      USE fldvar, only: d_p, ro_s, rop_s, theta_m, p_star
      use functions, only: wall_at, funlm
      USE param1, only: zero
      USE physprop, only: close_packed
      USE rdf, only: g_0
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! solids phase indices
      INTEGER, INTENT(IN) :: M, L

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: IJK
! index for storing solids-solids drag coefficients in the upper
! triangle of the matrix
      INTEGER :: LM
! particle diameters of phase M and phase L
      DOUBLE PRECISION :: D_pm, D_pl
! particle densities of phase M and phase L
      DOUBLE PRECISION :: RO_M, RO_L
! radial distribution function between phases M and L
      DOUBLE PRECISION :: G0_ML
! granular temperature of phase M and phase L
      DOUBLE PRECISION :: thetaM, thetaL
! solid-solid drag coefficient
      DOUBLE PRECISION :: lDss
! Sum of particle diameters
      DOUBLE PRECISION :: DPSUM
! Intermediate calculation
      DOUBLE PRECISION :: M_PM, M_PL, MPSUM, DPSUMo2
      DOUBLE PRECISION :: Ap_lm, Dp_lm, Bp_lm, R2p_lm
      DOUBLE PRECISION :: F_common_term
!---------------------------------------------------------------------//
      LM = FUNLM(L,M)

      DO IJK = ijkstart3, ijkend3
         IF (WALL_AT(IJK)) CYCLE
! Evaluate at all flow boundaries and fluid cellls. This is unlike
! calculation of the fluid-solid drag coefficient, which is only
! evaluated in fluid cells and pressure inflow cells

! setting aliases for easy reference
         D_PM = D_P(IJK,M)
         D_PL = D_P(IJK,L)
         RO_M = RO_S(IJK,M)
         RO_L = RO_S(IJK,L)
         ThetaM = Theta_M(IJK,M)
         ThetaL = Theta_M(IJK,L)
         G0_ML = G_0(IJK,L,M)

         DPSUM = D_PL + D_PM
         M_PM = (Pi/6.d0) * D_PM**3 *RO_M
         M_PL = (Pi/6.d0) * D_PL**3 *RO_L
         MPSUM = M_PM + M_PL
         DPSUMo2 = DPSUM/2.d0
         ldss = zero

! evaluating the solids-solids drag coefficient
         IF(ThetaM > ZERO .AND. ThetaL > ZERO) THEN
            Ap_lm = (M_PM*ThetaL+M_PL*ThetaM)/2.d0
            Bp_lm = (M_PM*M_PL*(ThetaL-ThetaM ))/(2.d0*MPSUM)
            Dp_lm = (M_PL*M_PM*(M_PM*ThetaM+M_PL*ThetaL ))/&
                (2.d0*MPSUM*MPSUM)

            R2p_lm = ( 1.d0/( 2.d0*Ap_lm**1.5 * Dp_lm*Dp_lm ) )+&
                     ( (3.d0*Bp_lm*Bp_lm)/( Ap_lm**2.5 * Dp_lm**3 ) )+&
                     ( (15.d0*Bp_lm**4)/( 2.d0*Ap_lm**3.5 * Dp_lm**4 ) )

            F_common_term = (DPSUMo2*DPSUMo2/4.d0)*(M_PM*M_PL/MPSUM)*&
                 G0_ML*(1.d0+C_E)*(M_PM*M_PL)**1.5

! Momentum source associated with relative velocity between solids
! phase m and solid phase l
! Factor of rop_sm and rop_sl included in calling routine
            ldss = F_common_term/(M_PM*M_PL)*DSQRT(PI)*R2p_lm*&
                        (ThetaM*ThetaL)**2
         ENDIF
         F_SS(IJK,LM) = lDss*ROP_S(IJK,M)*ROP_S(IJK,L)

! Gera: accounting for particle-particle drag due to enduring contact
! in a close-packed system.
         IF(CLOSE_PACKED(M) .AND. CLOSE_PACKED(L)) &
            F_SS(IJK,LM) = F_SS(IJK,LM) + &
               SEGREGATION_SLOPE_COEFFICIENT*P_star(IJK)

      ENDDO    ! end do (ijk=ijkstart3,ijkend3)

      RETURN
      END SUBROUTINE DRAG_SS_IA
