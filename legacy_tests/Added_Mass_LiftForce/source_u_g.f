!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_U_g                                              C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative.  The off-diagonal    C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 14-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_U_G(A_M, B_M, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE bc
      USE bodyforce
      USE compar
      USE constant
      USE cutcell
      USE drag
      USE fldvar
      USE geometry
      USE ghdtheory
      USE indices
      USE is
      USE mms
      USE physprop
      USE run, only: added_mass, drag_type_enum, hys, kt_type, m_am, model_b, odt, momentum_x_eq
      USE rxns
      USE scales
      USE tau_g
      USE toleranc
      USE visc_g
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IJKE, IPJK, IJKM, &
                 IPJKM, IMJK, IJMK, IPJMK, IJPK, IJKP
! Phase index
      INTEGER :: M, L, MM
! Internal surface
      INTEGER :: ISV
! Pressure at east cell
      DOUBLE PRECISION :: PgE
! Average volume fraction
      DOUBLE PRECISION :: EPGA
! Average density
      DOUBLE PRECISION :: ROPGA, ROGA
! Average viscosity
      DOUBLE PRECISION :: MUGA
! Average viscosity
      DOUBLE PRECISION :: EPMUGA
! Average W_g
      DOUBLE PRECISION :: Wge
! Average dW/Xdz
      DOUBLE PRECISION :: dWoXdz
! Source terms (Surface)
      DOUBLE PRECISION :: Sdp
! Source terms (Volumetric)
      DOUBLE PRECISION :: V0, Vpm, Vmt, Vbf, Vcf, Vtza
! Source terms (Volumetric) for GHD theory
      DOUBLE PRECISION :: Ghd_drag, avgRop
! Source terms for HYS drag relation
      DOUBLE PRECISION :: HYS_drag, avgDrag
! virtual (added) mass
      DOUBLE PRECISION :: ROP_MA, U_se, Usw, Vsw, Vse, Usn,&
                          Uss, Wsb, Wst, Wse, Usb, Ust
      DOUBLE PRECISION :: F_vir
! Lift force
      DOUBLE PRECISION :: F_Lift, Vgw, Vge, Vgc, U_slip, gradVg
! error message
      CHARACTER*80     LINE
!-----------------------------------------------

! Set reference phase to gas
      M = 0

      IF (.NOT.MOMENTUM_X_EQ(0)) RETURN


!$omp  parallel do default(shared)                                   &
!$omp  private(I, J, K, IJK, IJKE, IJKM, IPJK, IMJK, IPJKM,          &
!$omp          IJMK, IPJMK, IJPK, IJKP, EPGA, PGE, SDP, ROPGA,       &
!$omp           ROGA, ROP_MA, V0, ISV, MUGA, Vpm, Vmt, Vbf,          &
!$omp           U_se, Usw, Vsw, Vse, Usn, Uss, Wsb, Wst, Wse,        &
!$omp           Usb, Ust, F_vir, F_Lift, Vgw, Vge, Vgc, U_slip,      &
!$omp           gradVg, WGE, Vcf, EPMUGA, VTZA,                      &
!$omp           Ghd_drag, L, MM, avgRop, HYS_drag, avgDrag)

      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IJKE = EAST_OF(IJK)
         IJKM = KM_OF(IJK)
         IPJK = IP_OF(IJK)
         IMJK = IM_OF(IJK)
         IPJKM = IP_OF(IJKM)
         IJMK = JM_OF(IJK)
         IPJMK = IP_OF(IJMK)
         IJPK = JP_OF(IJK)
         IJKP = KP_OF(IJK)

         EPGA = AVG_X(EP_G(IJK),EP_G(IJKE),I)

! Impermeable internal surface
         IF (IP_AT_E(IJK)) THEN
            A_M(IJK,east,M) = ZERO
            A_M(IJK,west,M) = ZERO
            A_M(IJK,north,M) = ZERO
            A_M(IJK,south,M) = ZERO
            A_M(IJK,top,M) = ZERO
            A_M(IJK,bottom,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO

! Dilute flow
         ELSEIF (EPGA <= DIL_EP_S) THEN
            A_M(IJK,east,M) = ZERO
            A_M(IJK,west,M) = ZERO
            A_M(IJK,north,M) = ZERO
            A_M(IJK,south,M) = ZERO
            A_M(IJK,top,M) = ZERO
            A_M(IJK,bottom,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO
! set velocity equal to that of west or east cell if solids are present
! in those cells else set velocity equal to known value
            IF (EP_G(WEST_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,west,M) = ONE
            ELSE IF (EP_G(EAST_OF(IJK)) > DIL_EP_S) THEN
               A_M(IJK,east,M) = ONE
            ELSE
               B_M(IJK,M) = -U_G(IJK)
            ENDIF

! Cartesian grid implementation
         ELSEIF (BLOCKED_U_CELL_AT(IJK)) THEN
            A_M(IJK,east,M) = ZERO
            A_M(IJK,west,M) = ZERO
            A_M(IJK,north,M) = ZERO
            A_M(IJK,south,M) = ZERO
            A_M(IJK,top,M) = ZERO
            A_M(IJK,bottom,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = ZERO

! Normal case
         ELSE

! Surface forces
! Pressure term
            PGE = P_G(IJKE)
            IF (CYCLIC_X_PD) THEN
               IF (IMAP(I_OF(IJK)).EQ.IMAX1) PGE = P_G(IJKE) - DELP_X
            ENDIF
            IF (MODEL_B) THEN
               IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*(PGE - P_G(IJK))*AYZ(IJK)
               ELSE
                   SDP = -P_SCALE*(PGE * A_UPG_E(IJK) - P_G(IJK) * A_UPG_W(IJK) )
               ENDIF
            ELSE
               IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
                   SDP = -P_SCALE*EPGA*(PGE - P_G(IJK))*AYZ(IJK)
               ELSE
                   SDP = -P_SCALE*EPGA*(PGE * A_UPG_E(IJK) - P_G(IJK) * A_UPG_W(IJK) )
               ENDIF
            ENDIF

            IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
! Volumetric forces
               ROPGA = HALF * (VOL(IJK)*ROP_G(IJK) + &
                               VOL(IPJK)*ROP_G(IJKE))/VOL_U(IJK)
               ROGA  = HALF * (VOL(IJK)*RO_G(IJK) + &
                               VOL(IPJK)*RO_G(IJKE))/VOL_U(IJK)
! Previous time step
               V0 = HALF * (VOL(IJK)*ROP_GO(IJK) + &
                            VOL(IPJK)*ROP_GO(IJKE))*ODT/VOL_U(IJK)
! Added mass implicit transient term {Cv eps rop_g dU/dt}
               IF(Added_Mass) THEN
                  ROP_MA = AVG_X(ROP_g(IJK)*EP_s(IJK,M_AM),&
                                 ROP_g(IJKE)*EP_s(IJKE,M_AM),I)
                  V0 = V0 + Cv * ROP_MA * ODT
               ENDIF
            ELSE
! Volumetric forces
               ROPGA = (VOL(IJK)*ROP_G(IJK) + &
                        VOL(IPJK)*ROP_G(IJKE))/(VOL(IJK) + VOL(IPJK))
               ROGA  = (VOL(IJK)*RO_G(IJK)  + &
                        VOL(IPJK)*RO_G(IJKE) )/(VOL(IJK) + VOL(IPJK))
! Previous time step
               V0 = (VOL(IJK)*ROP_GO(IJK) + VOL(IPJK)*ROP_GO(IJKE))*&
                  ODT/(VOL(IJK) + VOL(IPJK))
! Added mass implicit transient term {Cv eps rop_g dU/dt}
               IF(Added_Mass) THEN
                  ROP_MA = (VOL(IJK)*ROP_g(IJK)*EP_s(IJK,M_AM) + &
                     VOL(IPJK)*ROP_g(IJKE)*EP_s(IJKE,M_AM) )/&
                     (VOL(IJK) + VOL(IPJK))
                  V0 = V0 + Cv * ROP_MA * ODT
               ENDIF
            ENDIF

! VIRTUAL MASS SECTION (explicit terms)
! adding transient term dvg/dt - dVs/dt to virtual mass term
            F_vir = ZERO
            F_lift = ZERO
            IF(Added_Mass.AND.(.NOT.CUT_U_TREATMENT_AT(IJK))) THEN
               F_vir = ( (U_s(IJK,M_AM) - U_sO(IJK,M_AM)) )*ODT*VOL_U(IJK)

! defining gas-particles velocity at momentum cell faces (or scalar cell center)
               Usw = AVG_X_E(U_S(IMJK,M_AM),U_s(IJK,M_AM),I)
               U_se = AVG_X_E(U_s(IJK,M_AM),U_s(IPJK,M_AM),IP1(I))
               Vsw = AVG_Y_N(V_s(IJMK,M_AM),V_s(IJK,M_AM))
               Vse = AVG_Y_N(V_s(IPJMK,M_AM),V_s(IPJK,M_AM))
               Uss = AVG_Y(U_s(IJMK,M_AM),U_s(IJK,M_AM),JM1(J))
               Usn = AVG_Y(U_s(IJK,M_AM),U_s(IJPK,M_AM),J)
               IF(DO_K) THEN
                  Wsb = AVG_Z_T(W_s(IJKM,M_AM),W_s(IJK,M_AM))
                  Wst = AVG_Z_T(W_s(IPJKM,M_AM),W_s(IPJK,M_AM))
                  Wse = AVG_X(Wsb,Wst,I)
                  Usb = AVG_Z(U_s(IJKM,M_AM),U_s(IJK,M_AM),KM1(K))
                  Ust = AVG_Z(U_s(IJK,M_AM),U_s(IJKP,M_AM),K)
                  F_vir = F_vir + Wse*OX_E(I) * (Ust - Usb) *AXY(IJK)
! centrifugal force
                  IF (CYLINDRICAL) F_vir = F_vir - Wse**2*OX_E(I)
              ENDIF
! adding convective terms (U dU/dx + V dU/dy + W dU/dz) to virtual mass
              F_vir = F_vir + U_s(IJK,M_AM)*(U_se - Usw)*AYZ(IJK) + &
                 AVG_X(Vsw,Vse,I) * (Usn - Uss)*AXZ(IJK)
              F_vir = F_vir * Cv * ROP_MA

! Lift force implemented for a special case
! i.e. this is not a general implementation.
              Vgw = AVG_Y_N(V_g(IJMK),V_g(IJK))
              Vge = AVG_Y_N(V_g(IPJMK),V_g(IPJK))
              Vgc = AVG_X(Vgw,Vge,I)
              U_slip = (Vgc - AVG_X(Vsw,Vse,I))
              gradVg = (Vge - Vgw)*AYZ(IJK)
              F_lift = U_slip * gradVg
              F_lift = F_lift * 0.288d0*ROP_MA ! Lift coefficient Cl = 0.288

            ENDIF

! pressure drop through porous media
            IF (SIP_AT_E(IJK)) THEN
               ISV = IS_ID_AT_E(IJK)
               MUGA = AVG_X(MU_G(IJK),MU_G(IJKE),I)
               VPM = MUGA/IS_PC(ISV,1)
               IF (IS_PC(ISV,2) /= ZERO) VPM = VPM + &
                  HALF*IS_PC(ISV,2)*ROPGA*ABS(U_G(IJK))
            ELSE
               VPM = ZERO
            ENDIF

! Interphase mass transfer
            IF(.NOT.CUT_U_TREATMENT_AT(IJK)) THEN
               VMT = HALF * (VOL(IJK)*SUM_R_G(IJK) + &
                             VOL(IPJK)*SUM_R_G(IJKE))/VOL_U(IJK)
            ELSE
               VMT = (VOL(IJK)*SUM_R_G(IJK) + VOL(IPJK)*SUM_R_G(IJKE))/&
                  (VOL(IJK) + VOL(IPJK))
            ENDIF

! Body force
            IF (MODEL_B) THEN
               VBF = ROGA*BFX_G(IJK)
            ELSE   ! Model A
               VBF = ROPGA*BFX_G(IJK)
            ENDIF

! Additional force for GHD from darg force sum(beta_ig * Joi/rhop_i)
            Ghd_drag = ZERO
            IF (TRIM(KT_TYPE) .EQ. 'GHD') THEN
               DO L = 1,SMAX
                  avgRop = AVG_X(ROP_S(IJK,L),ROP_S(IJKE,L),I)
                  if(avgRop > ZERO) Ghd_drag = Ghd_drag +&
                     AVG_X(F_GS(IJK,L),F_GS(IJKE,L),I) * JoiX(IJK,L) / avgRop
               ENDDO
            ENDIF

! Additional force for HYS drag force, do not use with mixture GHD theory
            avgDrag = ZERO
            HYS_drag = ZERO
            IF (DRAG_TYPE_ENUM .EQ. HYS .AND. TRIM(KT_TYPE) /= 'GHD') THEN
               DO MM=1,MMAX
                  DO L = 1,MMAX
                     IF (L /= MM) THEN
                        avgDrag = AVG_X(beta_ij(IJK,MM,L),beta_ij(IJKE,MM,L),I)
                        HYS_drag = HYS_drag + avgDrag * (U_g(ijk) - U_s(IJK,L))
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF

! Special terms for cylindrical coordinates
            IF (CYLINDRICAL) THEN
! centrifugal force
               WGE = AVG_X(HALF*(W_G(IJK)+W_G(IJKM)),&
                           HALF*(W_G(IPJK)+W_G(IPJKM)),I)
               VCF = ROPGA*WGE**2*OX_E(I)
! virtual mass contribution
               IF(Added_Mass) VCF = VCF + Cv*ROP_MA*WGE**2*OX_E(I)

! -(2mu/x)*(u/x) part of Tau_zz/X
               EPMUGA = AVG_X(MU_GT(IJK),MU_GT(IJKE),I)
               VTZA = 2.d0*EPMUGA*OX_E(I)*OX_E(I)
            ELSE
               VCF = ZERO
               VTZA = ZERO
            ENDIF

! Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,east,M)+A_M(IJK,west,M)+&
               A_M(IJK,north,M)+A_M(IJK,south,M)+A_M(IJK,top,M)+A_M(IJK,bottom,M)+&
               (V0+VPM+ZMAX(VMT)+VTZA)*VOL_U(IJK))
            B_M(IJK,M) = B_M(IJK,M) -(SDP + TAU_U_G(IJK) + &
               ( (V0+ZMAX((-VMT)))*U_GO(IJK) + VBF + &
               VCF + Ghd_drag + HYS_drag)*VOL_U(IJK) )

! adding explicit part of virtual mass force
            B_M(IJK,M) = B_M(IJK,M) - F_vir
! adding explicit part of the lift force
            B_M(IJK,M) = B_M(IJK,M) + F_lift

! MMS source term
            IF(USE_MMS) B_M(IJK,M) = &
               B_M(IJK,M) - MMS_U_G_SRC(IJK)*VOL_U(IJK)
         ENDIF   ! end branching on cell type (ip/dilute/block/else branches)
      ENDDO   ! end do loop over ijk
!$omp end parallel do

! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_U_G(A_M, B_M, IER)
! modifications for bc
      CALL SOURCE_U_G_BC (A_M, B_M, IER)
! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID) CALL CG_SOURCE_U_G_BC(A_M, B_M, IER)

      RETURN
      END SUBROUTINE SOURCE_U_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_U_g_BC                                           C
!  Purpose: Determine source terms for U_g momentum eq. The terms      C
!     appear in the center coefficient and RHS vector. The center      C
!     coefficient and source vector are negative. The off-diagonal     C
!     coefficients are positive.                                       C
!     The drag terms are excluded from the source at this stage.       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 15-MAY-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_U_G_BC(A_M, B_M, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run, only: k_epsilon
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE output
      USE compar
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Boundary condition
      INTEGER :: L
! Indices
      INTEGER ::  I,  J, K, IM, I1, I2, J1, J2, K1, K2, IJK,&
                  JM, KM, IJKW, IMJK, IP, IPJK
! Phase index
      INTEGER :: M
! Turbulent shear stress
      DOUBLE PRECISION  :: W_F_Slip
!-----------------------------------------------

! Set reference phase to gas
      M = 0


! Set the default boundary conditions
! The NS default setting is the where bc_type='dummy' or any default
! (i.e., bc_type=undefined) wall boundary regions are handled. Note that
! the east and west zy planes do not have to be explictly addressed for
! the u-momentum equation. In this direction the velocities are defined
! at the wall (due staggered grid). They are defined as zero for a
! no penetration condition (see zero_norm_vel subroutine and code under
! the ip_at_e branch in the above source routine).
! ---------------------------------------------------------------->>>
      IF (DO_K) THEN
! bottom xy plane
         K1 = 1
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IF (NS_WALL_AT(IJK)) THEN
! Setting the wall velocity to zero (set the boundary cell value equal
! and opposite to the adjacent fluid cell value)
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = -ONE
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ELSEIF (FS_WALL_AT(IJK)) THEN
! Setting the wall velocity equal to the adjacent fluid velocity (set
! the boundary cell value equal to adjacent fluid cell value)
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ONE
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO

! top xy plane
         K1 = KMAX2
         DO J1 = jmin3,jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IF (NS_WALL_AT(IJK)) THEN
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = -ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ELSEIF (FS_WALL_AT(IJK)) THEN
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO
      ENDIF   ! end if (do_k)

! south xz plane
      J1 = 1
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = -ONE
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ONE
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

! north xz plane
      J1 = JMAX2
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (NS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = -ONE
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ELSEIF (FS_WALL_AT(IJK)) THEN
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ONE
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO
! End setting the default boundary conditions
! ----------------------------------------------------------------<<<


! Setting user specified boundary conditions
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

! Setting wall boundary conditions
! ---------------------------------------------------------------->>>
            IF (BC_TYPE(L) == 'NO_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = ZERO
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN
                           A_M(IJK,north,M) = -ONE
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           A_M(IJK,south,M) = -ONE
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           A_M(IJK,top,M) = -ONE
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           A_M(IJK,bottom,M) = -ONE
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ELSEIF (BC_TYPE(L) == 'FREE_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = ZERO
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN
                           A_M(IJK,north,M) = ONE
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           A_M(IJK,south,M) = ONE
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           A_M(IJK,top,M) = ONE
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           A_M(IJK,bottom,M) = ONE
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ELSEIF (BC_TYPE(L) == 'PAR_SLIP_WALL' .AND. .NOT. K_Epsilon) THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
                        JM = JM1(J)
                        KM = KM1(K)
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = ZERO
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,north,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(J))
                              A_M(IJK,north,M) = -(HALF*BC_HW_G(L)-ODY_N(J))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,south,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,south,M) = -(HALF*BC_HW_G(L)-ODY_N(JM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODY_N(JM))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,top,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,0,M)=-(HALF*BC_HW_G(L)+ODZ_T(K)*OX_E(I))
                              A_M(IJK,top,M)=-(HALF*BC_HW_G(L)-ODZ_T(K)*OX_E(I))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           IF (BC_HW_G(L) == UNDEFINED) THEN
                              A_M(IJK,bottom,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_UW_G(L)
                           ELSE
                              A_M(IJK,bottom,M) = -(HALF*BC_HW_G(L)-ODZ_T(KM)*OX_E(I&
                                 ))
                              A_M(IJK,0,M) = -(HALF*BC_HW_G(L)+ODZ_T(KM)*OX_E(I&
                                 ))
                              B_M(IJK,M) = -BC_HW_G(L)*BC_UW_G(L)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

! Setting wall boundary conditions when K_EPSILON
            ELSEIF (BC_TYPE(L) == 'PAR_SLIP_WALL'   .OR.  &
                    BC_TYPE(L) == 'NO_SLIP_WALL'    .OR.  &
                    BC_TYPE(L) == 'FREE_SLIP_WALL'  .AND. &
                    K_Epsilon )THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IF (.NOT.WALL_AT(IJK)) CYCLE  ! skip redefined cells
                        JM = JM1(J)
                        KM = KM1(K)
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = ZERO
                        IF (FLUID_AT(NORTH_OF(IJK))) THEN
                           CALL Wall_Function(IJK,NORTH_OF(IJK),ODY_N(J),W_F_Slip)
                           A_M(IJK,north,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           CALL Wall_Function(IJK,SOUTH_OF(IJK),ODY_N(JM),W_F_Slip)
                           A_M(IJK,south,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           CALL Wall_Function(IJK,TOP_OF(IJK),ODZ_T(K)*OX_E(I),W_F_Slip)
                           A_M(IJK,top,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           CALL Wall_Function(IJK,BOTTOM_OF(IJK),ODZ_T(KM)*OX_E(I),W_F_Slip)
                           A_M(IJK,bottom,M) = W_F_Slip
                           A_M(IJK,0,M) = -ONE
                           IF (BC_TYPE(L) == 'PAR_SLIP_WALL') B_M(IJK,M) = -BC_UW_G(L)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
! end setting of wall boundary conditions
! ----------------------------------------------------------------<<<

! Setting p_inflow or p_outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L)=='P_INFLOW' .OR. BC_TYPE(L)=='P_OUTFLOW') THEN
               IF (BC_PLANE(L) == 'W') THEN
! if the fluid cell is on the west side of the outflow/inflow boundary
! then set the velocity in the boundary cell equal to the velocity of
! the adjacent fluid cell
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           A_M(IJK,east,M) = ZERO
                           A_M(IJK,west,M) = ONE
                           A_M(IJK,north,M) = ZERO
                           A_M(IJK,south,M) = ZERO
                           A_M(IJK,top,M) = ZERO
                           A_M(IJK,bottom,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
! end setting of p_inflow or p_otuflow flow boundary conditions
! ----------------------------------------------------------------<<<

! Setting outflow flow boundary conditions
! ---------------------------------------------------------------->>>
            ELSEIF (BC_TYPE(L) == 'OUTFLOW') THEN
               IF (BC_PLANE(L) == 'W') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           A_M(IJK,east,M) = ZERO
                           A_M(IJK,west,M) = ONE
                           A_M(IJK,north,M) = ZERO
                           A_M(IJK,south,M) = ZERO
                           A_M(IJK,top,M) = ZERO
                           A_M(IJK,bottom,M) = ZERO
                           A_M(IJK,0,M) = -ONE
                           B_M(IJK,M) = ZERO
                           IM = IM1(I)
                           IMJK = IM_OF(IJK)
                           A_M(IMJK,east,M) = ZERO
                           A_M(IMJK,west,M) = X_E(IM)/X_E(IM1(IM))
                           A_M(IMJK,north,M) = ZERO
                           A_M(IMJK,south,M) = ZERO
                           A_M(IMJK,top,M) = ZERO
                           A_M(IMJK,bottom,M) = ZERO
                           A_M(IMJK,0,M) = -ONE
                           B_M(IMJK,M) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ELSEIF (BC_PLANE(L) == 'E') THEN
                  I1 = BC_I_W(L)
                  I2 = BC_I_E(L)
                  J1 = BC_J_S(L)
                  J2 = BC_J_N(L)
                  K1 = BC_K_B(L)
                  K2 = BC_K_T(L)
                  DO K = K1, K2
                     DO J = J1, J2
                        DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                           IJK = FUNIJK(I,J,K)
                           IP = IP1(I)
                           IPJK = IP_OF(IJK)
                           A_M(IPJK,east,M) = X_E(IP)/X_E(I)
                           A_M(IPJK,west,M) = ZERO
                           A_M(IPJK,north,M) = ZERO
                           A_M(IPJK,south,M) = ZERO
                           A_M(IPJK,top,M) = ZERO
                           A_M(IPJK,bottom,M) = ZERO
                           A_M(IPJK,0,M) = -ONE
                           B_M(IPJK,M) = ZERO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
! end setting of outflow flow boundary conditions
! ----------------------------------------------------------------<<<

! Setting bc that are defined but not nsw, fsw, psw, p_inflow,
! p_outflow, or outflow (at this time, this section addresses
! mass_inflow and mass_outflow type boundaries)
! ---------------------------------------------------------------->>>
            ELSE
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
! setting the velocity in the boundary cell equal to what is known
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = -U_G(IJK)
                        IF (BC_PLANE(L) == 'W') THEN
! if the fluid cell is on the west side of the outflow/inflow boundary
! then set the velocity in the adjacent fluid cell equal to what is
! known in that cell
                           IJKW = WEST_OF(IJK)
                           A_M(IJKW,east,M) = ZERO
                           A_M(IJKW,west,M) = ZERO
                           A_M(IJKW,north,M) = ZERO
                           A_M(IJKW,south,M) = ZERO
                           A_M(IJKW,top,M) = ZERO
                           A_M(IJKW,bottom,M) = ZERO
                           A_M(IJKW,0,M) = -ONE
                           B_M(IJKW,M) = -U_G(IJKW)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF   ! end if/else (bc_type)
                    ! ns, fs, psw; else
                    ! p_inflow, p_outflow, or outflow; else
! end setting of 'else' flow boundary conditions
! (mass_inflow/mass_outflow)
! ----------------------------------------------------------------<<<

         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc

      RETURN
      END SUBROUTINE SOURCE_U_G_BC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: Wall_Function                                           C
!  Purpose: Calculate Slip velocity using wall functions               C
!                                                                      C
!  Author: S. Benyahia                                Date: MAY-13-04  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE Wall_Function(IJK1, IJK2, ODX_WF, W_F_Slip)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE physprop
      USE fldvar
      USE visc_g
      USE geometry
      USE indices
      USE bc
      USE compar
      USE turb
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! IJK indices for wall cell and fluid cell
      INTEGER :: IJK1, IJK2
! ODX_WF: 1/dx, and W_F_Slip: value of turb. shear stress at walls
      DOUBLE PRECISION ODX_WF, W_F_Slip
!-----------------------------------------------
! Local parameters
!-----------------------------------------------
! C_mu is constant in turbulent viscosity
      DOUBLE PRECISION, PARAMETER :: C_mu = 0.09D0
! Kappa is Von Karmen constant
      DOUBLE PRECISION, PARAMETER :: Kappa = 0.42D0
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!-----------------------------------------------

      IF(DABS(ODX_WF)>1.0D-5) THEN
! Avoid division by near-zero. This can occur when del_h is undefined
! for some bad cut-cells. Will probably need user-defined tolerance.

         W_F_Slip = (ONE - ONE/ODX_WF* RO_g(IJK2)*C_mu**0.25 * &
            SQRT(K_Turb_G(IJK2))/MU_gT(IJK2) * &
            Kappa/LOG(9.81D+0/(ODX_WF*2.D+0)*RO_g(IJK2)*C_mu**0.25*&
            SQRT(K_Turb_G(IJK2))/MU_g(IJK2)))
      ELSE
! Should it be set to another value in this case?
         W_F_Slip = ONE
      ENDIF


      RETURN
      END SUBROUTINE Wall_Function



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_U_G                                        C
!  Purpose: Adds point sources to the gas phase U-momentum equation.   C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_U_G(A_M, B_M, IER)

      use compar
      use constant
      use functions
      use geometry
      use indices
      use param
      use param1
      use physprop
      use ps

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV, M
      INTEGER :: lIE, lIW
! terms of bm expression
      DOUBLE PRECISION :: pSource
!-----------------------------------------------

! Set reference phase to gas
      M = 0

! Calculate the mass going into each IJK cell. This is done for each
! call in case the point source is time dependent.
      PS_LP: do PSV = 1, DIMENSION_PS
         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(abs(PS_U_g(PSV)) < small_number) cycle PS_LP

         if(PS_U_g(PSV) < 0.0d0) then
            lIW = PS_I_W(PSV) - 1
            lIE = PS_I_E(PSV) - 1
         else
            lIW = PS_I_W(PSV)
            lIE = PS_I_E(PSV)
         endif

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = lIW, lIE

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

            ijk = funijk(i,j,k)
            if(.NOT. fluid_at(ijk)) cycle

            pSource =  PS_MASSFLOW_G(PSV) * (VOL(IJK)/PS_VOLUME(PSV))

            B_M(IJK,M) = B_M(IJK,M) - pSource *                        &
               PS_U_g(PSV) * PS_VEL_MAG_g(PSV)

         enddo
         enddo
         enddo
      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_U_G
