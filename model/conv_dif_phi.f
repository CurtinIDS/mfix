!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_PHI                                            C
!  Purpose: Determine convection diffusion terms for a scalar eqn.     C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative;          C
!  See source_phi                                                      C
!                                                                      C
!  Diffusion at the flow boundaries is prevented by setting the        C
!  diffusion coefficients at boundary cells to zero and then using a   C
!  harmonic average to calculate the boundary diffusivity.  The value  C
!  diffusivities at the boundaries are checked in check_data_30.       C
!  Ensure that harmonic avergaing is used in this routine.             C
!  See source_phi                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_DIF_PHI(PHI, DIF, DISC, UF, VF, WF, &
                              Flux_E, Flux_N, Flux_T, M, A_M, B_M)


! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3, dimension_m
      USE run, only: def_cor
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Scalar
      DOUBLE PRECISION Phi(DIMENSION_3)
! Gamma -- diffusion coefficient
      DOUBLE PRECISION Dif(DIMENSION_3)
! Discretization index
      INTEGER, INTENT(IN) :: Disc
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: Uf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Vf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Wf(DIMENSION_3)
! Mass flux components
      DOUBLE PRECISION, INTENT(IN) :: Flux_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_T(DIMENSION_3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!---------------------------------------------------------------------//

! DEFERRED CORRECTION IS USED WITH THE SCALAR TRANSPORT EQN.
      IF(DEF_COR)THEN
         CALL CONV_DIF_PHI0(PHI, DIF, UF, VF, WF, &
            Flux_E, Flux_N, Flux_T, M, A_M)
         IF (DISC > 1) CALL CONV_DIF_PHI_DC(PHI, DISC, UF, VF, WF,&
                         Flux_E, Flux_N, Flux_T, M, B_M)
      ELSE

! DO NOT USE DEFERRED CORRECTION TO SOLVE THE SCALAR TRANSPORT EQN.
         IF (DISC == 0) THEN
            CALL CONV_DIF_PHI0(PHI, DIF, UF, VF, WF, &
               Flux_E, Flux_N, Flux_T, M, A_M)
         ELSE
            CALL CONV_DIF_PHI1(PHI, DIF, DISC, UF, VF, WF, &
               Flux_E, Flux_N, Flux_T, M, A_M)
         ENDIF
      ENDIF

      CALL DIF_PHI_IS (DIF, A_M, M)

      RETURN
      END SUBROUTINE CONV_DIF_PHI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  faces of a scalar cell. Note the fluxes are calculated at           C
!  all faces regardless of fuid_at condition of the west, south        C
!  or bottom cell.                                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_PHICELL_DIFF_TERMS(Dif, D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, IJK)

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_treatment_at, cut_cell_at

      USE functions, only: fluid_at
      USE functions, only: east_of, north_of, top_of
      USE functions, only: west_of, south_of, bottom_of
      USE functions, only: im_of, jm_of, km_of
      USE functions, only: ip_of, jp_of, kp_of

      USE fun_avg, only: avg_x_h, avg_y_h, avg_z_h

      USE geometry, only: odx_e, ody_n, odz_t
      USE geometry, only: do_k
      USE geometry, only: ox
      USE geometry, only: dx, dy, dz
      USE geometry, only: ayz, axz, axy

      USE indices, only: i_of, j_of, k_of
      USE indices, only: im1, jm1, km1

      USE param, only: dimension_3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Gamma -- diffusion coefficient
      DOUBLE PRECISION, INTENT(IN) :: Dif(DIMENSION_3)

! fluxes through faces of given ijk u-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: d_fe, d_fw
      DOUBLE PRECISION, INTENT(OUT) :: d_fn, d_fs
      DOUBLE PRECISION, INTENT(OUT) :: d_ft, d_fb
! ijk index
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: ipjk, ijpk, ijkp
      INTEGER :: i, j, k, im, jm, km
      INTEGER :: ijke, ijkw, ijkn, ijks, ijkt, ijkb

! area terms
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
!---------------------------------------------------------------------//

      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)


      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IM = IM1(I)
      JM = JM1(J)
      KM = KM1(K)

      IJKE = EAST_OF(IJK)
      IJKN = NORTH_OF(IJK)
      IJKS = SOUTH_OF(IJK)
      IJKW = WEST_OF(IJK)

      C_AE = ODX_E(I)*AYZ(IJK)
      C_AW = ODX_E(IM)*AYZ(IMJK)
      C_AN = ODY_N(J)*AXZ(IJK)
      C_AS = ODY_N(JM)*AXZ(IJMK)
      C_AT = OX(I)*ODZ_T(K)*AXY(IJK)
      C_AB = OX(I)*ODZ_T(KM)*AXY(IJKM)

      IF(CUT_TREATMENT_AT(IJK).AND.CUT_CELL_AT(IJK)) THEN
         IPJK = IP_OF(IJK)
         IJPK = JP_OF(IJK)
         IJKP = KP_OF(IJK)

         IF (.NOT.FLUID_AT(IPJK)) C_AE = ODX_E(I)*DY(J)*DZ(K)
         IF (.NOT.FLUID_AT(IMJK)) C_AW = ODX_E(IM)*DY(J)*DZ(K)
         IF (.NOT.FLUID_AT(IJPK)) C_AN = ODY_N(J)*DX(I)*DZ(K)
         IF (.NOT.FLUID_AT(IJMK)) C_AS = ODY_N(JM)*DX(I)*DZ(K)
         IF (.NOT.FLUID_AT(IJKP)) C_AT = OX(I)*ODZ_T(K)*DX(I)*DY(J)
         IF (.NOT.FLUID_AT(IJKM)) C_AB = OX(I)*ODZ_T(KM)*DX(I)*DY(J)
      ENDIF

! East face (i+1/2, j, k)
      D_Fe = AVG_X_H(DIF(IJK),DIF(IJKE),I)*C_AE

! West face (i-1/1, j, k)
      D_FW = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*C_AW

! North face (i, j+1/2, k)
      D_FN = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*C_AN

! South face (i, j-1/2, k)
      D_FS = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*C_AS

      IF (DO_K) THEN
         IJKT = TOP_OF(IJK)
         IJKB = BOTTOM_OF(IJK)

! Top face (i, j, k+1/2)
         D_FT = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*C_AT
! Bottom face (i, j, k-1/2)
         D_FB = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*C_AB
      ENDIF

      RETURN
      END SUBROUTINE GET_PHICELL_DIFF_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_PHI0                                           C
!  Purpose: Determine convection diffusion terms for Phi balance       C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative;          C
!  See source_phi                                                      C
!  Implement FOUP discretization                                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_DIF_PHI0(PHI, DIF, UF, VF, WF, &
                               Flux_E, Flux_N, Flux_T, M, A_M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE functions, only: fluid_at
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of
      USE geometry, only: do_k
      USE param
      USE param1, only: zero
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Scalar
      DOUBLE PRECISION, INTENT(IN) :: Phi(DIMENSION_3)
! Gamma -- diffusion coefficient
      DOUBLE PRECISION, INTENT(IN) :: Dif(DIMENSION_3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: Uf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Vf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Wf(DIMENSION_3)
! Mass flux components
      DOUBLE PRECISION, INTENT(IN) :: Flux_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_T(DIMENSION_3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
      INTEGER :: IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM
! Diffusion parameter
      DOUBLE PRECISION :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb
!---------------------------------------------------------------------//


!!!$omp      parallel do                                              &
!!!$omp&     private(IJK, IPJK, IJPK, IJKM, IMJK, IJMK, IJKM,         &
!!!$omp&             d_fe, df_w, df_n, df_s, df_t, df_b)
      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_PHICELL_DIFF_TERMS(dif, d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk)

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1/2, j, k)
            IF (UF(IJK) >= ZERO) THEN
               A_M(IJK,east,M) = D_Fe
               A_M(IPJK,west,M) = D_Fe + FLUX_E(IJK)
            ELSE
               A_M(IJK,east,M) = D_Fe - FLUX_E(IJK)
               A_M(IPJK,west,M) = D_Fe
            ENDIF
! West face (i-1/2, j, k)
            IF (.NOT.FLUID_AT(IMJK)) THEN
               IF (UF(IMJK) >= ZERO) THEN
                  A_M(IJK,west,M) = D_Fw + FLUX_E(IMJK)
               ELSE
                  A_M(IJK,west,M) = D_Fw
               ENDIF
            ENDIF


! North face (i, j+1/2, k)
            IF (VF(IJK) >= ZERO) THEN
               A_M(IJK,north,M) = D_Fn
               A_M(IJPK,south,M) = D_Fn + FLUX_N(IJK)
            ELSE
               A_M(IJK,north,M) = D_Fn - FLUX_N(IJK)
               A_M(IJPK,south,M) = D_Fn
            ENDIF
! South face (i, j-1/2, k)
            IF (.NOT.FLUID_AT(IJMK)) THEN
               IF (VF(IJMK) >= ZERO) THEN
                  A_M(IJK,south,M) = D_Fs + FLUX_N(IJMK)
               ELSE
                  A_M(IJK,south,M) = D_Fs
               ENDIF
            ENDIF


            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)

! Top face (i, j, k+1/2)
               IF (WF(IJK) >= ZERO) THEN
                  A_M(IJK,top,M) = D_FT
                  A_M(IJKP,bottom,M) = D_Ft + FLUX_T(IJK)
               ELSE
                  A_M(IJK,top,M) = D_Ft - FLUX_T(IJK)
                  A_M(IJKP,bottom,M) = D_Ft
               ENDIF

! Bottom face (i, j, k-1/2)

               IF (.NOT.FLUID_AT(IJKM)) THEN
                  IF (WF(IJKM) >= ZERO) THEN
                     A_M(IJK,bottom,M) = D_Fb + FLUX_T(IJKM)
                  ELSE
                     A_M(IJK,bottom,M) = D_Fb
                  ENDIF
               ENDIF
            ENDIF

         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_DIF_PHI0



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_PHI_DC                                         C
!  Purpose: Use deferred correction method to solve the scalar         C
!  transport equation. This method combines first order upwind and     C
!  a user specified higher order method                                C
!                                                                      C
!  Author: C. GUENTHER                                Date: 1-ARP-99   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_DIF_PHI_DC(PHI, DISC, UF, VF, WF, &
         Flux_E, Flux_N, Flux_T, M, B_M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE discretization, only: fpfoi_of

      USE function3, only: funijk3
      USE functions, only: fluid_at
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: do_k

      USE indices, only: i_of, j_of, k_of

      USE param, only: dimension_3, dimension_m, dimension_4
      USE param1, only: zero, one

      USE run, only: fpfoi
      USE sendrecv3, only: send_recv3

      USE xsi, only: calc_xsi
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Scalar
      DOUBLE PRECISION, INTENT(IN) :: Phi(DIMENSION_3)
! Discretization index
      INTEGER, INTENT(IN) :: Disc
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: Uf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Vf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Wf(DIMENSION_3)
! Mass flux components
      DOUBLE PRECISION, INTENT(IN) :: Flux_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_T(DIMENSION_3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)

! Dummy arguments
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
      INTEGER :: IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER :: IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! indication for shear
      INTEGER :: incr
! Deferred correction contribution from high order method
      DOUBLE PRECISION :: PHI_HO
! low order approximation
      DOUBLE PRECISION :: PHI_LO
! deferred correction contribution from each face
      DOUBLE PRECISION :: EAST_DC
      DOUBLE PRECISION :: WEST_DC
      DOUBLE PRECISION :: NORTH_DC
      DOUBLE PRECISION :: SOUTH_DC
      DOUBLE PRECISION :: TOP_DC
      DOUBLE PRECISION :: BOTTOM_DC
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TMP4

      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: XSI_e, XSI_n, XSI_t

!---------------------------------------------------------------------//

      allocate(tmp4(DIMENSION_4))

! Send recv the third ghost layer
      IF (FPFOI) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = PHI(IJK)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

! shear indicator:
      incr=0
      CALL CALC_XSI (DISC, PHI, UF, VF, WF, XSI_E, XSI_N, XSI_T, incr)

!!!$omp      parallel do                                             &
!!!$omp&     private(I, J, K, IJK,  IPJK, IJPK, IJKP,                &
!!!$omp&             IMJK, IJMK, IJKM  V_f,                          &
!!!$omp&             PHI_HO, PHI_LO, EAST_DC, WEST_DC, NORTH_DC,     &
!!!$omp&             SOUTH_DC, TOP_DC, BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK)) THEN

            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKP = KP_OF(IJK)
            IJKM = KM_OF(IJK)

! Third Ghost layer information
            IPPP  = IP_OF(IP_OF(IPJK))
            IPPP4 = funijk3(I_OF(IPPP), J_OF(IPPP), K_OF(IPPP))
            IMMM  = IM_OF(IM_OF(IMJK))
            IMMM4 = funijk3(I_OF(IMMM), J_OF(IMMM), K_OF(IMMM))
            JPPP  = JP_OF(JP_OF(IJPK))
            JPPP4 = funijk3(I_OF(JPPP), J_OF(JPPP), K_OF(JPPP))
            JMMM  = JM_OF(JM_OF(IJMK))
            JMMM4 = funijk3(I_OF(JMMM), J_OF(JMMM), K_OF(JMMM))
            KPPP  = KP_OF(KP_OF(IJKP))
            KPPP4 = funijk3(I_OF(KPPP), J_OF(KPPP), K_OF(KPPP))
            KMMM  = KM_OF(KM_OF(IJKM))
            KMMM4 = funijk3(I_OF(KMMM), J_OF(KMMM), K_OF(KMMM))


! East face (i+1/2, j, k)
            IF(UF(IJK)>= ZERO)THEN
               PHI_LO = PHI(IJK)
               IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IPJK), PHI(IJK), &
                                   PHI(IMJK), PHI(IM_OF(IMJK)))
            ELSE
               PHI_LO = PHI(IPJK)
               IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJK), PHI(IPJK), &
                                   PHI(IP_OF(IPJK)), TMP4(IPPP4))
            ENDIF
            IF (.NOT. FPFOI) PHI_HO = XSI_E(IJK)*PHI(IPJK)+&
                                      (1.0-XSI_E(IJK))*PHI(IJK)
            EAST_DC = FLUX_E(IJK)*(PHI_LO - PHI_HO)


! North face (i, j+1/2, k)
            IF(VF(IJK) >= ZERO)THEN
               PHI_LO = PHI(IJK)
               IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJPK), PHI(IJK), &
                                   PHI(IJMK), PHI(JM_OF(IJMK)))
            ELSE
               PHI_LO = PHI(IJPK)
               IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJPK), &
                                   PHI(JP_OF(IJPK)), TMP4(JPPP4))
            ENDIF
            IF (.NOT. FPFOI) PHI_HO = XSI_N(IJK)*PHI(IJPK)+&
                                      (1.0-XSI_N(IJK))*PHI(IJK)
            NORTH_DC = FLUX_N(IJK)*(PHI_LO - PHI_HO)


! Top face (i, j, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IF(WF(IJK) >= ZERO)THEN
                  PHI_LO = PHI(IJK)
                  IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJKP), PHI(IJK), &
                                      PHI(IJKM), PHI(KM_OF(IJKM)))
               ELSE
                  PHI_LO = PHI(IJKP)
                  IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJKP),  &
                                      PHI(KP_OF(IJKP)), TMP4(KPPP4))
               ENDIF
               IF (.NOT. FPFOI) PHI_HO = XSI_T(IJK)*PHI(IJKP)+&
                                         (1.0-XSI_T(IJK))*PHI(IJK)
               TOP_DC = FLUX_T(IJK)*(PHI_LO - PHI_HO)
            ELSE
               TOP_DC = ZERO
            ENDIF


! West face (i-1/2, j, k)
            IMJK = IM_OF(IJK)
            IF(UF(IMJK) >= ZERO)THEN
               PHI_LO = PHI(IMJK)
               IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJK), PHI(IMJK), &
                                   PHI(IM_OF(IMJK)), TMP4(IMMM4))
            ELSE
               PHI_LO = PHI(IJK)
               IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IMJK), PHI(IJK), &
                                   PHI(IPJK), PHI(IP_OF(IPJK)))
            ENDIF
            IF (.NOT. FPFOI) PHI_HO = XSI_E(IMJK)*PHI(IJK)+&
                                      (ONE-XSI_E(IMJK))*PHI(IMJK)
            WEST_DC = FLUX_E(IMJK)*(PHI_LO - PHI_HO)


! South face (i, j-1/2, k)
            IJMK = JM_OF(IJK)
            IF(VF(IJMK) >= ZERO)THEN
               PHI_LO = PHI(IJMK)
               IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJMK), &
                                   PHI(JM_OF(IJMK)), TMP4(JMMM4))
            ELSE
               PHI_LO = PHI(IJK)
               IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJMK), PHI(IJK), &
                                   PHI(IJPK), PHI(JP_OF(IJPK)))
            ENDIF
            IF (.NOT. FPFOI) PHI_HO = XSI_N(IJMK)*PHI(IJK)+&
                                      (ONE-XSI_N(IJMK))*PHI(IJMK)
            SOUTH_DC = FLUX_N(IJMK)*(PHI_LO - PHI_HO)


! Bottom face (i, j, k-1/2)
            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
               IF(WF(IJKM) >= ZERO)THEN
                  PHI_LO = PHI(IJKM)
                  IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJK), PHI(IJKM), &
                                      PHI(KM_OF(IJKM)), TMP4(KMMM4))
               ELSE
                  PHI_LO = PHI(IJK)
                  IF (FPFOI) PHI_HO = FPFOI_OF(PHI(IJKM), PHI(IJK), &
                                      PHI(IJKP), PHI(KP_OF(IJKP)))
               ENDIF
               IF (.NOT. FPFOI) PHI_HO = XSI_T(IJKM)*PHI(IJK)+&
                                         (1.0-XSI_T(IJKM))*PHI(IJKM)
               BOTTOM_DC = FLUX_T(IJKM)*(PHI_LO - PHI_HO)
            ELSE
               BOTTOM_DC = ZERO
            ENDIF


! CONTRIBUTION DUE TO DEFERRED CORRECTION
            B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-&
                         NORTH_DC+BOTTOM_DC-TOP_DC

         ENDIF   ! end if fluid_at
      ENDDO   ! end do ijk

      DEALLOCATE(tmp4)

      RETURN
      END SUBROUTINE CONV_DIF_PHI_DC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_Phi1                                           C
!  Purpose: Determine convection diffusion terms for scalar transport  C
!  equations. The off-diagonal coefficients calculated here must be    C
!  positive. The center coefficient and the source vector are          C
!  negative;                                                           C
!  See source_Phi                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 21-APR-97  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_DIF_PHI1(PHI, DIF, DISC, UF, VF, WF, &
                               Flux_E, Flux_N, Flux_T, M, A_M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE functions, only: fluid_at
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of
      USE geometry, only: do_k
      USE param
      USE param1, only: one
      USE xsi, only: calc_xsi
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Scalar
      DOUBLE PRECISION, INTENT(IN) :: Phi(DIMENSION_3)
! Gamma -- diffusion coefficient
      DOUBLE PRECISION, INTENT(IN) :: Dif(DIMENSION_3)
! Discretization index
      INTEGER, INTENT(IN) :: Disc
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: Uf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Vf(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Wf(DIMENSION_3)
! Mass flux components
      DOUBLE PRECISION, INTENT(IN) :: Flux_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: Flux_T(DIMENSION_3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
      INTEGER :: IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM
! indicator for shear
      INTEGER :: incr
! Diffusion parameter
      DOUBLE PRECISION :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb

      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: XSI_e, XSI_n, XSI_t

!---------------------------------------------------------------------//

! shear indicator
      incr=0
      CALL CALC_XSI (DISC, PHI, UF, VF, WF, XSI_E, XSI_N, XSI_T, incr)

!!!$omp      parallel do                                                 &
!!!$omp&     private(IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,            &
!!!$omp&             d_fe, d_fw, d_fn, d_fs, d_ft, d_fb)
      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK)) THEN
! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_PHICELL_DIFF_TERMS(dif, d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk)

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1/2, j, k)
            A_M(IJK,east,M) = D_Fe - XSI_E(IJK)*FLUX_E(IJK)
            A_M(IPJK,west,M) = D_Fe + (ONE - XSI_E(IJK))*FLUX_E(IJK)
!  West face (i-1/2, j, k)
            IF (.NOT.FLUID_AT(IMJK)) THEN
              A_M(IJK,west,M) = D_Fw + (ONE - XSI_E(IMJK))*FLUX_E(IMJK)
            ENDIF


! North face (i, j+1/2, k)
            A_M(IJK,north,M) = D_Fn - XSI_N(IJK)*FLUX_N(IJK)
            A_M(IJPK,south,M) = D_Fn + (ONE - XSI_N(IJK))*FLUX_N(IJK)
! South face (i, j-1/2, k)
            IF (.NOT.FLUID_AT(IJMK)) THEN
               A_M(IJK,south,M) = D_Fs + (ONE - XSI_N(IJMK))*FLUX_N(IJMK)
            ENDIF


            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)
! Top face (i, j, k+1/2)
               A_M(IJK,top,M) = D_Ft - XSI_T(IJK)*FLUX_T(IJK)
               A_M(IJKP,bottom,M)=D_Ft + (ONE-XSI_T(IJK))*FLUX_T(IJK)
! Bottom face (i, j, k-1/2)
               IF (.NOT.FLUID_AT(IJKM)) THEN
                  A_M(IJK,bottom,M) = D_Fb + (ONE - XSI_T(IJKM))*FLUX_T(IJKM)
               ENDIF
            ENDIF

         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_DIF_PHI1



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DIF_phi_IS                                              C
!  Purpose: Remove diffusive fluxes across internal surfaces.          C
!  (Make user defined internal surfaces non-conducting)                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DIF_PHI_IS(DIF, A_M, M)

! Modules
!---------------------------------------------------------------------//
      USE param
      USE geometry, only: do_k
      USE geometry, only: ody_n, odx_e, odz_t, ox
      USE geometry, only: axz, axy, ayz

      USE is, only: is_defined, is_plane
      USE is, only: is_i_w, is_i_e, is_j_s, is_j_n, is_k_t, is_k_b

      USE fun_avg, only: avg_x_h, avg_y_h, avg_z_h

      USE functions, only: funijk, ip_of, jp_of, kp_of
      USE functions, only: east_of, north_of, top_of

      USE compar, only: dead_cell_at
      USE compar, only: istart2, jstart2, kstart2
      USE compar, only: iend2, jend2, kend2
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Gamma -- diffusion coefficient
      DOUBLE PRECISION, INTENT(IN) :: Dif(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Solids phase
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! Internal surface
      INTEGER :: L
! Indices
      INTEGER :: I1, I2, J1, J2, K1, K2
      INTEGER :: I, J, K, IJK
      INTEGER :: IJKE, IJKN, IJKT, IPJK, IJPK, IJKP
! Diffusion parameter
      DOUBLE PRECISION :: D_f
!---------------------------------------------------------------------//


      DO L = 1, DIMENSION_IS
         IF (IS_DEFINED(L)) THEN
            I1 = IS_I_W(L)
            I2 = IS_I_E(L)
            J1 = IS_J_S(L)
            J2 = IS_J_N(L)
            K1 = IS_K_B(L)
            K2 = IS_K_T(L)

! Limit I1, I2 and all to local processor first ghost layer
            IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)
            IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)
            IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)
            IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)
            IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)
            IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)

            IF (IS_PLANE(L) == 'E') THEN
               DO K = K1, K2
               DO J = J1, J2
               DO I = I1, I2
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IJKE = EAST_OF(IJK)
                  IPJK = IP_OF(IJK)

                  D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*AYZ(IJK)
                  A_M(IJK,east,M) = A_M(IJK,east,M) - D_F
                  A_M(IPJK,west,M) = A_M(IPJK,west,M) - D_F
               ENDDO
               ENDDO
               ENDDO

            ELSEIF(IS_PLANE(L) == 'N') THEN
               DO K = K1, K2
               DO J = J1, J2
               DO I = I1, I2
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IJKN = NORTH_OF(IJK)
                  IJPK = JP_OF(IJK)

                  D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK)
                  A_M(IJK,north,M) = A_M(IJK,north,M) - D_F
                  A_M(IJPK,south,M) = A_M(IJPK,south,M) - D_F
               ENDDO
               ENDDO
               ENDDO

            ELSEIF(IS_PLANE(L) == 'T') THEN
               IF (DO_K) THEN
                  DO K = K1, K2
                  DO J = J1, J2
                  DO I = I1, I2
                     IJKT = TOP_OF(IJK)
                     IJKP = KP_OF(IJK)

                     D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*&
                        OX(I)*ODZ_T(K)*AXY(IJK)
                     A_M(IJK,top,M) = A_M(IJK,top,M) - D_F
                     A_M(IJKP,bottom,M) = A_M(IJKP,bottom,M) - D_F
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF   ! end if do_k
            ENDIF   ! endif/else is_plane
         ENDIF   ! endif is_defined
      ENDDO   ! end do dimension_is

      RETURN
      END SUBROUTINE DIF_PHI_IS
