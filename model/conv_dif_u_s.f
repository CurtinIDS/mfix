!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_DIF_U_s                                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be posotive.     C
!  The center coefficient and the source vector are negative.          C
!  See source_u_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-DEC-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_DIF_U_S(A_M, B_M)

! Modules
!---------------------------------------------------------------------//
      USE param, only: dimension_3, dimension_m
      USE run, only: kt_type_enum
      USE run, only: ghd_2007
      USE run, only: def_cor
      USE run, only: momentum_x_eq
      USE run, only: discretize
      USE physprop, only: mmax
      USE visc_s, only: mu_s
      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)

! Local Variables
!---------------------------------------------------------------------//
! Solids phase index
      INTEGER :: M
!---------------------------------------------------------------------//

      DO M = 1, MMAX
        IF(KT_TYPE_ENUM /= GHD_2007 .OR. &
           (KT_TYPE_ENUM == GHD_2007 .AND. M==MMAX)) THEN


          IF (MOMENTUM_X_EQ(M)) THEN

             IF(DEF_COR)THEN
! USE DEFERRED CORRECTION TO SOLVE U_S
                CALL STORE_A_U_S0 (A_M(1,-3,M), M)
                IF (DISCRETIZE(3) > 1) CALL STORE_A_U_SDC (M, B_M)

             ELSE
! DO NOT USE DEFERRED CORRECTION TO SOLVE FOR U_S
                IF (DISCRETIZE(3) == 0) THEN         ! 0 & 1 => FOUP
                   CALL STORE_A_U_S0 (A_M(1,-3,M), M)
                ELSE
                   CALL STORE_A_U_S1 (A_M(1,-3,M), M)
                ENDIF
             ENDIF

            CALL DIF_U_IS (MU_S(1,M), A_M, M)
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_DIF_U_S

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of velocity on the east, north,   C
!  and top face of a u-momentum cell                                   C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_SVTERMS(U, V, WW, M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE cutcell, only: cut_u_treatment_at
      USE cutcell, only: theta_ue, theta_ue_bar
      USE cutcell, only: theta_u_ne, theta_u_nw
      USE cutcell, only: theta_u_te, theta_u_tw
      USE cutcell, only: alpha_ue_c, alpha_un_c, alpha_ut_c

      USE fldvar, only: u_s, v_s, w_s

      USE fun_avg, only: avg_x_e, avg_x
      USE functions, only: ip_of
      USE geometry, only: do_k
      USE indices, only: i_of, ip1

      USE param, only: dimension_3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! phase index
      INTEGER, INTENT(IN) :: M
! velocity components
      DOUBLE PRECISION, INTENT(OUT) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: WW(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: IJK, I, IP, IPJK
! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//

!!!$omp parallel do private(IJK,I,IP,IPJK)
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         IP = IP1(I)
         IPJK = IP_OF(IJK)

         IF(CUT_U_TREATMENT_AT(IJK)) THEN

! East face (i+1, j, k)
            U(IJK) = (Theta_Ue_bar(IJK) * U_s(IJK,M) + &
                      Theta_Ue(IJK) * U_s(IPJK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK, M, 'U_MOMENTUM',&
               alpha_Ue_c(IJK), AW, HW, VELW)
            U(IJK) = U(IJK) * AW

! North face (i+1/2, j+1/2, k)
            V(IJK) = (Theta_U_nw(IJK) * V_s(IJK,M) + &
                      Theta_U_ne(IJK) * V_s(IPJK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK, M, 'U_MOMENTUM',&
               ALPHA_Un_c(IJK), AW, HW, VELW)
            V(IJK) = V(IJK) * AW

! Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN
               WW(IJK) = (Theta_U_tw(IJK) * W_s(IJK,M) + &
                          Theta_U_te(IJK) * W_s(IPJK,M))
               CALL GET_INTERPOLATION_TERMS_S(IJK, M, 'U_MOMENTUM',&
                  ALPHA_Ut_c(IJK), AW, HW, VELW)
               WW(IJK) = WW(IJK) * AW
            ENDIF

         ELSE
            U(IJK) = AVG_X_E(U_S(IJK,M),U_S(IPJK,M),IP)
            V(IJK) = AVG_X(V_S(IJK,M),V_S(IPJK,M),I)
            IF (DO_K) WW(IJK) = AVG_X(W_S(IJK,M),W_S(IPJK,M),I)
         ENDIF
      ENDDO   ! end do ijk

      RETURN
      END SUBROUTINE GET_UCELL_SVTERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the convective fluxes through the faces of a     C
!  u-momentum cell. Note the fluxes are calculated at all faces of     C
!  regardless of flow_at_e of condition of the west, south, or         C
!  bottom face.                                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_SCFLUX_TERMS(FLUX_E, FLUX_W, FLUX_N, &
         FLUX_S, FLUX_T, FLUX_B, IJK, M)

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_u_treatment_at
      USE cutcell, only: theta_ue, theta_ue_bar
      USE cutcell, only: theta_u_ne, theta_u_nw
      USE cutcell, only: theta_u_te, theta_u_tw
      USE cutcell, only: alpha_ue_c, alpha_un_c, alpha_ut_c
      USE functions, only: ip_of, im_of, jm_of, km_of
      USE geometry, only: do_k
      USE mflux, only: flux_se, flux_sn, flux_st

      USE param1, only: half
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! phase index
      INTEGER, INTENT(IN) :: M
! fluxes through faces of given ijk u-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: flux_e, flux_w
      DOUBLE PRECISION, INTENT(OUT) :: flux_n, flux_s
      DOUBLE PRECISION, INTENT(OUT) :: flux_t, flux_b
! ijk index
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: ipjk, ipjmk, ipjkm

! for cartesian grid
      DOUBLE PRECISION :: AW, HW, VELW
!---------------------------------------------------------------------//

! indices
      IPJK = IP_OF(IJK)
      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)
      IPJMK = IP_OF(IJMK)
      IPJKM = IP_OF(IJKM)

! First calculate the fluxes at the faces
      IF(CUT_U_TREATMENT_AT(IJK)) THEN
! east face: i+1, j, k
         Flux_e = (Theta_Ue_bar(IJK) * Flux_sE(IJK,M) + &
                   Theta_Ue(IJK) * Flux_sE(IPJK,M))
         CALL GET_INTERPOLATION_TERMS_S(IJK, M, 'U_MOMENTUM',&
            alpha_Ue_c(IJK), AW, HW, VELW)
         Flux_e = Flux_e * AW
! west face: i-1, j, k
         Flux_w = (Theta_Ue_bar(IMJK) * Flux_sE(IMJK,M) + &
                   Theta_Ue(IMJK) * Flux_sE(IJK,M))
         CALL GET_INTERPOLATION_TERMS_S(IJK, M, 'U_MOMENTUM',&
            alpha_Ue_c(IMJK), AW, HW, VELW)
         Flux_w = Flux_w * AW


! north face: i+1/2, j+1/2, k
         Flux_n = (Theta_U_nw(IJK) * Flux_sN(IJK,M) + &
                   Theta_U_ne(IJK) * Flux_sN(IPJK,M))
         CALL GET_INTERPOLATION_TERMS_S(IJK, M, 'U_MOMENTUM', &
            ALPHA_Un_c(IJK), AW, HW, VELW)
         Flux_n = Flux_n * AW
! south face: i+1/2, j-1/2, k
         Flux_s = (Theta_U_nw(IJMK) * Flux_sN(IJMK,M) + &
                   Theta_U_ne(IJMK) * Flux_sN(IPJMK,M))
         CALL GET_INTERPOLATION_TERMS_S(IJK, M, 'U_MOMENTUM',&
            ALPHA_Un_c(IJMK), AW, HW, VELW)
         Flux_s = Flux_s * AW

         IF (DO_K) THEN
! top face: i+1/2, j, k+1/2
            Flux_t = (Theta_U_tw(IJK) * Flux_sT(IJK,M) + &
                     Theta_U_te(IJK) * Flux_sT(IPJK,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK, M, 'U_MOMENTUM', &
               ALPHA_Ut_c(IJK), AW, HW, VELW)
            Flux_t = Flux_t * AW
! bottom face: i+1/2, j, k-1/2
            Flux_b = (Theta_U_tw(IJKM) * Flux_sT(IJKM,M) + &
                      Theta_U_te(IJKM) * Flux_sT(IPJKM,M))
            CALL GET_INTERPOLATION_TERMS_S(IJK, M, 'U_MOMENTUM',&
               ALPHA_Ut_c(IJKM), AW, HW, VELW)
            Flux_b = Flux_b * AW
         ENDIF
      ELSE
         Flux_e = HALF * (Flux_sE(IJK,M)  + Flux_sE(IPJK,M))
         Flux_w = HALF * (Flux_sE(IMJK,M) + Flux_sE(IJK,M))
         Flux_n = HALF * (Flux_sN(IJK,M)  + Flux_sN(IPJK,M))
         Flux_s = HALF * (Flux_sN(IJMK,M) + Flux_sN(IPJMK,M))

         IF (DO_K) THEN
            Flux_t = HALF * (Flux_sT(IJK,M)  + Flux_sT(IPJK,M))
            Flux_b = HALF * (Flux_sT(IJKM,M) + Flux_sT(IPJKM,M))
         ENDIF
      ENDIF   ! end if/else cut_u_treatment_at
      RETURN
      END SUBROUTINE GET_UCELL_SCFLUX_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Purpose: Calculate the components of diffusive flux through the     C
!  faces of a u-momentum cell. Note the fluxes are calculated at       C
!  all faces regardless of flow_at_e condition of the west, south      C
!  or bottom face.                                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GET_UCELL_SDIFF_TERMS(D_FE, D_FW, D_FN, D_FS, &
         D_FT, D_FB, IJK, M)

! Modules
!---------------------------------------------------------------------//
      USE cutcell, only: cut_u_treatment_at
      USE cutcell, only: oneodx_e_u, oneody_n_u, oneodz_t_u

      USE functions, only: wall_at
      USE functions, only: east_of, north_of, top_of
      USE functions, only: south_of, bottom_of
      USE functions, only: im_of, jm_of, km_of

      USE fun_avg, only: avg_x_h, avg_y_h, avg_z_h

      USE geometry, only: odx, ody_n, odz_t
      USE geometry, only: do_k
      USE geometry, only: ox_e
      USE geometry, only: ayz_u, axz_u, axy_u

      USE indices, only: i_of, j_of, k_of
      USE indices, only: ip1, jm1, km1

      USE visc_s, only: mu_s
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! phase index
      INTEGER, INTENT(IN) :: M
! diffusion through faces of given ijk u-momentum cell
      DOUBLE PRECISION, INTENT(OUT) :: d_fe, d_fw
      DOUBLE PRECISION, INTENT(OUT) :: d_fn, d_fs
      DOUBLE PRECISION, INTENT(OUT) :: d_ft, d_fb
! ijk index
      INTEGER, INTENT(IN) :: ijk

! Local variables
!---------------------------------------------------------------------//
! indices
      INTEGER :: imjk, ijmk, ijkm
      INTEGER :: i, j, k, ip, jm, km
      INTEGER :: ijkc, ijke, ijkn, ijkne, ijks, ijkse
      INTEGER :: ijkt, ijkte, ijkb, ijkbe

! length terms
      DOUBLE PRECISION :: C_AE, C_AW, C_AN, C_AS, C_AT, C_AB
!---------------------------------------------------------------------//

      IMJK = IM_OF(IJK)
      IJMK = JM_OF(IJK)
      IJKM = KM_OF(IJK)

      I = I_OF(IJK)
      J = J_OF(IJK)
      K = K_OF(IJK)
      IP = IP1(I)
      JM = JM1(J)
      KM = KM1(K)

      IJKE = EAST_OF(IJK)
      IF (WALL_AT(IJK)) THEN
         IJKC = IJKE
      ELSE
         IJKC = IJK
      ENDIF
      IJKN = NORTH_OF(IJK)
      IJKNE = EAST_OF(IJKN)
      IJKS = SOUTH_OF(IJK)
      IJKSE = EAST_OF(IJKS)

      IF(CUT_U_TREATMENT_AT(IJK)) THEN
         C_AE = ONEoDX_E_U(IJK)
         C_AW = ONEoDX_E_U(IMJK)
         C_AN = ONEoDY_N_U(IJK)
         C_AS = ONEoDY_N_U(IJMK)
         C_AT = ONEoDZ_T_U(IJK)
         C_AB = ONEoDZ_T_U(IJKM)
      ELSE
         C_AE = ODX(IP)
         C_AW = ODX(I)
         C_AN = ODY_N(J)
         C_AS = ODY_N(JM)
         C_AT = ODZ_T(K)
         C_AB = ODZ_T(KM)
      ENDIF

! East face (i+1, j, k)
      D_FE = MU_S(IJKE,M)*C_AE*AYZ_U(IJK)
! West face (i-1, j, k)
      D_FW = MU_S(IJKC,M)*C_AW*AYZ_U(IMJK)


! North face (i+1/2, j+1/2, k)
      D_FN = AVG_X_H(AVG_Y_H(MU_S(IJKC,M),MU_S(IJKN,M),J),&
                     AVG_Y_H(MU_S(IJKE,M),MU_S(IJKNE,M),J),I)*&
             C_AN*AXZ_U(IJK)
! South face (i+1/2, j-1/2, k)
      D_FS = AVG_X_H(AVG_Y_H(MU_S(IJKS,M),MU_S(IJKC,M),JM),&
                     AVG_Y_H(MU_S(IJKSE,M),MU_S(IJKE,M),JM),I)*&
             C_AS*AXZ_U(IJMK)

      IF (DO_K) THEN
         IJKT = TOP_OF(IJK)
         IJKTE = EAST_OF(IJKT)
         IJKB = BOTTOM_OF(IJK)
         IJKBE = EAST_OF(IJKB)

! Top face (i+1/2, j, k+1/2)
         D_FT = AVG_X_H(AVG_Z_H(MU_S(IJKC,M),MU_S(IJKT,M),K),&
                        AVG_Z_H(MU_S(IJKE,M),MU_S(IJKTE,M),K),I)*&
                OX_E(I)*C_AT*AXY_U(IJK)
! Bottom face (i+1/2, j, k-1/2)
         D_FB = AVG_X_H(AVG_Z_H(MU_S(IJKB,M),MU_S(IJKC,M),KM),&
                        AVG_Z_H(MU_S(IJKBE,M),MU_S(IJKE,M),KM),I)*&
                OX_E(I)*C_AB*AXY_U(IJKM)
      ENDIF
      RETURN
      END SUBROUTINE GET_UCELL_SDIFF_TERMS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_U_s0                                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs. C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative. See      C
!  source_u_s.                                                         C
!  Implement FOUP discretization.                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-APR-96  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_U_S0(A_U_S, M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE functions, only: flow_at_e
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: do_k

      USE param
      USE param1, only: zero
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_U_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_U_s(DIMENSION_3, -3:3, M:M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
      INTEGER :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
! Face mass flux
      DOUBLE PRECISION :: flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b
! Diffusion parameter
      DOUBLE PRECISION :: D_fe, d_fw, d_fn, d_fs, d_ft, d_fb

!---------------------------------------------------------------------//

!$omp     parallel do default(none)                                &
!$omp     private(IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,         &
!$omp             D_fe, d_fw, d_fn, d_fs, d_ft, d_fb,              &
!$omp             flux_e, flux_w, flux_n, flux_s, flux_t, flux_b)  &
!$omp     shared(ijkstart3, ijkend3, do_k, a_u_s, m)

      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_E(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_UCELL_SCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk, M)
            CALL GET_UCELL_SDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk, M)

            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1, j, k)
            IF (Flux_e >= ZERO) THEN
               A_U_S(IJK,east,M) = D_Fe
               A_U_S(IPJK,west,M) = D_Fe + Flux_e
            ELSE
               A_U_S(IJK,east,M) = D_Fe - Flux_e
               A_U_S(IPJK,west,M) = D_Fe
            ENDIF
! West face (i, j, k)
            IF (.NOT.FLOW_AT_E(IMJK)) THEN
               IF (Flux_w >= ZERO) THEN
                  A_U_S(IJK,west,M) = D_Fw + Flux_w
               ELSE
                  A_U_S(IJK,west,M) = D_Fw
               ENDIF
            ENDIF


! North face (i+1/2, j+1/2, k)
            IF (Flux_n >= ZERO) THEN
               A_U_S(IJK,north,M) = D_Fn
               A_U_S(IJPK,south,M) = D_Fn + Flux_n
            ELSE
               A_U_S(IJK,north,M) = D_Fn - Flux_n
               A_U_S(IJPK,south,M) = D_Fn
            ENDIF
! South face (i+1/2, j-1/2, k)
            IF (.NOT.FLOW_AT_E(IJMK)) THEN
               IF (Flux_s >= ZERO) THEN
                  A_U_S(IJK,south,M) = D_Fs + Flux_s
               ELSE
                  A_U_S(IJK,south,M) = D_Fs
               ENDIF
            ENDIF


            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)

! Top face (i+1/2, j, k+1/2)
               IF (Flux_t >= ZERO) THEN
                  A_U_S(IJK,top,M) = D_Ft
                  A_U_S(IJKP,bottom,M) = D_Ft + Flux_t
               ELSE
                  A_U_S(IJK,top,M) = D_Ft - Flux_t
                  A_U_S(IJKP,bottom,M) = D_Ft
               ENDIF
! Bottom face (i+1/2, j, k-1/2)
               IF (.NOT.FLOW_AT_E(IJKM)) THEN
                  IF (Flux_b >= ZERO) THEN
                     A_U_S(IJK,bottom,M) = D_Fb + Flux_b
                  ELSE
                     A_U_S(IJK,bottom,M) = D_Fb
                  ENDIF
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if (flow_at_e)
      ENDDO   !end do (ijk)
!$omp end parallel do


      RETURN
      END SUBROUTINE STORE_A_U_S0



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_U_sdc                                           C
!  Purpose: Use deferred correction method to solve the u-momentum     C
!  equation. This method combines first order upwind and a user        C
!  specified higher order method                                       C
!                                                                      C
!  Author: C. GUENTHER                                Date: 8-APR-99   C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_U_SDC(M, B_M)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3

      USE discretization, only: fpfoi_of

      USE fldvar, only: u_s

      USE function3, only: funijk3
      USE functions, only: flow_at_e
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: do_k

      USE indices, only: i_of, j_of, k_of

      USE param, only: dimension_3, dimension_m, dimension_4
      USE param1, only: zero

      USE run, only: discretize, fpfoi
      USE sendrecv3, only: send_recv3

      USE xsi, only: calc_xsi

      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! phase index
      INTEGER, INTENT(IN) :: M
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IMJK, IPJK, IJMK, IJPK, IJKM, IJKP
      INTEGER :: IJK4, IPPP, IPPP4, JPPP, JPPP4, KPPP, KPPP4
      INTEGER :: IMMM, IMMM4, JMMM, JMMM4, KMMM, KMMM4
! indication for shear
      INTEGER :: incr
! Deferred correction contribution from high order method
      DOUBLE PRECISION :: MOM_HO
! low order approximation
      DOUBLE PRECISION :: MOM_LO
! convection factor at each face
      DOUBLE PRECISION :: flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b
! deferred correction contribution from each face
      DOUBLE PRECISION :: EAST_DC
      DOUBLE PRECISION :: WEST_DC
      DOUBLE PRECISION :: NORTH_DC
      DOUBLE PRECISION :: SOUTH_DC
      DOUBLE PRECISION :: TOP_DC
      DOUBLE PRECISION :: BOTTOM_DC

! temporary use of global arrays:
! array1 (locally u)  - the x directional velocity
      DOUBLE PRECISION :: U(DIMENSION_3)
! array2 (locally v)  - the y directional velocity
      DOUBLE PRECISION :: V(DIMENSION_3)
! array3 (locally ww) - the z directional velocity
      DOUBLE PRECISION :: WW(DIMENSION_3)
!---------------------------------------------------------------------//
      DOUBLE PRECISION :: TMP4(DIMENSION_4)
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: XSI_e, XSI_n, XSI_t

      CALL GET_UCELL_SVTERMS(U, V, WW, M)

! Send recv the third ghost layer
      IF (FPFOI) THEN
         Do IJK = ijkstart3, ijkend3
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IJK4 = funijk3(I,J,K)
            TMP4(IJK4) = U_S(IJK,M)
         ENDDO
         CALL send_recv3(tmp4)
      ENDIF

! shear indicator:
      incr=1
      CALL CALC_XSI (DISCRETIZE(3), U_S(1,M), U, V, WW, XSI_E, XSI_N,&
                     XSI_T, incr)

!!!$omp      parallel do                                              &
!!!$omp&     private(IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,         &
!!!$omp&             flux_e, flux_w, flux_n, flux_s, flux_t, flux_b,  &
!!!$omp&             MOM_HO, MOM_LO, EAST_DC, WEST_DC, NORTH_DC,      &
!!!$omp&             SOUTH_DC, TOP_DC, BOTTOM_DC)
      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_E(IJK)) THEN

! Calculate convection fluxes through each of the faces
            CALL GET_UCELL_SCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk, M)

            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)

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


! East face (i+1, j, k)
            IF(U(IJK) >= ZERO)THEN
               MOM_LO = U_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IPJK,M), U_S(IJK,M), &
                                   U_S(IMJK,M), U_S(IM_OF(IMJK),M))
            ELSE
               MOM_LO = U_S(IPJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IPJK,M), &
                                   U_S(IP_OF(IPJK),M), TMP4(IPPP4))
            ENDIF

            IF (.NOT. FPFOI) MOM_HO = XSI_E(IJK)*U_S(IPJK,M)+ &
                                       (1.0-XSI_E(IJK))*U_S(IJK,M)
            EAST_DC = Flux_e *(MOM_LO - MOM_HO)
! West face (i, j, k)
            IF(U(IMJK) >= ZERO)THEN
               MOM_LO = U_S(IMJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IMJK,M), &
                                   U_S(IM_OF(IMJK),M), TMP4(IMMM4))
            ELSE
               MOM_LO = U_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IMJK,M), U_S(IJK,M), &
                          U_S(IPJK,M), U_S(IP_OF(IPJK),M))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_E(IMJK)*U_S(IJK,M)+ &
                                      (1.0-XSI_E(IMJK))*U_S(IMJK,M)
            WEST_DC = Flux_w * (MOM_LO - MOM_HO)


! North face (i+1/2, j+1/2, k)
            IF(V(IJK) >= ZERO)THEN
               MOM_LO = U_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJPK,M), U_S(IJK,M), &
                                   U_S(IJMK,M), U_S(JM_OF(IJMK),M))
            ELSE
               MOM_LO = U_S(IJPK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJPK,M), &
                                   U_S(JP_OF(IJPK),M), TMP4(JPPP4))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJK)*U_S(IJPK,M)+ &
                                      (1.0-XSI_N(IJK))*U_S(IJK,M)
            NORTH_DC = Flux_n * (MOM_LO - MOM_HO)

! South face (i+1/2, j-1/2, k)
            IF(V(IJMK) >= ZERO)THEN
               MOM_LO = U_S(IJMK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJMK,M), &
                                   U_S(JM_OF(IJMK),M), TMP4(JMMM4))
            ELSE
               MOM_LO = U_S(IJK,M)
               IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJMK,M), U_S(IJK,M), &
                                   U_S(IJPK,M), U_S(JP_OF(IJPK),M))
            ENDIF
            IF (.NOT. FPFOI) MOM_HO = XSI_N(IJMK)*U_S(IJK,M)+ &
                                      (1.0-XSI_N(IJMK))*U_S(IJMK,M)
            SOUTH_DC = Flux_s * (MOM_LO - MOM_HO)


            IF (DO_K) THEN
! Top face (i+1/2, j, k+1/2)
               IF(WW(IJK) >= ZERO)THEN
                  MOM_LO = U_S(IJK,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJKP,M), U_S(IJK,M), &
                                      U_S(IJKM,M), U_S(KM_OF(IJKM),M))
               ELSE
                  MOM_LO = U_S(IJKP,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJKP,M), &
                                      U_S(KP_OF(IJKP),M), TMP4(KPPP4))
               ENDIF
               IF (.NOT. FPFOI) MOM_HO = XSI_T(IJK)*U_S(IJKP,M)+ &
                                         (1.0-XSI_T(IJK))*U_S(IJK,M)
               TOP_DC = Flux_t *(MOM_LO - MOM_HO)

! Bottom face (i+1/2, j, k-1/2)
               IF(WW(IJK) >= ZERO)THEN
                  MOM_LO = U_S(IJKM,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJK,M), U_S(IJKM,M), &
                                      U_S(KM_OF(IJKM),M), TMP4(KMMM4))
               ELSE
                  MOM_LO = U_S(IJK,M)
                  IF (FPFOI) MOM_HO = FPFOI_OF(U_S(IJKM,M), U_S(IJK,M), &
                                      U_S(IJKP,M), U_S(KP_OF(IJKP),M))
               ENDIF
               IF (.NOT. FPFOI) MOM_HO = XSI_T(IJKM)*U_S(IJK,M)+ &
                                         (1.0-XSI_T(IJKM))*U_S(IJKM,M)
               BOTTOM_DC = Flux_b * (MOM_LO - MOM_HO)
            ELSE
               TOP_DC = ZERO
               BOTTOM_DC = ZERO
            ENDIF   ! end if do_k


! CONTRIBUTION DUE TO DEFERRED CORRECTION
            B_M(IJK,M) = B_M(IJK,M)+WEST_DC-EAST_DC+SOUTH_DC-NORTH_DC+&
                         BOTTOM_DC-TOP_DC

         ENDIF   ! end if flow_at_e
      ENDDO   ! end do ijk

      RETURN
      END SUBROUTINE STORE_A_U_SDC

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: STORE_A_U_s1                                            C
!  Purpose: Determine convection diffusion terms for U_s momentum eqs  C
!  The off-diagonal coefficients calculated here must be positive.     C
!  The center coefficient and the source vector are negative.          C
!  Implements higher order discretization.                             C
!  See source_u_s                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 20-MAR-97  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: To incorporate Cartesian grid modifications                C
!  Author: Jeff Dietiker                              Date: 01-Jul-09  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE STORE_A_U_S1(A_U_S, M)
! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE fldvar, only: u_s

      USE functions, only: flow_at_e
      USE functions, only: ip_of, jp_of, kp_of
      USE functions, only: im_of, jm_of, km_of

      USE geometry, only: do_k

      USE param
      USE param1, only: one

      USE run, only: discretize
      USE xsi, only: calc_xsi
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Solids phase
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_U_s
      DOUBLE PRECISION, INTENT(INOUT) :: A_U_s(DIMENSION_3, -3:3, M:M)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IPJK, IMJK, IJPK, IJMK, IJKP, IJKM
! indicator for shear
      INTEGER :: incr
! Diffusion parameter
      DOUBLE PRECISION :: d_fe, d_fw, d_fn, d_fs, d_ft, d_fb
! Face mass flux
      DOUBLE PRECISION :: Flux_e, flux_w, flux_n, flux_s
      DOUBLE PRECISION :: flux_t, flux_b

! temporary use of global arrays:
! array1 (locally u)  - the x directional velocity
      DOUBLE PRECISION :: U(DIMENSION_3)
! array2 (locally v)  - the y directional velocity
      DOUBLE PRECISION :: V(DIMENSION_3)
! array3 (locally ww) - the z directional velocity
      DOUBLE PRECISION :: WW(DIMENSION_3)
!---------------------------------------------------------------------//
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TMP4
      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: XSI_e, XSI_n, XSI_t

      CALL GET_UCELL_SVTERMS(U, V, WW, M)

! shear indicator:
      incr=1
      CALL CALC_XSI (DISCRETIZE(3), U_S(1,M), U, V, WW, XSI_E, XSI_N, &
                     XSI_T, incr)

!!!$omp      parallel do                                                 &
!!!$omp&     private(IJK, IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,            &
!!!$omp&             d_fe, d_fw, d_fn, d_fs, d_ft, d_fb,                 &
!!!$omp&             flux_e, flux_w, flux_n, flux_s, flux_t, flux_b)

      DO IJK = ijkstart3, ijkend3

         IF (FLOW_AT_E(IJK)) THEN

! Calculate convection-diffusion fluxes through each of the faces
            CALL GET_UCELL_SCFLUX_TERMS(flux_e, flux_w, flux_n, &
               flux_s, flux_t, flux_b, ijk, M)
            CALL GET_UCELL_SDIFF_TERMS(d_fe, d_fw, d_fn, d_fs, &
               d_ft, d_fb, ijk, M)

            IPJK = IP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJPK = JP_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1, j, k)
            A_U_S(IJK,east,M) = D_Fe - XSI_E(IJK) * Flux_e
            A_U_S(IPJK,west,M) = D_Fe + (ONE - XSI_E(IJK)) * Flux_e
! West face (i, j, k)
            IF (.NOT.FLOW_AT_E(IMJK)) THEN
               A_U_S(IJK,west,M) = D_Fw + (ONE - XSI_E(IMJK)) * Flux_w
            ENDIF


! North face (i+1/2, j+1/2, k)
            A_U_S(IJK,north,M) = D_Fn - XSI_N(IJK) * Flux_n
            A_U_S(IJPK,south,M) = D_Fn + (ONE - XSI_N(IJK)) * Flux_n
! South face (i+1/2, j-1/2, k)
            IF (.NOT.FLOW_AT_E(IJMK)) THEN
               A_U_S(IJK,south,M) = D_Fs + (ONE - XSI_N(IJMK)) * Flux_s
            ENDIF


! Top face (i+1/2, j, k+1/2)
            IF (DO_K) THEN
               IJKP = KP_OF(IJK)
               IJKM = KM_OF(IJK)
               A_U_S(IJK,top,M) = D_Ft - XSI_T(IJK) * Flux_t
               A_U_S(IJKP,bottom,M) = D_Ft + (ONE - XSI_T(IJK)) * Flux_t
! Bottom face (i+1/2, j, k-1/2)
               IF (.NOT.FLOW_AT_E(IJKM)) THEN
                  A_U_S(IJK,bottom,M) = D_Fb + (ONE - XSI_T(IJKM)) * Flux_b
               ENDIF
            ENDIF   ! end if (do_k)

         ENDIF   ! end if flow_at_e
      ENDDO   ! end do ijk

      RETURN
      END SUBROUTINE STORE_A_U_S1
