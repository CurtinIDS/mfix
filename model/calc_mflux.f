!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CALC_MFLUX                                              C
!  Purpose: Calculate the convection fluxes. Master routine.           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MFLUX()

! Modules
!---------------------------------------------------------------------//
      USE fldvar, only: u_g, v_g, w_g
      USE fldvar, only: u_s, v_s, w_s
      USE mflux, only: rop_ge, rop_gn, rop_gt
      USE mflux, only: rop_se, rop_sn, rop_st
      USE mflux, only: flux_ge, flux_gn, flux_gt
      USE mflux, only: flux_se, flux_sn, flux_st
      USE mflux, only: flux_gse, flux_gsn, flux_gst
      USE mflux, only: flux_sse, flux_ssn, flux_sst
      USE physprop, only: mmax
      USE run, only: added_mass, m_am
      IMPLICIT NONE

! Local variables
!---------------------------------------------------------------------//
! solids phase index
      INTEGER ::  M
!---------------------------------------------------------------------//

      IF(.NOT.Added_Mass) THEN
         CALL CALC_MFLUX0 (U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT, &
                           Flux_gE, Flux_gN, Flux_gT)
         DO M = 1, MMAX
           CALL CALC_MFLUX0 (U_s(1,M), V_s(1,M), W_s(1,M), &
                             ROP_sE(1,M), ROP_sN(1,M), ROP_sT(1,M), &
                             Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M))
         ENDDO

      ELSE
! New fluxes are defined for gas based on added mass
         CALL CALC_MFLUX_AM (U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT, &
                           ROP_sE(1,M_AM), ROP_sN(1,M_AM), ROP_sT(1,M_AM), &
                           Flux_gE, Flux_gN, Flux_gT, M_AM)

         CALL CALC_MFLUX0 (U_g, V_g, W_g, ROP_gE, ROP_gN, ROP_gT, &
                           Flux_gSE, Flux_gSN, Flux_gST)
         DO M = 1, MMAX
! New fluxes are defined for M = M_am where virtual mass force is added.
           IF(M==M_AM) THEN
              CALL CALC_MFLUX_AM (U_s(1,M), V_s(1,M), W_s(1,M), &
                             ROP_sE(1,M), ROP_sN(1,M), ROP_sT(1,M), &
                             ROP_gE, ROP_gN, ROP_gT, &
                             Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M),&
                             M_AM)
              CALL CALC_MFLUX0 (U_s(1,M), V_s(1,M), W_s(1,M), &
                             ROP_sE(1,M), ROP_sN(1,M), ROP_sT(1,M), &
                             Flux_sSE, Flux_sSN, Flux_sST)
           ELSE
              CALL CALC_MFLUX0 (U_s(1,M), V_s(1,M), W_s(1,M), &
                             ROP_sE(1,M), ROP_sN(1,M), ROP_sT(1,M), &
                             Flux_sE(1,M), Flux_sN(1,M), Flux_sT(1,M))
           ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CALC_MFLUX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine:: CALC_MFLUX0                                            C
!  Purpose: Calculate the convection fluxes.                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MFLUX0(U, V, W, ROP_E, ROP_N, ROP_T,&
                             Flux_E, Flux_N, Flux_T)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      USE functions, only: im_of, jm_of, km_of
      USE geometry, only: do_k
      USE geometry, only: ayz, axz, axy
      USE param, only: dimension_3
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)
! Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION, INTENT(IN) :: ROP_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: ROP_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: ROP_T(DIMENSION_3)
! Convective mass fluxes
      DOUBLE PRECISION, INTENT(OUT) :: Flux_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: Flux_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: Flux_T(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM
!---------------------------------------------------------------------//

!$omp  parallel do default(none) private( IJK, IMJK, IJMK, IJKM) &
!$omp              shared(ijkstart3, ijkend3, flux_e, flux_n, flux_t,&
!$omp                     rop_e, rop_n, rop_t, axz, axy, ayz, u, v, w, do_k)

      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK)) THEN

            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1/2, j, k)
            Flux_E(IJK) = ROP_E(IJK)*AYZ(IJK)*U(IJK)
! West face (i-1/2, j, k)
            IF (.NOT.FLUID_AT(IMJK)) THEN
               Flux_E(IMJK) = ROP_E(IMJK)*AYZ(IMJK)*U(IMJK)
            ENDIF

! North face (i, j+1/2, k)
            Flux_N(IJK) = ROP_N(IJK)*AXZ(IJK)*V(IJK)
! South face (i, j-1/2, k)
            IF (.NOT.FLUID_AT(IJMK)) THEN
              Flux_N(IJMK) = ROP_N(IJMK)*AXZ(IJMK)*V(IJMK)
            ENDIF

            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
! Top face (i, j, k+1/2)
               Flux_T(IJK) = ROP_T(IJK)*AXY(IJK)*W(IJK)
! Bottom face (i, j, k-1/2)
               IF (.NOT.FLUID_AT(IJKM)) THEN
                 Flux_T(IJKM) = ROP_T(IJKM)*AXY(IJKM)*W(IJKM)
               ENDIF
            ENDIF   ! end if do_k
         ENDIF   ! end if fluid_at
      ENDDO   ! end do ijk

      RETURN
      END SUBROUTINE CALC_MFLUX0

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!                                                                      C
!  Subroutine: CALC_MFLUX_AM                                           C
!  Purpose: Calculate the convection fluxes for case with added_mass.  C
!                                                                      C
!  Author: S. Benyahia                                 Date:           C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_MFLUX_AM(U, V, W, ROP_E, ROP_N, ROP_T, &
                             ROPa_E, ROPa_N, ROPa_T, Flux_E, &
                             Flux_N, Flux_T, M_AM)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      USE functions, only: im_of, jm_of, km_of
      USE fldvar, only: ro_s
      USE geometry, only: do_k
      USE geometry, only: ayz, axz, axy
      USE param, only: dimension_3
      USE param1, only: one
      USE physprop, only: cv
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)
! Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION, INTENT(IN) :: ROP_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: ROP_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: ROP_T(DIMENSION_3)
! Face value of density of added mass phase
      DOUBLE PRECISION, INTENT(IN) :: ROPa_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: ROPa_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: ROPa_T(DIMENSION_3)
! Convective mass fluxes
      DOUBLE PRECISION, INTENT(OUT) :: Flux_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: Flux_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: Flux_T(DIMENSION_3)
! phase index for phase of added_mass
      INTEGER, INTENT(IN) :: M_AM

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM
!---------------------------------------------------------------------//

!!!$omp  parallel do private( IJK, IMJK, IJMK, IJKM) &
!!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK)) THEN

            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1/2, j, k)
            Flux_E(IJK) = ROP_E(IJK)* &
                          (ONE + Cv*ROPa_E(IJK)/RO_S(IJK,M_AM))*&
                          AYZ(IJK)*U(IJK)
! West face (i-1/2, j, k)
            IF (.NOT.FLUID_AT(IMJK)) THEN
               Flux_E(IMJK) = ROP_E(IMJK)*&
                              (ONE + Cv*ROPa_E(IMJK)/RO_S(IJK,M_AM))*&
                              AYZ(IMJK)*U(IMJK)
            ENDIF

! North face (i, j+1/2, k)
            Flux_N(IJK) = ROP_N(IJK)*&
                          (ONE + Cv*ROPa_N(IJK)/RO_S(IJK,M_AM))*&
                          AXZ(IJK)*V(IJK)
! South face (i, j-1/2, k)
            IF (.NOT.FLUID_AT(IJMK)) THEN
              Flux_N(IJMK) = ROP_N(IJMK)*&
                             (ONE + Cv*ROPa_N(IJMK)/RO_S(IJK,M_AM))*&
                             AXZ(IJMK)*V(IJMK)
            ENDIF

            IF (DO_K) THEN
               IJKM = KM_OF(IJK)

! Top face (i, j, k+1/2)
               Flux_T(IJK) = ROP_T(IJK)*&
                             (ONE + Cv*ROPa_T(IJK)/RO_S(IJK,M_AM))*&
                             AXY(IJK)*W(IJK)
! Bottom face (i, j, k-1/2)
               IF (.NOT.FLUID_AT(IJKM)) THEN
                 Flux_T(IJKM) = ROP_T(IJKM)*&
                                (ONE + Cv*ROPa_T(IJKM)/RO_S(IJK,M_AM))*&
                                AXY(IJKM)*W(IJKM)
               ENDIF
            ENDIF   ! end if do_k
         ENDIF   ! end if fluid_at
      ENDDO    ! end do ijk

      RETURN
      END SUBROUTINE CALC_MFLUX_AM
