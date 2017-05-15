!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. Master routine.                                  C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP()

! Modules
!---------------------------------------------------------------------//
      USE fldvar, only: rop_g, u_g, v_g, w_g
      USE fldvar, only: rop_s, u_s, v_s, w_s
      USE mflux, only: rop_ge, rop_gn, rop_gt
      USE mflux, only: rop_se, rop_sn, rop_st
      USE physprop, only: mmax
      USE run, only: discretize
      IMPLICIT NONE

!---------------------------------------------------------------------//
! Local variables
!---------------------------------------------------------------------//
! solids phase index
      INTEGER :: M
!---------------------------------------------------------------------//

      IF (DISCRETIZE(1) == 0) THEN               ! 0 & 1 => first order upwinding
         CALL CONV_ROP0 (ROP_g, U_g, V_g, W_g, &
                         ROP_gE, ROP_gN, ROP_gT)
      ELSE
         CALL CONV_ROP1 (DISCRETIZE(1), ROP_g, U_g, V_g, W_g, &
                         ROP_gE, ROP_gN, ROP_gT)
      ENDIF

      IF (DISCRETIZE(2) == 0) THEN               ! 0 & 1 => first order upwinding
         DO M = 1, MMAX
            CALL CONV_ROP0 (ROP_s(1,M), U_s(1,M), V_s(1,M), W_s(1,M), &
                            ROP_sE(1,M), ROP_sN(1,M), ROP_sT(1,M))
        ENDDO
      ELSE
         DO M = 1, MMAX
            CALL CONV_ROP1 (DISCRETIZE(2), ROP_s(1,M), &
                            U_s(1,M), V_s(1,M), W_s(1,M), &
                            ROP_sE(1,M), ROP_sN(1,M), ROP_sT(1,M))
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CONV_ROP


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP                                                C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. FOU routine.                                     C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP0(ROP, U, V, W, ROP_E, ROP_N, ROP_T)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      USE functions, only: east_of, north_of, top_of
      USE functions, only: west_of, south_of, bottom_of
      USE functions, only: im_of, jm_of, km_of
      USE geometry, only: do_k
      USE param, only: dimension_3
      USE param1, only: zero
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! macroscopic density (rho_prime)
      DOUBLE PRECISION, INTENT(IN) :: ROP(DIMENSION_3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)
! Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_T(DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! Indices
      INTEGER :: IJK
      INTEGER :: IJKE, IJKN, IJKT
      INTEGER :: IJKW, IJKS, IJKB
      INTEGER :: IMJK, IJMK, IJKM
!---------------------------------------------------------------------//


!$omp  parallel do default(none) &
!$omp              private(IJK, IJKE, IJKN, IJKT, IJKW, &
!$omp                      IJKS, IJKB, IMJK, IJMK, IJKM) &
!$omp              shared(ijkstart3, ijkend3, u, v, w, do_k, rop, &
!$omp                     rop_e, rop_n, rop_t)
      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK)) THEN
            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKT = TOP_OF(IJK)

            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1/2, j, k)
            IF (U(IJK) >= ZERO) THEN
               ROP_E(IJK) = ROP(IJK)
            ELSE
               ROP_E(IJK) = ROP(IJKE)
            ENDIF
! West face (i-1/2, j, k)
            IF (.NOT.FLUID_AT(IMJK)) THEN
               IJKW = WEST_OF(IJK)
               IF (U(IMJK) >= ZERO) THEN
                  ROP_E(IMJK) = ROP(IJKW)
               ELSE
                  ROP_E(IMJK) = ROP(IJK)
               ENDIF
            ENDIF


! North face (i, j+1/2, k)
            IF (V(IJK) >= ZERO) THEN
               ROP_N(IJK) = ROP(IJK)
            ELSE
               ROP_N(IJK) = ROP(IJKN)
            ENDIF
! South face (i, j-1/2, k)
            IF (.NOT.FLUID_AT(IJMK)) THEN
               IJKS = SOUTH_OF(IJK)
               IF (V(IJMK) >= ZERO) THEN
                 ROP_N(IJMK) = ROP(IJKS)
               ELSE
                 ROP_N(IJMK) = ROP(IJK)
               ENDIF
            ENDIF


            IF (DO_K) THEN
               IJKM = KM_OF(IJK)
! Top face (i, j, k+1/2)
               IF (W(IJK) >= ZERO) THEN
                  ROP_T(IJK) = ROP(IJK)
               ELSE
                  ROP_T(IJK) = ROP(IJKT)
               ENDIF
! Bottom face (i, j, k-1/2)
               IF (.NOT.FLUID_AT(IJKM)) THEN
                  IJKB = BOTTOM_OF(IJK)
                  IF (W(IJKM) >= ZERO) THEN
                     ROP_T(IJKM) = ROP(IJKB)
                  ELSE
                     ROP_T(IJKM) = ROP(IJK)
                  ENDIF
               ENDIF
            ENDIF   ! end if do_k

         ENDIF   ! end if fluid_at
      ENDDO    ! end do ijk

      RETURN
      END SUBROUTINE CONV_ROP0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP1                                               C
!  Purpose: Calculate the face value of density used for calculating   C
!  convection fluxes. HR routine.  Here interpolate the face value of  C
!  density.                                                            C
!                                                                      C
!  Author: M. Syamlal                                 Date: 31-MAY-05  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CONV_ROP1(DISC, ROP, U, V, W, ROP_E, ROP_N, ROP_T)

! Modules
!---------------------------------------------------------------------//
      USE compar, only: ijkstart3, ijkend3
      USE functions, only: fluid_at
      USE functions, only: east_of, north_of, top_of
      USE functions, only: west_of, south_of, bottom_of
      USE functions, only: im_of, jm_of, km_of
      USE geometry, only: do_k
      USE param, only: dimension_3
      USE param1, only: one
      USE xsi, only: calc_xsi
      IMPLICIT NONE

! Dummy arguments
!---------------------------------------------------------------------//
! Discretization scheme
      INTEGER, INTENT(IN) :: DISC
! macroscopic density (rho_prime)
      DOUBLE PRECISION, INTENT(IN) :: ROP(DIMENSION_3)
! Velocity components
      DOUBLE PRECISION, INTENT(IN) :: U(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: V(DIMENSION_3)
      DOUBLE PRECISION, INTENT(IN) :: W(DIMENSION_3)
! Face value of density (for calculating convective fluxes)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_E(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_N(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: ROP_T(DIMENSION_3)
!
! Local variables
!---------------------------------------------------------------------//
      INTEGER :: IJK, IJKE, IJKN, IJKT
      INTEGER :: IJKW, IJKS, IJKB, IMJK, IJMK, IJKM
      Integer :: incr

      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: XSI_e, XSI_n, XSI_t

!---------------------------------------------------------------------//

! Calculate factors
      incr=0
      CALL CALC_XSI (DISC, ROP, U, V, W, XSI_E, XSI_N, XSI_T, incr)

!!!$omp  parallel do private(IJK, IJKE, IJKN, IJKT, IJKW, IJKS, IJKB, &
!!!$omp                      IMJK, IJMK, IJKM) &
!!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3

         IF (FLUID_AT(IJK)) THEN
            IJKE = EAST_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKT = TOP_OF(IJK)

            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)

! East face (i+1/2, j, k)
            ROP_E(IJK) = ((ONE-XSI_E(IJK))*ROP(IJK)+&
                         XSI_E(IJK)*ROP(IJKE))
! West face (i-1/2, j, k)
            IF (.NOT.FLUID_AT(IMJK)) THEN
               IJKW = WEST_OF(IJK)
               ROP_E(IMJK) = ((ONE - XSI_E(IMJK))*ROP(IJKW)+&
                             XSI_E(IMJK)*ROP(IJK))
            ENDIF


! North face (i, j+1/2, k)
            ROP_N(IJK) = ((ONE-XSI_N(IJK))*ROP(IJK)+&
                         XSI_N(IJK)*ROP(IJKN))
! South face (i, j-1/2, k)
            IF (.NOT.FLUID_AT(IJMK)) THEN
               IJKS = SOUTH_OF(IJK)
               ROP_N(IJMK) = ((ONE - XSI_N(IJMK))*ROP(IJKS)+&
                             XSI_N(IJMK)*ROP(IJK))
            ENDIF


            IF (DO_K) THEN
               IJKM = KM_OF(IJK)

! Top face (i, j, k+1/2)
               ROP_T(IJK) = ((ONE - XSI_T(IJK))*ROP(IJK)+&
                            XSI_T(IJK)*ROP(IJKT))
! Bottom face (i, j, k-1/2)
               IF (.NOT.FLUID_AT(IJKM)) THEN
                  IJKB = BOTTOM_OF(IJK)
                  ROP_T(IJKM) = ((ONE - XSI_T(IJKM))*ROP(IJKB)+&
                                XSI_T(IJKM)*ROP(IJK))
               ENDIF
            ENDIF   ! end if do_k

         ENDIF   ! end if fluid_at
      ENDDO    ! end do ijk

      RETURN
      END SUBROUTINE CONV_ROP1
