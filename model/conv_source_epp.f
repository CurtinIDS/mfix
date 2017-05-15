!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_SOURCE_EPp                                         C
!  Purpose: Determine convection & source terms for solids volume      C
!           fraction correction equation.  Master routine.             C
!                                                                      C
!  Notes: MCP must be defined to call this routine.                    C
!                                                                      C
!  Author: M. Syamlal                                 Date: 6-MAR-97   C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_SOURCE_EPP(A_M, B_M, B_mmax, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE fldvar
      USE geometry
      USE param
      USE param1
      USE run
      USE sendrecv
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax(DIMENSION_3, 0:DIMENSION_M)
! Lowest solids phase index of those solids phases that can
! close packed (M=MCP)
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------

      IF (DISCRETIZE(2) == 0) THEN               ! 0 & 1 => first order upwinding
         CALL CONV_SOURCE_EPP0 (A_M, B_M, B_MMAX, M)
      ELSE
         CALL CONV_SOURCE_EPP1 (A_M, B_M, B_MMAX, M)
      ENDIF

      RETURN
      END SUBROUTINE CONV_SOURCE_EPP

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_SOURCE_EPp0                                        C
!  Purpose: Determine convection & source terms for solids volume      C
!           fraction correction equation.  First order upwinding.      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-SEP-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_SOURCE_EPP0(A_M, B_M, B_MMAX, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE pgcor
      USE physprop
      USE pscor
      USE run
      USE rxns
      USE sendrecv
      USE solids_pressure
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax(DIMENSION_3, 0:DIMENSION_M)
! Lowest solids phase index of those solids phases that can
! close packed (M=MCP)
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM, IJKE, IJKW, IJKN, IJKS
      INTEGER :: IJKB, IJKT
! dPodEP_s(EP_s(IJK, M)); the derivative of the plastic or friction
! pressure with respect to ep_s
      DOUBLE PRECISION :: K_P
! Mass source
      DOUBLE PRECISION :: Src
! error message
      CHARACTER(LEN=80) :: LINE(1)
! terms of bm expression
      DOUBLE PRECISION :: bma, bme, bmw, bmn, bms, bmt, bmb, bmr
!-----------------------------------------------

!!$omp    parallel do                                            &
!!$omp&   private(I, J, K, IJK, IPJK, IJPK, IJKP,                &
!!$omp&           IMJK, IJMK, IJKM,                              &
!!$omp&           IJKE, IJKW, IJKN, IJKS, IJKT, IJKB,            &
!!$omp&           K_P, SRC, bma, bme, bmw, bmn, bms, bmt, bmb, bmr )
      DO IJK = ijkstart3, ijkend3
! Determine whether IJK falls within 1 ghost layer........
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF (FLUID_AT(IJK)) THEN
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

            IJKE = EAST_OF(IJK)
            IJKW = WEST_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKS = SOUTH_OF(IJK)
            IJKT = TOP_OF(IJK)
            IJKB = BOTTOM_OF(IJK)

! initializing
            A_M(IJK,0,0) = ZERO
            B_M(IJK,0) = ZERO
            K_P = K_CP(IJK)

! Calculate convection-diffusion fluxes through each of the faces

! East face (i+1/2, j, k)
            IF (U_S(IJK,M) < ZERO) THEN
               A_M(IJK,east,0) = (ROP_S(IJKE,M)*E_E(IJK)*K_CP(IJKE)-&
                 RO_S(IJKE,M)*U_S(IJK,M))*AYZ(IJK)
               A_M(IJK,0,0) = A_M(IJK,0,0)+&
                  ROP_S(IJKE,M)*E_E(IJK)*K_P*AYZ(IJK)
               bme = (-ROP_S(IJKE,M)*U_S(IJK,M))*AYZ(IJK)
               B_M(IJK,0) = B_M(IJK,0) +  bme
            ELSE
               A_M(IJK,east,0) = (ROP_S(IJK,M)*&
                  E_E(IJK)*K_CP(IJKE))*AYZ(IJK)
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*&
                  E_E(IJK)*K_P+RO_S(IJK,M)*U_S(IJK,M))*AYZ(IJK)
               bme = (-ROP_S(IJK,M)*U_S(IJK,M))*AYZ(IJK)
               B_M(IJK,0) = B_M(IJK,0) +  bme
            ENDIF


! West face (i-1/2, j, k)
            IF (U_S(IMJK,M) > ZERO) THEN
               A_M(IJK,west,0) = (ROP_S(IJKW,M)*E_E(IMJK)*K_CP(IJKW)+&
                  RO_S(IJKW,M)*U_S(IMJK,M))*AYZ(IMJK)
               A_M(IJK,0,0) = A_M(IJK,0,0) + &
                  ROP_S(IJKW,M)*E_E(IMJK)*K_P*AYZ(IMJK)
               bmw = (ROP_S(IJKW,M)*U_S(IMJK,M))*AYZ(IMJK)
               B_M(IJK,0) = B_M(IJK,0) + bmw
            ELSE
               A_M(IJK,west,0) = (ROP_S(IJK,M)*&
                  E_E(IMJK)*K_CP(IJKW))*AYZ(IMJK)
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*&
                  E_E(IMJK)*K_P-RO_S(IJK,M)*U_S(IMJK,M))*AYZ(IMJK)
               bmw = (ROP_S(IJK,M)*U_S(IMJK,M))*AYZ(IMJK)
               B_M(IJK,0) = B_M(IJK,0) + bmw
            ENDIF


! North face (i, j+1/2, k)
            IF (V_S(IJK,M) < ZERO) THEN
               A_M(IJK,north,0) = (ROP_S(IJKN,M)*E_N(IJK)*K_CP(IJKN)-&
                  RO_S(IJKN,M)*V_S(IJK,M))*AXZ(IJK)
               A_M(IJK,0,0)=A_M(IJK,0,0)+&
                  ROP_S(IJKN,M)*E_N(IJK)*K_P*AXZ(IJK)
               bmn = (-ROP_S(IJKN,M)*V_S(IJK,M))*AXZ(IJK)
               B_M(IJK,0) = B_M(IJK,0) + bmn
            ELSE
               A_M(IJK,north,0) = (ROP_S(IJK,M)*&
                  E_N(IJK)*K_CP(IJKN))*AXZ(IJK)
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*&
                  E_N(IJK)*K_P+RO_S(IJK,M)*V_S(IJK,M))*AXZ(IJK)
               bmn = (-ROP_S(IJK,M)*V_S(IJK,M))*AXZ(IJK)
               B_M(IJK,0) = B_M(IJK,0) + bmn
            ENDIF


! South face (i, j-1/2, k)
            IF (V_S(IJMK,M) > ZERO) THEN
               A_M(IJK,south,0) = (ROP_S(IJKS,M)*E_N(IJMK)*K_CP(IJKS)+&
                  RO_S(IJKS,M)*V_S(IJMK,M))*AXZ(IJMK)
               A_M(IJK,0,0) = A_M(IJK,0,0) + &
                  ROP_S(IJKS,M)*E_N(IJMK)*K_P*AXZ(IJMK)
               bms = (ROP_S(IJKS,M)*V_S(IJMK,M))*AXZ(IJMK)
               B_M(IJK,0) = B_M(IJK,0) + bms
            ELSE
               A_M(IJK,south,0) = (ROP_S(IJK,M)*&
                  E_N(IJMK)*K_CP(IJKS))*AXZ(IJMK)
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*&
                  E_N(IJMK)*K_P-RO_S(IJK,M)*V_S(IJMK,M))*AXZ(IJMK)
               bms = (ROP_S(IJK,M)*V_S(IJMK,M))*AXZ(IJMK)
               B_M(IJK,0) = B_M(IJK,0) + bms
            ENDIF


            IF (DO_K) THEN
! Top face (i, j, k+1/2)
               IF (W_S(IJK,M) < ZERO) THEN
                  A_M(IJK,top,0) = (ROP_S(IJKT,M)*E_T(IJK)*K_CP(IJKT)-&
                     RO_S(IJKT,M)*W_S(IJK,M))*AXY(IJK)
                  A_M(IJK,0,0) = A_M(IJK,0,0) + &
                     ROP_S(IJKT,M)*E_T(IJK)*K_P*AXY(IJK)
                  bmt = (-ROP_S(IJKT,M)*W_S(IJK,M))*AXY(IJK)
                  B_M(IJK,0)=B_M(IJK,0) + bmt
               ELSE
                  A_M(IJK,top,0) = (ROP_S(IJK,M)*&
                     E_T(IJK)*K_CP(IJKT))*AXY(IJK)
                  A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*&
                     E_T(IJK)*K_P+RO_S(IJK,M)*W_S(IJK,M))*AXY(IJK)
                  bmt = (-ROP_S(IJK,M)*W_S(IJK,M))*AXY(IJK)
                  B_M(IJK,0) = B_M(IJK,0) + bmt
               ENDIF


! Bottom face (i, j, k-1/2)
               IF (W_S(IJKM,M) > ZERO) THEN
                  A_M(IJK,bottom,0) = (ROP_S(IJKB,M)*E_T(IJKM)*K_CP(IJKB)+&
                     RO_S(IJKB,M)*W_S(IJKM,M))*AXY(IJKM)
                  A_M(IJK,0,0) = A_M(IJK,0,0) + &
                     ROP_S(IJKB,M)*E_T(IJKM)*K_P*AXY(IJKM)
                  bmb = (ROP_S(IJKB,M)*W_S(IJKM,M))*AXY(IJKM)
                  B_M(IJK,0) = B_M(IJK,0) + bmb
               ELSE
                  A_M(IJK,bottom,0) = (ROP_S(IJK,M)*&
                     E_T(IJKM)*K_CP(IJKB))*AXY(IJKM)
                  A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_S(IJK,M)*&
                     E_T(IJKM)*K_P-RO_S(IJK,M)*W_S(IJKM,M))*AXY(IJKM)
                  bmb = (ROP_S(IJK,M)*W_S(IJKM,M))*AXY(IJKM)
                  B_M(IJK,0)=B_M(IJK,0) + bmb
               ENDIF

            ELSE   ! not(DO_K) branch
               bmt = zero
               bmb = zero
            ENDIF   ! end if/else (do_K)

            IF (ROP_S(IJK,M)>ZERO .AND. SUM_R_S(IJK,M)<ZERO) THEN
               SRC = VOL(IJK)*(-SUM_R_S(IJK,M))/ROP_S(IJK,M)
            ELSE
               SRC = ZERO
            ENDIF

            A_M(IJK,0,0) = -(A_M(IJK,0,0)+VOL(IJK)*ODT*RO_S(IJK,M)+SRC*RO_S(IJK,M))

            bma = (ROP_S(IJK,M)-ROP_SO(IJK,M))*VOL(IJK)*ODT
            bmr = SUM_R_S(IJK,M)*VOL(IJK)
            B_M(IJK,0) = -(B_M(IJK,0) - bma + bmr)
            B_MMAX(IJK,0) = max(abs(bma), abs(bme), abs(bmw), abs(bmn), abs(bms), abs(bmt), abs(bmb), abs(bmr) )

            IF ((-A_M(IJK,0,0)) < SMALL_NUMBER) THEN
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN
                  A_M(IJK,0,0) = -ONE            ! Equation is undefined.
                  B_M(IJK,0) = ZERO              ! Use existing value
               ELSE
!!$omp             critical
                  WRITE (LINE, '(A,I6,A,G12.5)') &
                     'Error: At IJK = ', IJK, &
                     ' A = 0 and b = ', B_M(IJK,0)
                  CALL WRITE_ERROR ('CONV_SOURCE_EPp0', LINE, 1)
!!$omp             end critical
               ENDIF
            ENDIF
         ELSE
            A_M(IJK,0,0) = -ONE
            B_M(IJK,0) = ZERO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_SOURCE_EPP0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_SOURCE_EPp1                                        C
!  Purpose: Determine convection & source terms for solids volume      C
!           fraction correction equation.  Higher order scheme.        C
!                                                                      C
!  Author: M. Syamlal                                 Date: 24-SEP-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_SOURCE_EPP1(A_M, B_M, B_MMAX, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
      USE constant
      USE fldvar
      USE functions
      USE geometry
      USE indices
      USE parallel
      USE param
      USE param1
      USE pgcor
      USE physprop
      USE pscor
      USE run
      USE rxns
      USE sendrecv
      USE solids_pressure
      USE vshear
      USE xsi
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax(DIMENSION_3, 0:DIMENSION_M)
! Lowest solids phase index of those solids phases that can
! close packed (M=MCP)
      INTEGER, INTENT(IN) :: M

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM, IJKE, IJKW, IJKN, IJKS
      INTEGER :: IJKB, IJKT
! loezos : used for including shearing
      INTEGER :: incr
! dPodEP_s(EP_s(IJK, M)); the derivative of the plastic or friction
! pressure with respect to ep_s
      DOUBLE PRECISION :: K_P
! Mass source
      DOUBLE PRECISION :: Src
! face value of ROP_s
      DOUBLE PRECISION :: ROP_sf
! terms of bm expression
      DOUBLE PRECISION :: bma, bme, bmw, bmn, bms, bmt, bmb, bmr
! error message
      CHARACTER(LEN=80) :: LINE(1)
! temporary use of global arrays:
! xsi_array: convection weighting factors
      DOUBLE PRECISION :: XSI_e(DIMENSION_3), &
                          XSI_n(DIMENSION_3),&
                          XSI_t(DIMENSION_3)
!-----------------------------------------------

! Loezos:
      incr=0

! Calculate convection factors
      CALL CALC_XSI (DISCRETIZE(2), ROP_S(1,M), U_S(1,M), V_S(1,M), W_S(1,M), &
         XSI_E, XSI_N, XSI_T,incr)

! Loezos:
! update to true velocity
      IF (SHEAR) THEN
!!$omp parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN
               V_S(IJK,m)=V_s(IJK,m)+VSH(IJK)
            ENDIF
         ENDDO
      ENDIF

!!$omp parallel do                                                      &
!!$omp&   private(I, J, K, IJK, IPJK, IJPK, IJKP,                &
!!$omp&           IMJK, IJMK, IJKM, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB, &
!!$omp&           K_P,ROP_SF,SRC, bma, bme, bmw, bmn, bms, bmt, bmb, bmr )
      DO IJK = ijkstart3, ijkend3
! Determine if IJK falls within 1 ghost layer........
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF (FLUID_AT(IJK)) THEN
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

            IJKE = EAST_OF(IJK)
            IJKW = WEST_OF(IJK)
            IJKN = NORTH_OF(IJK)
            IJKS = SOUTH_OF(IJK)
            IJKT = TOP_OF(IJK)
            IJKB = BOTTOM_OF(IJK)

! initializing
            A_M(IJK,0,0) = ZERO
            B_M(IJK,0) = ZERO
            K_P = K_CP(IJK)

! Calculate convection-diffusion fluxes through each of the faces

! East face (i+1/2, j, k)
            ROP_SF = ROP_S(IJKE,M)*XSI_E(IJK) + ROP_S(IJK,M)*(ONE - XSI_E(IJK))
            A_M(IJK,east,0) = (ROP_SF*E_E(IJK)*K_CP(IJKE)-RO_S(IJK,M)*U_S(IJK,M)*XSI_E&
               (IJK))*AYZ(IJK)
            A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_E(IJK)*K_P+RO_S(IJK,M)*U_S(IJK,&
               M)*(ONE-XSI_E(IJK)))*AYZ(IJK)
            bme = (-ROP_SF*U_S(IJK,M))*AYZ(IJK)
            B_M(IJK,0) = B_M(IJK,0) + bme

! West face (i-1/2, j, k)
            ROP_SF=ROP_S(IJK,M)*XSI_E(IMJK)+ROP_S(IJKW,M)*(ONE-XSI_E(IMJK))
            A_M(IJK,west,0) = (ROP_SF*E_E(IMJK)*K_CP(IJKW)+RO_S(IJK,M)*U_S(IMJK,M)*(&
               ONE-XSI_E(IMJK)))*AYZ(IMJK)
            A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_E(IMJK)*K_P-RO_S(IJK,M)*U_S(&
               IMJK,M)*XSI_E(IMJK))*AYZ(IMJK)
            bmw = (ROP_SF*U_S(IMJK,M))*AYZ(IMJK)
            B_M(IJK,0) = B_M(IJK,0) +  bmw


! North face (i, j+1/2, k)
            ROP_SF = ROP_S(IJKN,M)*XSI_N(IJK) + ROP_S(IJK,M)*(ONE - XSI_N(IJK))
            A_M(IJK,north,0) = (ROP_SF*E_N(IJK)*K_CP(IJKN)-RO_S(IJK,M)*V_S(IJK,M)*XSI_N&
               (IJK))*AXZ(IJK)
            A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_N(IJK)*K_P+RO_S(IJK,M)*V_S(IJK,&
               M)*(ONE-XSI_N(IJK)))*AXZ(IJK)
            bmn = (-ROP_SF*V_S(IJK,M))*AXZ(IJK)
            B_M(IJK,0) = B_M(IJK,0) + bmn


! South face (i, j-1/2, k)
            ROP_SF=ROP_S(IJK,M)*XSI_N(IJMK)+ROP_S(IJKS,M)*(ONE-XSI_N(IJMK))
            A_M(IJK,south,0) = (ROP_SF*E_N(IJMK)*K_CP(IJKS)+RO_S(IJK,M)*V_S(IJMK,M)*(&
               ONE-XSI_N(IJMK)))*AXZ(IJMK)
            A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_N(IJMK)*K_P-RO_S(IJK,M)*V_S(&
               IJMK,M)*XSI_N(IJMK))*AXZ(IJMK)
            bms = (ROP_SF*V_S(IJMK,M))*AXZ(IJMK)
            B_M(IJK,0) = B_M(IJK,0) + bms


            IF (DO_K) THEN
! Top face (i, j, k+1/2)
               ROP_SF=ROP_S(IJKT,M)*XSI_T(IJK)+ROP_S(IJK,M)*(ONE-XSI_T(IJK))
               A_M(IJK,top,0) = (ROP_SF*E_T(IJK)*K_CP(IJKT)-RO_S(IJK,M)*W_S(IJK,M)*&
                  XSI_T(IJK))*AXY(IJK)
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_T(IJK)*K_P+RO_S(IJK,M)*W_S(&
                  IJK,M)*(ONE-XSI_T(IJK)))*AXY(IJK)
               bmt = (-ROP_SF*W_S(IJK,M))*AXY(IJK)
               B_M(IJK,0) = B_M(IJK,0) + bmt

! Bottom face (i, j, k-1/2)
               ROP_SF = ROP_S(IJK,M)*XSI_T(IJKM) + ROP_S(IJKB,M)*(ONE - XSI_T(&
                  IJKM))
               A_M(IJK,bottom,0) = (ROP_SF*E_T(IJKM)*K_CP(IJKB)+RO_S(IJK,M)*W_S(IJKM,M)*&
                  (ONE-XSI_T(IJKM)))*AXY(IJKM)
               A_M(IJK,0,0) = A_M(IJK,0,0) + (ROP_SF*E_T(IJKM)*K_P-RO_S(IJK,M)*W_S(&
                  IJKM,M)*XSI_T(IJKM))*AXY(IJKM)
               bmb = (ROP_SF*W_S(IJKM,M))*AXY(IJKM)
               B_M(IJK,0) = B_M(IJK,0) + bmb
            ELSE   ! not(do_k) branch
              bmt = zero
              bmb = zero
           ENDIF   ! end if/else (do_k)

            IF (ROP_S(IJK,M)>ZERO .AND. SUM_R_S(IJK,M)<ZERO) THEN
               SRC = VOL(IJK)*(-SUM_R_S(IJK,M))/ROP_S(IJK,M)
            ELSE
               SRC = ZERO
            ENDIF

            A_M(IJK,0,0) = -(A_M(IJK,0,0)+VOL(IJK)*ODT*RO_S(IJK,M)+SRC*RO_S(IJK,M))

            bma = (ROP_S(IJK,M)-ROP_SO(IJK,M))*VOL(IJK)*ODT
            bmr = SUM_R_S(IJK,M)*VOL(IJK)
            B_M(IJK,0) = -(B_M(IJK,0)- bma + bmr)
            B_MMAX(IJK,0) = max(abs(bma), abs(bme), abs(bmw), abs(bmn), abs(bms), abs(bmt), abs(bmb), abs(bmr) ) !

            IF (ABS(A_M(IJK,0,0)) < SMALL_NUMBER) THEN
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN
                  A_M(IJK,0,0) = -ONE            ! Equation is undefined.
                  B_M(IJK,0) = ZERO              ! Use existing value
               ELSE
!!$omp             critical
                  WRITE (LINE(1), '(A,I6,A,G12.5)') &
                     'Error: At IJK = ', IJK, &
                     ' A = 0 and b = ', B_M(IJK,0)
! Having problem to compile this statement on SGI
                  CALL WRITE_ERROR ('CONV_SOURCE_EPp1', LINE, 1)
!!$omp             end critical
               ENDIF
            ENDIF
         ELSE
            A_M(IJK,0,0) = -ONE
            B_M(IJK,0) = ZERO
         ENDIF
      ENDDO

! Loezos:
      IF (SHEAR) THEN
!!$omp parallel do private(IJK)
         DO IJK = ijkstart3, ijkend3
            IF (FLUID_AT(IJK)) THEN
               V_S(IJK,m)=V_s(IJK,m)-VSH(IJK)
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CONV_SOURCE_EPP1

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_EPP                                        C
!  Purpose: Adds point sources to the solids volume fraction           C
!           correction equation.                                       C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative. See          C
!         conv_Pp_g                                                    C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_EPP(B_M, B_MMAX, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar
      use constant
      use geometry
      use indices
      use physprop
      use ps
      use pscor
      use run
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! maximum term in b_m expression
      DOUBLE PRECISION, INTENT(INOUT) :: B_mmax(DIMENSION_3, 0:DIMENSION_M)
! Lowest solids phase index of those solids phases that can
! close packed (M=MCP)
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, I, J, K
      INTEGER :: PSV

! terms of bm expression
      DOUBLE PRECISION :: pSource

!-----------------------------------------------

      PS_LP: do PSV = 1, DIMENSION_PS

         if(.NOT.PS_DEFINED(PSV)) cycle PS_LP
         if(PS_MASSFLOW_S(PSV,M) < small_number) cycle PS_LP

         do k = PS_K_B(PSV), PS_K_T(PSV)
         do j = PS_J_S(PSV), PS_J_N(PSV)
         do i = PS_I_W(PSV), PS_I_E(PSV)

            if(.NOT.IS_ON_myPE_plus2layers(I,J,K)) cycle
            IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
            ijk = funijk(i,j,k)
            if(fluid_at(ijk)) then

               pSource =  PS_MASSFLOW_S(PSV,M) * &
                  (VOL(IJK)/PS_VOLUME(PSV))

               B_M(IJK,0) = B_M(IJK,0) - pSource
               B_MMAX(IJK,0) = max(abs(B_MMAX(IJK,0)), abs(B_M(IJK,0)))
            endif
         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_EPP
