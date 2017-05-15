!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP_s                                              C
!  Purpose: Determine convection terms for solids continuity           C
!           equation.  Master routine.                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-MAR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_ROP_S(A_M, B_M, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE run
      USE compar
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Solids phase index
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------

      IF (DISCRETIZE(2) == 0) THEN               ! 0 & 1 => first order upwinding
         CALL CONV_ROP_S0 (A_M, M)
      ELSE
         CALL CONV_ROP_S1 (A_M, M)
      ENDIF

      RETURN
      END SUBROUTINE CONV_ROP_S


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP_s0                                             C
!  Purpose: Determine convection terms for solids continuity equation. C
!                                                                      C
!  Notes: The off-diagonal coefficients calculated here must be        C
!         positive. The center coefficient and the source vector       C
!         are negative. -- First order upwinding                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 3-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_ROP_S0(A_M, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE run
      USE parallel
      USE physprop
      USE geometry
      USE indices
      USE pgcor
      USE pscor
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Solids phase index
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM
!-----------------------------------------------

! Calculate convection-diffusion fluxes through each of the faces

!!$omp  parallel do private( I, J, K, IJK, IPJK, IJPK, IJKP,  &
!!$omp&  IMJK, IJMK, IJKM) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (PHASE_4_P_G(IJK)/=M .AND. PHASE_4_P_S(IJK)/=M) THEN
! if either phase_4_p_g or phase_4_p_s are assigned the index of the
! current solids phase then this section is skipped. currently,
! phase_4_p_g and phase_4_p_s are always assigned to 0 and undefined,
! respectively.  Hence this section should always be entered.

! it was previously possible for phase_4_p_g to be assigned the index of
! the gas phase in some cells, and the index of a solids phase in other
! cells if that solids phase has close_packed=F and was in higher
! concentrations than the gas. in such a case this branch would become
! skipped, while the corresponding gas continuity would become
! activated. moreover, that solids phase's volume fraction would be
! corrected rather than the gas phases void fraction in calc_vol_fr.

! if phase_4_p_s were to be assigned the index of a solids phase that
! can close pack if the current cell exhibited close packing, then
! this branch would become skipped.  moreover, that solids phase's
! volume fraction would then be corected based on the value of
! maximum packing in that cell and the sum of all other solids that
! can close pack in calc_vol_fr.
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

! East face (i+1/2, j, k)
            A_M(IJK,east,M) = ZMAX((-U_S(IJK,M)))*AYZ(IJK)
            A_M(IPJK,west,M) = ZMAX(U_S(IJK,M))*AYZ(IJK)

! North face (i, j+1/2, k)
            A_M(IJK,north,M) = ZMAX((-V_S(IJK,M)))*AXZ(IJK)
            A_M(IJPK,south,M) = ZMAX(V_S(IJK,M))*AXZ(IJK)

! Top face (i, j, k+1/2)
            IF (DO_K) THEN
               A_M(IJK,top,M) = ZMAX((-W_S(IJK,M)))*AXY(IJK)
               A_M(IJKP,bottom,M) = ZMAX(W_S(IJK,M))*AXY(IJK)
            ENDIF

! Modify west (i-1/2,j,k), south (i j-1/2,k) and bottom (i,j,k-1/2)
! faces if the neighboring west, south, bottom cells have
! phase_4_p_g of m or phase_4_p_s of m.
            IF (PHASE_4_P_G(IMJK)==M .OR. PHASE_4_P_S(IMJK)==M) &
               A_M(IJK,west,M) = ZMAX(U_S(IMJK,M))*AYZ(IMJK)
            IF (PHASE_4_P_G(IJMK)==M .OR. PHASE_4_P_S(IJMK)==M) &
               A_M(IJK,south,M) = ZMAX(V_S(IJMK,M))*AXZ(IJMK)

            IF (DO_K) THEN
               IF (PHASE_4_P_G(IJKM)==M .OR. PHASE_4_P_S(IJKM)==M) &
                  A_M(IJK,bottom,M) = ZMAX(W_S(IJKM,M))*AXY(IJKM)
            ENDIF
         ENDIF
      END DO

      RETURN
      END SUBROUTINE CONV_ROP_S0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP_s1                                             C
!  Purpose: Determine convection terms for solids continuity equation. C
!  Notes: The off-diagonal coefficients calculated here must be        C
!         positive. The center coefficient and the source vector       C
!         are negative. -- Higher order methods                        C
!         See conv_rop_s0 for additional details.                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 18-MAR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_ROP_S1(A_M, M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE compar
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
      USE xsi
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Solids phase index
      INTEGER, INTENT(IN) :: M
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM
! loezos
      INTEGER :: incr
! temporary use of global arrays:
! xsi_array: convection weighting factors
      DOUBLE PRECISION :: XSI_e(DIMENSION_3), XSI_n(DIMENSION_3),&
                          XSI_t(DIMENSION_3)
!-----------------------------------------------

! Loezos:
      incr=0

! Calculate convection factors
      CALL CALC_XSI(DISCRETIZE(2), ROP_S(1,M), U_S(1,M), V_S(1,M),&
         W_S(1,M), XSI_E, XSI_N, XSI_T, incr)

! Calculate convection-diffusion fluxes through each of the faces

!!$omp  parallel do private( I, J, K, IJK, IPJK, IJPK, IJKP,  &
!!$omp&  IMJK, IJMK, IJKM) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (PHASE_4_P_G(IJK)/=M .AND. PHASE_4_P_S(IJK)/=M) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IPJK = IP_OF(IJK)
            IJPK = JP_OF(IJK)
            IJKP = KP_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

! East face (i+1/2, j, k)
            A_M(IJK,east,M) = -XSI_E(IJK)*U_S(IJK,M)*AYZ(IJK)
            A_M(IPJK,west,M) = (ONE - XSI_E(IJK))*U_S(IJK,M)*AYZ(IJK)

! North face (i, j+1/2, k)
            A_M(IJK,north,M) = -XSI_N(IJK)*V_S(IJK,M)*AXZ(IJK)
            A_M(IJPK,south,M) = (ONE - XSI_N(IJK))*V_S(IJK,M)*AXZ(IJK)

! Top face (i, j, k+1/2)
            IF (DO_K) THEN
               A_M(IJK,top,M) = -XSI_T(IJK)*W_S(IJK,M)*AXY(IJK)
               A_M(IJKP,bottom,M) = (ONE - XSI_T(IJK))*W_S(IJK,M)*AXY(IJK)
            ENDIF

            IF (PHASE_4_P_G(IMJK)==M .OR. PHASE_4_P_S(IMJK)==M) A_M(IJK,west,M) = &
               (ONE - XSI_E(IMJK))*U_S(IMJK,M)*AYZ(IMJK)
            IF (PHASE_4_P_G(IJMK)==M .OR. PHASE_4_P_S(IJMK)==M) A_M(IJK,south,M) = &
               (ONE - XSI_N(IJMK))*V_S(IJMK,M)*AXZ(IJMK)
            IF (DO_K) THEN
               IF (PHASE_4_P_G(IJKM)==M .OR. PHASE_4_P_S(IJKM)==M) A_M(IJK,bottom,M)&
                   = (ONE - XSI_T(IJKM))*W_S(IJKM,M)*AXY(IJKM)
            ENDIF
         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE CONV_ROP_S1

