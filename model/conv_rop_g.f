!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP_g                                              C
!  Purpose: Determine convection terms for gas continuity              C
!           equation.  Master routine.                                 C
!                                                                      C
!  Author: M. Syamlal                                 Date: 19-MAR-97  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_ROP_G(A_M, B_M)

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
!-----------------------------------------------

      IF (DISCRETIZE(1) == 0) THEN               ! 0 & 1 => first order upwinding
         CALL CONV_ROP_G0 (A_M)
      ELSE
         CALL CONV_ROP_G1 (A_M)
      ENDIF


      RETURN
      END SUBROUTINE CONV_ROP_G


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP_g0                                             C
!  Purpose: Determine convection terms for gas continuity equation     C
!                                                                      C
!  Notes: The off-diagonal coefficients calculated here must be        C
!         positive. The center coefficient and the source vector       C
!         are negative.  -- First order upwinding                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 2-JUL-96   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_ROP_G0(A_M)

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
      USE compar
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IPJK, IJPK, IJKP
      INTEGER :: IMJK, IJMK, IJKM
!-----------------------------------------------


! Calculate convection fluxes through each of the faces

!!$omp  parallel do private( J, K, IJK, IPJK, IJPK, IJKP, &
!!$omp&  IMJK, IJMK, IJKM) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (PHASE_4_P_G(IJK) /= 0) THEN
! if phase_4_p_g is assigned the index of a solids phase then this
! section is entered. currently, phase_4_p_g is always assigned to
! the gas phase index 0 and so this section will never be entered.

! it was previously possible for phase_4_p_g to be assigned the index of
! the gas phase in some cells, and the index of a solids phase in other
! cells if that solids phase has close_packed=F and was in higher
! concentrations than the gas. in such a case this branch would become
! activated, while the corresponding solids continuity would become
! skipped. moreover, that solids phase's volume fraction would be
! corrected rather than the gas phases void fraction in calc_vol_fr.
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
            A_M(IJK,east,0) = ZMAX((-U_G(IJK)))*AYZ(IJK)
            A_M(IPJK,west,0) = ZMAX(U_G(IJK))*AYZ(IJK)

! North face (i, j+1/2, k)
            A_M(IJK,north,0) = ZMAX((-V_G(IJK)))*AXZ(IJK)
            A_M(IJPK,south,0) = ZMAX(V_G(IJK))*AXZ(IJK)

! Top face (i, j, k+1/2)
            IF (DO_K) THEN
               A_M(IJK,top,0) = ZMAX((-W_G(IJK)))*AXY(IJK)
               A_M(IJKP,bottom,0) = ZMAX(W_G(IJK))*AXY(IJK)
            ENDIF

! Modify west (i-1/2,j,k), south (i j-1/2,k) and bottom (i,j,k-1/2)
! faces if the neighboring west, south, bottom cells have
! phase_4_p_g of 0
            IF(PHASE_4_P_G(IMJK)==0) A_M(IJK,west,0)=ZMAX(U_G(IMJK))*AYZ(IMJK)
            IF(PHASE_4_P_G(IJMK)==0) A_M(IJK,south,0)=ZMAX(V_G(IJMK))*AXZ(IJMK)
            IF (DO_K) THEN
               IF(PHASE_4_P_G(IJKM)==0) A_M(IJK,bottom,0)=ZMAX(W_G(IJKM))*AXY(IJKM)
            ENDIF

         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_ROP_G0


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CONV_ROP_g1                                             C
!  Purpose: Determine convection terms for gas continuity equation.    C
!                                                                      C
!  Notes: The off-diagonal coefficients calculated here must be        C
!         positive. The center coefficient and the source vector       C
!         are negative.  -- Higher order schemes                       C
!                                                                      C
!  Author: M. Syamlal                                 Date:19-MAR-97   C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CONV_ROP_G1(A_M)

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
      USE run
      USE xsi
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
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


! loezos
      incr=0

! Calculate convection factors
      CALL CALC_XSI (DISCRETIZE(1), ROP_G, U_G, V_G, W_G, XSI_E, &
         XSI_N, XSI_T,incr)

! Calculate convection fluxes through each of the faces

!!$omp  parallel do private( J, K, IJK, IPJK, IJPK, IJKP,  &
!!$omp&  IMJK, IJMK, IJKM) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (PHASE_4_P_G(IJK) /= 0) THEN
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
            A_M(IJK,east,0) = -XSI_E(IJK)*U_G(IJK)*AYZ(IJK)
            A_M(IPJK,west,0) = (ONE - XSI_E(IJK))*U_G(IJK)*AYZ(IJK)

! North face (i, j+1/2, k)
            A_M(IJK,north,0) = -XSI_N(IJK)*V_G(IJK)*AXZ(IJK)
            A_M(IJPK,south,0) = (ONE - XSI_N(IJK))*V_G(IJK)*AXZ(IJK)

! Top face (i, j, k+1/2)
            IF (DO_K) THEN
               A_M(IJK,top,0) = -XSI_T(IJK)*W_G(IJK)*AXY(IJK)
               A_M(IJKP,bottom,0) = (ONE - XSI_T(IJK))*W_G(IJK)*AXY(IJK)
            ENDIF

! Modify west (i-1/2,j,k), south (i j-1/2,k) and bottom (i,j,k-1/2)
! faces if the neighboring west, south, bottom cells have
! phase_4_p_g of 0
            IF (PHASE_4_P_G(IMJK) == 0) A_M(IJK,west,0) = (ONE - XSI_E(IMJK))*U_G(&
               IMJK)*AYZ(IMJK)
            IF (PHASE_4_P_G(IJMK) == 0) A_M(IJK,south,0) = (ONE - XSI_N(IJMK))*V_G(&
               IJMK)*AXZ(IJMK)
            IF (DO_K) THEN
               IF (PHASE_4_P_G(IJKM) == 0) A_M(IJK,bottom,0) = (ONE - XSI_T(IJKM))*&
                  W_G(IJKM)*AXY(IJKM)
            ENDIF

         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CONV_ROP_G1


