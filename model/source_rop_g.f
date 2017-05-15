!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_ROP_g                                            C
!  Purpose: Determine source terms for continuity equation.            C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative.              C
!         See conv_rop_g for additional details.                       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 2 -JUL-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_ROP_G(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE rxns
      USE run
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
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! DEL dot V
      DOUBLE PRECISION :: DEL_V
! Mass source
      DOUBLE PRECISION :: Src
! Indices
      INTEGER :: I, J, K, IJK, IMJK, IJMK, IJKM
! error message
      CHARACTER(LEN=80) :: LINE

!-----------------------------------------------

!!$omp  parallel do private( I, J, K, IJK, IMJK, IJMK, IJKM,  DEL_V, &
!!$omp&  Src, LINE) &
!!$omp&  schedule(static)
      DO IJK = ijkstart3, ijkend3
         IF (FLUID_AT(IJK) .AND. PHASE_4_P_G(IJK)/=0) THEN
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

            DEL_V = U_G(IJK)*AYZ(IJK) - U_G(IMJK)*AYZ(IMJK) + &
                    V_G(IJK)*AXZ(IJK) - V_G(IJMK)*AXZ(IJMK) + &
                    W_G(IJK)*AXY(IJK) - W_G(IJKM)*AXY(IJKM)

            IF (ROP_G(IJK) > ZERO) THEN
               SRC = VOL(IJK)*ZMAX((-SUM_R_G(IJK)))/ROP_G(IJK)
            ELSE
               SRC = ZERO
            ENDIF

            A_M(IJK,0,0) = -(A_M(IJK,east,0)+A_M(IJK,west,0)+A_M(IJK,north,0)+&
                             A_M(IJK,south,0)+A_M(IJK,top,0)+A_M(IJK,bottom,0)+&
                             VOL(IJK)*ODT+ZMAX(DEL_V)+SRC)
            B_M(IJK,0) = -(ROP_GO(IJK)*VOL(IJK)*ODT+&
                           ZMAX((-DEL_V))*ROP_G(IJK)+&
                           ZMAX(SUM_R_G(IJK))*VOL(IJK))

            IF (ABS(A_M(IJK,0,0)) < SMALL_NUMBER) THEN
               IF (ABS(B_M(IJK,0)) < SMALL_NUMBER) THEN
                  A_M(IJK,0,0) = -ONE
                  B_M(IJK,0) = ZERO
               ELSE
!!$omp             critical
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', 0, ' A = 0 and b = ', B_M(IJK,0)
                  CALL WRITE_ERROR ('SOURCE_ROP_g', LINE, 1)
!!$omp             end critical
               ENDIF
            ENDIF
         ELSE
! set the value of rop_g in all wall and flow boundary cells to what is
! known for that cell
            A_M(IJK,east,0) = ZERO
            A_M(IJK,west,0) = ZERO
            A_M(IJK,north,0) = ZERO
            A_M(IJK,south,0) = ZERO
            A_M(IJK,top,0) = ZERO
            A_M(IJK,bottom,0) = ZERO
            A_M(IJK,0,0) = -ONE
            B_M(IJK,0) = -ROP_G(IJK)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE SOURCE_ROP_G


