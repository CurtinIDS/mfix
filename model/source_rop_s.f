!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOURCE_ROP_s                                            C
!  Purpose: Determine source terms for solids continuity equation.     C
!                                                                      C
!  Notes: The off-diagonal coefficients are positive. The center       C
!         coefficient and the source vector are negative.              C
!         See conv_rop_s0 for additional details.                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 3 -JUL-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOURCE_ROP_S(A_M, B_M, M)

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
      USE pscor
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
! solids phase index
      INTEGER, INTENT(IN) :: M
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
         IF (FLUID_AT(IJK) .AND. PHASE_4_P_G(IJK)/=M .AND. &
             PHASE_4_P_S(IJK)/=M) THEN

            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            IMJK = IM_OF(IJK)
            IJMK = JM_OF(IJK)
            IJKM = KM_OF(IJK)

            DEL_V = U_S(IJK,M)*AYZ(IJK) - U_S(IMJK,M)*AYZ(IMJK) +&
                    V_S(IJK,M)*AXZ(IJK) - V_S(IJMK,M)*AXZ(IJMK) + &
                    W_S(IJK,M)*AXY(IJK) - W_S(IJKM,M)*AXY(IJKM)

            IF (ROP_S(IJK,M) > ZERO) THEN
               SRC = VOL(IJK)*ZMAX((-SUM_R_S(IJK,M)))/ROP_S(IJK,M)
            ELSE
               SRC = ZERO
            ENDIF

            A_M(IJK,0,M) = -(A_M(IJK,east,M)+A_M(IJK,west,M)+A_M(IJK,north,M)+&
                             A_M(IJK,south,M)+A_M(IJK,top,M)+A_M(IJK,bottom,M)+&
                             VOL(IJK)*ODT+ZMAX(DEL_V)+SRC)
            B_M(IJK,M) = -(ROP_SO(IJK,M)*VOL(IJK)*ODT+&
                           ZMAX((-DEL_V))*ROP_S(IJK,M)+&
                           ZMAX(SUM_R_S(IJK,M))*VOL(IJK))

            IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN
               IF (ABS(B_M(IJK,M)) < SMALL_NUMBER) THEN
                  A_M(IJK,0,M) = -ONE            ! Equation is undefined.
                  B_M(IJK,M) = -ROP_S(IJK,M)     ! Use existing value
               ELSE
!!$omp             critical
                  WRITE (LINE, '(A,I6,A,I1,A,G12.5)') 'Error: At IJK = ', IJK, &
                     ' M = ', M, ' A = 0 and b = ', B_M(IJK,M)
                  CALL WRITE_ERROR ('SOURCE_ROP_s', LINE, 1)
!!$omp             end critical
               ENDIF
            ENDIF
         ELSE
! set the value of rop_s in all wall and flow boundary cells to what is
! known for that cell
            A_M(IJK,east,M) = ZERO
            A_M(IJK,west,M) = ZERO
            A_M(IJK,north,M) = ZERO
            A_M(IJK,south,M) = ZERO
            A_M(IJK,top,M) = ZERO
            A_M(IJK,bottom,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = -ROP_S(IJK,M)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE SOURCE_ROP_S

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: POINT_SOURCE_ROP_S                                      C
!  Purpose: Adds point sources to the solids continuity equation.      C
!                                                                      C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE POINT_SOURCE_ROP_S(B_M, M)

      use compar
      use constant
      use geometry
      use indices
      use param, only: dimension_m, dimension_3
      use param1, only: small_number
      use physprop
      use ps
      use run
      use functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Solids phase index.
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

               B_M(IJK,M) = B_M(IJK,M) - pSource
            endif
         enddo
         enddo
         enddo

      enddo PS_LP

      RETURN
      END SUBROUTINE POINT_SOURCE_ROP_S
