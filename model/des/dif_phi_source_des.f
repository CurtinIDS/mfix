!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SOURCE_phi                                              !
!  Author: M. Syamlal                                 Date: 30-APR-97  !
!                                                                      !
!  Purpose: Determine source terms for phi eq. The terms               !
!     appear in the center coefficient and RHS vector. The center      !
!     coefficient and source vector are negative. The off-diagonal     !
!     coefficients are positive. S_p must be positive.                 !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DIF_PHI_SOURCE_DES(PHI, M, A_M, B_M, lDT)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE scales
      USE physprop
      USE fldvar
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE compar
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Source term on RHS
      DOUBLE PRECISION :: S_C
! Phi
      DOUBLE PRECISION, INTENT(IN) :: Phi(DIMENSION_3)
! phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
! Diffusion time-step
      DOUBLE PRECISION, INTENT(IN) :: lDT
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK

      DOUBLE PRECISION :: APO
      DOUBLE PRECISION :: lOoDT
!-----------------------------------------------


      lOoDT = 1.0d0/lDT

      DO IJK = IJKSTART3, IJKEND3

         IF (FLUID_AT(IJK)) THEN

            APO = VOL(IJK)*lOoDT

! Collect the terms
            A_M(IJK,0,M) = -(A_M(IJK,east,M)+A_M(IJK,west,M)+&
                             A_M(IJK,north,M)+A_M(IJK,south,M)+&
                             A_M(IJK,top,M)+A_M(IJK,bottom,M)+ APO)

            S_C = APO*PHI(IJK)

            IF(B_M(IJK,M) < S_C .OR. PHI(IJK) == ZERO) THEN
               B_M(IJK,M) = B_M(IJK,M) - S_C
            ELSE
               A_M(IJK,0,M) = A_M(IJK,0,M) - B_M(IJK,M)/PHI(IJK)
               B_M(IJK,M) = -S_C
            ENDIF

         ENDIF
      ENDDO


      RETURN
      END SUBROUTINE DIF_PHI_SOURCE_DES
