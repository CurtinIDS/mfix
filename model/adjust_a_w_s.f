!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: ADJUST_A_W_s                                            C
!  Purpose: Handle the special case of the center coefficient in       C
!           W_s momentum eq. becoming zero.                            C
!                                                                      C
!  Author: M. Syamlal                                 Date:  2-AUG-96  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ADJUST_A_W_S(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE physprop
      USE geometry
      USE run
      USE indices
      USE compar
      USE sendrecv
      USE fun_avg
      USE functions
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKT, IJKM
! Phase index
      INTEGER :: M
!-----------------------------------------------

      DO M = 1, MMAX
         IF (DRAG_TYPE_ENUM == GHD_2007 .AND. M /= MMAX) CYCLE
         IF (MOMENTUM_Z_EQ(M)) THEN

!!$omp  parallel do private(IJK,IJKT,IJKM)
            DO IJK = ijkstart3, ijkend3
               IF (ABS(A_M(IJK,0,M)) < SMALL_NUMBER) THEN
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  IF (B_M(IJK,M) < ZERO) THEN
                     IJKT = TOP_OF(IJK)
                     IF (ROP_S(IJKT,M)*AXY_W(IJK) > SMALL_NUMBER) THEN
                        B_M(IJK,M) = SQRT((-B_M(IJK,M)/(ROP_S(IJKT,M)*&
                           AVG_Z_T(ONE,ZERO)*AXY_W(IJK))))
                     ELSE
                        B_M(IJK,M) = ZERO
                     ENDIF
                  ELSE IF (B_M(IJK,M) > ZERO) THEN
                     IJKM = KM_OF(IJK)
                     IF (ROP_S(IJK,M)*AXY_W(IJKM) > SMALL_NUMBER) THEN
                        B_M(IJK,M) = SQRT(B_M(IJK,M)/(ROP_S(IJK,M)*&
                           AVG_Z_T(ZERO,ONE)*AXY_W(IJKM)))
                     ELSE
                        B_M(IJK,M) = ZERO
                     ENDIF
                  ENDIF
               ENDIF    ! end if (abs(a_m(ijk,0,m))<small_number)
            ENDDO    ! end do loop (ijk=ijkstart3,ijkend3)

         ENDIF   ! end if (momentum_z_eq(m))
      ENDDO   ! end do loop (m=1,mmax)

      RETURN
      END SUBROUTINE ADJUST_A_W_S
