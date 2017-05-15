!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLID_DRAG                                              C
!  Purpose: Accounting for the equal and opposite drag force on a      C
!  continuous solids phase due to discrete solid particles by          C
!  introducing the drag as a source term. Face centered.               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SOLID_DRAG_U(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE physprop
      USE indices
      USE compar
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION :: B_M(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKE, I
! Solids phase indices
      INTEGER :: M
!-----------------------------------------------

! currently no difference between interpolated and non-interpolated
! implementation of solid-solid drag

      DO M = 1, MMAX
         DO IJK = IJKSTART3, IJKEND3
            IF(FLUID_AT(IJK)) THEN

               I = I_OF(IJK)
               IJKE = EAST_OF(IJK)

               A_M(IJK,0,M) = A_M(IJK,0,M) - VOL_U(IJK) *              &
                  AVG_X(SDRAG_AM(IJK,M), SDRAG_AM(IJKE,M), I)

               B_M(IJK,M) = B_M(IJK,M) - VOL_U(IJK) *                  &
                  AVG_X(SDRAG_BM(IJK,1,M), SDRAG_BM(IJKE,1,M), I)

               ENDIF   ! end if (fluid_at(ijk))
         ENDDO   ! end do (ijk=ijkstart3,ijkend3)
      ENDDO   ! end do (cm=1,mmax)


      RETURN
      END SUBROUTINE SOLID_DRAG_U

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLID_DRAG                                              C
!  Purpose: Accounting for the equal and opposite drag force on a      C
!  continuous solids phase due to discrete solid particles by          C
!  introducing the drag as a source term. Face centered.               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLID_DRAG_V(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE physprop
      USE indices
      USE compar
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION :: B_M(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKN, J
! Solids phase indices
      INTEGER :: M

      DO M = 1, MMAX
         DO IJK = IJKSTART3, IJKEND3
            IF(FLUID_AT(IJK)) THEN
                  J = J_OF(IJK)
                  IJKN = NORTH_OF(IJK)

                  A_M(IJK,0,M) = A_M(IJK,0,M) - VOL_V(IJK) *           &
                     AVG_Y(SDRAG_AM(IJK,M), SDRAG_AM(IJKN,M), J)

                  B_M(IJK,M) = B_M(IJK,M) - VOL_V(IJK) *               &
                     AVG_Y(SDRAG_BM(IJK,2,M), SDRAG_BM(IJKN,2,M), J)

            ENDIF   ! end if (fluid_at(ijk))
         ENDDO   ! end do (ijk=ijkstart3,ijkend3)
      ENDDO   ! end do (cm=1,mmax)


      RETURN
      END SUBROUTINE SOLID_DRAG_V


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLID_DRAG_W                                            C
!  Purpose: Accounting for the equal and opposite drag force on a      C
!  continuous solids phase due to discrete solid particles by          C
!  introducing the drag as a source term. Face centered.               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SOLID_DRAG_W(A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE physprop
      USE indices
      USE compar
      USE discretelement
      USE fun_avg
      USE functions

      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m
      DOUBLE PRECISION :: A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION :: B_M(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IJKT, K
! Solids phase indices
      INTEGER :: M


      DO M = 1, MMAX
         DO IJK = IJKSTART3, IJKEND3
            IF(FLUID_AT(IJK)) THEN
               IJKT = TOP_OF(IJK)
               K = K_OF(IJK)

               A_M(IJK,0,M) = A_M(IJK,0,M) + VOL_W(IJK) *              &
                  AVG_Z(SDRAG_AM(IJK,M), SDRAG_AM(IJKT,M), K)
               B_M(IJK,M) = B_M(IJK,M) + VOL_W(IJK) *                  &
                  AVG_Z(SDRAG_BM(IJK,3,M), SDRAG_BM(IJKT,3,M), K)

            ENDIF   ! end if (fluid_at(ijk))
         ENDDO   ! end do (ijk=ijkstart3,ijkend3)
      ENDDO   ! end do (cm=1,mmax)


      RETURN
      END SUBROUTINE SOLID_DRAG_W
