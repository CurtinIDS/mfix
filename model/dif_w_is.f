!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DIF_W_IS                                                C
!  Purpose: Remove diffusive fluxes across internal surfaces.          C
!  (Make user defined internal surfaces non-conducting)                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DIF_W_IS(DIF, A_M, M)

! Modules
!---------------------------------------------------------------------//
      USE param
      USE geometry, only: odx_e, ody_n, axz_w, ayz_w

      USE is, only: is_defined, is_plane
      USE is, only: is_i_w, is_i_e, is_j_s, is_j_n, is_k_t, is_k_b

      USE fun_avg, only: avg_x_h, avg_y_h, avg_z_h

      USE functions, only: funijk, ip_of, jp_of
      USE functions, only: east_of, north_of, top_of
      USE functions, only: is_on_mype_plus2layers

      USE compar, only: dead_cell_at
      IMPLICIT NONE
! Dummy arguments
!---------------------------------------------------------------------//
! Gamma -- diffusion coefficient
      DOUBLE PRECISION, INTENT(IN) :: Dif(DIMENSION_3)
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Solids phase
      INTEGER, INTENT(IN) :: M

! Local variables
!---------------------------------------------------------------------//
! Diffusion parameter
      DOUBLE PRECISION :: D_f
! Internal surface
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK
      INTEGER :: IJKE, IJKN, IJKT, IJPK, IPJK, IJKTN, IJKTE
!---------------------------------------------------------------------//


      DO L = 1, DIMENSION_IS
         IF (IS_DEFINED(L)) THEN
            I1 = IS_I_W(L)
            I2 = IS_I_E(L)
            J1 = IS_J_S(L)
            J2 = IS_J_N(L)
            K1 = IS_K_B(L)
            K2 = IS_K_T(L)

            IF (IS_PLANE(L) == 'N') THEN
               DO K = K1, K2
               DO J = J1, J2
               DO I = I1, I2
                  IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IJKT = TOP_OF(IJK)
                  IJKN = NORTH_OF(IJK)
                  IJKTN = TOP_OF(IJKN)
                  IJPK = JP_OF(IJK)

                  D_F = AVG_Z_H(AVG_Y_H(DIF(IJK),DIF(IJKN),J),&
                                AVG_Y_H(DIF(IJKT),DIF(IJKTN),J),K)*&
                        ODY_N(J)*AXZ_W(IJK)

                  A_M(IJK,north,M) = A_M(IJK,north,M) - D_F
                  A_M(IJPK,south,M) = A_M(IJPK,south,M) - D_F
               ENDDO
               ENDDO
               ENDDO

            ELSEIF (IS_PLANE(L) == 'E') THEN
               DO K = K1, K2
               DO J = J1, J2
               DO I = I1, I2
                  IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IJKE = EAST_OF(IJK)
                  IJKT = TOP_OF(IJK)
                  IJKTE = EAST_OF(IJKT)
                  IPJK = IP_OF(IJK)

                  D_F = AVG_Z_H(AVG_X_H(DIF(IJK),DIF(IJKE),I),&
                                AVG_X_H(DIF(IJKT),DIF(IJKTE),I),K)*&
                        ODX_E(I)*AYZ_W(IJK)

                  A_M(IJK,east,M) = A_M(IJK,east,M) - D_F
                  A_M(IPJK,west,M) = A_M(IPJK,west,M) - D_F
               ENDDO
               ENDDO
               ENDDO
            ENDIF

         ENDIF   ! end if is_defined
      ENDDO   ! end do dimension_is
      RETURN
      END SUBROUTINE DIF_W_IS
