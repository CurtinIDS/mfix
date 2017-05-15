!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: DIF_U_IS                                                C
!  Purpose: Remove diffusive fluxes across internal surfaces.          C
!  (Make user defined internal surfaces non-conducting)                C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DIF_U_IS(DIF, A_M, M)

! Modules
!---------------------------------------------------------------------//
      USE param
      USE geometry, only: do_k
      USE geometry, only: ody_n, ox_e, odz_t, axz_u, axy_u

      USE is, only: is_defined, is_plane
      USE is, only: is_i_w, is_i_e, is_j_s, is_j_n, is_k_t, is_k_b

      USE fun_avg, only: avg_x_h, avg_y_h, avg_z_h

      USE functions, only: funijk, jp_of, kp_of
      USE functions, only: east_of, north_of, top_of
      USE functions, only: is_on_mype_plus2layers

      USE compar, only: dead_cell_at
      USE compar, only: istart2, jstart2, kstart2
      USE compar, only: iend2, jend2, kend2
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
! Internal surface
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK
      INTEGER :: IJKE, IJKN, IJKT, IJPK, IJKP, IJKNE, IJKTE
! Diffusion parameter
      DOUBLE PRECISION :: D_f
!---------------------------------------------------------------------//

      DO L = 1, DIMENSION_IS
         IF (IS_DEFINED(L)) THEN
            I1 = IS_I_W(L)
            I2 = IS_I_E(L)
            J1 = IS_J_S(L)
            J2 = IS_J_N(L)
            K1 = IS_K_B(L)
            K2 = IS_K_T(L)

! Limit I1, I2 and all to local processor first ghost layer
            IF(I1.LE.IEND2)   I1 = MAX(I1, ISTART2)
            IF(J1.LE.JEND2)   J1 = MAX(J1, JSTART2)
            IF(K1.LE.KEND2)   K1 = MAX(K1, KSTART2)
            IF(I2.GE.ISTART2) I2 = MIN(I2, IEND2)
            IF(J2.GE.JSTART2) J2 = MIN(J2, JEND2)
            IF(K2.GE.KSTART2) K2 = MIN(K2, KEND2)

            IF (IS_PLANE(L) == 'N') THEN
               DO K = K1, K2
               DO J = J1, J2
               DO I = I1, I2
                  IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                  IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                  IJK = FUNIJK(I,J,K)
                  IJKE = EAST_OF(IJK)
                  IJKN = NORTH_OF(IJK)
                  IJKNE = EAST_OF(IJKN)
                  IJPK = JP_OF(IJK)

                  D_F = AVG_X_H(AVG_Y_H(DIF(IJK),DIF(IJKN),J),&
                                AVG_Y_H(DIF(IJKE),DIF(IJKNE),J),I)*&
                        ODY_N(J)*AXZ_U(IJK)

                  A_M(IJK,north,M) = A_M(IJK,north,M) - D_F
                  A_M(IJPK,south,M) = A_M(IJPK,south,M) - D_F
               ENDDO
               ENDDO
               ENDDO

            ELSEIF (IS_PLANE(L) == 'T') THEN
               IF (DO_K) THEN
                  DO K = K1, K2
                  DO J = J1, J2
                  DO I = I1, I2
                     IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                     IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                     IJK = FUNIJK(I,J,K)
                     IJKE = EAST_OF(IJK)
                     IJKT = TOP_OF(IJK)
                     IJKTE = EAST_OF(IJKT)
                     IJKP = KP_OF(IJK)

                     D_F = AVG_X_H(AVG_Z_H(DIF(IJK),DIF(IJKT),K),&
                                   AVG_Z_H(DIF(IJKE),DIF(IJKTE),K),I)*&
                           OX_E(I)*ODZ_T(K)*AXY_U(IJK)

                     A_M(IJK,top,M) = A_M(IJK,top,M) - D_F
                     A_M(IJKP,bottom,M) = A_M(IJKP,bottom,M) - D_F
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF
            ENDIF

         ENDIF   ! end if is_defined
      ENDDO   ! end do dimension_is


      RETURN
      END SUBROUTINE DIF_U_IS

