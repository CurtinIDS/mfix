!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DIFFUSE_MEAN_FIELDS                                     !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DIF_PHI_DES(M, DIF, A_M, B_M)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE toleranc
      USE run
      USE geometry
      USE compar
      USE sendrecv
      USE indices
      USE fun_avg
      USE functions
      USE cutcell

      IMPLICIT NONE

! Phase index
      INTEGER, INTENT(IN) :: M

!  Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)

! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)

      DOUBLE PRECISION, INTENT(IN) :: DIF(DIMENSION_3)

! Fluid Cell indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IMJK, IM, IJKW, IPJK, IJKE
      INTEGER :: IJMK, JM, IJKS, IJPK, IJKN
      INTEGER :: IJKM, KM, IJKB, IJKP, IJKT
!
! Diffusion parameter
      DOUBLE PRECISION :: D_f


!  Calculate convection-diffusion fluxes through each of the faces
      DO IJK = ijkstart3, ijkend3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF(.NOT.FLUID_AT(IJK)) CYCLE

         IPJK = IP_OF(IJK)
         IJPK = JP_OF(IJK)

         IJKE = EAST_OF(IJK)
         IJKN = NORTH_OF(IJK)

! East face (i+1/2, j, k)
         D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*AYZ(IJK)
         IF(CUT_TREATMENT_AT(IJK)) THEN
            IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IPJK))) THEN
               D_F = AVG_X_H(DIF(IJK),DIF(IJKE),I)*ODX_E(I)*DY(J)*DZ(K)
            ENDIF
         ENDIF

         A_M(IJK,east,M) = D_F
         A_M(IPJK,west,M) = D_F

! West face (i-1/2, j, k)
         IMJK = IM_OF(IJK)
         IF (.NOT.FLUID_AT(IMJK)) THEN
            IM = IM1(I)
            IJKW = WEST_OF(IJK)
            D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*AYZ(IMJK)
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IMJK))) THEN
                  D_F = AVG_X_H(DIF(IJKW),DIF(IJK),IM)*ODX_E(IM)*DY(J)*DZ(K)
               ENDIF
            ENDIF
            A_M(IJK,west,M) = D_F
         ENDIF


! North face (i, j+1/2, k)
         D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*AXZ(IJK)
         IF(CUT_TREATMENT_AT(IJK)) THEN
            IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IJPK))) THEN
               D_F = AVG_Y_H(DIF(IJK),DIF(IJKN),J)*ODY_N(J)*DX(I)*DZ(K)
            ENDIF
         ENDIF
         A_M(IJK,north,M) = D_F
         A_M(IJPK,south,M) = D_F

! South face (i, j-1/2, k)
         IJMK = JM_OF(IJK)
         IF (.NOT.FLUID_AT(IJMK)) THEN
            JM = JM1(J)
            IJKS = SOUTH_OF(IJK)
            D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*AXZ(IJMK)
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IJMK))) THEN
                  D_F = AVG_Y_H(DIF(IJKS),DIF(IJK),JM)*ODY_N(JM)*DX(I)*DZ(K)
               ENDIF
            ENDIF
            A_M(IJK,south,M) = D_F
         ENDIF


! Top face (i, j, k+1/2)
         IF (DO_K) THEN
            IJKP = KP_OF(IJK)
            IJKT = TOP_OF(IJK)
            D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*AXY(IJK)
            IF(CUT_TREATMENT_AT(IJK)) THEN
               IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IJKP))) THEN
                  D_F = AVG_Z_H(DIF(IJK),DIF(IJKT),K)*OX(I)*ODZ_T(K)*DX(I)*DY(J)
               ENDIF
            ENDIF
            A_M(IJK,top,M) = D_F
            A_M(IJKP,bottom,M) = D_F


! Bottom face (i, j, k-1/2)
            IJKM = KM_OF(IJK)
            IF (.NOT.FLUID_AT(IJKM)) THEN
               KM = KM1(K)
               IJKB = BOTTOM_OF(IJK)
               D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX(I)*ODZ_T(KM)*AXY(IJKM)
               IF(CUT_TREATMENT_AT(IJK)) THEN
                  IF(CUT_CELL_AT(IJK).AND.(.NOT.FLUID_AT(IJKM))) THEN
                     D_F = AVG_Z_H(DIF(IJKB),DIF(IJK),KM)*OX(I)*ODZ_T(KM)*DX(I)*DY(J)
                  ENDIF
               ENDIF
               A_M(IJK,bottom,M) = D_F
            ENDIF
         ENDIF
      END DO

      CALL DIF_PHI_IS(DIF, A_M, M)


      RETURN
      END SUBROUTINE DIF_PHI_DES
