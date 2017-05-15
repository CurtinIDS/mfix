!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_PG_GRAD                                            !
!  Purpose: Calculate cell centered pressure force exerted on the      !
!           particles in the cell by the gas/fluid phase               !
!           Note that P_force is evaluated as -dp/dx                   !
!                                                                      !
!  Notes: This pressure force only needs to be calculated once during  !
!         the DEM loop (at the beginning) since the gas/fluid phase    !
!         is essentially static at that point (i.e., gas field is not  !
!         updated during DEM loop                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_GRAD_DES(PHI, DEL_PHI)

! Flag for cut-cells
      use cutcell, only: CARTESIAN_GRID

      use param, only: DIMENSION_3

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: DEL_PHI(3,DIMENSION_3)
!......................................................................!

      IF(CARTESIAN_GRID) THEN
         CALL CALC_GRAD_DES_CG(PHI, DEL_PHI)
      ELSE
         CALL CALC_GRAD_DES_STD(PHI, DEL_PHI)
      ENDIF

      RETURN
      END SUBROUTINE CALC_GRAD_DES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_GRAD_DES_STD                                       !
!                                                                      !
!  Purpose: Calculate cell centered gradient (DEL_PHI) of scalar PHI.  !
!           This routine calculates the average scalar value at the    !
!           cell faces to calculate the graident across the cell.      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_GRAD_DES_STD(PHI, DEL_PHI)

! Modules
!-----------------------------------------------

      use geometry, only: DO_K
      use geometry, only: oDX, oDY, oDZ
      USE geometry, only: IMIN1, JMIN1, KMIN1
      USE geometry, only: IMAX1, JMAX1, KMAX1
      USE indices, only: I_OF, J_OF, K_OF
      USE compar, only: IJKSTART3, IJKEND3

      use functions, only: FLUID_AT
      use functions, only: IM_OF, JM_OF, KM_OF
      use functions, only: IP_OF, JP_OF, KP_OF
      use fun_avg, only: AVG_X, AVG_Y, AVG_Z

      USE param, only: DIMENSION_3
      USE param1, only: ZERO

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: DEL_PHI(3,DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! general i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IPJK, IJPK, IJKP, IMJK, IJMK, IJKM
!......................................................................!


      DO IJK = IJKSTART3, IJKEND3
         DEL_PHI(:,IJK) = ZERO
         IF(.NOT.FLUID_AT(IJK)) CYCLE

         I = I_OF(IJK)
         IMJK = IM_OF(IJK)
         IPJK = IP_OF(IJK)

         IF(FLUID_AT(IMJK) .AND. FLUID_AT(IPJK)) THEN
            DEL_PHI(1,IJK) = oDX(I)*(AVG_X(PHI(IJK),PHI(IPJK),I) -     &
               AVG_X(PHI(IMJK),PHI(IJK),I-1))
         ELSEIF(FLUID_AT(IPJK)) THEN
            DEL_PHI(1,IJK) = 2.0d0*oDX(I) *                            &
               (AVG_X(PHI(IJK),PHI(IPJK),I) -  PHI(IJK))
         ELSEIF(FLUID_AT(IMJK)) THEN
            DEL_PHI(1,IJK) = 2.0d0*oDX(I) *                            &
               (PHI(IJK) - AVG_X(PHI(IMJK), PHI(IJK), I-1))
         ELSE
            DEL_PHI(1,IJK) = ZERO
         ENDIF


         J = J_OF(IJK)
         IJMK = JM_OF(IJK)
         IJPK = JP_OF(IJK)

         IF(FLUID_AT(IJMK) .AND. FLUID_AT(IJPK)) THEN
            DEL_PHI(2,IJK) = oDY(J)*(AVG_Y(PHI(IJK),PHI(IJPK),J) -     &
               AVG_Y(PHI(IJMK),PHI(IJK),J-1))
         ELSEIF(FLUID_AT(IJPK)) THEN
            DEL_PHI(2,IJK) = 2.0d0*oDY(J) *                            &
               (AVG_Y(PHI(IJK),PHI(IJPK),J) - PHI(IJK))
         ELSEIF(FLUID_AT(IJMK)) THEN
            DEL_PHI(2,IJK) = 2.0d0*oDY(J) *                            &
               (PHI(IJK)- AVG_Y(PHI(IJMK),PHI(IJK),J-1))
         ELSE
            DEL_PHI(2,IJK) = ZERO
         ENDIF

         IF(DO_K) THEN

            K = K_OF(IJK)
            IJKM = KM_OF(IJK)
            IJKP = KP_OF(IJK)

            IF(FLUID_AT(IJKM) .AND. FLUID_AT(IJKP)) THEN
               DEL_PHI(3,IJK) = oDZ(K)*(AVG_Z(PHI(IJK),PHI(IJKP),K) -  &
                  AVG_Z(PHI(IJKM),PHI(IJK),K-1))
            ELSEIF(FLUID_AT(IJKP)) THEN
               DEL_PHI(3,IJK) = 2.0d0*oDZ(K) *                         &
                  (AVG_Z(PHI(IJK),PHI(IJKP),K) - PHI(IJK))
            ELSEIF(FLUID_AT(IJKM)) THEN
               DEL_PHI(3,IJK) = 2.0d0*oDZ(K) *                         &
                  (PHI(IJK) - AVG_Z(PHI(IJKM),PHI(IJK),K-1))
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_GRAD_DES_STD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CALC_GRAD_DES_CG                                        !
!                                                                      !
!  Purpose: Calculate cell centered gradient (DEL_PHI) of scalar PHI.  !
!           This routine calculates the graidents at the cell face     !
!           and averages the face values to the cell center.           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CALC_GRAD_DES_CG(PHI, DEL_PHI)

! Modules
!---------------------------------------------------------------------//
      use geometry, only: DO_K
      use geometry, only: DX, DY, DZ
      use indices, only: I_OF, J_OF, K_OF
      use compar, only: IJKSTART3, IJKEND3

      use functions, only: FLUID_AT
      use functions, only: IM_OF, JM_OF, KM_OF
      use functions, only: IP_OF, JP_OF, KP_OF

      use param, only: DIMENSION_3
      use param1, only: ZERO

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)
      DOUBLE PRECISION, INTENT(OUT) :: DEL_PHI(3,DIMENSION_3)

! Local variables
!---------------------------------------------------------------------//
! General i, j, k indices
      INTEGER :: I, J, K, IJK
      INTEGER :: IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
! Counter
      DOUBLE PRECISION :: dLC
!......................................................................!


      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         DEL_PHI(:,IJK) = ZERO

         IF(.NOT.FLUID_AT(IJK)) CYCLE

         dLC = 0.0d0
         IJKE = IP_OF(IJK)
         IJKW = IM_OF(IJK)

         IF(FLUID_AT(IJKE)) THEN
            DEL_PHI(1,IJK) = DEL_PHI(1,IJK) +                          &
               2.0d0*(PHI(IJKE) - PHI(IJK))/(DX(I) + DX(I_OF(IJKE)))
            dLC = dLC + 1.0d0
         ENDIF
         IF(FLUID_AT(IJKW)) THEN
            DEL_PHI(1,IJK) = DEL_PHI(1,IJK) +                          &
               2.0d0*(PHI(IJK) - PHI(IJKW))/(DX(I) + DX(I_OF(IJKW)))
            dLC = dLC + 1.0d0
         ENDIF
         DEL_PHI(1,IJK) = DEL_PHI(1,IJK)/max(1.0d0,dLC)


         dLC = 0.0d0
         IJKN = JP_OF(IJK)
         IJKS = JM_OF(IJK)

         IF(FLUID_AT(IJKN)) THEN
            DEL_PHI(2,IJK) = DEL_PHI(2,IJK) +                          &
               2.0d0*(PHI(IJKN) - PHI(IJK))/(DY(J) + DY(J_OF(IJKN)))
            dLC = dLC + 1.0d0
         ENDIF

         IF(FLUID_AT(IJKS)) THEN
            DEL_PHI(2,IJK) = DEL_PHI(2,IJK) +                          &
               2.d0*(PHI(IJK) - PHI(IJKS))/(DY(J) + DY(J_OF(IJKS)))
            dLC = dLC + 1.0d0
         ENDIF
         DEL_PHI(2,IJK) = DEL_PHI(2,IJK)/max(1.0d0, dLC)

         IF(DO_K) THEN

            dLC = 0.0d0
            IJKT = KP_OF(IJK)
            IJKB = KM_OF(IJK)

            IF(FLUID_AT(IJKT)) THEN
               DEL_PHI(3,IJK) = DEL_PHI(3,IJK) +                       &
                  2.0d0*(PHI(IJKT) - PHI(IJK))/(DZ(K) + DZ(K_OF(IJKT)))
               dLC = dLC + 1.0d0
            ENDIF
            IF(FLUID_AT(IJKB)) THEN
               DEL_PHI(3,IJK) = DEL_PHI(3,IJK) +                       &
                  2.0d0*(PHI(IJK) - PHI(IJKB))/(DZ(K) + DZ(K_OF(IJKB)))
               dLC = dLC + 1.0d0
            ENDIF
            DEL_PHI(3,IJK) = DEL_PHI(3,IJK)/max(1.0d0, dLC)
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CALC_GRAD_DES_CG
