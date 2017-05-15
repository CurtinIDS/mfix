!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG_W                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_U(A_M, B_M, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use particle_filter, only: DES_INTERP_SCHEME_ENUM, DES_INTERP_GARG
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Coefficient at cell corners added to the gas momentum A matix.
      use discretelement, only: DRAG_AM, F_GDS
! Coefficient at cell corners added to gas momentum B vector.
      use discretelement, only: DRAG_BM
! Volume of X-momentum cell
      use geometry, only: VOL_U
! Flag to calculate Z direction
      use geometry, only: DO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and solids phase..
      use param, only: DIMENSION_3, DIMENSION_M
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! The I, J, and K values that comprise an IJK
      use indices, only: I_OF, J_OF, K_OF
! Flag: Flow at East face of fluid cell
      use functions, only: FLOW_AT_E
! IJK of cell to east.
      use functions, only: EAST_OF
! IJK function for I,J,K that includes mapped indices.
      use compar, only: FUNIJK_MAP_C
! Function for averaging to a scalar cell's east face.
      use fun_avg, only: AVG_X
! Domain index bounds.
      use compar, only: ISTART2, JSTART2, KSTART2
      use compar, only: IEND2, JEND2, KEND2

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3,-3:3,0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER


! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K, IJK, IJMK, IJKM, IJMKM, IJKE
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION :: tmp_A, tmp_B
! Averaging factor
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

! Average the interpoalted drag force from the cell corners to the cell face.
      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG)THEN

         AVG_FACTOR = merge(0.25d0, 0.5d0, DO_K)

!!$omp parallel do schedule(guided,50) default(none)                    &
!!$omp shared(IJKSTART3, IJKEND3, FUNIJK_MAP_C, I_OF, J_OF,             &
!!$omp    K_OF, DO_K, AVG_FACTOR, DRAG_AM, DRAG_BM, A_M, B_M, VOL_U)    &
!!$omp private(IJK, I, J, K, IJMK, IJKM, IJMKM, tmp_A, tmp_B)
         DO IJK = IJKSTART3, IJKEND3

            IF(FLOW_AT_E(IJK)) THEN

               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)

               IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
               IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
               IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

               IJMK = FUNIJK_MAP_C(I, J-1, K)

               tmp_A = -AVG_FACTOR*(DRAG_AM(IJK) + DRAG_AM(IJMK))
               tmp_B = -AVG_FACTOR*(DRAG_BM(IJK,1) + DRAG_BM(IJMK,1))

               IF(DO_K) THEN
                  IJKM = FUNIJK_MAP_C(I, J, K-1)
                  IJMKM = FUNIJK_MAP_C(I, J-1, K-1)
                  tmp_A = tmp_A - AVG_FACTOR*                          &
                     (DRAG_AM(IJKM) + DRAG_AM(IJMKM))
                  tmp_B = tmp_B - AVG_FACTOR*                          &
                     (DRAG_BM(IJKM,1) + DRAG_BM(IJMKM,1))
               ENDIF

               A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A*VOL_U(IJK)
               B_M(IJK,0) = B_M(IJK,0) + tmp_B*VOL_U(IJK)
            ENDIF

         ENDDO
!!$omp end parallel do

      ELSE

!$omp parallel do default(none) &
!$omp private(IJK, I, IJKE, tmp_A, tmp_B) &
!$omp shared(IJKSTART3,IJKEND3,I_OF, F_GDS, DRAG_BM, A_M, B_M, VOL_U)
         DO IJK = IJKSTART3, IJKEND3
            IF(FLOW_AT_E(IJK)) THEN
               I = I_OF(IJK)
               IJKE = EAST_OF(IJK)

               tmp_A = AVG_X(F_GDS(IJK), F_GDS(IJKE), I)
               tmp_B = AVG_X(DRAG_BM(IJK,1), DRAG_BM(IJKE,1), I)

               A_M(IJK,0,0) = A_M(IJK,0,0) - VOL_U(IJK) * tmp_A
               B_M(IJK,0) = B_M(IJK,0) - VOL_U(IJK) * tmp_B
            ENDIF
         ENDDO
!$omp end parallel do
      ENDIF

      END SUBROUTINE GAS_DRAG_U

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG_V                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_V(A_M, B_M, IER)


! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use particle_filter, only: DES_INTERP_SCHEME_ENUM, DES_INTERP_GARG
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Coefficient at cell corners added to the gas momentum A matix.
      use discretelement, only: DRAG_AM, F_GDS
! Coefficient at cell corners added to gas momentum B vector.
      use discretelement, only: DRAG_BM
! Volume of Y-momentum cell
      use geometry, only: VOL_V
! Flag to calculate Z direction
      use geometry, only: DO_K

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and solids phase..
      use param, only: DIMENSION_3, DIMENSION_M
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! The I, J, and K values that comprise an IJK
      use indices, only: I_OF, J_OF, K_OF
! Flag: Flow at North face of cell
      use functions, only: FLOW_AT_N
! IJK of cell to north.
      use functions, only: NORTH_OF
! IJK function for I,J,K that includes mapped indices.
      use compar, only: FUNIJK_MAP_C
! Function for averaging to a scalar cell's north face.
      use fun_avg, only: AVG_Y
! Domain index bounds.
      use compar, only: ISTART2, JSTART2, KSTART2
      use compar, only: IEND2, JEND2, KEND2


      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K, IJK, IMJK, IJKM, IMJKM, IJKN
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG)THEN

         AVG_FACTOR = merge(0.25d0, 0.5d0, DO_K)

!!$omp parallel do schedule (guided,50) default(none)                   &
!!$omp shared(IJKSTART3, IJKEND3, FUNIJK_MAP_C, I_OF, J_OF,             &
!!$omp    K_OF, DO_K, AVG_FACTOR, DRAG_AM, DRAG_BM, A_M, B_M, VOL_V)    &
!!$omp private(IJK, I, J, K, IMJK, IJKM, IMJKM, tmp_A, tmp_B)
         DO IJK = IJKSTART3, IJKEND3
            IF(FLOW_AT_N(IJK)) THEN

               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)

               IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
               IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
               IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

               IMJK = FUNIJK_MAP_C(I-1,J,K)

               tmp_A = -AVG_FACTOR*(DRAG_AM(IJK) + DRAG_AM(IMJK))
               tmp_B = -AVG_FACTOR*(DRAG_BM(IJK,2) + DRAG_BM(IMJK,2))

               IF(DO_K) THEN

                  IJKM = FUNIJK_MAP_C(I,J,K-1)
                  IMJKM = FUNIJK_MAP_C(I-1,J,K-1)

                  tmp_A = tmp_A - AVG_FACTOR*                          &
                     (DRAG_AM(IJKM) + DRAG_AM(IMJKM))
                  tmp_B = tmp_B - AVG_FACTOR*                          &
                     (DRAG_BM(IJKM,2) + DRAG_BM(IMJKM,2))
               ENDIF

               A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A*VOL_V(IJK)
               B_M(IJK,0) = B_M(IJK,0) + tmp_B*VOL_V(IJK)

            ENDIF
         ENDDO
!!$omp end parallel do

      ELSE

!$omp parallel do default(none) &
!$omp private(IJK, J, IJKN, tmp_a, tmp_B) &
!$omp shared(IJKSTART3, IJKEND3, J_OF, F_GDS, DRAG_AM, DRAG_BM, A_M, B_M, VOL_V)
         DO IJK = IJKSTART3, IJKEND3
            IF(FLOW_AT_N(IJK)) THEN
               J = J_OF(IJK)
               IJKN = NORTH_OF(IJK)

               tmp_A = AVG_Y(F_GDS(IJK), F_GDS(IJKN), J)
               tmp_B = AVG_Y(DRAG_BM(IJK,2), DRAG_BM(IJKN,2), J)

               A_M(IJK,0,0) = A_M(IJK,0,0) - VOL_V(IJK) * tmp_A
               B_M(IJK,0) = B_M(IJK,0) - VOL_V(IJK) * tmp_B
            ENDIF
         ENDDO
!$omp end parallel do
      ENDIF

      END SUBROUTINE GAS_DRAG_V

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: GAS_DRAG_W                                              !
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  !
!                                                                      !
!  Purpose: Account for the equal and opposite drag force on the gas   !
!           phase due to particles by introducing the drag as a        !
!           source term.  Face centered.                               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GAS_DRAG_W(A_M, B_M, IER)

! Global Variables:
!---------------------------------------------------------------------//
! Runtime Flag: Interpolate DES values
      use particle_filter, only: DES_INTERP_SCHEME_ENUM, DES_INTERP_GARG
! Flag: Gas sees the effect of particles in gas/solids flows.
      use discretelement, only: DES_ONEWAY_COUPLED
! Coefficient at cell corners added to the gas momentum A matix.
      use discretelement, only: DRAG_AM, F_GDS
! Coefficient at cell corners added to gas momentum B vector.
      use discretelement, only: DRAG_BM
! Volume of Z-momentum cell
      use geometry, only: VOL_W

! Global Parameters:
!---------------------------------------------------------------------//
! Size of IJK arrays and solids phase..
      use param, only: DIMENSION_3, DIMENSION_M
! Fluid grid loop bounds.
      use compar, only: IJKStart3, IJKEnd3
! The I, J, and K values that comprise an IJK
      use indices, only: I_OF, J_OF, K_OF
! Flag: Flow at top face of fluid cell
      use functions, only: FLOW_AT_T
! IJK of cell to top.
      use functions, only: TOP_OF
! IJK function for I,J,K that includes mapped indices.
      use compar, only: FUNIJK_MAP_C
! Function for averaging to a scalar cell's north face.
      use fun_avg, only: AVG_Z
! Domain index bounds.
      use compar, only: ISTART2, JSTART2, KSTART2
      use compar, only: IEND2, JEND2, KEND2

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_M(DIMENSION_3,-3:3,0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_M(DIMENSION_3, 0:DIMENSION_M)
! Error index
      INTEGER, INTENT(INOUT) :: IER

! Local Variables:
!---------------------------------------------------------------------//
! Grid cell indices
      INTEGER :: I, J, K, IJK, IMJK, IJMK, IMJMK, IJKT
! temporary variables for matrix A_M and vector B_M
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor
! (=0.25 in 3D and =0.5 in 2D)
      DOUBLE PRECISION :: AVG_FACTOR
!......................................................................!

! Initialize error flag.
      IER = 0

! Skip this routine if the gas/solids are only one-way coupled.
      IF(DES_ONEWAY_COUPLED) RETURN

      IF(DES_INTERP_SCHEME_ENUM == DES_INTERP_GARG)THEN

         AVG_FACTOR = 0.25d0

!!$omp parallel do schedule (guided,50) default(none)                   &
!!$omp shared(IJKSTART3, IJKEND3, FUNIJK_MAP_C, I_OF, J_OF,             &
!!$omp    K_OF, AVG_FACTOR, DRAG_AM, DRAG_BM, A_M, B_M, VOL_W)          &
!!$omp private(IJK, I, J, K, IMJK, IJMK, IMJMK, tmp_A, tmp_B)
         DO IJK = IJKSTART3, IJKEND3
            IF(FLOW_AT_T(IJK)) THEN

               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)

               IF (I.LT.ISTART2 .OR. I.GT.IEND2) CYCLE
               IF (J.LT.JSTART2 .OR. J.GT.JEND2) CYCLE
               IF (K.LT.KSTART2 .OR. K.GT.KEND2) CYCLE

               IMJK = FUNIJK_MAP_C(I-1,J,K)
               IJMK = FUNIJK_MAP_C(I,J-1,K)
               IMJMK = FUNIJK_MAP_C(I-1,J-1,K)

               tmp_A = -AVG_FACTOR*(DRAG_AM(IJK) + DRAG_AM(IMJK) +     &
                  DRAG_AM(IJMK) + DRAG_AM(IMJMK))

               tmp_B = -AVG_FACTOR*(DRAG_BM(IJK,3) + DRAG_BM(IMJK,3) + &
                  DRAG_BM(IJMK,3) + DRAG_BM(IMJMK,3))

               A_M(IJK,0,0) = A_M(IJK,0,0) + tmp_A*VOL_W(IJK)
               B_M(IJK,0) = B_M(IJK,0) + tmp_B*VOL_W(IJK)

            ENDIF
         ENDDO
!!$omp end parallel do

      ELSE

!$omp parallel do default(none) &
!$omp private(IJK, K, IJKT, tmp_A, tmp_B) &
!$omp shared(IJKSTART3,IJKEND3,K_OF, F_GDS, DRAG_BM, A_M, B_M, VOL_W)
         DO IJK = IJKSTART3, IJKEND3
            IF(FLOW_AT_T(IJK)) THEN
               K = K_OF(IJK)
               IJKT = TOP_OF(IJK)

               tmp_A = AVG_Z(F_GDS(IJK), F_GDS(IJKT), K)
               tmp_B = AVG_Z(DRAG_BM(IJK,3), DRAG_BM(IJKT,3), K)

               A_M(IJK,0,0) = A_M(IJK,0,0) - VOL_W(IJK) * tmp_A
               B_M(IJK,0) = B_M(IJK,0) - VOL_W(IJK) * tmp_B
            ENDIF
         ENDDO
!$omp end parallel do

      ENDIF

      RETURN
      END SUBROUTINE GAS_DRAG_W
