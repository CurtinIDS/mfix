!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: DIFFUSE_MEAN_FIELDS                                     !
!  Author: J.Musser                                   Date: 11-NOV-14  !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE DIF_PHI_BC_DES(PHI, M, A_M, B_M)

      USE param
      USE param1
      USE geometry
      USE indices
      USE bc
      USE compar
      USE cutcell, only : CARTESIAN_GRID, CG_SAFE_MODE
      USE fun_avg
      USE functions

      IMPLICIT NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------


      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)
! Phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Boundary condition index
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK, &
                 IM, JM, KM
!-----------------------------------------------

! Set up the default walls (i.e., bc_type='dummy' or undefined/default
! boundaries) as non-conducting...
! ---------------------------------------------------------------->>>
! when setting up default walls do not use cutcells to avoid conflict

      IF(.NOT.CARTESIAN_GRID) THEN

! west yz plane
         I1 = imin2
         DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(IP_OF(IJK),west,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,:,M) = ZERO
               A_M(IJK,east,M) = ONE
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
         ENDDO

! east yz plane
         I1 = IMAX2
         DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(IM_OF(IJK),east,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,:,M) = ZERO
               A_M(IJK,west,M) = ONE
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
         ENDDO

! south xz plane
         J1 = 1
         DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(JP_OF(IJK),south,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,:,M) = ZERO
               A_M(IJK,north,M) = ONE
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
         ENDDO

! north xz plane
         J1 = JMAX2
         DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(JM_OF(IJK),north,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,:,M) = ZERO
               A_M(IJK,south,M) = ONE
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
         ENDDO

         IF (DO_K) THEN
! bottom xy plane
            K1 = 1
            DO J1 = jmin3, jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)

               IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
                  A_M(KP_OF(IJK),bottom,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value (set
! the boundary cell value equal to adjacent fluid cell value)
                  A_M(IJK,:,M) = ZERO
                  A_M(IJK,top,M) = ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
            ENDDO

! top xy plane
            K1 = KMAX2
            DO J1 = jmin3, jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
                  A_M(KM_OF(IJK),top,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
                  A_M(IJK,:,M) = ZERO
                  A_M(IJK,bottom,M) = ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
            ENDDO
         ENDIF

      ENDIF

! Set user defined wall boundary conditions
! ---------------------------------------------------------------->>>
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN

            I1 = BC_I_W(L)
            I2 = BC_I_E(L)
            J1 = BC_J_S(L)
            J2 = BC_J_N(L)
            K1 = BC_K_B(L)
            K2 = BC_K_T(L)

            DO K = K1, K2
            DO J = J1, J2
            DO I = I1, I2

               IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells

               IJK = FUNIJK(I,J,K)


! Set the boundary cell equal to the known value in that
! cell
               A_M(IJK,:,M) = ZERO
               B_M(IJK,M) = ZERO
! second modify the matrix equation according to the user specified
! boundary condition
               IF (FLUID_AT(EAST_OF(IJK))) THEN
                  A_M(IJK,0,M) = -ODX_E(I)
                  A_M(IJK,east,M) =  ODX_E(I)

               ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                  IM = IM1(I)
                  A_M(IJK,west,M) =  ODX_E(IM)
                  A_M(IJK,0,M) = -ODX_E(IM)

               ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                  A_M(IJK,0,M) = -ODY_N(J)
                  A_M(IJK,north,M) =  ODY_N(J)

               ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                  JM = JM1(J)
                  A_M(IJK,south,M) =  ODY_N(JM)
                  A_M(IJK,0,M) = -ODY_N(JM)

               ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                  A_M(IJK,0,M) = -OX(I)*ODZ_T(K)
                  A_M(IJK,top,M) =  OX(I)*ODZ_T(K)

               ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                  KM = KM1(K)
                  A_M(IJK,bottom,M) =  OX(I)*ODZ_T(KM)
                  A_M(IJK,0,M) = -OX(I)*ODZ_T(KM)

               ENDIF
            ENDDO
            ENDDO
            ENDDO
         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc

! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID .AND. .NOT.(CG_SAFE_MODE(1)==1)) &
         CALL DIF_PHI_BC_DES_CG(PHI, 0, A_M, B_M)

      RETURN
      END SUBROUTINE DIF_PHI_BC_DES


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: BC_PHI_CG                                               C
!  Purpose: Modify boundary conditions for cartesian grid cut-cell     C
!           implementation                                             C
!                                                                      C
!  Author: Jeff Dietiker                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE DIF_PHI_BC_DES_CG(PHI, M, A_M, B_M)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE geometry
      USE indices
      USE bc
      USE compar
      USE cutcell
      USE fun_avg
      USE functions

      IMPLICIT NONE

!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! The field variable being solved for:
!     e.g., T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
      DOUBLE PRECISION, INTENT(IN) :: PHI(DIMENSION_3)

! Phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: I, J, K, IJK, IM, JM, KM

      LOGICAL :: ALONG_GLOBAL_GHOST_LAYER

!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         ALONG_GLOBAL_GHOST_LAYER = (I < IMIN1) .OR. ( I>IMAX1 )       &
            .OR. (J < JMIN1) .OR. (J > JMAX1)

         IF(DO_K) ALONG_GLOBAL_GHOST_LAYER = ALONG_GLOBAL_GHOST_LAYER  &
            .OR. (K < KMIN1) .OR. (K > KMAX1)

! Setting the value in the boundary cell equal to what is known
         IF(BLOCKED_CELL_AT(IJK)) THEN
            A_M(IJK,:,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = -PHI(IJK)
         ENDIF

         IF(BLOCKED_CELL_AT(IJK).OR.ALONG_GLOBAL_GHOST_LAYER) THEN


            IF (CUT_CELL_AT(IP_OF(IJK))) THEN
               IF( BC_ID(IP_OF(IJK)) > 0 ) THEN
                  A_M(IJK,0,M) = -ODX_E(I)
                  A_M(IJK,east,M) =  ODX_E(I)
                  B_M(IJK,M) =  ZERO
               ENDIF

            ELSEIF (CUT_CELL_AT(IM_OF(IJK))) THEN
               IF(BC_ID(IM_OF(IJK)) > 0 ) THEN
                  IM = IM1(I)
                  A_M(IJK,west,M) =  ODX_E(IM)
                  A_M(IJK,0,M) = -ODX_E(IM)
                  B_M(IJK,M) = ZERO
               ENDIF

            ELSEIF (CUT_CELL_AT(JP_OF(IJK))) THEN
               IF(BC_ID(JP_OF(IJK)) > 0 ) THEN
                  A_M(IJK,0,M) = -ODY_N(J)
                  A_M(IJK,north,M) =  ODY_N(J)
                  B_M(IJK,M) = ZERO
               ENDIF

            ELSEIF (CUT_CELL_AT(JM_OF(IJK))) THEN
               IF(BC_ID(JM_OF(IJK)) > 0 ) THEN
                  JM = JM1(J)
                  A_M(IJK,south,M) =  ODY_N(JM)
                  A_M(IJK,0,M) = -ODY_N(JM)
                  B_M(IJK,M) = ZERO
               ENDIF

            ELSEIF (CUT_CELL_AT(KP_OF(IJK))) THEN
               IF(BC_ID(KP_OF(IJK)) > 0 ) THEN
                  A_M(IJK,0,M) = -OX(I)*ODZ_T(K)
                  A_M(IJK,top,M) =  OX(I)*ODZ_T(K)
                  B_M(IJK,M)  = ZERO
               ENDIF

            ELSEIF (CUT_CELL_AT(KM_OF(IJK))) THEN
               IF(BC_ID(KM_OF(IJK)) > 0 ) THEN
                  KM = KM1(K)
                  A_M(IJK,bottom,M) =  OX(I)*ODZ_T(KM)
                  A_M(IJK,0,M) = -OX(I)*ODZ_T(KM)
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      RETURN

      END SUBROUTINE DIF_PHI_BC_DES_CG


