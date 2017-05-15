!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: BC_phi                                                  C
!  Purpose: Set up the phi boundary conditions                         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 30-APR-97  C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE BC_PHI(VAR, BC_PHIF, BC_PHIW, BC_HW_PHI, &
                        BC_C_PHI, M, A_M, B_M)


! Modules
!--------------------------------------------------------------------//
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

! Dummy arguments
!--------------------------------------------------------------------//
! The field variable being solved for:
!     e.g., T_g, T_s, x_g, x_s, Theta_m, scalar, K_Turb_G,
!     e_Turb_G
      DOUBLE PRECISION, INTENT(IN) :: VAR(DIMENSION_3)
! Boundary conditions specifications
! bc_phif = flow boundary value
! bc_phiw = wall boundary value
! bc_hw_phi = transfer coefficient
!      = 0 value means specified flux (neumann type)
!      = undefined value means specified wall value (dirichlet type)
!      = other value means mixed type
! bc_C_phi = transfer flux
      DOUBLE PRECISION, INTENT(IN) :: BC_phif(DIMENSION_BC), &
                                      BC_Phiw(DIMENSION_BC), &
                                      BC_hw_Phi(DIMENSION_BC), &
                                      BC_C_Phi(DIMENSION_BC)
! Phase index
      INTEGER, INTENT(IN) :: M
! Septadiagonal matrix A_m
      DOUBLE PRECISION, INTENT(INOUT) :: A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m
      DOUBLE PRECISION, INTENT(INOUT) :: B_m(DIMENSION_3, 0:DIMENSION_M)

! Local variables
!--------------------------------------------------------------------//
! Boundary condition index
      INTEGER :: L
! Indices
      INTEGER :: I, J, K, I1, I2, J1, J2, K1, K2, IJK, &
                 IM, JM, KM
!--------------------------------------------------------------------//

! Set up the default walls (i.e., bc_type='dummy' or undefined/default
! boundaries) as non-conducting...
! ---------------------------------------------------------------->>>
      IF(.NOT.CARTESIAN_GRID) THEN
! when setting up default walls do not use cutcells to avoid conflict

      IF (DO_K) THEN
! bottom xy plane
         K1 = 1
!!$omp    parallel do private(IJK, J1, I1)
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
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ONE
                  A_M(IJK,bottom,M) = ZERO
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO

! top xy plane
         K1 = KMAX2
!!$omp    parallel do private(IJK, J1, I1)
         DO J1 = jmin3, jmax3
            DO I1 = imin3, imax3
               IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
               IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
               IJK = FUNIJK(I1,J1,K1)
               IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
                  A_M(KM_OF(IJK),top,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
                  A_M(IJK,east,M) = ZERO
                  A_M(IJK,west,M) = ZERO
                  A_M(IJK,north,M) = ZERO
                  A_M(IJK,south,M) = ZERO
                  A_M(IJK,top,M) = ZERO
                  A_M(IJK,bottom,M) = ONE
                  A_M(IJK,0,M) = -ONE
                  B_M(IJK,M) = ZERO
               ENDIF
            ENDDO
         ENDDO
      ENDIF

! south xz plane
      J1 = 1
!!$omp    parallel do private(IJK, K1, I1)
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(JP_OF(IJK),south,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ONE
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

! north xz plane
      J1 = JMAX2
!!$omp    parallel do private(IJK, K1, I1)
      DO K1 = kmin3, kmax3
         DO I1 = imin3, imax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(JM_OF(IJK),north,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ONE
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

! west yz plane
      I1 = imin2
!!$omp    parallel do private(IJK, K1, J1)
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(IP_OF(IJK),west,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,east,M) = ONE
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

! east yz plane
      I1 = IMAX2
!!$omp    parallel do private(IJK, K1, J1)
      DO K1 = kmin3, kmax3
         DO J1 = jmin3, jmax3
            IF (.NOT.IS_ON_myPE_plus2layers(I1,J1,K1)) CYCLE
            IF (DEAD_CELL_AT(I1,J1,K1)) CYCLE  ! skip dead cells
            IJK = FUNIJK(I1,J1,K1)
            IF (DEFAULT_WALL_AT(IJK)) THEN
! Cutting the neighbor link between fluid cell and wall cell
               A_M(IM_OF(IJK),east,M) = ZERO
! Setting the wall value equal to the adjacent fluid cell value
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ONE
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = ZERO
            ENDIF
         ENDDO
      ENDDO

      ENDIF   !(.NOT.CARTESIAN_GRID)

! End setting the default boundary conditions
! ----------------------------------------------------------------<<<


! Set user defined wall boundary conditions
! ---------------------------------------------------------------->>>
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN
            IF (BC_TYPE_ENUM(L)==NO_SLIP_WALL .OR. &
                BC_TYPE_ENUM(L)==FREE_SLIP_WALL .OR. &
                BC_TYPE_ENUM(L)==PAR_SLIP_WALL) THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
!!$omp    parallel do private(IJK, K, J, I, IM, JM, KM)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
                        IM = IM1(I)
                        JM = JM1(J)
                        KM = KM1(K)
! first set the boundary cell value equal to the known value in that
! cell
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = VAR(IJK)
! second modify the matrix equation according to the user specified
! boundary condition
                        IF (FLUID_AT(EAST_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
! specified wall value (i.e., dirichlet type boundary)
                              A_M(IJK,east,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
! if bc_hw__phi=0 then specified flux boundary (i.e., neumann type
! boundary) otherwise a mixed type boundary
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODX_E(I))
                              A_M(IJK,east,M) = -(HALF*BC_HW_PHI(L)-ODX_E(I))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))
                           ENDIF
                        ELSEIF (FLUID_AT(WEST_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,west,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
                              A_M(IJK,west,M) = -(HALF*BC_HW_PHI(L)-ODX_E(IM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODX_E(IM))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))
                           ENDIF
                        ELSEIF (FLUID_AT(NORTH_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,north,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODY_N(J))
                              A_M(IJK,north,M) = -(HALF*BC_HW_PHI(L)-ODY_N(J))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))
                           ENDIF
                        ELSEIF (FLUID_AT(SOUTH_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,south,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
                              A_M(IJK,south,M) = -(HALF*BC_HW_PHI(L)-ODY_N(JM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODY_N(JM))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))
                           ENDIF
                        ELSEIF (FLUID_AT(TOP_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,top,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
                              A_M(IJK,0,M)=-(HALF*BC_HW_PHI(L)+OX(I)*ODZ_T(K))
                              A_M(IJK,top,M)=-(HALF*BC_HW_PHI(L)-OX(I)*ODZ_T(K))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))
                           ENDIF
                        ELSEIF (FLUID_AT(BOTTOM_OF(IJK))) THEN
                           IF (BC_HW_PHI(L) == UNDEFINED) THEN
                              A_M(IJK,bottom,M) = -HALF
                              A_M(IJK,0,M) = -HALF
                              B_M(IJK,M) = -BC_PHIW(L)
                           ELSE
                              A_M(IJK,bottom,M) = -(HALF*BC_HW_PHI(L)-&
                                               OX(I)*ODZ_T(KM))
                              A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+&
                                               OX(I)*ODZ_T(KM))
                              B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+&
                                             BC_C_PHI(L))
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
           ENDIF   ! end if (ns, fs, psw)
         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc
! end setting of wall boundary conditions
! ----------------------------------------------------------------<<<


! Set user defined boundary conditions for non-wall cells
! Setting p_inflow, p_outflow, mass_outflow or outflow flow boundary
! conditions
! ---------------------------------------------------------------->>>
      DO L = 1, DIMENSION_BC
         IF (BC_DEFINED(L)) THEN
            IF (BC_TYPE_ENUM(L)==P_OUTFLOW .OR. &
                BC_TYPE_ENUM(L)==MASS_OUTFLOW .OR. &
                BC_TYPE_ENUM(L)==OUTFLOW) THEN
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
!!$omp    parallel do private(IJK, K, J, I)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                       IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                       IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
! first set the flow boundary cell value equal to zero
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
                        B_M(IJK,M) = ZERO
! now set the flow boundary cell value equal to the adjacent fluid
! cell value
                        SELECT CASE (TRIM(BC_PLANE(L)))
                        CASE ('E')
! fluid cell on the east side
                           A_M(IJK,east,M) = ONE
                        CASE ('W')
! fluid cell on the west side
                           A_M(IJK,west,M) = ONE
                        CASE ('N')
                           A_M(IJK,north,M) = ONE
                        CASE ('S')
                           A_M(IJK,south,M) = ONE
                        CASE ('T')
                           A_M(IJK,top,M) = ONE
                        CASE ('B')
                           A_M(IJK,bottom,M) = ONE
                        END SELECT
                     ENDDO
                  ENDDO
               ENDDO
! end setting p_outflow, mass_outflow or outflow flow boundary
! conditions
! ----------------------------------------------------------------<<<

            ELSEIF(BC_TYPE_ENUM(L)==P_INFLOW .OR. &
                   BC_TYPE_ENUM(L)==MASS_INFLOW) THEN

! Setting bc that are defined but not nsw, fsw, psw, p_outflow,
! mass_outflow or outflow (at this time, this section addresses
! p_inflow and mass_inflow type boundaries)
! ----------------------------------------------------------------<<<
               I1 = BC_I_W(L)
               I2 = BC_I_E(L)
               J1 = BC_J_S(L)
               J2 = BC_J_N(L)
               K1 = BC_K_B(L)
               K2 = BC_K_T(L)
!!$omp    parallel do private(IJK, K, J, I)
               DO K = K1, K2
                  DO J = J1, J2
                     DO I = I1, I2
                        IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                        IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
                        IJK = FUNIJK(I,J,K)
! setting the value in the boundary cell equal to what is known
                        A_M(IJK,east,M) = ZERO
                        A_M(IJK,west,M) = ZERO
                        A_M(IJK,north,M) = ZERO
                        A_M(IJK,south,M) = ZERO
                        A_M(IJK,top,M) = ZERO
                        A_M(IJK,bottom,M) = ZERO
                        A_M(IJK,0,M) = -ONE
!                        B_M(IJK,M) = -BC_PHIF(L)  !does not allow the profile to be changed, e.g., from usr1
                        B_M(IJK,M) = -VAR(IJK)
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF   ! end if/else (bc_type)
! end setting of p_inflow or mass_inflow boundary conditions
! ----------------------------------------------------------------<<<

         ENDIF   ! end if (bc_defined)
      ENDDO   ! end L do loop over dimension_bc



! modifications for cartesian grid implementation
      IF(CARTESIAN_GRID .AND. .NOT.(CG_SAFE_MODE(1)==1)) &
         CALL BC_PHI_CG(VAR, BC_PHIF, BC_PHIW, BC_HW_PHI, &
                        BC_C_PHI, M, A_M, B_M)

      RETURN
      END SUBROUTINE BC_PHI


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: BC_PHI_CG                                               C
!  Purpose: Modify boundary conditions for cartesian grid cut-cell     C
!           implementation                                             C
!                                                                      C
!  Author: Jeff Dietiker                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE BC_PHI_CG(VAR, BC_PHIF, BC_PHIW, BC_HW_PHI, &
                        BC_C_PHI, M, A_M, B_M)

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
      DOUBLE PRECISION, INTENT(IN) :: VAR(DIMENSION_3)
! Boundary conditions specifications
! bc_phif = flow boundary value
! bc_phiw = wall boundary value
! bc_hw_phi = transfer coefficient
!      = 0 value means specified flux (neumann type)
!      = undefined value means specified wall value (dirichlet type)
!      = other value means mixed type
! bc_C_phi = transfer flux
      DOUBLE PRECISION, INTENT(IN) :: BC_phif(DIMENSION_BC), &
                                      BC_Phiw(DIMENSION_BC), &
                                      BC_hw_Phi(DIMENSION_BC), &
                                      BC_C_Phi(DIMENSION_BC)
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
! Boundary condition index
      INTEGER :: L
! Boundary identifiers
      INTEGER :: BCV
      INTEGER :: BCT

      LOGICAL :: ALONG_GLOBAL_GHOST_LAYER

!-----------------------------------------------

      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         ALONG_GLOBAL_GHOST_LAYER = (I<IMIN1).OR.(I>IMAX1).OR.(J<JMIN1).OR.(J>JMAX1)

         IF(DO_K) ALONG_GLOBAL_GHOST_LAYER = ALONG_GLOBAL_GHOST_LAYER.OR.(K<KMIN1).OR.(K>KMAX1)

         IF(BLOCKED_CELL_AT(IJK)) THEN
! setting the value in the boundary cell equal to what is known
            A_M(IJK,east,M) = ZERO
            A_M(IJK,west,M) = ZERO
            A_M(IJK,north,M) = ZERO
            A_M(IJK,south,M) = ZERO
            A_M(IJK,top,M) = ZERO
            A_M(IJK,bottom,M) = ZERO
            A_M(IJK,0,M) = -ONE
            B_M(IJK,M) = -VAR(IJK)
         ENDIF


         IF(BLOCKED_CELL_AT(IJK).OR.ALONG_GLOBAL_GHOST_LAYER.AND.(.NOT.FLOW_AT(IJK))) THEN
            IM = IM1(I)
            JM = JM1(J)
            KM = KM1(K)

            IF (CUT_CELL_AT(IP_OF(IJK))) THEN
               BCV = BC_ID(IP_OF(IJK))
               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               IF (BCT==CG_NSW.OR.BCT==CG_FSW.OR.BCT==CG_PSW) THEN
                  L = BCV
                  IF (BC_HW_PHI(L) == UNDEFINED) THEN
                     A_M(IJK,east,M) = -HALF
                     A_M(IJK,0,M) = -HALF
                     B_M(IJK,M) = -BC_PHIW(L)
                  ELSE
                     A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODX_E(I))
                     A_M(IJK,east,M) = -(HALF*BC_HW_PHI(L)-ODX_E(I))
                     B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L))
                  ENDIF
               ENDIF

            ELSEIF (CUT_CELL_AT(IM_OF(IJK))) THEN
               BCV = BC_ID(IM_OF(IJK))
               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               IF (BCT==CG_NSW.OR.BCT==CG_FSW.OR.BCT==CG_PSW) THEN
                  L = BCV
                  IF (BC_HW_PHI(L) == UNDEFINED) THEN
                     A_M(IJK,west,M) = -HALF
                     A_M(IJK,0,M) = -HALF
                     B_M(IJK,M) = -BC_PHIW(L)
                  ELSE
                     A_M(IJK,west,M) = -(HALF*BC_HW_PHI(L)-ODX_E(IM))
                     A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODX_E(IM))
                     B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L))
                  ENDIF
               ENDIF

            ELSEIF (CUT_CELL_AT(JP_OF(IJK))) THEN
               BCV = BC_ID(JP_OF(IJK))
               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               IF (BCT==CG_NSW.OR.BCT==CG_FSW.OR.BCT==CG_PSW) THEN
                  L = BCV
                  IF (BC_HW_PHI(L) == UNDEFINED) THEN
                     A_M(IJK,north,M) = -HALF
                     A_M(IJK,0,M) = -HALF
                     B_M(IJK,M) = -BC_PHIW(L)
                  ELSE
                     A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODY_N(J))
                     A_M(IJK,north,M) = -(HALF*BC_HW_PHI(L)-ODY_N(J))
                     B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L))
                  ENDIF
               ENDIF

            ELSEIF (CUT_CELL_AT(JM_OF(IJK))) THEN
               BCV = BC_ID(JM_OF(IJK))
               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               IF (BCT==CG_NSW.OR.BCT==CG_FSW.OR.BCT==CG_PSW) THEN
                  L = BCV
                  IF (BC_HW_PHI(L) == UNDEFINED) THEN
                     A_M(IJK,south,M) = -HALF
                     A_M(IJK,0,M) = -HALF
                     B_M(IJK,M) = -BC_PHIW(L)
                  ELSE
                     A_M(IJK,south,M) = -(HALF*BC_HW_PHI(L)-ODY_N(JM))
                     A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+ODY_N(JM))
                     B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L))
                  ENDIF
               ENDIF

            ELSEIF (CUT_CELL_AT(KP_OF(IJK))) THEN
               BCV = BC_ID(KP_OF(IJK))
               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               IF (BCT==CG_NSW.OR.BCT==CG_FSW.OR.BCT==CG_PSW) THEN
                  L = BCV
                  IF (BC_HW_PHI(L) == UNDEFINED) THEN
                     A_M(IJK,top,M) = -HALF
                     A_M(IJK,0,M) = -HALF
                     B_M(IJK,M) = -BC_PHIW(L)
                  ELSE
                     A_M(IJK,0,M)=-(HALF*BC_HW_PHI(L)+OX(I)*ODZ_T(K))
                     A_M(IJK,top,M)=-(HALF*BC_HW_PHI(L)-OX(I)*ODZ_T(K))
                     B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L))
                  ENDIF
               ENDIF

            ELSEIF (CUT_CELL_AT(KM_OF(IJK))) THEN
               BCV = BC_ID(KM_OF(IJK))
               IF(BCV > 0 ) THEN
                  BCT = BC_TYPE_ENUM(BCV)
               ELSE
                  BCT = NONE
               ENDIF

               IF (BCT==CG_NSW.OR.BCT==CG_FSW.OR.BCT==CG_PSW) THEN
                  L = BCV
                  IF (BC_HW_PHI(L) == UNDEFINED) THEN
                     A_M(IJK,bottom,M) = -HALF
                     A_M(IJK,0,M) = -HALF
                     B_M(IJK,M) = -BC_PHIW(L)
                  ELSE
                     A_M(IJK,bottom,M) = -(HALF*BC_HW_PHI(L)-OX(I)*ODZ_T(KM))
                     A_M(IJK,0,M) = -(HALF*BC_HW_PHI(L)+OX(I)*ODZ_T(KM))
                     B_M(IJK,M) = -(BC_HW_PHI(L)*BC_PHIW(L)+BC_C_PHI(L))
                  ENDIF
               ENDIF  ! BCT

            ENDIF   ! end if/else cut cell next to blocked cell

         ENDIF   ! end if (blocked_cell_at(ijk))

         IF(CUT_CELL_AT(IJK)) THEN
            BCV = BC_ID(IJK)
            IF(BCV > 0 ) THEN
               BCT = BC_TYPE_ENUM(BCV)
            ELSE
               BCT = NONE
            ENDIF

            IF (BCT==CG_MI.OR.BCT==CG_PO) THEN
               L = BCV
               A_M(IJK,east,M) = ZERO
               A_M(IJK,west,M) = ZERO
               A_M(IJK,north,M) = ZERO
               A_M(IJK,south,M) = ZERO
               A_M(IJK,top,M) = ZERO
               A_M(IJK,bottom,M) = ZERO
               A_M(IJK,0,M) = -ONE
               B_M(IJK,M) = -BC_PHIF(L)
            ENDIF
         ENDIF

      ENDDO   ! end do (ijk=ijkstart3,ijkend3)


      RETURN
      END SUBROUTINE BC_PHI_CG
