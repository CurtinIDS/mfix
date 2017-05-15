!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BOUNDARY_CONDITIONS                               !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Check boundary condition specifications                    !
!     - convert physical locations to i, j, k's (GET_FLOW_BC)          !
!     - compute area of boundary surfaces (GET_BC_AREA)                !
!     - convert mass and volumetric flows to velocities (FLOW_TO_VEL)  !
!     - check specification of physical quantities                     !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BOUNDARY_CONDITIONS

! Global Variables:
!---------------------------------------------------------------------//
! Total number of (actual) continuum solids.
      use physprop, only: SMAX, MMAX
! Total number of discrete solids.
      use discretelement, only: DES_MMAX
! Flag: BC dimensions or Type is specified
      use bc, only: BC_DEFINED
! Use specified BC type
      use bc
! User specified BC solids bulk density
      use bc, only: BC_ROP_s
! Solids volume fraction at BC
      use bc, only: BC_EP_s
      use bc, only: BC_EP_g
! Run-time flag for DEM solids
      use run, only: DEM_SOLIDS
! Run-time flag for PIC solids
      use run, only: PIC_SOLIDS

! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constants
      use param1, only: ZERO, ONE, UNDEFINED
! Maximum number of BCs
      use param, only: DIMENSION_BC
! Maximum number of disperse phases
      use param, only: DIM_M

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! Loop counter for BCs
      INTEGER :: BCV
! Total number of solids phases (continuum + discrete)
      INTEGER :: MMAX_TOT
! Flag to skip checks on indexed solid phase.
      LOGICAL :: SKIP(1:DIM_M)
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BOUNDARY_CONDITIONS")

! Determine which BCs are DEFINED
      CALL CHECK_BC_GEOMETRY

! Total number of solids. (this won't work for GHD/hybrid)
      MMAX_TOT = SMAX + DES_MMAX

! Loop over each defined BC and check the user data.
      DO BCV = 1, DIMENSION_BC

         IF (BC_DEFINED(BCV)) THEN

! Determine which solids phases are present.
            SKIP=(BC_ROP_S(BCV,:)==UNDEFINED.OR.BC_ROP_S(BCV,:)==ZERO) &
               .AND.(BC_EP_S(BCV,:)==UNDEFINED.OR.BC_EP_S(BCV,:)==ZERO)

            IF(MMAX_TOT == 1 .AND. BC_EP_g(BCV)/=ONE) SKIP(1) = .FALSE.

            SELECT CASE (BC_TYPE_ENUM(BCV))

            CASE (MASS_INFLOW)
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)
               CALL CHECK_BC_MASS_INFLOW(MMAX_TOT, SKIP, BCV)
               CALL CHECK_BC_INFLOW(MMAX_TOT,SKIP,BCV)

            CASE (P_INFLOW)
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)
               CALL CHECK_BC_P_INFLOW(MMAX_TOT, SKIP, BCV)
               CALL CHECK_BC_INFLOW(MMAX_TOT, SKIP, BCV)
               CALL CHECK_BC_OUTFLOW(MMAX_TOT, BCV)

            CASE (OUTFLOW)
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)
               CALL CHECK_BC_OUTFLOW(MMAX_TOT, BCV)

            CASE (MASS_OUTFLOW)
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)
               CALL CHECK_BC_MASS_OUTFLOW(MMAX_TOT, BCV)
               CALL CHECK_BC_OUTFLOW(MMAX_TOT, BCV)

            CASE (P_OUTFLOW)
               CALL CHECK_BC_GEOMETRY_FLOW(BCV)
               CALL CHECK_BC_P_OUTFLOW(MMAX_TOT, BCV)
               CALL CHECK_BC_OUTFLOW(MMAX_TOT, BCV)

            CASE (FREE_SLIP_WALL)
               CALL CHECK_BC_GEOMETRY_WALL(BCV)
               CALL CHECK_BC_WALLS(MMAX_TOT, SKIP, BCV)

            CASE (NO_SLIP_WALL)
               CALL CHECK_BC_GEOMETRY_WALL(BCV)
               CALL CHECK_BC_WALLS(MMAX_TOT, SKIP, BCV)

            CASE (PAR_SLIP_WALL)
               CALL CHECK_BC_GEOMETRY_WALL(BCV)
               CALL CHECK_BC_WALLS(MMAX_TOT, SKIP, BCV)

            END SELECT

! Check whether BC values are specified for undefined BC locations
         ELSEIF(BC_TYPE_ENUM(BCV) /= DUMMY .AND.                          &
            .NOT.IS_CG(BC_TYPE_ENUM(BCV))) THEN

            CALL CHECK_BC_RANGE(BCV)

         ENDIF
      ENDDO

      MMAX_TOT = MMAX+DES_MMAX
! Additional checks needed for DEM boundaries
      IF(DEM_SOLIDS) CALL CHECK_BC_DEM(MMAX_TOT)
! Additional checks needed for PIC inflow/outflow boundaries
      IF(PIC_SOLIDS) CALL CHECK_BC_PIC(MMAX_TOT)

! Cleanup and exit.
      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_BOUNDARY_CONDITIONS



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_RANGE                                          !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Verify that data was not given for undefined BC regions.   !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_RANGE(BCV)

! Global Variables:
!---------------------------------------------------------------------//
! Gas phase BC varaibles
      use bc, only: BC_EP_g, BC_T_g, BC_X_g, BC_P_g
      use bc, only: BC_U_g, BC_V_g, BC_W_g
! Solids phase BC variables.
      USE bc, only: BC_EP_s, BC_ROP_s, BC_T_s, BC_X_s
      use bc, only: BC_U_s, BC_V_s, BC_W_s
! Scalar equation BC variables.
      USE bc, only: BC_SCALAR


! Global Parameters:
!---------------------------------------------------------------------//
! Parameter constant for unspecified values.
      use param1, only: UNDEFINED
! Maximum number of disperse phases.
      use param, only: DIM_M
! Maximum number of species gas/solids
      use param, only: DIMENSION_N_G, DIMENSION_N_S
! Maximum number of scalar equations.
      use param, only: DIM_SCALAR


! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Dummy Arguments:
!---------------------------------------------------------------------/
! Boundary condition index.
      INTEGER, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Generic loop varaibles.
      INTEGER :: M, N
!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_BC_RANGE")


! Check gas phase variables.
      IF(BC_U_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_U_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF(BC_V_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_V_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF (BC_W_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_W_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF (BC_EP_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_EP_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF (BC_P_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_P_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF
      IF (BC_T_G(BCV) /= UNDEFINED) THEN
         WRITE(ERR_MSG,1100) trim(iVar('BC_T_g',BCV))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      DO N = 1, DIMENSION_N_G
         IF(BC_X_G(BCV,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_X_g',BCV,N))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO

! Check solids phase variables.
      DO M = 1, DIM_M
         IF(BC_ROP_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_ROP_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_EP_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_EP_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_U_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_U_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_V_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_V_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         IF(BC_W_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_W_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
         IF(BC_T_S(BCV,M) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_T_s',BCV,M))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

         DO N = 1, DIMENSION_N_S
            IF(BC_X_S(BCV,M,N) /= UNDEFINED) THEN
               WRITE(ERR_MSG,1100) trim(iVar('BC_X_s',BCV,M,N))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      ENDDO

! Check scalar equation variables.
      DO N = 1, DIM_SCALAR
         IF(BC_Scalar(BCV,N) /= UNDEFINED) THEN
            WRITE(ERR_MSG,1100) trim(iVar('BC_Scalar',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDDO


      CALL FINL_ERR_MSG


      RETURN

 1100 FORMAT('Error 1100:',A,' specified for an undefined BC location')

      END SUBROUTINE CHECK_BC_RANGE
