!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: SET_BC_FLOW                                             !
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
      SUBROUTINE SET_BC_FLOW

! Global Variables:
!---------------------------------------------------------------------//
! Total number of (actual) continuum solids.
      use physprop, only: SMAX
! Total number of discrete solids.
      use discretelement, only: DES_MMAX
! Flag: BC dimensions or Type is specified
      use bc, only: BC_DEFINED
! Use specified BC type
      use bc
! User specifed BC solids bulk density
      use bc, only: BC_ROP_s
! Solids volume fraction at BC
      use bc, only: BC_EP_s
      use bc, only: BC_EP_g

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
      CALL INIT_ERR_MSG("SET_BC_FLOW")

! Total number of solids.
      MMAX_TOT = SMAX + DES_MMAX

! Loop over each defined BC and check the user data.
      DO BCV = 1, DIMENSION_BC

         IF(.NOT.BC_DEFINED(BCV)) CYCLE

! Determine which solids phases are present.
         SKIP=(BC_ROP_S(BCV,:)==UNDEFINED.OR.BC_ROP_S(BCV,:)==ZERO) &
            .AND.(BC_EP_S(BCV,:)==UNDEFINED.OR.BC_EP_S(BCV,:)==ZERO)

         IF(MMAX_TOT == 1 .AND. BC_EP_g(BCV)/=ONE) SKIP(1) = .FALSE.

         SELECT CASE (BC_TYPE_ENUM(BCV))

         CASE (MASS_INFLOW)
            CALL FLOW_TO_VEL_NEW(.TRUE., MMAX_TOT, SKIP, BCV)
            CALL CHECK_BC_VEL_INFLOW(MMAX_TOT, SKIP, BCV)

         CASE (MASS_OUTFLOW)
            CALL FLOW_TO_VEL_NEW(.TRUE., MMAX_TOT, SKIP, BCV)
            CALL CHECK_BC_VEL_OUTFLOW(MMAX_TOT, SKIP, BCV)
         END SELECT
      ENDDO

! Cleanup and exit.
      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE SET_BC_FLOW



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_VEL_INFLOW                                      !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
! Comments:                                                            !
!     The velocities at the inflow face are fixed and the momentum     !
!     equations are not solved in the inflow cells. Since the flow is  !
!     into the domain all other scalars that are used need to be       !
!     specified (e.g., mass fractions, void fraction, etc.,)           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_VEL_INFLOW(M_TOT, SKIP, BCV)

      USE param, only: DIM_M
      USE param1, only: ZERO
      USE param1, only: UNDEFINED

      use geometry, only: NO_I
      use geometry, only: NO_J
      use geometry, only: NO_K

      use bc

      use error_manager

      IMPLICIT NONE


      INTEGER, INTENT(in) :: BCV
      INTEGER, INTENT(in) :: M_TOT

      LOGICAL, INTENT(in) :: SKIP(DIM_M)

! loop/variable indices
      INTEGER :: M

      CALL INIT_ERR_MSG("CHECK_BC_VEL_INFLOW")


! Check that gas phase velocities are defined.
      IF(BC_U_G(BCV) == UNDEFINED) THEN
         IF(NO_I) THEN
            BC_U_G(BCV) = ZERO
         ELSE
            WRITE(ERR_MSG,1000) trim(iVar('BC_U_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      IF (BC_V_G(BCV) == UNDEFINED) THEN
         IF (NO_J) THEN
            BC_V_G(BCV) = ZERO
         ELSE
            WRITE(ERR_MSG,1000) trim(iVar('BC_V_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      IF(BC_W_G(BCV) == UNDEFINED) THEN
         IF (NO_K) THEN
            BC_W_G(BCV) = ZERO
         ELSE
            WRITE(ERR_MSG,1000) trim(iVar('BC_W_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

! Check that solids phase velocities are defined.
      DO M = 1, M_TOT
         IF(BC_U_S(BCV,M) == UNDEFINED) THEN
            IF(SKIP(M) .OR. NO_I) THEN
               BC_U_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_U_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_V_S(BCV,M) == UNDEFINED) THEN
            IF(SKIP(M) .OR. NO_J) THEN
               BC_V_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_V_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_W_S(BCV,M) == UNDEFINED) THEN
            IF(SKIP(M) .OR. NO_K) THEN
               BC_W_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_W_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDDO

! Check that gas phase velocities are consistent.
      SELECT CASE (BC_PLANE(BCV))

      CASE ('W')
         IF(BC_U_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_U_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_U_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('E')
         IF(BC_U_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_U_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_U_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M)), '>'
              CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('S')
         IF(BC_V_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_V_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_V_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('N')
         IF(BC_V_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_V_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_V_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('B')
         IF(BC_W_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_W_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_W_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('T')
         IF(BC_W_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_W_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_W_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      END SELECT

 1300 FORMAT('Error 1300: Invalid flow direction. ',A,' should be ',   &
         A,' zero. ',/'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_VEL_INFLOW

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_VEL_OUTFLOW                                     !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Provided a detailed error message when the sum of volume    !
!                                                                      !
! Comments:                                                            !
!     The velocities at the outflow face are fixed and the momentum    !
!     equations are not solved in the outflow cells. Since the flow    !
!     is out of the domain none of the other scalars should need to    !
!     be specified (e.g., mass fractions, void fraction, etc.,).       !
!     Such values will become defined according to their adjacent      !
!     fluid cell                                                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_VEL_OUTFLOW(M_TOT, SKIP, BCV)

      USE param
      USE param1
      USE geometry
      USE fldvar
      USE physprop
      USE run
      USE bc
      USE indices
      USE funits
      USE scalars
      USE compar
      USE sendrecv
      USE discretelement
      USE mfix_pic
      USE cutcell

      use error_manager

      IMPLICIT NONE

! loop/variable indices
      INTEGER, intent(in) :: BCV
      INTEGER, intent(in) :: M_TOT
      LOGICAL, intent(in) :: SKIP(DIM_M)

! Loop variable
      INTEGER :: M

      CALL INIT_ERR_MSG("CHECK_BC_VEL_OUTFLOW")

! Check that gas phase velocities are defined.
      IF(BC_U_G(BCV) == UNDEFINED) THEN
         IF(NO_I) THEN
            BC_U_G(BCV) = ZERO
         ELSE
            WRITE(ERR_MSG,1000) trim(iVar('BC_U_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      IF (BC_V_G(BCV) == UNDEFINED) THEN
         IF (NO_J) THEN
            BC_V_G(BCV) = ZERO
         ELSE
            WRITE(ERR_MSG,1000) trim(iVar('BC_V_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      IF(BC_W_G(BCV) == UNDEFINED) THEN
         IF (NO_K) THEN
            BC_W_G(BCV) = ZERO
         ELSE
            WRITE(ERR_MSG,1000) trim(iVar('BC_W_g',BCV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

! Check that solids phase velocities are defined.
      DO M = 1, M_TOT
         IF(BC_U_S(BCV,M) == UNDEFINED) THEN
            IF(SKIP(M) .OR. NO_I) THEN
               BC_U_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_U_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_V_S(BCV,M) == UNDEFINED) THEN
            IF(SKIP(M) .OR. NO_J) THEN
               BC_V_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_V_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_W_S(BCV,M) == UNDEFINED) THEN
            IF(SKIP(M) .OR. NO_K) THEN
               BC_W_S(BCV,M) = ZERO
            ELSE
               WRITE(ERR_MSG,1000) trim(iVar('BC_W_s',BCV,M))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF
      ENDDO


! Check that gas phase velocities are consistent.
      SELECT CASE (BC_PLANE(BCV))

      CASE ('W')
         IF(BC_U_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_U_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_U_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('E')
         IF(BC_U_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_U_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_U_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_U_s',BCV,M)), '<'
              CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('S')
         IF(BC_V_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_V_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_V_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('N')
         IF(BC_V_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_V_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_V_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_V_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('B')
         IF(BC_W_G(BCV) < ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_W_g',BCV)), '>'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_W_S(BCV,M) < ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M)), '>'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      CASE('T')
         IF(BC_W_G(BCV) > ZERO) THEN
            WRITE(ERR_MSG,1300) trim(iVar('BC_W_g',BCV)), '<'
            CALL FLUSH_ERR_MSG
         ENDIF
         DO M = 1, M_TOT
            IF(BC_W_S(BCV,M) > ZERO) THEN
               WRITE(ERR_MSG, 1300) trim(iVar('BC_W_s',BCV,M)), '<'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDDO

      END SELECT

 1300 FORMAT('Error 1300: Invalid flow direction. ',A,' should be ',   &
         A,' zero. ',/'Please correct the mfix.dat file.')
      CALL FINL_ERR_MSG


      RETURN


 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

      END SUBROUTINE CHECK_BC_VEL_OUTFLOW
