!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: CHECK_BC_GEOMETRY                                        !
! Author: J.Musser                                    Date: 01-Mar-14  !
!                                                                      !
! Purpose: Determine if BCs are "DEFINED" and that they contain the    !
! minimum amount of geometry data.                                     !
!                                                                      !
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
      SUBROUTINE CHECK_BC_GEOMETRY

! Global Variables:
!---------------------------------------------------------------------//
! Flag: BC contains geometric data and/or specified type
      use bc, only: BC_DEFINED
! User specified BC
      use bc
! User specified: BC geometry
      use bc, only: BC_X_e, BC_X_w, BC_I_e, BC_I_w
      use bc, only: BC_Y_n, BC_Y_s, BC_J_n, BC_J_s
      use bc, only: BC_Z_t, BC_Z_b, BC_K_t, BC_K_b
! User specified: System geometry
      use geometry, only: NO_I, XLENGTH
      use geometry, only: NO_J, YLENGTH
      use geometry, only: NO_K, ZLENGTH
! Flag: Reinitializing field conditions.
      use run, only: REINITIALIZING

! Global Parameters:
!---------------------------------------------------------------------//
! The max number of BCs.
      use param, only: DIMENSION_BC
! Parameter constants
      use param1, only: ZERO, UNDEFINED, UNDEFINED_I, UNDEFINED_C

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Local Variables:
!---------------------------------------------------------------------//
! loop/variable indices
      INTEGER :: BCV, I
! Error flag
      LOGICAL :: RECOGNIZED_BC_TYPE
!......................................................................!

! Skip this routine if reinitializing as BC locations cannot be changed.
      IF(REINITIALIZING) RETURN

      CALL INIT_ERR_MSG("CHECK_BC_GEOMETRY")

      BC_TYPE_ENUM = BLANK

      L50: DO BCV = 1, DIMENSION_BC

         BC_DEFINED(BCV) = .FALSE.
         IF(BC_X_W(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE.
         IF(BC_X_E(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE.
         IF(BC_Y_S(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE.
         IF(BC_Y_N(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE.
         IF(BC_Z_B(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE.
         IF(BC_Z_T(BCV) /= UNDEFINED)   BC_DEFINED(BCV) = .TRUE.
         IF(BC_I_W(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE.
         IF(BC_I_E(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE.
         IF(BC_J_S(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE.
         IF(BC_J_N(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE.
         IF(BC_K_B(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE.
         IF(BC_K_T(BCV) /= UNDEFINED_I) BC_DEFINED(BCV) = .TRUE.
         IF(BC_TYPE(BCV) == 'CG_NSW')   BC_DEFINED(BCV) = .TRUE.
         IF(BC_TYPE(BCV) == 'CG_FSW')   BC_DEFINED(BCV) = .TRUE.
         IF(BC_TYPE(BCV) == 'CG_PSW')   BC_DEFINED(BCV) = .TRUE.
         IF(BC_TYPE(BCV) == 'CG_MI')    BC_DEFINED(BCV) = .TRUE.
         IF(BC_TYPE(BCV) == 'CG_PO')    BC_DEFINED(BCV) = .TRUE.

         IF (BC_TYPE(BCV) == 'DUMMY') BC_DEFINED(BCV) = .FALSE.

         IF(BC_TYPE(BCV)/=UNDEFINED_C .AND. BC_TYPE(BCV)/='DUMMY')THEN

            RECOGNIZED_BC_TYPE = .FALSE.
            DO I = 1, DIM_BCTYPE
                VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
                IF(VALID_BC_TYPE(I) == trim(BC_TYPE(BCV))) THEN
                   BC_TYPE_ENUM(BCV) = VALID_BC_TYPE_ENUM(I)
                   RECOGNIZED_BC_TYPE = .TRUE.
                   EXIT
                ENDIF
            ENDDO

            IF(.NOT.RECOGNIZED_BC_TYPE) THEN
               WRITE(ERR_MSG, 1100) trim(iVar('BC_TYPE',BCV)), &
                  BC_TYPE_ENUM(BCV), VALID_BC_TYPE
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(.NOT.BC_DEFINED(BCV)) CYCLE
         IF(IS_CG(BC_TYPE_ENUM(BCV))) CYCLE

         IF(BC_X_W(BCV)==UNDEFINED .AND. BC_I_W(BCV)==UNDEFINED_I) THEN
            IF(NO_I) THEN
               BC_X_W(BCV) = ZERO
            ELSE
               WRITE(ERR_MSG,1101) BCV, 'BC_X_w and BC_I_w'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_X_E(BCV)==UNDEFINED .AND. BC_I_E(BCV)==UNDEFINED_I) THEN
            IF(NO_I) THEN
               BC_X_E(BCV) = XLENGTH
            ELSE
               WRITE(ERR_MSG, 1101) BCV, 'BC_X_e and BC_I_e'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_Y_S(BCV)==UNDEFINED .AND. BC_J_S(BCV)==UNDEFINED_I) THEN
            IF(NO_J) THEN
               BC_Y_S(BCV) = ZERO
            ELSE
               WRITE(ERR_MSG, 1101) BCV, 'BC_Y_s and BC_J_s'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_Y_N(BCV)==UNDEFINED .AND. BC_J_N(BCV)==UNDEFINED_I) THEN
            IF(NO_J) THEN
               BC_Y_N(BCV) = YLENGTH
            ELSE
               WRITE(ERR_MSG, 1101) BCV, 'BC_Y_n and BC_J_n'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_Z_B(BCV)==UNDEFINED .AND. BC_K_B(BCV)==UNDEFINED_I) THEN
            IF(NO_K) THEN
               BC_Z_B(BCV) = ZERO
            ELSE
               WRITE(ERR_MSG, 1101) BCV, 'BC_Z_b and BC_K_b'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

         IF(BC_Z_T(BCV)==UNDEFINED .AND. BC_K_T(BCV)==UNDEFINED_I) THEN
            IF(NO_K) THEN
               BC_Z_T(BCV) = ZLENGTH
            ELSE
               WRITE(ERR_MSG, 1101) BCV, 'BC_Z_t and BC_K_t'
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF
         ENDIF

 1101 FORMAT('Error 1101: Boundary condition ',I3,' is ill-defined.',/ &
         A,' are not specified.',/'Please correct the mfix.dat file.')

! Swap BC aliases for the "full name" complement.
         DO I = 1, DIM_BCTYPE
            VALID_BC_TYPE(I) = TRIM(VALID_BC_TYPE(I))
            IF(VALID_BC_TYPE(I) == BC_TYPE(BCV)) THEN
               IF(MOD(I,2) == 0) THEN
                  BC_TYPE_ENUM(BCV) = VALID_BC_TYPE_ENUM(I-1)
               ELSE
                  BC_TYPE_ENUM(BCV) = VALID_BC_TYPE_ENUM(I)
               ENDIF
               CYCLE  L50
            ENDIF
         ENDDO

         WRITE(ERR_MSG, 1100) trim(iVar('BC_TYPE',BCV)),               &
            BC_TYPE_ENUM(BCV), VALID_BC_TYPE
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ENDDO L50   ! end loop over (bcv=1,dimension_bc)

      CALL FINL_ERR_MSG

      RETURN


 1100 FORMAT('Error 1100: Illegal entry: ',A,' = ',A,/'Valid entries:',&
         ' ',10(/5X,A,2x,A),/5X,A)

      END SUBROUTINE CHECK_BC_GEOMETRY



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_GEOMETRY_WALL                                  !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for walls BC's         !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_GEOMETRY_WALL(BCV)

! Global Variables:
!---------------------------------------------------------------------//
! Boundary condition locations and corresponding grid index
      use bc, only: BC_X_w, BC_X_e, BC_I_w, BC_I_e
      use bc, only: BC_Y_s, BC_Y_n, BC_J_s, BC_J_n
      use bc, only: BC_Z_b, BC_Z_t, BC_K_b, BC_K_t
! Basic grid information
      use geometry, only: NO_I, XLENGTH, DX, IMAX, IMAX2, XMIN
      use geometry, only: NO_J, YLENGTH, DY, JMAX, JMAX2
      use geometry, only: NO_K, ZLENGTH, DZ, KMAX, KMAX2
! Flag: Reinitializing field conditions.
      use run, only: REINITIALIZING
! Function to compare two values
      use toleranc, only: COMPARE

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, UNDEFINED_I, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      IMPLICIT NONE


! Dummy Arguments:
!---------------------------------------------------------------------//
! Index of boundary condition.
      INTEGER, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Calculated indices of the wall boundary
      INTEGER :: I_w , I_e , J_s , J_n , K_b , K_t
! Integer error flag
      INTEGER :: IER
!......................................................................!


! Skip this routine if reinitializing as BC locations cannot be changed.
      IF(REINITIALIZING) RETURN

      CALL INIT_ERR_MSG("CHECK_BC_GEOMETRY_WALL")

      IF(BC_X_W(BCV)/=UNDEFINED .AND. BC_X_E(BCV)/=UNDEFINED) THEN

! setting indices to 1 if there is no variation in the i (x) direction
         IF (NO_I) THEN
            I_W = 1
            I_E = 1
         ELSE
            CALL CALC_CELL (XMIN, BC_X_W(BCV), DX, IMAX, I_W)
            I_W = I_W + 1
            CALL CALC_CELL (XMIN, BC_X_E(BCV), DX, IMAX, I_E)
! BC along zy plane, checking if far west or far east of domain
            IF(BC_X_W(BCV) == BC_X_E(BCV)) THEN
               IF(COMPARE(BC_X_W(BCV),XMIN)) THEN
                  I_W = 1
                  I_E = 1
               ELSEIF(COMPARE(BC_X_W(BCV),XMIN+XLENGTH)) THEN
                  I_W = IMAX2
                  I_E = IMAX2
               ENDIF
            ENDIF
         ENDIF

! checking/setting corresponding i indices according to specified x
! coordinates
         IF(BC_I_W(BCV)/=UNDEFINED_I .OR. BC_I_E(BCV)/=UNDEFINED_I) THEN
            CALL LOCATION_CHECK (BC_I_W(BCV), I_W, BCV, 'BC - west')
            CALL LOCATION_CHECK (BC_I_E(BCV), I_E, BCV, 'BC - east')
         ELSE
            BC_I_W(BCV) = I_W
            BC_I_E(BCV) = I_E
         ENDIF
      ENDIF


      IF(BC_Y_S(BCV)/=UNDEFINED .AND. BC_Y_N(BCV)/=UNDEFINED) THEN
! setting indices to 1 if there is no variation in the j (y) direction
         IF(NO_J) THEN
            J_S = 1
            J_N = 1
         ELSE
            CALL CALC_CELL (ZERO, BC_Y_S(BCV), DY, JMAX, J_S)
            J_S = J_S + 1
            CALL CALC_CELL (ZERO, BC_Y_N(BCV), DY, JMAX, J_N)
! BC along xz plane, checking if far south or far north of domain
            IF(BC_Y_S(BCV) == BC_Y_N(BCV)) THEN
               IF(COMPARE(BC_Y_S(BCV),ZERO)) THEN
                  J_S = 1
                  J_N = 1
               ELSE IF (COMPARE(BC_Y_S(BCV),YLENGTH)) THEN
                  J_S = JMAX2
                  J_N = JMAX2
               ENDIF
            ENDIF
         ENDIF
! checking/setting corresponding j indices according to specified y
! coordinates
         IF(BC_J_S(BCV)/=UNDEFINED_I .OR. BC_J_N(BCV)/=UNDEFINED_I) THEN
            CALL LOCATION_CHECK (BC_J_S(BCV), J_S, BCV, 'BC - south')
            CALL LOCATION_CHECK (BC_J_N(BCV), J_N, BCV, 'BC - north')
         ELSE
            BC_J_S(BCV) = J_S
            BC_J_N(BCV) = J_N
         ENDIF
      ENDIF

      IF(BC_Z_B(BCV)/=UNDEFINED .AND. BC_Z_T(BCV)/=UNDEFINED) THEN
! setting indices to 1 if there is no variation in the k (z) direction
         IF(NO_K)THEN
            K_B = 1
            K_T = 1
         ELSE
            CALL CALC_CELL (ZERO, BC_Z_B(BCV), DZ, KMAX, K_B)
            K_B = K_B + 1
            CALL CALC_CELL (ZERO, BC_Z_T(BCV), DZ, KMAX, K_T)
! BC along xy plane, checking if far bottom or far top of domain
            IF(BC_Z_B(BCV) == BC_Z_T(BCV)) THEN
               IF(COMPARE(BC_Z_B(BCV),ZERO)) THEN
                  K_B = 1
                  K_T = 1
               ELSEIF(COMPARE(BC_Z_B(BCV),ZLENGTH)) THEN
                  K_B = KMAX2
                  K_T = KMAX2
               ENDIF
            ENDIF
         ENDIF
! checking/setting corresponding j indices according to specified y
! coordinates
         IF(BC_K_B(BCV)/=UNDEFINED_I .OR.BC_K_T(BCV)/=UNDEFINED_I) THEN
            CALL LOCATION_CHECK (BC_K_B(BCV), K_B, BCV, 'BC - bottom')
            CALL LOCATION_CHECK (BC_K_T(BCV), K_T, BCV, 'BC - top')
         ELSE
            BC_K_B(BCV) = K_B
            BC_K_T(BCV) = K_T
         ENDIF
      ENDIF


! CHECK FOR VALID VALUES
      IER = 0
      IF (BC_K_B(BCV)<1 .OR. BC_K_B(BCV)>KMAX2) IER = 1
      IF (BC_J_S(BCV)<1 .OR. BC_J_S(BCV)>JMAX2) IER = 1
      IF (BC_I_W(BCV)<1 .OR. BC_I_W(BCV)>IMAX2) IER = 1
      IF (BC_K_T(BCV)<1 .OR. BC_K_T(BCV)>KMAX2) IER = 1
      IF (BC_J_N(BCV)<1 .OR. BC_J_N(BCV)>JMAX2) IER = 1
      IF (BC_I_E(BCV)<1 .OR. BC_I_E(BCV)>IMAX2) IER = 1
      IF (BC_K_B(BCV) > BC_K_T(BCV)) IER = 1
      IF (BC_J_S(BCV) > BC_J_N(BCV)) IER = 1
      IF (BC_I_W(BCV) > BC_I_E(BCV)) IER = 1

      IF(IER /= 0)THEN
         WRITE(ERR_MSG,1100) BCV,                                      &
            'X', BC_X_W(BCV), BC_X_E(BCV),'I',BC_I_W(BCV),BC_I_E(BCV), &
            'Y', BC_Y_S(BCV), BC_Y_N(BCV),'J',BC_J_S(BCV),BC_J_N(BCV), &
            'Z', BC_Z_B(BCV), BC_Z_T(BCV),'K',BC_K_B(BCV),BC_K_T(BCV)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Invalid location specified for BC ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),/         &
         'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_BC_GEOMETRY_WALL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: CHECK_BC_GEOMETRY_FLOW                                  !
!  Author: P. Nicoletti                               Date: 10-DEC-91  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for flow BC's. Also    !
!           set value of bc_plane for flow BC's.                       !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_BC_GEOMETRY_FLOW(BCV)

! Global Variables:
!---------------------------------------------------------------------//
! Boundary condition locations and corresponding grid index
      use bc, only: BC_X_w, BC_X_e, BC_I_w, BC_I_e
      use bc, only: BC_Y_s, BC_Y_n, BC_J_s, BC_J_n
      use bc, only: BC_Z_b, BC_Z_t, BC_K_b, BC_K_t
! Basic grid information
      use geometry, only: NO_I, DX, IMAX, IMAX2, XMIN
      use geometry, only: NO_J, DY, JMAX, JMAX2
      use geometry, only: NO_K, DZ, KMAX, KMAX2
! Flag: Reinitializing field conditions.
      use run, only: REINITIALIZING

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ZERO, UNDEFINED_I, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE


! Dummy Arguments:
!---------------------------------------------------------------------//
! Index of boundary condition.
      INTEGER, INTENT(in) :: BCV

! Local Variables:
!---------------------------------------------------------------------//
! Calculated indices of the wall boundary
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t
! Indices for error checking
      INTEGER :: IER

! surface indictors:
! a value of T indicates that the defined boundary region does not
! vary in indicated coordinate direction. that is, if bc_x_w is
! equal to bc_x_e then the boundary region must be in the yz plane
      LOGICAL :: X_CONSTANT, Y_CONSTANT, Z_CONSTANT
!......................................................................!

! Skip this routine if reinitializing as BC locations cannot be changed.
      IF(REINITIALIZING) RETURN

      CALL INIT_ERR_MSG("CHECK_BC_GEOMETRY_FLOW")

      X_CONSTANT = .TRUE.
      Y_CONSTANT = .TRUE.
      Z_CONSTANT = .TRUE.

      IF (BC_X_W(BCV)/=UNDEFINED .AND. BC_X_E(BCV)/=UNDEFINED) THEN
         CALL CALC_CELL (XMIN, BC_X_W(BCV), DX, IMAX, I_W)
         CALL CALC_CELL (XMIN, BC_X_E(BCV), DX, IMAX, I_E)
         IF (BC_X_W(BCV) /= BC_X_E(BCV)) THEN
            X_CONSTANT = .FALSE.
            I_W = I_W + 1
            IF(BC_I_W(BCV)/=UNDEFINED_I.OR.BC_I_E(BCV)/=UNDEFINED_I)THEN
               CALL LOCATION_CHECK (BC_I_W(BCV), I_W, BCV, 'BC - west')
               CALL LOCATION_CHECK (BC_I_E(BCV), I_E, BCV, 'BC - east')
            ENDIF
         ENDIF
         BC_I_W(BCV) = I_W
         BC_I_E(BCV) = I_E
      ELSE
         IF(BC_I_W(BCV) /= UNDEFINED_I) &
            CALL CALC_LOC (XMIN,DX,BC_I_W(BCV),BC_X_W(BCV))
         IF(BC_I_E(BCV) /= UNDEFINED_I) &
            CALL CALC_LOC (XMIN,DX,BC_I_E(BCV),BC_X_E(BCV))
         IF(BC_X_W(BCV) /= BC_X_E(BCV)) X_CONSTANT = .FALSE.
      ENDIF

! If there is no variation in the I direction set indices to 1
      IF(NO_I) THEN
         BC_I_W(BCV) = 1
         BC_I_E(BCV) = 1
      ENDIF

      IF (BC_Y_S(BCV)/=UNDEFINED .AND. BC_Y_N(BCV)/=UNDEFINED) THEN
         CALL CALC_CELL (ZERO, BC_Y_S(BCV), DY, JMAX, J_S)
         CALL CALC_CELL (ZERO, BC_Y_N(BCV), DY, JMAX, J_N)
         IF(BC_Y_S(BCV) /= BC_Y_N(BCV)) THEN
            Y_CONSTANT = .FALSE.
            J_S = J_S + 1
            IF(BC_J_S(BCV)/=UNDEFINED_I.OR.BC_J_N(BCV)/=UNDEFINED_I)THEN
               CALL LOCATION_CHECK (BC_J_S(BCV), J_S, BCV, 'BC - south')
               CALL LOCATION_CHECK (BC_J_N(BCV), J_N, BCV, 'BC - north')
            ENDIF
         ENDIF
         BC_J_S(BCV) = J_S
         BC_J_N(BCV) = J_N
      ELSE
         IF(BC_J_S(BCV) /= UNDEFINED_I) &
            CALL CALC_LOC (ZERO,DY,BC_J_S(BCV),BC_Y_S(BCV))
         IF(BC_J_N(BCV) /= UNDEFINED_I) &
            CALL CALC_LOC (ZERO,DY,BC_J_N(BCV),BC_Y_N(BCV))
         IF (BC_Y_S(BCV) /= BC_Y_N(BCV)) Y_CONSTANT = .FALSE.
      ENDIF

! If there is no variation in the J direction set indices to 1
      IF(NO_J) THEN
         BC_J_S(BCV) = 1
         BC_J_N(BCV) = 1
      ENDIF

      IF(BC_Z_B(BCV)/=UNDEFINED .AND. BC_Z_T(BCV)/=UNDEFINED) THEN
         CALL CALC_CELL (ZERO, BC_Z_B(BCV), DZ, KMAX, K_B)
         CALL CALC_CELL (ZERO, BC_Z_T(BCV), DZ, KMAX, K_T)
         IF(BC_Z_B(BCV) /= BC_Z_T(BCV)) THEN
            Z_CONSTANT = .FALSE.
            K_B = K_B + 1
            IF(BC_K_B(BCV)/=UNDEFINED_I.OR.BC_K_T(BCV)/=UNDEFINED_I)THEN
               CALL LOCATION_CHECK (BC_K_B(BCV), K_B, BCV, 'BC - bottom')
               CALL LOCATION_CHECK (BC_K_T(BCV), K_T, BCV, 'BC - top')
            ENDIF
         ENDIF
         BC_K_B(BCV) = K_B
         BC_K_T(BCV) = K_T
      ELSE
         IF(BC_K_B(BCV) /= UNDEFINED_I) &
            CALL CALC_LOC (ZERO,DZ,BC_K_B(BCV),BC_Z_B(BCV))
         IF(BC_K_T(BCV) /= UNDEFINED_I) &
            CALL CALC_LOC (ZERO,DZ,BC_K_T(BCV),BC_Z_T(BCV))
         IF(BC_Z_B(BCV) /= BC_Z_T(BCV)) Z_CONSTANT = .FALSE.
      ENDIF

! If there is no variation in the K direction set indices to 1
      IF(NO_K) THEN
         BC_K_B(BCV) = 1
         BC_K_T(BCV) = 1
      ENDIF

! Check whether the boundary is a plane parallel to one of the three
! coordinate planes
      IF(BC_X_W(BCV)/=UNDEFINED .AND. BC_Y_S(BCV)/=UNDEFINED .AND. &
         BC_Z_B(BCV)/=UNDEFINED) CALL CHECK_PLANE (X_CONSTANT, &
         Y_CONSTANT, Z_CONSTANT, BCV, 'BC')


! CHECK FOR VALID VALUES
      IER = 0
      IF(BC_I_W(BCV)<1 .OR. BC_I_W(BCV)>IMAX2) IER = 1
      IF(BC_I_E(BCV)<1 .OR. BC_I_E(BCV)>IMAX2) IER = 1
      IF(BC_J_S(BCV)<1 .OR. BC_J_S(BCV)>JMAX2) IER = 1
      IF(BC_J_N(BCV)<1 .OR. BC_J_N(BCV)>JMAX2) IER = 1
      IF(BC_K_B(BCV)<1 .OR. BC_K_B(BCV)>KMAX2) IER = 1
      IF(BC_K_T(BCV)<1 .OR. BC_K_T(BCV)>KMAX2) IER = 1
      IF(BC_K_B(BCV) > BC_K_T(BCV)) IER = 1
      IF(BC_J_S(BCV) > BC_J_N(BCV)) IER = 1
      IF(BC_I_W(BCV) > BC_I_E(BCV)) IER = 1

      IF(IER /= 0)THEN
         WRITE(ERR_MSG,1100) BCV,                                      &
            'X', BC_X_W(BCV), BC_X_E(BCV),'I',BC_I_W(BCV),BC_I_E(BCV), &
            'Y', BC_Y_S(BCV), BC_Y_N(BCV),'J',BC_J_S(BCV),BC_J_N(BCV), &
            'Z', BC_Z_B(BCV), BC_Z_T(BCV),'K',BC_K_B(BCV),BC_K_T(BCV)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Invalid location specified for BC ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),/         &
         'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_BC_GEOMETRY_FLOW
