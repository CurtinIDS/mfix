!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: GET_IS                                                  !
!  Author: M. Syamlal                                 Date: 21-OCT-92  !
!                                                                      !
!  Purpose: Find and validate i, j, k locations for IS's               !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE GET_IS(ISV)

      USE param
      USE param1
      USE geometry
      USE is
      USE indices
      USE funits
      USE compar

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Loop/variable indices
      INTEGER, INTENT(in) :: ISV

! Local Variables:
!---------------------------------------------------------------------//
! Error flag.
      INTEGER :: IER
! Calculated indices of the wall boundary
      INTEGER :: I_w , I_e , J_s , J_n , K_b , K_t
! Surface indictors
      LOGICAL :: X_CONSTANT, Y_CONSTANT, Z_CONSTANT
!......................................................................!

      CALL INIT_ERR_MSG('GET_IS')

      X_CONSTANT = .TRUE.
      Y_CONSTANT = .TRUE.
      Z_CONSTANT = .TRUE.

      IF(IS_X_W(ISV)/=UNDEFINED .AND. IS_X_E(ISV)/=UNDEFINED) THEN
         CALL CALC_CELL(XMIN, IS_X_W(ISV), DX, IMAX, I_W)
         CALL CALC_CELL(XMIN, IS_X_E(ISV), DX, IMAX, I_E)
         IF (IS_X_W(ISV) /= IS_X_E(ISV)) THEN
            X_CONSTANT = .FALSE.
            I_W = I_W + 1
            IF(IS_I_W(ISV)/=UNDEFINED_I .OR.                           &
               IS_I_E(ISV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK(IS_I_W(ISV), I_W, ISV, 'IS - west')
               CALL LOCATION_CHECK(IS_I_E(ISV), I_E, ISV, 'IS - east')
            ENDIF
         ENDIF
         IS_I_W(ISV) = I_W
         IS_I_E(ISV) = I_E
      ELSE
         IF(IS_I_W(ISV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(XMIN, DX, IS_I_W(ISV), IS_X_W(ISV))
         IF(IS_I_E(ISV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(XMIN, DX, IS_I_E(ISV), IS_X_E(ISV))
         IF (IS_X_W(ISV) /= IS_X_E(ISV)) X_CONSTANT = .FALSE.
      ENDIF

!  If there is no variation in the I direction set indices to 1
      IF (NO_I) THEN
         IS_I_W(ISV) = 1
         IS_I_E(ISV) = 1
      ENDIF
!
      IF (IS_Y_S(ISV)/=UNDEFINED .AND. IS_Y_N(ISV)/=UNDEFINED) THEN
         CALL CALC_CELL(ZERO, IS_Y_S(ISV), DY, JMAX, J_S)
         CALL CALC_CELL(ZERO, IS_Y_N(ISV), DY, JMAX, J_N)
         IF (IS_Y_S(ISV) /= IS_Y_N(ISV)) THEN
            Y_CONSTANT = .FALSE.
            J_S = J_S + 1
            IF(IS_J_S(ISV)/=UNDEFINED_I .OR.                           &
               IS_J_N(ISV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK(IS_J_S(ISV), J_S, ISV, 'IS - south')
               CALL LOCATION_CHECK(IS_J_N(ISV), J_N, ISV, 'IS - north')
            ENDIF
         ENDIF
         IS_J_S(ISV) = J_S
         IS_J_N(ISV) = J_N
      ELSE
         IF(IS_J_S(ISV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(ZERO, DY, IS_J_S(ISV), IS_Y_S(ISV))
         IF(IS_J_N(ISV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(ZERO, DY, IS_J_N(ISV), IS_Y_N(ISV))
         IF (IS_Y_S(ISV) /= IS_Y_N(ISV)) Y_CONSTANT = .FALSE.
      ENDIF

! If there is no variation in the J direction set indices to 1
      IF (NO_J) THEN
         IS_J_S(ISV) = 1
         IS_J_N(ISV) = 1
      ENDIF

      IF (IS_Z_B(ISV)/=UNDEFINED .AND. IS_Z_T(ISV)/=UNDEFINED) THEN
         CALL CALC_CELL(ZERO, IS_Z_B(ISV), DZ, KMAX, K_B)
         CALL CALC_CELL(ZERO, IS_Z_T(ISV), DZ, KMAX, K_T)
         IF (IS_Z_B(ISV) /= IS_Z_T(ISV)) THEN
            Z_CONSTANT = .FALSE.
            K_B = K_B + 1
            IF (IS_K_B(ISV)/=UNDEFINED_I .OR.                          &
               IS_K_T(ISV)/=UNDEFINED_I) THEN
               CALL LOCATION_CHECK(IS_K_B(ISV), K_B, ISV, 'IS - bottom')
               CALL LOCATION_CHECK(IS_K_T(ISV), K_T, ISV, 'IS - top')
            ENDIF
         ENDIF
         IS_K_B(ISV) = K_B
         IS_K_T(ISV) = K_T
      ELSE
         IF(IS_K_B(ISV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(ZERO, DZ, IS_K_B(ISV), IS_Z_B(ISV))
         IF(IS_K_T(ISV) /= UNDEFINED_I)                                &
            CALL CALC_LOC(ZERO, DZ, IS_K_T(ISV), IS_Z_T(ISV))
         IF (IS_Z_B(ISV) /= IS_Z_T(ISV)) Z_CONSTANT = .FALSE.
      ENDIF

!  If there is no variation in the K direction set indices to 1
      IF (NO_K) THEN
         IS_K_B(ISV) = 1
         IS_K_T(ISV) = 1
      ENDIF

!  Check whether the boundary is a plane parallel to one of the three
!  coordinate planes, else check whether a direction is specified by IS_TYPE
      IF(X_CONSTANT .OR. Y_CONSTANT .OR. Z_CONSTANT) THEN
         IF(IS_X_W(ISV)/=UNDEFINED .AND. IS_Y_S(ISV)/=UNDEFINED .AND.  &
            IS_Z_B(ISV)/=UNDEFINED) CALL CHECK_PLANE(X_CONSTANT,       &
            Y_CONSTANT, Z_CONSTANT,ISV, 'IS')
      ELSE
         SELECT CASE(IS_TYPE(ISV)(1:1))
         CASE('X','Y','Z') ! Do Nothing
         CASE DEFAULT
            WRITE(ERR_MSG, 1100) trim(iVar('IS_TYPE',ISV)),            &
               trim(IS_TYPE(ISV))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         END SELECT
      ENDIF

 1100 FORMAT('Error 1100: ',A,' = ',A,' is specified as a volume',/    &
         'and therefore should have should have a directional prefix;',&
         /'Valid prefix values: X_, Y_, Z_',/'Please correct the ',/   &
         'mfix.dat file.')

! CHECK FOR VALID VALUES
      IER = 0
      IF(IS_I_W(ISV)<1 .OR. IS_I_W(ISV)>IMAX2) IER = 1
      IF(IS_I_E(ISV)<1 .OR. IS_I_E(ISV)>IMAX2) IER = 1
      IF(IS_J_S(ISV)<1 .OR. IS_J_S(ISV)>JMAX2) IER = 1
      IF(IS_J_N(ISV)<1 .OR. IS_J_N(ISV)>JMAX2) IER = 1
      IF(IS_K_B(ISV)<1 .OR. IS_K_B(ISV)>KMAX2) IER = 1
      IF(IS_K_T(ISV)<1 .OR. IS_K_T(ISV)>KMAX2) IER = 1
      IF(IS_K_B(ISV) > IS_K_T(ISV)) IER = 1
      IF(IS_J_S(ISV) > IS_J_N(ISV)) IER = 1
      IF(IS_I_W(ISV) > IS_I_E(ISV)) IER = 1

      IF(IER /= 0)THEN
         WRITE(ERR_MSG,1101) ISV,                                      &
            'X', IS_X_W(ISV), IS_X_E(ISV),'I',IS_I_W(ISV),IS_I_E(ISV), &
            'Y', IS_Y_S(ISV), IS_Y_N(ISV),'J',IS_J_S(ISV),IS_J_N(ISV), &
            'Z', IS_Z_B(ISV), IS_Z_T(ISV),'K',IS_K_B(ISV),IS_K_T(ISV)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1101 FORMAT('Error 1101: Invalid location specified for IS ',I3,'.',  &
         3(/3x,A1,': ',g12.5,',',g12.5,8x,A1,': ',I8,',',I8),/         &
         'Please correct the mfix.dat file.')

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GET_IS
