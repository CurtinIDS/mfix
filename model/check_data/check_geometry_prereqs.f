!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_GEOMETRY_PREREQS                                  !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GEOMETRY_PREREQS



! Global Variables:
!---------------------------------------------------------------------//
! Domain partitions in various directions.
      use geometry, only: IMAX, NO_I, XMIN
      use geometry, only: JMAX, NO_J
      use geometry, only: KMAX, NO_K, DZ, ZLENGTH

! Runtime flag specifying 2D simulations
!      use geometry, only: NO_K

      use geometry, only: COORDINATES, CYLINDRICAL
      use geometry, only: CYCLIC_X, CYCLIC_X_PD
!      use geometry, only: COORDINATES

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: ONE, ZERO, UNDEFINED_I, UNDEFINED

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      use toleranc

      implicit none

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GEOMETRY_PREREQS")

! Verify that the domain decomposition was specified.
      IF(IMAX == UNDEFINED_I .OR. JMAX == UNDEFINED_I .OR.             &
         (.NOT.NO_K .AND. KMAX == UNDEFINED_I) ) THEN
         WRITE(ERR_MSG,1000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1000 FORMAT('Error 1000: IMAX or JMAX or KMAX not specified in ',     &
          'mfix.dat')

! If no variation in a direction is considered, the number of cells in
! that direction should be 1
      IF(NO_I) THEN
         WRITE(ERR_MSG, 1100) 'I','I','east and west'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(NO_J) THEN
         WRITE(ERR_MSG, 1100) 'J','J','north and south'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Illegal geometry: NO_',A1,' is disabled. ',  &
         'The same functionality',/'is achieved with one cell (',A1,   &
         'MAX=1) and making the ',A,' walls',/'free-slip. Please ',    &
         'correct the mfix.dat file.')

      IF (XMIN < ZERO) THEN
         WRITE(ERR_MSG, 1101) 'XMIN'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1101 FORMAT('Error 1101: Illegal geometry: ',A,' cannot be less ',    &
         'than zero.',/'Please correct the mfix.dat file.')


      SELECT CASE(trim(COORDINATES))
      CASE ('CYLINDRICAL')
         CYLINDRICAL = .TRUE.
         IF(CYCLIC_X .OR. CYCLIC_X_PD) THEN
            WRITE(ERR_MSG, 1102)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1102 FORMAT('Error 1102: X-axis cannot be CYCLIC in cylindrical ',    &
         'coordinates',/'Please correct the mfix.dat file.')

      CASE ('CARTESIAN')
         CYLINDRICAL = .FALSE.

      CASE DEFAULT
         WRITE(ERR_MSG, 1103)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1103 FORMAT('Error 1103: Unknown COORDINATES specified. Please ',     &
         'correct the ',/'mfix.dat file.')

      END SELECT


      IF(NO_K) THEN
         IF(KMAX == UNDEFINED_I) THEN
            KMAX = 1
         ELSEIF(KMAX /= 1) THEN
            WRITE(ERR_MSG, 1110) 'KMAX','NO_K'
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1110 FORMAT('Error 1110: Illegal geometry: ',A,' must remain ',       &
         'UNDEFINED_I or 1 when',/A,' is TRUE. Please correct the ',   &
         'mfix.dat file.')

         IF(DZ(1)==UNDEFINED) THEN
            IF(ZLENGTH==UNDEFINED) THEN
               IF(CYLINDRICAL) THEN
                  DZ(1) = 8.*ATAN(ONE)
                  ZLENGTH = 8.*ATAN(ONE)
               ELSE
                  DZ(1) = ONE
                  ZLENGTH = ONE
               ENDIF
            ELSE
               DZ(1) = ZLENGTH
            ENDIF
         ELSE
            IF(ZLENGTH==UNDEFINED) THEN
               ZLENGTH = DZ(1)
            ELSE
               IF(.NOT.COMPARE(ZLENGTH,DZ(1)))THEN
                  WRITE(ERR_MSG, 1111) 'DZ(1) and ZLENGTH'
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF

  1111 FORMAT('Error 1111: Illegal geometry: ',A,' are not equal.',/   &
           'Please correct the mfix.dat file.')

            ENDIF
         ENDIF
      ENDIF

      CALL FINL_ERR_MSG

      RETURN



      END SUBROUTINE CHECK_GEOMETRY_PREREQS
