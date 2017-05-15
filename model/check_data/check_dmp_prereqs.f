!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_DMP_PREREQS                                       !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_DMP_PREREQS


! Global Variables:
!---------------------------------------------------------------------//
! Number of ranks.
      use compar, only: numPEs
! DMP grid partitioning data:
      use compar, only: NODESI  ! Partitions along x-axis
      use compar, only: NODESJ  ! Partitions along y-axis
      use compar, only: NODESK  ! Partitions along z-axis
! Domain partitions in various directions.
      use geometry, only: IMAX
      use geometry, only: JMAX
      use geometry, only: KMAX

! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

! Local Variables:
!---------------------------------------------------------------------//
! The approximate number of domain partitions that will be assigned
! to each process for DMP runs. (eg IMAX/NODESI)
      INTEGER :: LAYERS

! Local Parameters:
!---------------------------------------------------------------------//
! The minimum number of computational cell layers required.
      INTEGER, PARAMETER :: DMP_MIN = 3

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_DMP_PREREQS")

! Verify that DMP partitioning information is provided given if
! there is more than one rank.
      IF( numPEs > 1 ) then
         IF(NODESI .EQ. UNDEFINED_I .AND.                              &
            NODESJ .EQ. UNDEFINED_I .AND.                              &
            NODESK .EQ. UNDEFINED_I) THEN
            WRITE(ERR_MSG,1000)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

! Initialize NODE values if undefined. If this is a DMP run, then a
! warning message is passed to the user.
      IF (NODESI .EQ. UNDEFINED_I) THEN
         IF(numPEs > 1) THEN
            WRITE(ERR_MSG,1001)'I','I'
            CALL FLUSH_ERR_MSG
         ENDIF
         NODESI = 1
! Verify that the DMP partition is appropriate for the domain.
      ELSEIF(NODESI > 1) THEN
         LAYERS=int(IMAX/NODESI)
         IF(LAYERS < DMP_MIN) THEN
            WRITE(ERR_MSG,1002) 'X', DMP_MIN, 'I', 'I', LAYERS
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      IF (NODESJ .EQ. UNDEFINED_I) THEN
         IF(numPEs > 1) THEN
            WRITE(ERR_MSG,1001)'J','J'
            CALL FLUSH_ERR_MSG
         ENDIF
         NODESJ = 1

! Verify that the DMP partition is appropriate for the domain.
      ELSEIF(NODESJ > 1) THEN
         LAYERS=int(JMAX/NODESJ)
         IF(LAYERS < DMP_MIN) THEN
            WRITE(ERR_MSG,1002) 'Y', DMP_MIN, 'J', 'J', LAYERS
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      IF (NODESK .EQ. UNDEFINED_I) THEN
         IF(numPEs > 1) THEN
            WRITE(ERR_MSG,1001)'K','K'
            CALL FLUSH_ERR_MSG
         ENDIF
         NODESK = 1
! Verify that the DMP partition is appropriate for the domain.
      ELSEIF(NODESK > 1) THEN
         LAYERS=int(KMAX/NODESK)
         IF(LAYERS < DMP_MIN) THEN
            WRITE(ERR_MSG,1002) 'Z', DMP_MIN, 'K', 'K', LAYERS
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

! Verify that the number of requested processes (munPEs) matches
! the domain decomposition (nodesi * nodesj * nodesk).
      IF(numPEs .NE. (NODESi*NODESj*NODESk)) THEN
         WRITE(ERR_MSG,1003) numPEs, (NODESi*NODESj*NODESk)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


      CALL FINL_ERR_MSG

      RETURN


 1000 FORMAT('Error 1000: No DMP grid partitioning data provided in ', &
         'mfix.dat.',/'NODESI, NODESJ, and NODESK are all undefined.',/&
         'Refer to the users manual for required input and make the ', &
         'necessary',/'corrections to the input data file.')

 1001 FORMAT('Warning 1001: Setting NODES',A1,' to default: ',         &
         'NODES',A1,'=1.')

 1002 FORMAT('Error 1002:  Too many DMP partitions specified for ',    &
         A1,' axis.',/'There must be at least ',I2,' computational ',  &
         'cells per DMP parition.',/' >>> Computational Cells/DMP ',   &
         'Partition = int(',A1,'MAX/NODES',A1,') = ',I2,/'Refer to ',  &
         'the users manual for required input and make the necessary',/&
         'corrections to the input data file.')

 1003 FORMAT('Error 1003: The number of requested processors is ',     &
         'inconsistent',/'with the domain decomposition, (NODESi * ',  &
         'NODESj * NODESk).',/' These numbers must match.',2/,         &
         '  Number of requested processes: ',I8,/,                     &
         '  Domain decomposition : ',I8,/)


      END SUBROUTINE CHECK_DMP_PREREQS
