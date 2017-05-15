!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_GEOMETRY                                          !
!  Purpose: Check the distributed parallel namelist variables.         !
!                                                                      !
!  Author: P. Nicoletti                               Date: 14-DEC-99  !
!  Reviewer: J.Musser                                 Date: 16-Jan-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GEOMETRY(SHIFT)

         use check_data_cg, only: get_dxyz_from_control_points

! Global Variables:
!---------------------------------------------------------------------//
! Domain partitions in various directions.
      use geometry, only: DX, XLENGTH
      use geometry, only: DY, YLENGTH
      use geometry, only: DZ, ZLENGTH

      use geometry, only: NO_I, IMIN1, IMAX, IMAX1, IMAX3
      use geometry, only: NO_J, JMIN1, JMAX, JMAX1, JMAX3
      use geometry, only: NO_K, KMIN1, KMAX, KMAX1, KMAX3

! Runtime flag specifying 2D simulations
!      use geometry, only: NO_K

      use geometry, only: CYLINDRICAL
      use geometry, only: CYCLIC_X, CYCLIC_X_PD
      use geometry, only: CYCLIC_Y, CYCLIC_Y_PD
      use geometry, only: CYCLIC_Z, CYCLIC_Z_PD
!      use geometry, only: COORDINATES

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none


      LOGICAL, intent(IN) :: SHIFT
      LOGICAL, external :: COMPARE

! Local Variables:
!---------------------------------------------------------------------//


! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GEOMETRY")

      CALL GET_DXYZ_FROM_CONTROL_POINTS

      CALL CHECK_AXIS(IMAX, IMAX3, XLENGTH, DX, 'X', 'I', NO_I, SHIFT)
      CALL CHECK_AXIS(JMAX, JMAX3, YLENGTH, DY, 'Y', 'J', NO_J, SHIFT)
      CALL CHECK_AXIS(KMAX, KMAX3, ZLENGTH, DZ, 'Z', 'K', NO_K, SHIFT)

      IF(SHIFT) CALL SHIFT_DXYZ

!  Ensure that the cell sizes across cyclic boundaries are comparable
      IF(CYCLIC_X .OR. CYCLIC_X_PD) THEN
         IF(DX(IMIN1) /= DX(IMAX1)) THEN
            WRITE(ERR_MSG,1100) 'DX(IMIN1)',DX(IMIN1),'DX(IMAX1)',DX(IMAX1)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      IF(CYCLIC_Y .OR. CYCLIC_Y_PD) THEN
         IF(DY(JMIN1) /= DY(JMAX1)) THEN
            WRITE(ERR_MSG,1100) 'DY(JMIN1)',DY(JMIN1),'DY(JMAX1)',DY(JMAX1)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

      IF(CYCLIC_Z .OR. CYCLIC_Z_PD .OR. CYLINDRICAL) THEN
         IF (DZ(KMIN1) /= DZ(KMAX1)) THEN
            WRITE(ERR_MSG,1100) 'DZ(KMIN1)',DZ(KMIN1),'DZ(KMAX1)',DZ(KMAX1)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF
      ENDIF

 1100 FORMAT('Error 1100: Cells adjacent to cyclic boundaries must ',  &
         'be of same size:',/2X,A,' = ',G12.5,/2x,A,' = ',G12.5,/      &
         'Please correct the mfix.dat file.')


      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_GEOMETRY


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_GEOMETRY_DES                                      !
!  Author: Pradeep Gopalakrishnan                     Date:    Nov-11  !
!                                                                      !
!  Purpose: Checks the des grid input parameters.                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_GEOMETRY_DES

! Global Variables:
!---------------------------------------------------------------------//
! Domain partition for DEM background mesh.
      use discretelement, only: DESGRIDSEARCH_IMAX
      use discretelement, only: DESGRIDSEARCH_JMAX
      use discretelement, only: DESGRIDSEARCH_KMAX
! Domain size specified by the user.
      use geometry, only: XLENGTH, YLENGTH, ZLENGTH, NO_K
! Maximum particle size.
      use discretelement, only: MAX_RADIUS


! Global Parameters:
!---------------------------------------------------------------------//
      use param1, only: UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager


      implicit none


! Local Variables:
!---------------------------------------------------------------------//
! Maximum particle diameter.
      DOUBLE PRECISION :: MAX_DIAM
! Calculated cell dimension based on particle size
      DOUBLE PRECISION :: WIDTH
!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_GEOMETRY_DES")

! Calculate the max particle diameter and cell width.
      MAX_DIAM = 2.0d0*MAX_RADIUS
      WIDTH = 3.0d0*(max_diam)

! Calculate and/or verify the grid in the X-axial direction.
      IF(DESGRIDSEARCH_IMAX == UNDEFINED_I) THEN
         DESGRIDSEARCH_IMAX = max(int(XLENGTH/WIDTH), 1)
      ELSEIF((XLENGTH/dble(DESGRIDSEARCH_IMAX)) < MAX_DIAM) THEN
         WRITE(ERR_MSG, 1100) 'X', MAX_DIAM,                           &
            XLENGTH/dble(DESGRIDSEARCH_IMAX)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Calculate and/or verify the grid in the Y-axial direction.
      IF(DESGRIDSEARCH_JMAX == UNDEFINED_I) THEN
         DESGRIDSEARCH_JMAX = max(int(YLENGTH/WIDTH), 1)
      ELSEIF((YLENGTH/dble(DESGRIDSEARCH_JMAX)) < MAX_DIAM) THEN
         WRITE(ERR_MSG, 1100) 'Y', MAX_DIAM,                           &
            YLENGTH/dble(DESGRIDSEARCH_JMAX)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

! Calculate and/or verify the grid in the Z-axial direction.
      IF(NO_K) THEN
         DESGRIDSEARCH_KMAX = 1
      ELSEIF(DESGRIDSEARCH_KMAX == UNDEFINED_I) THEN
         DESGRIDSEARCH_KMAX = max(int(ZLENGTH/WIDTH), 1)
      ELSEIF((ZLENGTH/dble(DESGRIDSEARCH_KMAX)) < MAX_DIAM) THEN
         WRITE(ERR_MSG, 1100) 'Z', MAX_DIAM,                           &
            ZLENGTH/dble(DESGRIDSEARCH_KMAX)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      CALL FINL_ERR_MSG

 1100 FORMAT('Error 1100: The des search grid is too fine in the ',A1, &
         '-direction. The',/'maximum particle diameter is larger than',&
         ' the cell width:',/2x,'MAX DIAM:   ',g12.5,/2x,'CELL ',      &
         'WIDTH: ',g12.5,/'Decrease the values for DESGRIDSEARCH in ', &
         'the mfix.dat file.')

      RETURN
      END SUBROUTINE CHECK_GEOMETRY_DES
