!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  SUBROUTINE: CHECK_SOLIDS_MPPIC                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 02-FEB-14  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_SOLIDS_MPPIC


! Global Variables:
!---------------------------------------------------------------------//
! Domain partitions in various directions.
!      use geometry, only: IMAX
!      use geometry, only: JMAX
!      use geometry, only: KMAX
! Runtime flag specifying 2D simulations
!      use geometry, only: NO_K


      USE param1
      USE geometry
      USE funits
      USE discretelement
      USE constant
      USE physprop
      USE fldvar
      USE toleranc
      USE mfix_pic
      USE cutcell
      USE functions

      USE mpi_utility


! Global Parameters:
!---------------------------------------------------------------------//
!      use param1, only: UNDEFINED_I

! Use the error manager for posting error messages.
!---------------------------------------------------------------------//
      use error_manager

      implicit none

!-----------------------------------------------

! Initialize the error manager.
      CALL INIT_ERR_MSG("CHECK_SOLIDS_MPPIC")




      IF(MPPIC_COEFF_EN1 == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) 'MPPIC_COEFF_EN1'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(MPPIC_COEFF_EN1 > ONE .OR.                                &
         MPPIC_COEFF_EN1 < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'MPPIC_COEFF_EN1',                       &
            trim(iVal(MPPIC_COEFF_EN1))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(MPPIC_COEFF_EN2 == UNDEFINED) THEN
         WRITE(ERR_MSG, 1000) 'MPPIC_COEFF_EN2'
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

      ELSEIF(MPPIC_COEFF_EN2 > ONE .OR.                                &
         MPPIC_COEFF_EN2 < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'MPPIC_COEFF_EN2',                       &
            trim(iVal(MPPIC_COEFF_EN2))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(MPPIC_COEFF_EN_WALL > ONE .OR.                                &
         MPPIC_COEFF_EN_WALL < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'MPPIC_COEFF_EN_WALL',                   &
            trim(iVal(MPPIC_COEFF_EN_WALL))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

      IF(MPPIC_COEFF_ET_WALL > ONE .OR.                                &
         MPPIC_COEFF_ET_WALL < ZERO) THEN
         WRITE(ERR_MSG, 1001) 'MPPIC_COEFF_ET_WALL',                   &
            trim(iVal(MPPIC_COEFF_ET_WALL))
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF


 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
         'correct the mfix.dat file.')

 1001 FORMAT('Error 1001: Illegal or unknown input: ',A,' = ',A,/   &
         'Please correct the mfix.dat file.')


      CALL FINL_ERR_MSG

      RETURN

      END SUBROUTINE CHECK_SOLIDS_MPPIC
