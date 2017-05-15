!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_OUT3                                             !
!  Author: M. Syamlal                                 Date: 10-JAN-92  !
!                                                                      !
!  Purpose: To write cpu and wall time used by the code                !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_OUT3(CPU, WALL, IO)

      use error_manager
      use run, only: get_tunit

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: CPU
      DOUBLE PRECISION, INTENT(INOUT) :: WALL
      DOUBLE PRECISION, INTENT(INOUT) :: IO
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------


      CHARACTER(len=4) :: UNIT_CPU
      CHARACTER(len=4) :: UNIT_WALL
      CHARACTER(len=4) :: UNIT_IO


      WRITE(ERR_MSG, "(2/1x,70('*'))")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      CALL GET_TUNIT(CPU, UNIT_CPU)
      WRITE(ERR_MSG, 1000) 'CPU', trim(iVal(CPU)), UNIT_CPU
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      CALL GET_TUNIT(IO, UNIT_IO)
      WRITE(ERR_MSG, 1000) 'CPU IO', trim(iVal(IO)), UNIT_IO
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      CALL GET_TUNIT(WALL, UNIT_WALL)
      WRITE(ERR_MSG, 1000) 'wall time', trim(iVal(WALL)), UNIT_WALL
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

      WRITE(ERR_MSG, "(1x,70('*'))")
      CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)

 1000 FORMAT(' Total ',A,' used = ',A,1x,A)

      RETURN
      END SUBROUTINE WRITE_OUT3
