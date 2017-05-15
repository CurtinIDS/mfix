!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_USR0                                             C
!  Purpose: Write initial part of user-defined output                  C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_USR0

      use compar, only: myPE, PE_IO

      IMPLICIT NONE

      IF(myPE /= PE_IO) RETURN

      CALL WRITE_DAT_HEADER('POST_POS.dat','Pos')
      CALL WRITE_DAT_HEADER('POST_VEL.dat','Vel')

      RETURN

      CONTAINS

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_DAT_HEADER(FNAME, VAR)

      use discretelement, only: DES_ONEWAY_COUPLED
      use discretelement, only: DES_EXPLICITLY_COUPLED
      use particle_filter, only: DES_INTERP_ON
      use particle_filter, only: DES_INTERP_SCHEME
      use particle_filter, only: DES_INTERP_MEAN_FIELDS
      use particle_filter, only: DES_INTERP_WIDTH
      use particle_filter, only: DES_DIFFUSE_MEAN_FIELDS
      use particle_filter, only: DES_DIFFUSE_WIDTH

      use run, only: DESCRIPTION

      use param1, only: UNDEFINED

      IMPLICIT NONE

      CHARACTER(len=*) :: FNAME
      CHARACTER(len=*) :: VAR

! logical used for testing is the data file already exists
      LOGICAL :: EXISTS
! file unit for heat transfer data
      INTEGER, PARAMETER :: fUNIT = 2030

      INQUIRE(FILE=FNAME,EXIST=EXISTS)
      IF (.NOT.EXISTS) THEN
         OPEN(UNIT=fUNIT,FILE=FNAME,STATUS='NEW')
         WRITE(fUNIT, 1000) trim(DESCRIPTION)
      ELSE
         OPEN(UNIT=fUNIT,FILE=FNAME,POSITION="APPEND",STATUS='OLD')
      ENDIF

      WRITE(fUNIT, 1100) DES_ONEWAY_COUPLED

      WRITE(fUNIT, 1110, ADVANCE='NO') DES_INTERP_ON
      IF(DES_INTERP_ON) THEN
         WRITE(fUNIT, 1111, ADVANCE='YES') DES_INTERP_SCHEME
      ELSE
         WRITE(fUNIT, *) '   '
      ENDIF

      WRITE(fUNIT, 1120, ADVANCE='NO') DES_INTERP_MEAN_FIELDS
      IF(DES_INTERP_WIDTH /= UNDEFINED) THEN
         WRITE(fUNIT, 1121, ADVANCE='YES') DES_INTERP_WIDTH
      ELSE
         WRITE(fUNIT, *) '   '
      ENDIF

      WRITE(fUNIT, 1130, ADVANCE='NO') DES_DIFFUSE_MEAN_FIELDS
      IF(DES_DIFFUSE_MEAN_FIELDS) THEN
         WRITE(fUNIT, 1131, ADVANCE='YES') DES_DIFFUSE_WIDTH
      ELSE
         WRITE(fUNIT, *) '   '
      ENDIF

      WRITE(fUNIT, 1140, ADVANCE='YES') DES_EXPLICITLY_COUPLED

      WRITE(fUNIT, 1250) VAR, VAR

 1000 FORMAT(2/,25x,A)

 1100 FORMAT(2/,7x,'DES_ONEWAY_COUPLED =',6x,L1)

 1110 FORMAT(7x,'DES_INTERP_ON =',11x,L1)
 1111 FORMAT(5x,'DES_INTERP_SCHEME = ',A)

 1120 FORMAT(7x,'DES_INTERP_MEAN_FIELDS =',2x,L1)
 1121 FORMAT(5x,'DES_INTERP_WIDTH =  ',F9.6)

 1130 FORMAT(7x,'DES_DIFFUSE_MEAN_FIELDS =',1x,L1)
 1131 FORMAT(5x,'DES_DIFFUSE_WIDTH = ',F7.2)

 1140 FORMAT(7x,'DES_EXPLICITLY_COUPLED =',2x,L1)

 1250 FORMAT(/10X,'Time',17X,A,13X,A,'_MFIX',9X,'%REL DIFF')

      CLOSE(fUNIT)
      RETURN
      END SUBROUTINE WRITE_DAT_HEADER

      END SUBROUTINE WRITE_USR0
