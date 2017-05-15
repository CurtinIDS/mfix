!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: WRITE_USR0                                             !
!  Purpose: Write initial part of user-defined output                  !
!                                                                      !
!  Author:                                            Date: dd-mmm-yy  !
!  Reviewer:                                          Date: dd-mmm-yy  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE WRITE_USR0

      use compar, only: myPE, PE_IO

      IMPLICIT NONE

      IF(myPE /= PE_IO) RETURN

      CALL WRITE_DAT_HEADER('POST_POS1.dat','Pos')
      CALL WRITE_DAT_HEADER('POST_POS2.dat','Pos')

      RETURN

      CONTAINS

!----------------------------------------------------------------------!
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE WRITE_DAT_HEADER(FNAME, VAR)

      use run, only: DESCRIPTION
      use discretelement, only: DES_INTG_METHOD
      use discretelement, only: KN, KN_W
      use discretelement, only: DES_EN_INPUT, DES_EN_WALL_INPUT 

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

      WRITE(fUNIT, 1100) trim(adjustl(DES_INTG_METHOD))

      WRITE(fUNIT, 1110) KN, KN_W
      WRITE(fUNIT, 1120) DES_EN_WALL_INPUT(1), DES_EN_WALL_INPUT(1) 

      WRITE(fUNIT, 1200) VAR, VAR

 1000 FORMAT(2/,25x,A)

 1100 FORMAT(2/,7x,'Time Stepping Scheme: ',A)

 1110 FORMAT(/7x,'Normal collision spring sonstant. (N/m)',/&
         10x,'KN = ',T30,G12.4,/&
         10x,'KN_W = ',T30,G12.4)

 1120 FORMAT(/7x,'Restitution coefficient. (1)',/&
         10x,'DES_EN_INPUT = ',T30,G12.4,/&
         10x,'DES_EN_WALL_INPUT = ',T30,G12.4)

 1200 FORMAT(2/11X,'Time',17X,A,12X,A,'_MFIX',10X,'%REL ERR')


      CLOSE(fUNIT)
      RETURN
      END SUBROUTINE WRITE_DAT_HEADER

      END SUBROUTINE WRITE_USR0
