!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_PROGRESS_BAR                                     C
!  Purpose: Displays a progress bar on the screen                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 30-JAN-09  C
!  Reviewer:                                          Date: **-***-**  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_PROGRESS_BAR(I,I_MAX,JUSTIFICATION)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------

      USE compar
      USE exit, only: mfix_exit
      USE fldvar
      USE funits
      USE mpi_utility
      USE param
      USE param1
      USE physprop
      USE progress_bar
      USE run
      USE rxns
      USE scalars
      IMPLICIT NONE

      INTEGER :: I,I_MAX,ISKIP,PROGRESS
      INTEGER :: P,P1,P2
      CHARACTER (LEN=9) :: TEXT
      CHARACTER (LEN=4) :: BAR_STATUS
      CHARACTER (LEN=BAR_WIDTH) :: PROGRESSBAR
      CHARACTER (LEN=1) :: JUSTIFICATION
      DOUBLE PRECISION :: PERCENT

      IF(.NOT.PRINT_PROGRESS_BAR)  RETURN

      IF(myPE /= PE_IO) RETURN

      ISKIP = INT(BAR_RESOLUTION * 0.01 *FLOAT(I_MAX))

      IF((MOD(I,ISKIP)/=0).AND.(I/=I_MAX)) RETURN

      CALL ERASE_PROGRESS_BAR(BAR_WIDTH,BAR_STATUS,JUSTIFICATION)

      BAR_STATUS =''
      IF((BAR_WIDTH<10).OR.(BAR_WIDTH>80)) RETURN

      PERCENT  = FLOAT(I)/FLOAT(I_MAX) * 100.0
      PROGRESS = INT(0.01*PERCENT * BAR_WIDTH)  + 1

      WRITE(TEXT,10) PERCENT
10    FORMAT(' ',F5.1,' % ')

      DO P = 1, PROGRESS
         PROGRESSBAR(P:P)= BAR_CHAR
      ENDDO

      DO P = PROGRESS+1,BAR_WIDTH
         PROGRESSBAR(P:P)= ' '
      ENDDO

      SELECT CASE(JUSTIFICATION)
         CASE('L')

            WRITE(*,15,ADVANCE='NO')TEXT,'|',PROGRESSBAR,'|'

         CASE('C')

            P1 = BAR_WIDTH / 2 - 3
            P2 = BAR_WIDTH / 2 + 5

            PROGRESSBAR(P1:P2)= TEXT

            WRITE(*,20,ADVANCE='NO')'|',PROGRESSBAR,'|'

         CASE('R')

            WRITE(*,15,ADVANCE='NO')'|',PROGRESSBAR,'|',TEXT

         CASE('N')

            WRITE(*,20,ADVANCE='NO')'|',PROGRESSBAR,'|'

         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: WRITE_PROGRESS_BAR.'
            WRITE(*,*)'INCORRECT JUSTIFICATION DESCRIPTOR:',JUSTIFICATION
            WRITE(*,*)'ACCEPTABLE VALUES ARE:'
            WRITE(*,*)'L : LEFT JUSTIFICATION'
            WRITE(*,*)'C : CENTER JUSTIFICATION'
            WRITE(*,*)'R : RIGHT JUSTIFICATION'
            WRITE(*,*)'N : NO TEXT'
            call mfix_exit(myPE)
      END SELECT

      IF(PERCENT>=100.0) THEN
         BAR_STATUS ='DONE'
         WRITE(*,*)
      ENDIF

15    FORMAT(A,A,A,A)
20    FORMAT(A,A,A)


      RETURN
      END SUBROUTINE WRITE_PROGRESS_BAR

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ERASE_PROGRESS_BAR                                     C
!  Purpose: Erases a progress bar on the screen                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 30-JAN-09  C
!  Reviewer:                                          Date: **-***-**  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE ERASE_PROGRESS_BAR(BAR_WIDTH,BAR_STATUS,JUSTIFICATION)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE exit, only: mfix_exit
      USE fldvar
      USE funits
      USE mpi_utility
      USE parallel
      USE param
      USE param1
      USE physprop
      USE run
      USE rxns
      USE scalars
      USE sendrecv
      IMPLICIT NONE

      INTEGER :: I,BAR_WIDTH,NERASE
      CHARACTER (LEN=4) :: BAR_STATUS
      CHARACTER (LEN=1) :: JUSTIFICATION

      IF(myPE /= PE_IO) RETURN

      IF(BAR_STATUS=='DONE') THEN
         BAR_STATUS = ''
         RETURN
      ENDIF

      IF((BAR_WIDTH<10).OR.(BAR_WIDTH>80)) RETURN


      SELECT CASE(JUSTIFICATION)
         CASE('L')
            NERASE = BAR_WIDTH + 11
         CASE('C')
            NERASE = BAR_WIDTH + 2
         CASE('R')
            NERASE = BAR_WIDTH + 11
         CASE('N')
            NERASE = BAR_WIDTH + 2
         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: WRITE_PROGRESS_BAR.'
            WRITE(*,*)'INCORRECT JUSTIFICATION DESCRIPTOR:',JUSTIFICATION
            WRITE(*,*)'ACCEPTABLE VALUES ARE:'
            WRITE(*,*)'L : LEFT JUSTIFICATION'
            WRITE(*,*)'C : CENTER JUSTIFICATION'
            WRITE(*,*)'R : RIGHT JUSTIFICATION'
            WRITE(*,*)'N : NO TEXT'
            call mfix_exit(myPE)
      END SELECT

      DO I = 1,NERASE
         WRITE(*,10,ADVANCE='NO')CHAR(8)
      ENDDO
10    FORMAT(A)

      RETURN
      END SUBROUTINE ERASE_PROGRESS_BAR


