!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: READ_SPX0(READ_SPX)                                    C
!  Purpose: Read the initial records (REAL)                            C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR  C
!                        ID_MINUTE, ID_SECOND                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC, VERSION                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE READ_SPX0(READ_SPX)
!
!
      USE funits
      USE machine
      USE param
      USE param1
      USE post3d
      USE run
      IMPLICIT NONE
!
      logical read_spx(*)
!                loop counter
      INTEGER    LC

!
!                file version ID
      CHARACTER(LEN=512) :: VERSION
!
      DO 100 LC = 1,N_SPX
        IF (READ_SPX(LC) .AND. SPX_OPEN(LC)) THEN
           READ (UNIT_SPX+LC,REC=1) VERSION
           READ(VERSION(6:512),*)VERSION_NUMBER
           READ (UNIT_SPX+LC,REC=2)RUN_NAME,ID_MONTH,ID_DAY,ID_YEAR,&
                ID_HOUR,ID_MINUTE,ID_SECOND
!
!  The first field contains the pointer to the next record.
!  The second field contains the number of records written each time step
!
           READ (UNIT_SPX+LC,REC=3) LAST_REC(LC),NUM_REC(LC)
        ENDIF
100   CONTINUE
      RETURN
      END
