!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_FILE                                              C
!  Purpose: open a file                                                C
!                                                                      C
!  Author: P. Nicoletti                               Date: 12-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OPEN_FILE(FILENAME, NB, IUNIT, EXT, FULL_NAME,        &
         OPEN_STAT, OPEN_ACCESS, OPEN_FORM, IRECL, IER)

      use cdist
      use compar

      IMPLICIT NONE

! Dummy Arguments
!---------------------------------------------------------------------//
! FILENAME (without extension)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
! File extension.
      CHARACTER(LEN=*), INTENT(IN) :: EXT
! FILENAME + EXTENSION
      CHARACTER(LEN=*), INTENT(INOUT) :: FULL_NAME
! File status (NEW, OLD, UNKNOWN)
      CHARACTER(LEN=*), INTENT(IN) :: OPEN_STAT
! File access method ('SEQUENTIAL', 'DIRECT')
      CHARACTER(LEN=*), INTENT(IN) :: OPEN_ACCESS
! Open form ('FORMATTED' or 'UNFORMATTED')
      CHARACTER(LEN=*)   OPEN_FORM
! Index to first blank character in FILENAME
      INTEGER, INTENT(IN) :: NB
! Unit number to open
      INTEGER, INTENT(IN) :: IUNIT
! Record length
      INTEGER, INTENT(IN) :: IRECL
! Integer Error index:
! 000 - no error
! 100 - NEW run with existing files in directory
! 101 - OLD run missing RES and/or SPx files
! 102 - Unknown OPEN_STAT
      INTEGER, INTENT(OUT) :: IER

! Local Variables
!---------------------------------------------------------------------//
! Logical used to store result of file INQUIRE
      LOGICAL :: FILE_EXISTS

! Logicals that determine if files should be index.
      LOGICAL :: RES_IDX  ! Index RES files
      LOGICAL :: SPX_IDX  ! Index SPx files
      LOGICAL :: USE_IDX  ! Use the IDX value

! Initialize the error flag.
      IER = 0

! Conditions for indexing the RES files for distributed IO.
      RES_IDX = (myPE .NE. PE_IO) .OR. (.NOT.bStart_with_one_RES)
! Conditions for indexing the SPX files for distributed IO.
      SPX_IDX = .TRUE.

! Flag for indexing files.
      USE_IDX = bDist_IO .AND. (                                       &
         (SPX_IDX .AND. (EXT(2:3) .EQ. 'SP')) .OR.                     &
         (RES_IDX .AND. (EXT(2:4) .EQ. 'RES')))

! Construct the file name.
      FULL_NAME = ''
      IF(USE_IDX)THEN
         WRITE(FULL_NAME,1000) FILENAME(1:NB-1), myPE, EXT(1:4)
      ELSE
         WRITE(FULL_NAME,1001) FILENAME(1:NB-1), EXT(1:4)
      ENDIF

! Check to see if the file already exists in the run directory.
      INQUIRE(FILE=trim(FULL_NAME),EXIST=FILE_EXISTS)

! NEW files should not be in the run directory.
      IF(FILE_EXISTS .AND. (OPEN_STAT == 'NEW')) THEN
         IER = 100; RETURN
! OLD files must be in the run directory.
      ELSEIF(.NOT. FILE_EXISTS .AND. OPEN_STAT .EQ. 'OLD') THEN
         IER = 101; RETURN
      ENDIF

! Open direct access files.
      IF (OPEN_ACCESS == 'DIRECT') THEN
         OPEN(CONVERT='BIG_ENDIAN',UNIT=IUNIT, FILE=trim(FULL_NAME), STATUS=OPEN_STAT,     &
            RECL=IRECL, ACCESS=OPEN_ACCESS, FORM=OPEN_FORM, IOSTAT=IER)
      ELSE
! No matter the status passed to the routine, the file is created as
! NEW if it doesn't exist in the run directory.
         IF(.NOT.FILE_EXISTS) THEN
            OPEN(CONVERT='BIG_ENDIAN',UNIT=IUNIT, FILE=trim(FULL_NAME), STATUS='NEW',       &
               ACCESS=OPEN_ACCESS, FORM=OPEN_FORM, IOSTAT=IER)
         ELSEIF(OPEN_STAT == 'REPLACE') THEN
            OPEN(CONVERT='BIG_ENDIAN',UNIT=IUNIT, FILE=trim(FULL_NAME), STATUS=OPEN_STAT,   &
               ACCESS=OPEN_ACCESS, FORM=OPEN_FORM, IOSTAT=IER)
         ELSEIF(OPEN_STAT == 'APPEND' .OR. OPEN_STAT == 'UNKNOWN') THEN
            OPEN(CONVERT='BIG_ENDIAN',UNIT=IUNIT, FILE=trim(FULL_NAME), STATUS='UNKNOWN',   &
               ACCESS=OPEN_ACCESS, FORM=OPEN_FORM, POSITION='APPEND',  &
               IOSTAT=IER)
         ELSE
            IER = 102
         ENDIF
      ENDIF

      RETURN

 1000 FORMAT(A,'_',I5.5,A4)
 1001 FORMAT(A,A4)

      END SUBROUTINE OPEN_FILE
