!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARSE_LINE(LINE, LMAX, RXN_FLAG, READ_FLAG)            C
!  Purpose: Parse input line                                           C
!                                                                      C
!  Author: M. Syamlal                                 Date: 27-JUN-97  C
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
!
      SUBROUTINE PARSE_LINE(LINE, LMAX, RXN_FLAG, READ_FLAG)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE des_rxns
      USE param
      USE param1
      USE parse
      USE rxns

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

! Line to be parsed.
      CHARACTER(len=*), INTENT(IN) :: LINE

! Length of of LINE.
      INTEGER, INTENT(IN) :: LMAX

! Indicate whether currently reading chemical reaction data. The
! namelist read is skipped when reading a chemical reaction.
      LOGICAL, INTENT(OUT) :: RXN_FLAG

! Indicate whether to do a namelist read on the line. A namelist read
! is still preformed when an arithmetic operation is found.
      LOGICAL, INTENT(OUT) :: READ_FLAG

! Start and end locations for the search parameters.
      INTEGER LSTART, LEND


! The string is empty. No need to parse.
      IF (LMAX == 0) THEN
         READ_FLAG = .FALSE.
         RETURN
      ENDIF

! Check to see if the input line contains '@('. If this string is found,
! then the line contains information to parsed. This string indicates
! that one of two actions need to occur"
! 1) there is an expression to evalute; @(6.0/2.0) = 3.0
! 2) this is the start of a reaction block; @(RXNS...
      LSTART = INDEX(LINE,START_STR)

! If the returned index is not zero, the input line contain the string.
      IF (LSTART /= 0) THEN

! Look for the ending parenthesis. If none exists, flag the error and
! exit MFiX.
         LEND = LSTART - 1 + INDEX(LINE(LSTART:LMAX),END_STR)
         IF (LEND <= LSTART) THEN
            WRITE (*, 1000) myPE,LINE(LSTART:LMAX)
            CALL MFIX_EXIT(myPE)
         ENDIF

! Check to see if this is the end of a reaction block.
         IF (END_RXN(LINE(LSTART:LEND),LEND-LSTART)) THEN
! This is the end of a reaction block, but either no reaction block
! initializer '@(RXNS)' preceded it, or the preceding reaction block
! was already closed by another '@(END) statement.
            IF(.NOT.RXN_FLAG) THEN
               WRITE (*, 1010)
               CALL MFiX_EXIT(0)
            ENDIF

! Set flags indicating that no additional rate information will be
! processed.
            RXN_FLAG = .FALSE.
            READ_FLAG = .FALSE.
            CALL END_PARSE_RXN()
            RETURN
         ENDIF

! Check to see if this is the start of a reaction block.
         IF (START_DES_RXN(LINE(LSTART:LEND),LEND-LSTART)) THEN
            DES_RXN = .TRUE.
            RXN_FLAG = .TRUE.
            READ_FLAG = .FALSE.
! Initialize logicals for parsing reaction data.
            CALL INIT_PARSE_DES_RXN()
            RETURN
         ELSEIF(START_RXN(LINE(LSTART:LEND),LEND-LSTART)) THEN
            TFM_RXN = .TRUE.
            RXN_FLAG = .TRUE.
            READ_FLAG = .FALSE.
! Initialize logicals for parsing reaction data.
            CALL INIT_PARSE_RXN()
            RETURN
         ENDIF
      ENDIF ! IF (LSTART /= 0) THEN


      IF(TFM_RXN) THEN
         CALL PARSE_RXN (LINE, NO_OF_RXNS, RXN_NAME, RXN_CHEM_EQ,      &
            usrDH, usrfDH)
         READ_FLAG = .FALSE.
         RETURN
      ELSEIF(DES_RXN) THEN
         CALL PARSE_RXN (LINE, NO_OF_DES_RXNS, DES_RXN_NAME,           &
            DES_RXN_CHEM_EQ, DES_usrDH, DES_usrfDH)
         READ_FLAG = .FALSE.
         RETURN
      ENDIF
!
      LSTART = INDEX(LINE,START_STR)             !Arithmetic processing ?
      IF (LSTART /= 0) CALL PARSE_ARITH (LINE, LMAX)
      READ_FLAG = .TRUE.
!
      RETURN

 1000 FORMAT(//1X,70('*')/' (PE ',I6,'): From: PARSE_LINE',/&
         ' Message: An evaluation statement "@(" was found in ',&
         'the input line,',/' but no ending parenthesis was located:',/&
         ' INPUT: ',A,/1X,70('*')//)

 1010 FORMAT(/1X,70('*')/': From: PARSE_LINE',/&
         ' Error: END keyword before a start keyword in line: ',       &
          /1X,A,/1X,70('*')/)


      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: START_RXN(LINE, LMAX)                                !
!                                                                      !
!  Purpose: Returns a value of TRUE if this is the start of a reaction !
!           block. Otherwise, the return value is FALSE.               !
!                                                                      !
!  Author: M. Syamlal                                 Date: 27-JUN-97  !
!                                                                      !
!  Revision Number: 1                                                  !
!  Author: J. Musser                                  Date: 13-SPT-12  !
!  Reviewer:                                          Date: dd-mmm-yy  !
!                                                                      !
!  Literature/Document References: None                                !
!                                                                      !
!  Variables referenced: RXN_BLK - string indicating a reaction block  !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION START_RXN (LINE, LMAX)

! Input line containing an '@(' statment.
      CHARACTER(len=*), INTENT(IN) :: LINE
! Length of of LINE.
      INTEGER LMAX

! Check to see if the line contains 'RXNS'
      IF (INDEX(LINE(1:LMAX),RXN_BLK) == 0) THEN
! 'RXNS' was not found. This is not the start of a reaction block.
         START_RXN = .FALSE.
      ELSE
! 'RXNS' was found. This is the start of a reaction block.
         START_RXN = .TRUE.
      ENDIF

      RETURN
      END FUNCTION START_RXN

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: START_RXN(LINE, LMAX)                                !
!                                                                      !
!  Purpose: Returns a value of TRUE if this is the start of a reaction !
!           block. Otherwise, the return value is FALSE.               !
!                                                                      !
!  Author: J. Musser                                  Date: 31-Oct-12  !
!                                                                      !
!  Reviewer:                                          Date: dd-mmm-yy  !
!                                                                      !
!  Literature/Document References: None                                !
!                                                                      !
!  Variables referenced: RXN_BLK - string indicating a reaction block  !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION START_DES_RXN (LINE, LMAX)

! Input line containing an '@(' statment.
      CHARACTER(len=*), INTENT(IN) :: LINE
! Length of of LINE.
      INTEGER LMAX

! Check to see if the line contains 'RXNS'
      IF (INDEX(LINE(1:LMAX),DES_RXN_BLK) == 0) THEN
! 'RXNS' was not found. This is not the start of a reaction block.
         START_DES_RXN = .FALSE.
      ELSE
! 'RXNS' was found. This is the start of a reaction block.
         START_DES_RXN = .TRUE.
      ENDIF

      RETURN
      END FUNCTION START_DES_RXN

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: END_RXN(LINE, LMAX)                                  !
!                                                                      !
!  Purpose: Check for the end of rxn block                             !
!                                                                      !
!  Author: M. Syamlal                                 Date: 27-JUN-97  !
!                                                                      !
!  Revision Number: 1                                                  !
!  Purpose: Add additional comments.                                   !
!  Author:                                            Date: dd-mmm-yy  !
!  Reviewer:                                          Date: dd-mmm-yy  !
!                                                                      !
!  Literature/Document References: None                                !
!                                                                      !
!  Variables referenced:                                               !
!                                                                      !
!   - END_BLK - string indicating the end of a reaction block.         !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION END_RXN (LINE, LMAX)

! Input line containing an '@(' statment.
      CHARACTER(len=*), INTENT(IN) :: LINE
! Length of of LINE.
      INTEGER LMAX

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:LMAX),END_BLK) == 0) THEN
! 'END' was not found. This is not the end of a reaction block.
         END_RXN = .FALSE.
      ELSE
! 'END' was found. This is the end of a reaction block.
         END_RXN = .TRUE.
      ENDIF
!
      RETURN
      END FUNCTION END_RXN


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: INIT_PARSE_RXN()                                     !
!                                                                      !
!  Purpose: Initialize variables for the reaction parser.              !
!                                                                      !
!  Author: J. Musser                                  Date: 14-SPT-12  !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified:                                                 !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE INIT_PARSE_RXN()

! Allocate the necessary storage arrays for chemical reaction data
! read from the data file. These arrays are 'allocatable' so that after
! processing in CHECK_DATA_09, they can be deallocated as they are no
! longer necessary.
!-----------------------------------------------------------------------
! Reaction Names: Allocate/Initialize
      IF(Allocated( RXN_NAME )) GoTo 100
      Allocate( RXN_NAME( DIMENSION_RXN ))
      RXN_NAME(:) = ''
! Chemical Equations: Allocate/Initialize
      IF(Allocated( RXN_CHEM_EQ )) GoTo 100
      Allocate( RXN_CHEM_EQ( DIMENSION_RXN ))
      RXN_CHEM_EQ(:) = ''
! User defined heat of reaction: Allocate/Initialize
      IF(Allocated( usrDH )) GoTo 100
      Allocate( usrDH( DIMENSION_RXN ))
      usrDH(:) = UNDEFINED
! User defined heat of reaction partitions: Allocate/Initialize
      IF(Allocated( usrfDH )) GoTo 100
      Allocate( usrfDH( DIMENSION_RXN, 0:DIM_M ))
      usrfDH(:,:) = UNDEFINED
! Logical indicating that the code is in the middle of parsing a
! reaction construct.
      IN_CONSTRUCT = .FALSE.
! Number of reactions found in data file.
      NO_OF_RXNS = 0

! Flag indicating that the chemical equation is specified over
! multiple lines.
      MORE_ChemEq = .FALSE.

      RETURN

 100  WRITE(*,1001)
      WRITE(*,1000)
      CALL MFiX_EXIT(0)

 1001 FORMAT(/1X,70('*')/' From: PARSE_LINE --> INIT_PARSE_RXN',/      &
         ' Error 1001: More than one reaction block has been located!',&
         ' A data file',/' can only contain one reaction block',       &
         ' [@(RXNS)...@(END)].'/)

 1000 FORMAT(' Please refer to the Readme file for chemical equation', &
         ' input formats',/' and correct the data file.',/1X,70('*')/)

      END SUBROUTINE INIT_PARSE_RXN


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: INIT_PARSE_DES_RXN()                                 !
!                                                                      !
!  Purpose: Initialize variables for the DES reaction parser.          !
!                                                                      !
!  Author: J. Musser                                  Date: 14-SPT-12  !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified:                                                 !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE INIT_PARSE_DES_RXN()

! Allocate the necessary storage arrays for chemical reaction data
! read from the data file. These arrays are 'allocatable' so that after
! processing in CHECK_DATA_09, they can be deallocated as they are no
! longer necessary.
!-----------------------------------------------------------------------
! Reaction Names: Allocate/Initialize
      IF(Allocated( DES_RXN_NAME )) GoTo 100
      Allocate( DES_RXN_NAME( DIMENSION_RXN ))
      DES_RXN_NAME(:) = ''
! Chemical Equations: Allocate/Initialize
      IF(Allocated( DES_RXN_CHEM_EQ )) GoTo 100
      Allocate( DES_RXN_CHEM_EQ( DIMENSION_RXN ))
      DES_RXN_CHEM_EQ(:) = ''
! User defined heat of reaction: Allocate/Initialize
      IF(Allocated( DES_usrDH )) GoTo 100
      Allocate( DES_usrDH( DIMENSION_RXN ))
      DES_usrDH(:) = UNDEFINED
! User defined heat of reaction partitions: Allocate/Initialize
      IF(Allocated( DES_usrfDH )) GoTo 100
      Allocate( DES_usrfDH( DIMENSION_RXN, 0:DIM_M ))
      DES_usrfDH(:,:) = UNDEFINED
! Logical indicating that the code is in the middle of parsing a
! reaction construct.
      IN_DES_CONSTRUCT = .FALSE.
! Number of reactions found in data file.
      NO_OF_DES_RXNS = 0

! Flag indicating that the chemical equation is specified over
! multiple lines.
      MORE_ChemEq = .FALSE.

      RETURN

 100  WRITE(*,1001)
      WRITE(*,1000)
      CALL MFiX_EXIT(0)

 1001 FORMAT(/1X,70('*')/' From: PARSE_LINE --> INIT_PARSE_DES_RXN',/  &
         ' Error 1001: More than one DES reaction block has been',     &
         ' located! A data',/' file can only contain one reaction',    &
         ' block [@(DES_RXNS)...@(END)].'/)

 1000 FORMAT(' Please refer to the Readme file for chemical equation', &
         ' input formats',/' and correct the data file.',/1X,70('*')/)

      END SUBROUTINE INIT_PARSE_DES_RXN


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: END_PARSE_RXN()                                      !
!                                                                      !
!  Purpose: Initialize variables for the reaction parser.              !
!                                                                      !
!  Author: J. Musser                                  Date: 14-SPT-12  !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified:                                                 !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE END_PARSE_RXN()

      READING_RXN = .FALSE.
      READING_RATE = .FALSE.
      DES_RXN = .FALSE.
      TFM_RXN = .FALSE.

      RETURN
      END SUBROUTINE END_PARSE_RXN


      END SUBROUTINE PARSE_LINE



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARSE_ARITH(LINE, LMAX)                                C
!  Purpose: Complete arithmetic operations and expand the line         C
!                                                                      C
!  Author: M. Syamlal                                 Date: 10-AUG-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
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
!
      SUBROUTINE PARSE_ARITH(LINE, LMAX)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE param
      USE param1
      USE parse
      USE utilities, ONLY: seek_end
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!C
!                      The part of LINE containing input
      INTEGER LMAX
!                      Input line with arithmetic operations.  Out put
!                      line with completed arithmetic statements.
!
      CHARACTER LINE*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                      Value of pi
      DOUBLE PRECISION PI
!
!                      Cumulative value and sub value
      DOUBLE PRECISION VALUE, SUB_VALUE
!
!                      Start and end locations for the arithmetic operation
      INTEGER          LSTART, LEND
!
!                      Length of arithmetic operation string
      INTEGER          LENGTH
!
!                      22 - LENGTH
      INTEGER          LDIF
!
!                      Locations in SUB_STR, and LINE
      INTEGER          LSUB, L
!
!                      Operator symbol (Legal values: *, /)
      CHARACTER(LEN=1) :: OPERATION
!
!                      Substring taken from LINE
      CHARACTER(LEN=80) :: SUB_STR
!
!-----------------------------------------------
!
!
      PI = 4.0D0*ATAN(ONE)
!
!  Search for arithmetic operation
!
   10 CONTINUE
      LMAX = SEEK_END(LINE,LEN(LINE))
!
      LSTART = INDEX(LINE,START_STR)
!
      IF (LSTART == 0) RETURN
!
      LEND = LSTART - 1 + INDEX(LINE(LSTART:LMAX),END_STR)
      IF (LEND <= LSTART) THEN
         WRITE (*, 1000) myPE,LINE(LSTART:LMAX)
         CALL MFIX_EXIT(myPE)
      ENDIF
!
!    Do the arithmetic
!
      VALUE = ONE
      OPERATION = '*'
      LSUB = 1
      DO L = LSTART + 2, LEND
         IF (LINE(L:L)=='*' .OR. LINE(L:L)=='/' .OR. LINE(L:L)==END_STR) THEN
            IF (LSUB == 1) THEN
               WRITE (*, 1015) myPE,LINE(LSTART:LEND)
               CALL MFIX_EXIT(myPE)
            ENDIF
            IF (SUB_STR(1:LSUB-1) == 'PI') THEN
               SUB_VALUE = PI
            ELSE
               READ (SUB_STR(1:LSUB-1), *, ERR=900) SUB_VALUE
            ENDIF
            IF (OPERATION == '*') THEN
               VALUE = VALUE*SUB_VALUE
            ELSE IF (OPERATION == '/') THEN
               VALUE = VALUE/SUB_VALUE
            ENDIF
            LSUB = 1
            OPERATION = LINE(L:L)
         ELSE IF (LINE(L:L) == ' ') THEN
         ELSE
            SUB_STR(LSUB:LSUB) = LINE(L:L)
            LSUB = LSUB + 1
         ENDIF
      END DO
      LENGTH = LEND - LSTART + 1
      IF (LENGTH > 22) THEN
         DO L = LSTART + 22, LEND
            LINE(L:L) = ' '
         END DO
      ELSE IF (LENGTH < 22) THEN
         LMAX = SEEK_END(LINE,LEN(LINE))
         LDIF = 22 - LENGTH
         IF (LMAX + LDIF > LEN(LINE)) THEN
            WRITE (*, 1020) myPE,LINE(1:80)
            CALL MFIX_EXIT(myPE)
         ENDIF
         DO L = LMAX, LEND + 1, -1
            LINE(L+LDIF:L+LDIF) = LINE(L:L)
         END DO
      ENDIF
!
!  Transfer the value to LINE
!
      WRITE (SUB_STR, '(G22.15)') VALUE
      L = LSTART
      DO LSUB = 1, 22
         LINE(L:L) = SUB_STR(LSUB:LSUB)
         L = L + 1
      END DO
      GO TO 10
!
  900 CONTINUE
      WRITE (*, 1010) myPE, SUB_STR(1:LSUB-1)
      CALL MFIX_EXIT(myPE)
 1000 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: No ending ) found in the input line: ',/9X,A,/1X,70('*')/)
 1010 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: Error reading the input string: ',/9X,A,/1X,70('*')/)
 1015 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: Invalid operator in the input string: ',/9X,A,/1X,70('*')/)
 1020 FORMAT(/1X,70('*')//'(PE ',I6,'): From: PARSE_ARITH',/&
         ' Message: Too many arithmetic operations in the line: ',/1X,A,/1X,70(&
         '*')/)
      END SUBROUTINE PARSE_ARITH

!// Comments on the modifications for DMP version implementation
!// 001 Include header file and common declarations for parallelization
!//PAR_I/O added myPE stamp in output
