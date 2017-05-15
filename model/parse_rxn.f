!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: PARSE_RXN(LINE, LMAX)                                  C
!  Purpose: Parse input line                                           C
!                                                                      C
!  Author: P. Nicoletti                               Date: 30-JUN-97  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: This routine was complete rewritten as part of the effort  C
!  to simplify reaction inputs in MFiX.
!  Author: J. Musser                                  Date: 01-Oct-12  C
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
      SUBROUTINE PARSE_RXN(LINE, lNoOfRxns, lName, lChemEq, lDH, lFDH)

!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE compar
      USE funits
      USE param
      USE param1
      USE parse

      IMPLICIT NONE

! Input line from mfix.dat.
      CHARACTER(len=*), INTENT(IN) :: LINE
! Array of reaction names.
      INTEGER, INTENT(INOUT) :: lNoOfRxns
! Array of reaction names.
      CHARACTER(len=*), INTENT(INOUT), DIMENSION(DIMENSION_RXN) :: lName
! Array of Chemical reaction equations.
      CHARACTER(len=*), INTENT(INOUT), DIMENSION(DIMENSION_RXN) :: lChemEq
! Array of User defined heat of reactions.
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(DIMENSION_RXN) :: lDH
! Array of User defined heat of reaction phase partitions.
      DOUBLE PRECISION, INTENT(INOUT), DIMENSION(DIMENSION_RXN, 0:DIM_M) :: lFDH


!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      CHARACTER, PARAMETER :: CT_BEG = '{'
      CHARACTER, PARAMETER :: CT_END = '}'

! Positions of braces {...}
      INTEGER bIDX, eIDX
! Reaction Name
      CHARACTER(LEN=128) :: INPUT
! Index of reaction.
      INTEGER IDX

! Copy line to input for processing.
      INPUT = TRIM(ADJUSTL(LINE))

! Look for the start and end of a reaction construct by checking for
! left and right braces.
      bIDX = INDEX(INPUT,CT_BEG)
      eIDX = INDEX(INPUT,CT_END)

! If not already inside/reading from a reaction construct, check to see
! if this is the start of a construct.
      IF(.NOT.IN_CONSTRUCT) THEN

! An identifier for the end of a construct was found.
         IF(eIDX .GT. 0) THEN
            IF(eIDX .GT. bIDX) THEN
! The reaction construct is specified in a single line.
! rxn_name { A + B --> AB }
               IDX = getReactionIndex(lNoOfRxns, 'NEW')
! Pull off the reaction construct name.
               CALL getName(INPUT,(bIDX-1), lNAME(IDX))
! Process the rest of the line.
               IF(isFracDH(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE(*, 1002) 'FracDH', trim(adjustl(INPUT))
                  WRITE(*, 1000)
                  CALL MFIX_EXIT(myPE)
               ELSEIF(isDH(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (*, 1002) 'DH', trim(adjustl(INPUT))
                  WRITE(*, 1000)
                  CALL MFIX_EXIT(myPE)
               ELSEIF(.NOT.isChemEq(INPUT(bIDX+1:eIDX-1)))THEN
                  WRITE (*, 1003) trim(adjustl(INPUT))
                  WRITE(*, 1000)
                  CALL MFIX_EXIT(myPE)
               ENDIF
               CALL getChemEq(INPUT(bIDX+1:eIDX-1), lChemEq(IDX))
            ELSE
! The format given in the deck file is incorrect. Brace mismatch.
               WRITE(*, 1001) trim(adjustl(LINE))
               WRITE(*, 1000)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSE
! This is the start of a reaction construct.
            IF(bIDX .GT. 0) THEN
! Get the reaction index.
               IDX = getReactionIndex(lNoOfRxns, 'NEW')
! Extract the reaction name.
               CALL getName(INPUT,(bIDX-1), lNAME(IDX))
! Process any data.
               IF(LEN_TRIM(ADJUSTL(INPUT(bIDX+1:eIDX-1))) .GT. 0)      &
                 CALL readConstruct(INPUT(bIDX+1:eIDX-1),              &
                    lChemEq(IDX), lDH(IDX), lFDH(IDX,:))
               IN_CONSTRUCT = .TRUE.
            ELSE
! Format Error.
               WRITE(*, 1004) trim(adjustl(INPUT))
               WRITE(*, 1000)
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
      ELSE

         IF(bIDX .GT. 0) THEN
! Format Error.
            WRITE(*, 1005) trim(adjustl(INPUT))
            WRITE(*, 1000)
            CALL MFIX_EXIT(myPE)
! This is the last line of the reaction construct which may or may not
! contain additional data.
         ELSEIF(eIDX .GT. 0) THEN
           IDX = getReactionIndex(lNoOfRxns)
           CALL readConstruct(INPUT(bIDX+1:eIDX-1), lChemEq(IDX),      &
              lDH(IDX), lFDH(IDX,:))
            IN_CONSTRUCT = .FALSE.

! Reading from somewhere inside of a reaction construct.
         ELSE
           IDX = getReactionIndex(lNoOfRxns)
           CALL readConstruct(INPUT(bIDX+1:), lChemEq(IDX),            &
              lDH(IDX), lFDH(IDX,:))
         ENDIF
      ENDIF

      RETURN

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1001: Mismatch of braces "{...}" in reaction ',       &
         ' construct.',//' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1002: Input format error in reaction construct.',     &
         ' Opening and',/' closing braces were found on the same line',&
         ' along with the',/' keyword ',A,'.',/' Single line',         &
         ' constructs can only contain a chemical equation.',//        &
         ' INPUT: ',A,//                                               &
         ' Example 1: RXN_NAME { chem_eq = "A + B --> AB" }',//        &
         ' Example 2: RXN_NAME {',/14X,'chem_eq = "A + B --> AB"',/14X,&
         'DH = 2.5d4',/14X,'fracDH(0) = 1.0',/12X,'}')

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1003: Input format error in reaction construct.',     &
         ' Opening and',/' closing braces were found on the same line',&
         ' and chem_eq was NOT found.',/' Single line constructs can', &
         ' only contain a chemical equation.',//' INPUT: ',A,//        &
         ' Example 1: RXN_NAME { chem_eq = "A + B --> AB" }',//        &
         ' Example 2: RXN_NAME {',/14X,'chem_eq = "A + B --> AB"',/14X,&
         'DH = 2.5d4',/14X,'fracDH(0) = 1.0',/12X,'}')

 1004 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1004: Data within the reaction block was identified', &
         ' outside of a',/' reaction construct. ',//' INPUT: ',A)

 1005 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1005: The start of a new reaction construct was',     &
         ' found before the',/' closing of the previous construct.',// &
         ' INPUT: ',A)

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: getReactionIndex()                                   !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      INTEGER FUNCTION getReactionIndex(lNOR, STAT)

      use rxns

      IMPLICIT NONE
! Number of reactions.
      INTEGER, INTENT(INOUT) :: lNOR
! Status
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: STAT

      IF(.NOT.PRESENT(STAT)) THEN
         getReactionIndex = lNOR

      ELSE
         IF(STAT == 'NEW') THEN
! Increment the number of reactions processed from the data file and
! return the new value as the index of the reactoin being processed.
            lNOR = lNOR + 1
            getReactionIndex = lNOR
         ELSE
            WRITE(*,*) ' Unknown status'
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      RETURN
      END FUNCTION getReactionIndex



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: readConstruct(IN, ChemEq, uDH, uFDH)                 !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE readConstruct(IN, ChemEq, uDH, uFDH)

      IMPLICIT NONE

! Input string being parsed.
      CHARACTER(len=*), INTENT(IN) :: IN
! Chemical equation.
      CHARACTER(len=*), INTENT(OUT) :: ChemEq
! User defined heat of reaction.
      DOUBLE PRECISION, INTENT(OUT) :: uDH
! User defined splitting of heat of reaction
      DOUBLE PRECISION, INTENT(OUT) :: uFDH(0:DIM_M)

! The input line contains no additional data.
      IF(LEN_TRIM(ADJUSTL(IN)) == 0) RETURN

! The input contains chemical equation data.
      IF(MORE_ChemEq .OR. isChemEq(IN)) THEN
         CALL getChemEq(IN, ChemEq)
! The input contains heat of reaction parsing data.
      ELSEIF(isFracDH(IN)) THEN
         CALL getFracDH(IN, uFDH(:))
! The input contains heat of reaction data.
      ELSEIF(isDH(IN)) THEN
         CALL getDH(IN, uDH)
! The entry doesn't match any of the keywords.
      ELSE
! Unidentified keyword.
         WRITE(*, 1001) trim(adjustl(IN))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF

      RETURN

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> readConstruct',/       &
         ' Error 1001: Unidentified keyword in reaction construct.'//, &
         ' INPUT: ',A)

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

      END SUBROUTINE readConstruct



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isChemEq(INPUT)                                      !
!                                                                      !
!  Purpose: Checks if the line contains the chemical Eq.               !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isChemEq(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"CHEM_EQ") == 0) THEN
! 'CHEM_EQ' was not found. This line does not contains a chemical eq.
         isChemEq = .FALSE.
      ELSE
! 'CHEM_EQ' was found. This line contains all or part of a chemical eq.
         isChemEq = .TRUE.
      ENDIF

      END FUNCTION isChemEq


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isDH(INPUT)                                          !
!                                                                      !
!  Purpose: Checks if the line contains user defined heat of reaction. !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isDH(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"DH") == 0) THEN
! 'DH' was not found. This line does not contains a heat of reaction.
         isDH = .FALSE.
      ELSE
! 'DH' was found. This line contains the heat of reaction
         isDH = .TRUE.
      ENDIF

      END FUNCTION isDH

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: isDH(INPUT)                                          !
!                                                                      !
!  Purpose: Checks if the line contains user defined heat of reaction. !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      LOGICAL FUNCTION isFracDH(INPUT)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT

! Check to see if the line contains 'END'
      IF (INDEX(LINE(1:),"FRACDH") == 0) THEN
! 'FRACDH' was not found.
         isFracDH = .FALSE.
      ELSE
! 'FRACDH' was found.
         isFracDH = .TRUE.
      ENDIF

      END FUNCTION isFracDH

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: get_ChemEq(INPUT, lNAME, IER)                      !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getName(INPUT, rPOS, lNAME)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! End of search location for reaction name.
      INTEGER, INTENT(IN) :: rPOS
! Name of reaction pulled from input.
      CHARACTER(LEN=32) , INTENT(OUT) :: lNAME

      INTEGER NAME_LEN

! Initialize the return value.
      lNAME = ''
! Verify that the name is not too long. This should be caught by
! preprocessing of the data file. However, if the user changed the
! reaction name after compiling (an error check for later) this check
! prevents and overflow error.
      NAME_LEN = len_trim(adjustl(INPUT(1:rPOS)))
      IF(NAME_LEN .GT. 32) THEN
         WRITE(*, 1001) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
! Verify that the name was not deleted after compiling.
! prevents and overflow error.
      ELSEIF(NAME_LEN .EQ. 0) THEN
         WRITE(*, 1002) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ELSE
         lNAME = trim(adjustl(INPUT(1:rPOS)))
      ENDIF

! There shouldn't be any crazy characters at this point because the
! code should fail to compile if the reaction names are not defined
! or contain invalid characters.

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> get_ChemEq',/          &
         ' Error 1001: Reaction name too long! Reaaction names are',   &
         ' limited to 32',/' characters.',//' Reaction Name: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> get_ChemEq',/          &
         ' Error 1002: Unable to determine reaction name.',//          &
         ' INPUT: ',A)

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

      RETURN
      END SUBROUTINE getName


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getDH(INPUT, lDH)                                  !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getDH(INPUT, lDH)

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Name of reaction pulled from input.
      DOUBLE PRECISION, INTENT(OUT) :: lDH

      INTEGER lQ
      INTEGER lLMAX
! read/write output status
      INTEGER IOS

      lLMAX = LEN_TRIM(INPUT)

      IF(INDEX(INPUT,"DH") .EQ. 0) THEN
         WRITE (*, 1100) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF

      lQ = INDEX(INPUT(:lLMAX),'=')

      IF(lQ .EQ. 0) THEN
         WRITE (*, 1001) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Convert the entrying into an double precision value.
      READ(INPUT(lQ+1:),*,IOSTAT=IOS) lDH
      IF(IOS .NE. 0) THEN
         WRITE(*, 1002) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF


 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getDH',/               &
         ' Error 1001: Input format error for DH.',//' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getDH',/               &
         ' Error 1002: Unable to determine DH value from input.',/     &
         ' Cannot convert specified value to double precision value.',/&
         /' INPUT: ',A)

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1105: DH was initially located within the input line',&
         /' however its location cannot be determined.',&
         ' INPUT: ',A)

      RETURN
      END SUBROUTINE getDH


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getfracDH(INPUT, lChemEq)                          !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getFracDH(INPUT, lFracDH)

      USE param
      USE param1

      IMPLICIT NONE

! Input line.
      CHARACTER(len=*), INTENT(IN) :: INPUT
! Name of reaction pulled from input.
      DOUBLE PRECISION, INTENT(OUT) :: lFracDH(0:DIM_M)


      INTEGER POS, lP, rP, lQ
      INTEGER lLMAX
! read/write output status
      INTEGER IOS
! Phase Index
      INTEGER pIDX

      lLMAX = LEN_TRIM(INPUT)
      POS = INDEX(INPUT,"FRACDH")

      IF(POS == 0) THEN
         WRITE (*, 1100) trim(adjustl(INPUT))
         WRITE(*, 1100)
         CALL MFIX_EXIT(myPE)
      ENDIF

      lP = INDEX(INPUT(:lLMAX),'(')
      rP = INDEX(INPUT(:lLMAX),')')
      lQ = INDEX(INPUT(:lLMAX),'=')

      IF(lP .EQ. rP .AND. lP .EQ. ZERO) THEN
         WRITE(*, 1001) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ELSEIF(lP .GE. rP) THEN
         WRITE(*, 1002) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ELSEIF(rP .GE. lQ) THEN
         WRITE(*, 1002) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF
! Convert the entrying into an integer value.
      READ(INPUT(lP+1:rP-1),*,IOSTAT=IOS) pIDX
      IF(IOS .NE. 0) THEN
         WRITE(*, 1003) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ELSEIF(pIDX .LT. 0 .OR. pIDX .GT. DIM_M)THEN
         WRITE(*, 1004) trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF

! Convert the entrying into an double precision value.
      READ(INPUT(lQ+1:),*,IOSTAT=IOS) lFracDH(pIDX)
      IF(IOS .NE. 0) THEN
         WRITE(*, 1005)trim(adjustl(INPUT))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF

      RETURN

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1001: Unable to determine phase association for',     &
         ' fracDH. When',/' specifying heat of reaction (DH), the',    &
         ' fraction of DH assigned to',/' each phase must be',         &
         ' given explicitly.',//' Example: fracDH(0) = 0.25  ! 25% of',&
         ' DH is assigned to gas phase',/'          fracDH(1) = 0.75 ',&
         ' ! 75% of DH is assigned to solids phase 1',//' Note:',      &
         ' fracDH(0) + fracDH(1) + ... + frachDH(M) == 1.0',//         &
         ' INPUT: ',A)

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1002: Input format error for fracDH.',//              &
         ' INPUT: ',A)

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1003: Unable to determine phase index for fracDH.',//&
         ' INPUT: ',A)

 1004 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1004: Phase index for fracDH exceeds DIM_M!',//       &
         ' INPUT: ',A)

 1005 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getFracDH',/           &
         ' Error 1005: Unable to determine fracDH value from input.',/ &
         ' Cannot convert specified value to double precision value.',/&
         /' INPUT: ',A)


 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

 1100 FORMAT(//1X,70('*')/' From: PARSE_RXN',/                         &
         ' Error 1105: fracDH was initially located within the',       &
         ' input line,',/' however its location cannot be determined.',&
         ' INPUT: ',A)



      END SUBROUTINE getFracDH

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine name: getChemEq(INPUT, lChemEq)                          !
!                                                                      !
!  Purpose: Extract the reaction name from a construct.                !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE getChemEq(IN, lChemEq)

      IMPLICIT NONE

! Input line.
      CHARACTER(len=*), INTENT(IN) :: IN
! Name of reaction pulled from input.
      CHARACTER(len=*), INTENT(OUT) :: lChemEq

! read/write output status
      INTEGER IOS

      INTEGER POS, lPOS, rPOS, ldP, lsP, aPOS
      INTEGER lLMAX

! Chemical equations start with the keyword: CHEM_EQ. If this is not a
! continuation of the previous line, search for the keyword and
! flag an error if not found.
      IF(.NOT.MORE_ChemEq) THEN
         lLMAX = LEN_TRIM(IN)
         POS = INDEX(IN,"CHEM_EQ")

         IF(POS == 0) THEN
            WRITE (*, 1105) 'Chem_Eq'
            CALL MFIX_EXIT(myPE)
         ENDIF
! Initialize
         lChemEq = ''
! Update POS to skip over the keyword: CHEM_EQ
         POS = POS+7
      ELSE
         POS = 1
      ENDIF

! Search for quote marks bounding the chemical equation.
      ldP = POS + INDEX(IN(POS:),'"')  ! double quote "
      lsP = POS + INDEX(IN(POS:),"'")  ! single quote '

      IF(ldP .GT. POS .AND. lsP .EQ. POS) THEN
! The chemical equation is bounded by double quotes
         lPOS = ldP
! Search for the second quote mark.
         rPOS = lPOS + INDEX(IN(lPOS+1:),'"')
      ELSEIF(ldP .EQ. POS .AND. lsP .GT. POS) THEN
! The chemical equation is bounded by single quotes
         lPOS = lsP
! Search for the second quote mark.
         rPOS = lPOS + INDEX(IN(lPOS+1:),"'")
      ELSE
! Different errors are thrown depending if this is a continuation
! (MORE_ChemEq) or the start of a chemical equation.
         IF(.NOT.MORE_ChemEq) THEN
            WRITE(*, 1001) trim(adjustl(IN))
            WRITE(*, 1000)
         ELSE
            IF(isFracDH(IN)) THEN
              WRITE(*, 1002) 'Keyword fracDH was found inside',        &
                 ' the chemical equation!', trim(adjustl(IN))
              WRITE(*, 1000)
            ELSEIF(isDH(IN)) THEN
              WRITE(*, 1002) 'Keyword DH was found inside',            &
                 ' the chemical equation!', trim(adjustl(IN))
              WRITE(*, 1000)
            ELSE
              WRITE(*, 1002) 'Unbalanced or missing parentheses', '',  &
                 trim(adjustl(IN))
              WRITE(*, 1000)
            ENDIF
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF

! Mismatch/Unbalanced parentheses
      IF(lPOS .EQ. rPOS) THEN
! Different errors are thrown depending if this is a continuation
! (MORE_ChemEq) or the start of a chemical equation.
         IF(.NOT.MORE_ChemEq) THEN
            WRITE(*, 1001) trim(adjustl(IN))
            WRITE(*, 1000)
         ELSE
           WRITE(*, 1002) 'Unbalanced or missing parentheses', '',  &
              trim(adjustl(IN))
           WRITE(*, 1000)
         ENDIF
         CALL MFIX_EXIT(myPE)
      ENDIF

! Search for an ampersand.
      aPOS = lPOS + INDEX(IN(lPOS+1:),'&')
! An ampersand was found.
      IF(aPOS .GT. lPOS) THEN
         MORE_ChemEq = .TRUE.
! The ampersand should be further to the right than the last quote mark.
         IF(aPOS .LE. rPOS) THEN
            WRITE(*, 1003) trim(adjustl(IN))
            WRITE(*, 1000)
            CALL MFIX_EXIT(myPE)
         ENDIF
      ELSE
         MORE_ChemEq = .FALSE.
      ENDIF

! Store the chemical equation.
      WRITE(lChemEq,"(A,1X,A)",IOSTAT=IOS) trim(lChemEq), &
         trim(adjustl(IN(lPOS:rPOS-1)))
      IF(IOS .NE. 0) THEN
         WRITE(*, 1004) trim(lChemEq), trim(adjustl(IN))
         WRITE(*, 1000)
         CALL MFIX_EXIT(myPE)
      ENDIF

      RETURN

 1001 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1001: Unbalanced or missing parentheses for chem_eq.',&
         //' INPUT: ',A,//' Example 1: RXN_NAME { chem_eq = ',         &
         '"A + B --> AB" }',//' Example 2: RXN_NAME {',/14X,           &
         'chem_eq = "A + B --> AB"',/14X, 'DH = 2.5d4',/14X,           &
         'fracDH(0) = 1.0',/12X,'}')

 1002 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1002: Chemical equation continuation input error.',   &
         //'  > ',2A//' INPUT: ',A)

 1003 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1003: Input format error for chem_eq. An amperand',   &
         ' (&)',/' was located within the parentheses.',//' INPUT: ',A)

 1004 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1004: Unable to process chemical equation input.',/   &
         ' A possible error is variable overflow as the total length', &
         ' is limited',/' to 512 characters.',//' lChemEq: ',A,//      &
         ' INPUT: ',A)

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

 1105 FORMAT(//1X,70('*')/' From: PARSE_RXN --> getChemEq',/           &
         ' Error 1105: chem_eq was initially located within the',      &
         ' input line,',/' however its location cannot be determined.',&
         ' INPUT: ',A)

      END SUBROUTINE getChemEq

      END SUBROUTINE PARSE_RXN
