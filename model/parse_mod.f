MODULE parse

   USE compar
   USE exit, only: mfix_exit
   USE funits
   USE param
   USE param1

      IMPLICIT NONE

! Strings indicating arithmetic operation and reaction blocks.
      CHARACTER(LEN=2), PARAMETER :: START_STR = '@('  ! start
      CHARACTER(LEN=1), PARAMETER :: END_STR = ')'     ! end

! Strings indicating reaction blocks.
      CHARACTER(LEN=4), PARAMETER :: RXN_BLK     = 'RXNS'      ! start block
      CHARACTER(LEN=8), PARAMETER :: DES_RXN_BLK = 'DES_RXNS'  ! start block
      CHARACTER(LEN=3), PARAMETER :: END_BLK     = 'END'       ! end block

      LOGICAL READING_RXN
      LOGICAL READING_RATE

      LOGICAL DES_RXN
      LOGICAL TFM_RXN

! Logical indicating that the start of a reaction construct has
! been identified.
      LOGICAL IN_CONSTRUCT

! Logical indicating that the chemical equation spans additional lines.
      LOGICAL MORE_ChemEq

! Reaction names
      CHARACTER(len=32),  DIMENSION(:),   ALLOCATABLE :: RXN_NAME
! chemical Equations
      CHARACTER(len=512), DIMENSION(:),   ALLOCATABLE :: RXN_CHEM_EQ
! User defined heat of reaction
      DOUBLE PRECISION,   DIMENSION(:),   ALLOCATABLE :: usrDH
! User defined heat of reaction partitions.
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: usrfDH


! Logical indicating that the start of a reaction construct has
! been identified.
      LOGICAL IN_DES_CONSTRUCT

! Reaction names
      CHARACTER(len=32),  DIMENSION(:),   ALLOCATABLE :: DES_RXN_NAME
! chemical Equations
      CHARACTER(len=512), DIMENSION(:),   ALLOCATABLE :: DES_RXN_CHEM_EQ
! User defined heat of reaction
      DOUBLE PRECISION,   DIMENSION(:),   ALLOCATABLE :: DES_usrDH
! User defined heat of reaction partitions.
      DOUBLE PRECISION,   DIMENSION(:,:), ALLOCATABLE :: DES_usrfDH

      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: setReaction                                          !
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
      SUBROUTINE setReaction(RxN, lNg, lSAg, lM, lNs, lSAs, lDH, lfDH)

      use rxn_com
      use toleranc

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(INOUT) :: RxN
! Number of gas speices
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lM
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs
! User defined heat of reaction.
      DOUBLE PRECISION, INTENT(IN) :: lDH
! User defined heat of reaction partition.
      DOUBLE PRECISION, DIMENSION(0:DIM_M), INTENT(IN) :: lfDH

! Local Variables:
!---------------------------------------------------------------------//
! Alias, phase, species, stoich coeff :: Reactants    Products
      CHARACTER(LEN=32), DIMENSION(50)     :: rAlias    ,  pAlias
      DOUBLE PRECISION, DIMENSION(50) :: rCoeff    ,  pCoeff

! Number of products and reactants
      INTEGER rNo, pNo
! Positions in ChemEq delineating the reactants and products.
      INTEGER rEnd, pStart
! Loop counters
      INTEGER L, LL, M, lN

! Sum of user specified heat of reaction partitions. If fracDH is set
! by the user, they must sum to one over all phases.
      DOUBLE PRECISION sumFDH
! Local storage for chemical equations, left-adjusted and trimmed
      CHARACTER(LEN=512) lChemEq
! Local storage for reaction name, left-adjusted and trimmed
      CHARACTER(LEN=32)  lName
! Logical indicating the reaction is skipped.
      LOGICAL Skip

      LOGICAL pMap(0:lM)

      INTEGER nSpecies, nPhases

      LOGICAL blankAlias(0:(DIM_N_g + lM*DIM_N_s))

! Initialize local reaction name and chemical equation variables.
      lName = trim(adjustl(RxN%Name))
      lChemEq = trim(adjustl(RxN%ChemEq))

      RxN%Classification = "Undefined"
      RxN%Calc_DH = .TRUE.
      RxN%nSpecies = 0
      RxN%nPhases = 0

! Verify that the reactants are separated by --> or = signs. If the
! chemical equation is NONE, the reaction is skipped.
      CALL checkSplit(lName, lChemEq, rEnd, pStart, Skip)
      IF(Skip) THEN
         RxN%nSpecies = 0
         RETURN
      ENDIF
! Set the flag to calculate heat of reaction.
      RxN%Calc_DH = .TRUE.
      IF(lDH /= UNDEFINED) RxN%Calc_DH = .FALSE.

! Pull off the reactants from the chemical equations.
      CALL splitEntries(lName, lChemEq, 1, rEnd, rNo, rAlias, rCoeff)
! Pull off the products from the chemical equations.
      CALL splitEntries(lName, lChemEq, pStart, len_trim(lChemEq),     &
         pNo, pAlias, pCoeff)

      nSpecies = rNo + pNo
      RxN%nSpecies = nSpecies
      Allocate( RxN%Species( nSpecies ))

      CALL checkBlankAliases(lNg, lSAg, lM, lNs, lSAs, blankAlias)

! Check that species in the chemical equation match a species alias
! in one of the phases.
      CALL mapAliases(lName, lChemEq, lNg, lSAg, lM, lNs, lSAs, rNo,   &
         rAlias, rCoeff, -ONE, 0, blankAlias, RxN)

! Check that species in the chemical equation match a species alias
! in one of the phases.
      CALL mapAliases(lName, lChemEq, lNg, lSAg, lM, lNs, lSAs, pNo,   &
         pAlias, pCoeff, ONE, rNo, blankAlias, RxN)

! All the relevant data has been collected at this point. Build the
! reaction block data structure.
      L = max(1,lM)
      LL = (L * (L-1)/2)
      Allocate( RxN%rPhase( LL+L ))


! Initialize local map and global values
      pMap(:) = .FALSE.
      nPhases = 0
      DO lN = 1, nSpecies
         M = RxN%Species(lN)%pMap

         RxN%Species(lN)%mXfr = M
         RxN%Species(lN)%xXfr = ZERO

         IF(.NOT.pMap(M)) THEN
            pMap(M) = .TRUE.
            nPhases = nPhases + 1
         ENDIF
      ENDDO
      RxN%nPhases = nPhases

! Initialize sum of heat of reaction partitions.
      sumFDH = ZERO
! The user specified the heat of reaction.
      IF(.NOT.RxN%Calc_DH) THEN
! Allocate and initialize the heat of reaction storage array.
         Allocate( RxN%HoR( 0:lM ))
         RxN%HoR(:) = ZERO
         DO M=0,lM
! The phase is referenced by the reaction and heat of reaction is
! allocated (in part or fully) this this phase.
            IF(pMap(M) .AND. lFDH(M) .NE. UNDEFINED) THEN
! Store the heat of reaction.
               RxN%HoR(M) = lFDH(M) * lDH
               sumFDH = sumFDH + lFDH(M)
! The phase is not referenced by the reaction, but the heat of reaction
! is allocated (in part or fully) to this phase. Flag error and exit.
            ELSEIF(.NOT.pMap(M) .AND. lFDH(M) .NE. UNDEFINED) THEN
               IF(myPE == PE_IO) THEN
                  write(*,1000) trim(lName)
                  write(*,1001)
                  write(*,1010)
! Write message to log file.
                  write(UNIT_LOG,1000) trim(lName)
                  write(UNIT_LOG,1001)
                  write(UNIT_LOG,1010)
               ENDIF
! Exit MFIX
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
! Logical check: No partition was assigned to an undefined phase.
         DO M=lM+1,DIM_M
            IF(.NOT.RxN%Calc_DH .AND. lFDH(M) .NE. UNDEFINED) THEN
               IF(myPE == PE_IO) THEN
                  write(*,1000) trim(lName)
                  write(*,1001)
                  write(*,1010)
! Write message to log file.
                  write(UNIT_LOG,1000) trim(lName)
                  write(UNIT_LOG,1001)
                  write(UNIT_LOG,1010)
! Exit MFIX
               ENDIF
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDDO
      ENDIF

! Verify that the heat of reaction partitions sum to one.
      IF(.NOT.RxN%Calc_DH .AND. .NOT. COMPARE(sumFDH, ONE)) THEN
         IF(myPE == PE_IO) THEN
            write(*,1002) trim(lName)
            write(*,1010)
! Write message to log file.
            write(UNIT_LOG,1002) trim(lName)
            write(UNIT_LOG,1010)
! Exit MFIX
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      RETURN

 1000 FORMAT(/1X,70('*')/' From: From: setReaction:',/                 &
         ' Message: Heat of reaction is proportioned to a phase not',  &
         ' referenced',/' by the chemical equation for reaction ',A,'.')

 1001 FORMAT(/' If this is a catalytic reaction, reference one of the',&
         ' species of the',/' catalyst phase within the chemical',     &
         ' equation with a stoichiometric',/' coefficient of zero.'/)

 1002 FORMAT(/1X,70('*')/' From: From: setReaction:',/                 &
         ' Message: The heat of reaction partitions (fracDH) to all',  &
         ' phases do',/' not sum to one for reaction ',A,'.')

 1010 FORMAT(' Please refer to the Readme file for chemical equation', &
         ' input formats',/' and correct the data file.',/1X,70('*')/)


      END SUBROUTINE setReaction


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: checkSplit                                           !
!                                                                      !
!  Purpose: Determine the location of reactatns and products within    !
!  the chemical equation. If the entry is NONE, flag that the reaction !
!  is to be skipped for further processing.                            !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE checkSplit(lName, lChemEq, lrEnd, lpStart, lSkip)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Chemical reaction name.
      CHARACTER(len=*), INTENT(IN) :: lName
! Chemical Equation from deck file.
      CHARACTER(len=*), INTENT(IN) :: lChemEq
! Position specifying the end of the reactants in lChemEq
      INTEGER, INTENT(OUT) :: lrEnd
! Position specifying the start of the products in lChemEq
      INTEGER, INTENT(OUT) :: lpStart
! If the chemical equation is NONE, the the skip flag is set to avoid
! further processing.
      LOGICAL, INTENT(OUT) :: lSkip

! Local Variables:
!---------------------------------------------------------------------//
! Position of the head (>) and tail (-) of an arror (-->)
      INTEGER hArr, tArr
! Postion of the head and tail of equal signs (=, ==, ===, ...)
      INTEGER hEqs, tEqs
! Position of the head/tail of a reverse arrow (<--)
      INTEGER hRArr, tRArr
! A flag generated to point out the location of the entry error.
      CHARACTER(LEN=512) FLAG
      FLAG = ''

! If the chemical equation is set to 'none', then the reaction is
! skipped for the simulation.
      lSkip = .FALSE.
      IF(INDEX(lChemEq,'NONE') > 0) THEN
         lSkip = .TRUE.
         lrEnd = UNDEFINED_I
         lpStart = UNDEFINED_I
         RETURN
      ENDIF

! Search for > (part of -->) and search for < (part of <--)
      tArr = INDEX(lChemEq,'-', BACK=.FALSE.)
      hArr = INDEX(lChemEq,">", BACK=.TRUE.)
! Search for the first and last instances of equal signs.
      tEqs = INDEX(lChemEq,"=", BACK=.FALSE.)
      hEqs = INDEX(lChemEq,"=", BACK=.TRUE.)
! Search for < (as part of <-- or <-->). Illegal chemical equation.
      hRArr = INDEX(lChemEq,"<", BACK=.FALSE.)
      tRArr = INDEX(lChemEq,"-", BACK=.TRUE.)

! An illegal arrow was found! Flag error and exit.
      IF(hRArr > 0) THEN
         IF(myPE == PE_IO) THEN
! Construct the error flag.
            IF(hArr > 0) THEN
               FLAG = setFlag(20, hRArr, hArr)
            ELSEIF(tRArr > 0) THEN
               FLAG = setFlag(20, hRArr, tRArr)
            ELSE
               FLAG = setFlag(20, hRArr)
            ENDIF
! Write the message to the screen.
            write(*,1000) trim(lName)
            write(*,1002)'Illegal'
            write(*,1010) trim(lChemEq), trim(Flag)
            write(*,1001)
! Write message to log file.
            write(UNIT_LOG,1000) trim(lName)
            write(UNIT_LOG,1002)'Illegal'
            write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
            write(UNIT_LOG,1001)
! Exit MFIX
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF
! If there are more than one operator, flag error and exit.
      IF(hArr /= 0 .AND. hEqs /= 0) THEN
         IF(myPE == PE_IO) THEN
! Construct the error flag.
            FLAG = setFlag(20, hArr, hEqs)
! Write the message to the screen.
            write(*,1000) trim(lName)
            write(*,1002)'Too many'
            write(*,1010) trim(lChemEq), trim(Flag)
            write(*,1001)
! Write message to log file.
            write(UNIT_LOG,1000) trim(lName)
            write(UNIT_LOG,1002)'Too many'
            write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
            write(UNIT_LOG,1001)
! Exit MFIX
            CALL MFIX_EXIT(myPE)
         ENDIF
! If there is no operator (--> or =), flag error and exit.
      ELSEIF(hArr == 0 .AND. hEqs == 0) THEN
         IF(myPE == PE_IO) THEN
            write(*,1000) trim(lName)
            write(*,1002) 'No'
            write(*,1011) trim(lChemEq)
            write(*,1001)
! Write message to log file.
            write(UNIT_LOG,1000) trim(lName)
            write(UNIT_LOG,1002) 'No'
            write(UNIT_LOG,1011) trim(lChemEq)
            write(UNIT_LOG,1001)
! Exit MFIX
            CALL MFIX_EXIT(myPE)
         ENDIF
! The head of an arrow was found.
      ELSEIF(hArr /= 0) THEN
! Verify that a tail was found.
         IF(tArr == 0) THEN
! Construct the error flag.
            FLAG = setFlag(20, hArr)
            IF(myPE == PE_IO) THEN
               write(*,1000) trim(lName)
               write(*,1003) 'Missing the tail; -->'
               write(*,1010) trim(lChemEq), trim(Flag)
               write(*,1001)
! Write message to log file.
               write(UNIT_LOG,1000) trim(lName)
               write(UNIT_LOG,1003) 'Missing the tail; -->'
               write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
               write(UNIT_LOG,1001)
! Exit MFIX
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSEIF(tArr > hArr) THEN
            IF(myPE == PE_IO) THEN
               FLAG = setFlag(20, hArr, INDEX(lChemEq,'-',BACK=.TRUE.))
               write(*,1000) trim(lName)
               write(*,1003) 'Arror head preceeds the tail; -->'
               write(*,1010) trim(lChemEq), trim(Flag)
               write(*,1001)
! Write message to log file.
               write(UNIT_LOG,1000) trim(lName)
               write(UNIT_LOG,1003) 'Arror head preceeds the tail; -->'
               write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
               write(UNIT_LOG,1001)
! Exit MFIX
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSE
! An arror was used to seperate reactants and products. Send back the
! ending index of reactants and the starting index for products.
            lrEnd = tArr - 1
            lpStart = hArr + 1
         ENDIF
! Equals sign(s) were used to specify the reaction. Send back the ending
! index of reactants and the starting index for products.
      ELSEIF(hEqs /= 0) THEN
         lrEnd = tEqs - 1
         lpStart = hEqs + 1
! Fatal Error. One of the above checks should have caught any problems
! and sent out an error message.
      ELSE
         IF(myPE == PE_IO) THEN
            write(*,1000) trim(lName)
            write(*,1004)
            write(*,1011) trim(lChemEq)
            write(*,1001)
! Write message to log file.
            write(UNIT_LOG,1000) trim(lName)
            write(UNIT_LOG,1004)
            write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
            write(UNIT_LOG,1001)
! Exit MFIX
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

      RETURN

 1000 FORMAT(/1X,70('*')/' From: From: setReaction --> checkSplit',/   &
         ' Message: Error in determining the reactants and products',  &
         ' in the',/' chemical equation for reaction ',A,'.')

 1001 FORMAT(' Please refer to the Readme file for chemical equation', &
         ' input formats',/' and correct the data file.',/1X,70('*')/)

 1002 FORMAT(/1X,A,' operators were found!')

 1003 FORMAT(' Incorrect operator format! ',A)

 1004 FORMAT(' FATAL ERROR: All logical checks failed.')

 1010 FORMAT(/' Chemical Equation: ',A,/1X, A/)

 1011 FORMAT(/' Chemical Equation: ',A,/)


      END SUBROUTINE checkSplit


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: splitEntries()                                       !
!                                                                      !
!  Purpose: Takes a string of either reactants or products and splits  !
!  the string into individual species/stoichiometric entries.          !
!                                                                      !
!  A call to splitAliasAndCoeff is made to further split the entries   !
!  into species aliases and matching stoichiometric coefficients.      !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE splitEntries(lName, lChemEq, lStart, lEnd, lNo,       &
         lAlias, lCoeff)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Chemical reaction name.
      CHARACTER(len=*), INTENT(IN) :: lName
! Chemical Equation.
      CHARACTER(len=*), INTENT(IN) :: lChemEq
! Starting position for substring analysis.
      INTEGER, INTENT(IN) :: lStart
! Ending position for substring analysis.
      INTEGER, INTENT(IN) :: lEnd
! The number of individual species found in lSpecies.
      INTEGER, INTENT(OUT) :: lNo
! Species Aliases from the chemical equation.
      CHARACTER(LEN=32), DIMENSION(50), INTENT(OUT) :: lAlias
! Stoichiometric coefficient pulled from the chemical equation.
      DOUBLE PRECISION, DIMENSION(50), INTENT(OUT) :: lCoeff

! Local Variables:
!---------------------------------------------------------------------//
! Flag indicating that there are more entries to process.
      LOGICAL MORE
! Starting position for left-to-right search.
      INTEGER lPOS
! Position of plus sign charachter found in search.
      INTEGER rPOS

! Initialize storage variables.
      lNo = 0
      lAlias(:) = ''
      lCoeff(:) = UNDEFINED

! Initialiase local variables.
      lPOS = lStart
      MORE = .TRUE.
! Loop through the string, splitting entries separated by a plus sign.
      DO WHILE(MORE)
! Increment the species counter.
         lNo = lNo + 1
! Locate the plus sign. (Left to right)
         rPOS = (lPOS-1) + INDEX(lChemEq(lPOS:lEnd),"+", BACK=.FALSE.)
! A plus sign was found.
         IF(rPOS .GT. lPOS) THEN
! Extract the entry and split it into the species alias and
! stoichiometric coefficient.
            CALL splitAliasAndCoeff(lName, lChemEq, lPOS, rPOS-1,      &
               lAlias(lNo), lCoeff(lNo))
! Indicate that there are more entries to process.
            MORE = .TRUE.
! No plus sign was found. This is the last entry.
         ELSE
! Extract the entry and split it into the species alias and
! stoichiometric coefficient.
            CALL splitAliasAndCoeff(lName, lChemEq, lPOS, lEnd,        &
               lAlias(lNo), lCoeff(lNo))
! Indicate that there are no more entries to process.
            MORE = .FALSE.
         ENDIF
! Move past the found plus sign for next search.
         lPOS = rPOS + 1
      ENDDO

      RETURN
      END SUBROUTINE splitEntries

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: splitAliasAndCoeff()                                 !
!                                                                      !
!  Purpose: Take a string containing a species alias and stoichio-     !
!  metric coefficient and splits them into their respective parts.     !
!                                                                      !
!  If no numerical coefficient is found, it is set to one.             !
!                                                                      !
!  If present, asterisks (*) are assumed to seperate numerical         !
!  coefficients and the species alias. If more than one asterisk is    !
!  found, and error is reported and MFIX_EXIT is called.               !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE splitAliasAndCoeff(lName, lChemEq, lStart, lEnd,      &
         lAlias, lCoeff)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Chemical reaction name.
      CHARACTER(len=*), INTENT(IN) :: lName
! Chemical Equation.
      CHARACTER(len=*), INTENT(IN) :: lChemEq
! Starting position for substring analysis.
      INTEGER, INTENT(IN) :: lStart
! Ending position for substring analysis.
      INTEGER, INTENT(IN) :: lEnd
! Species Aliases from the chemical equation.
      CHARACTER(LEN=32), INTENT(OUT) :: lAlias
! Stoichiometric coefficient pulled from the chemical equation.
      DOUBLE PRECISION, INTENT(OUT) :: lCoeff

! Local Variables:
!---------------------------------------------------------------------//
! Flag
      LOGICAL MATCH
      INTEGER nPOS

      INTEGER L, N, IOS, aPOS, a2POS

      CHARACTER(LEN=12), PARAMETER :: Numbers = '.0123456789'

! A flag generated to point out the location of the entry error.
      CHARACTER(LEN=512) FLAG
      FLAG = ''

! Locate the first asterisk (if any). Search left-to-right.
      aPOS = INDEX(lChemEq(lStart:lEnd),"*", BACK=.FALSE.)
! An asterisk was found.
      IF(aPOS .GT. ZERO) THEN
! Make sure that there isn't another asterisk further down the string.
         a2POS = INDEX(lChemEq(lStart:lEnd),"*", BACK=.TRUE.)
         IF(aPOS /= a2POS) THEN
            IF(myPE == PE_IO) THEN
! Construct the error flag.
               FLAG = setFlag(20, lStart+aPOS, lStart+a2POS)
! Write the message to the screen.
               write(*,1000) trim(lName)
               write(*,1002)'Too many'
               write(*,1010) trim(lChemEq), trim(Flag)
               write(*,1001)
! Write message to log file.
               write(UNIT_LOG,1000) trim(lName)
               write(UNIT_LOG,1002)'Too many'
               write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
               write(UNIT_LOG,1001)
! Exit MFIX
               CALL MFIX_EXIT(myPE)
            ENDIF
         ELSE
! Store left-of-asterisk as the coefficient. If an error occurs in
! converting the string to double precision, flag the problem and
! call MFIX_EXIT.
            READ(lChemEq(lStart:(lStart+aPOS-2)),*,IOSTAT=IOS) lCoeff
            IF(IOS .NE. 0 .AND. myPE == PE_IO) THEN
! Construct the error flag.
               FLAG = setFlag(20, lStart + int(aPOS/2))
! Write the message to the screen.
               write(*,1000) trim(lName)
               write(*,1010) trim(lChemEq), trim(Flag)
               write(*,1001)
! Write message to log file.
               write(UNIT_LOG,1000) trim(lName)
               write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
               write(UNIT_LOG,1001)
! Exit MFIX
               CALL MFIX_EXIT(myPE)
            ENDIF
! Store right-of-asterisk as the species alias.
            WRITE(lAlias,"(A)") &
               trim(adjustl(lChemEq((lStart+aPOS):lEnd)))
         ENDIF
! If no asterisk was found, search for numbers and spaces.
      ELSE
! Initialize the position of last consecutive number.
         nPOS = 0
! In a left-to-right search, check if the characters in the entry are
! numbers or punctuation.
         DO L=lStart,lEnd
            MATCH = .FALSE.
            DO N=1,12
               IF(lChemEq(L:L) /= Numbers(N:N)) CYCLE
! Note the position of the number.
               nPOS = L
! Flag that a match was made.
               MATCH = .TRUE.
            ENDDO
! If no match, assume the end of the coefficient was found.
            IF(.NOT.MATCH) EXIT
         ENDDO
! If no numbers or punctuation was found, assumed the stoichiometric
! coefficient is one.
         IF(trim(lChemEq(lStart:nPOS)) =='') THEN
            lCoeff = 1.0d0
         ELSE
! If leading numbers were found, store as the stoich-coeff.
            READ(lChemEq(lStart:nPOS),*,IOSTAT=IOS) lCoeff
! Report any problems in converting the string to double precision.
            IF(IOS .NE. 0 .AND. myPE == PE_IO) THEN
! Construct the error flag.
               FLAG = setFlag(20, &
                  lStart+int(len_trim(lChemEq(lStart:nPOS))/2))
! Write the message to the screen.
               write(*,1000) trim(lName)
               write(*,1010) trim(lChemEq), trim(Flag)
               write(*,1001)
! Write message to log file.
               write(UNIT_LOG,1000) trim(lName)
               write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
               write(UNIT_LOG,1001)
! Exit MFIX
               CALL MFIX_EXIT(myPE)
            ENDIF
         ENDIF
! Store right-of-coefficient as the species alias.
         READ(lChemEq(nPOS+1:lEnd),*,IOSTAT=IOS) lAlias
      ENDIF
! Quick check to make sure that the species alias is not empty.
      IF(LEN_TRIM(lAlias) == 0 .AND. myPE == PE_IO) THEN
! Construct the error flag.
         FLAG = setFlag(20, lStart + int(lEnd/2))
! Write the message to the screen.
         write(*,1003) trim(lName)
         write(*,1010) trim(lChemEq), trim(Flag)
         write(*,1001)
! Write message to log file.
         write(UNIT_LOG,1003) trim(lName)
         write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
         write(UNIT_LOG,1001)
! Exit MFIX
         CALL MFIX_EXIT(myPE)
      ENDIF

      RETURN

 1000 FORMAT(/1X,70('*')/' From: From: setReaction -->',               &
         ' splitAliasAndCoeff',/' Message: Error determining the',     &
         ' stoichiometric coefficient in the',/' chemical equation',   &
         ' for reaction ',A,'.')


 1001 FORMAT(' Please refer to the Readme file for chemical equation', &
         ' input formats',/' and correct the data file.',/1X,70('*')/)

 1002 FORMAT(/1X,A,' operators were found!')

 1003 FORMAT(/1X,70('*')/' From: From: setReaction -->',               &
         ' splitAliasAndCoeff',/' Message: Error determining the',     &
         ' speices in the chemical equation for',/' reaction ',A,'.'/)

 1010 FORMAT(/' Chemical Equation: ',A,/1X, A/)

 1011 FORMAT(/' Chemical Equation: ',A,/)

      END SUBROUTINE splitAliasAndCoeff

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: checkBlankAliases()                                  !
!                                                                      !
!  Purpose: Take a string containing a species alias and stoichio-     !
!  metric coefficient and splits them into their respective parts.     !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE checkBlankAliases(lNg, lSAg, lM, lNs, lSAs, lBA)

      IMPLICIT NONE

! Number of gas speices
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lM
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs

      LOGICAL, INTENT(OUT) :: lBA(0:(DIM_N_g + lM*DIM_N_s))

      INTEGER M, N

! Loop counter for continuum and discrete species
      INTEGER Nsp

! Initialize counters
      Nsp = 0

      lBA(0) = .FALSE.
      DO N = 1, lNg
         Nsp = Nsp + 1
         lBA(Nsp) = .FALSE.
         IF(len_trim(lSAg(N)) == 0) THEN
            lBA(Nsp) = .TRUE.
            lBA(0) = .TRUE.
         ENDIF
      ENDDO

! Compare aliaes between solids phases
      DO M = 1, lM
         DO N = 1, lNs(M)
            Nsp = Nsp + 1
            lBA(Nsp) = .FALSE.
            IF(len_trim(lSAs(M,N)) == 0) THEN
               lBA(Nsp) = .TRUE.
               lBA(0) = .TRUE.
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE checkBlankAliases

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: mapAliases()                                         !
!                                                                      !
!  Purpose: Take a string containing a species alias and stoichio-     !
!  metric coefficient and splits them into their respective parts.     !
!                                                                      !
!  If no numerical coefficient is found, it is set to one.             !
!                                                                      !
!  If present, asterisks (*) are assumed to seperate numerical         !
!  coefficients and the species alias. If more than one asterisk is    !
!  found, and error is reported and MFIX_EXIT is called.               !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE mapAliases(lName, lChemEq, lNg, lSAg, lM, lNs, lSAs,  &
         lNo, lAlias, lCoeff, lSgn, lStart, lBA, lRxN)

      use rxn_com

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Chemical reaction name.
      CHARACTER(len=*), INTENT(IN) :: lName
! Chemical Equation.
      CHARACTER(len=*), INTENT(IN) :: lChemEq
! Number of gas speices
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lM
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs
! Number of products (or reactants)
      INTEGER, INTENT(IN) :: lNo
! Species Alaises pulled from the chemical equation.
      CHARACTER(LEN=32), DIMENSION(50), INTENT(IN) :: lAlias

      DOUBLE PRECISION, DIMENSION(50), INTENT(IN) :: lCoeff

      DOUBLE PRECISION, INTENT(IN) :: lSgn

      INTEGER, INTENT(IN) :: lStart

      LOGICAL, INTENT(IN) :: lBA(0:(DIM_N_g + lM*DIM_N_s))

! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(INOUT) :: lRxN


! Local Variables:
!---------------------------------------------------------------------//
! Loop counters.
      INTEGER L, M, NN
! A flag generated to point out the location of the entry error.
      CHARACTER(LEN=512) FLAG
! Location in string to locate error.
      INTEGER lPOS, rPOS

! Initialize the error flag.
      FLAG = ''

! Loop over the number of (reactants/products)
      ALOOP : DO L=1, lNo
! Compare entry with gas phase species.
         DO NN = 1, lNg
            IF( checkMatch(lSAg(NN), lAlias(L))) THEN
               lRxN%Species(lStart + L)%pMap = 0
               lRxN%Species(lStart + L)%sMap = NN
               lRxN%Species(lStart + L)%Coeff = lSgn * lCoeff(L)
               CYCLE ALOOP
            ENDIF
         ENDDO
! Compare entry with solids phase species.
         DO M = 1, lM
            DO NN = 1, lNs(M)
               IF(checkMatch(lSAs(M,NN),lAlias(L))) THEN
                  lRxN%Species(lStart + L)%pMap = M
                  lRxN%Species(lStart + L)%sMap = NN
                  lRxN%Species(lStart + L)%Coeff = lSgn * lCoeff(L)
                  CYCLE ALOOP
               ENDIF
            ENDDO
         ENDDO
! No matching species was located. Flag an error and exit.
         IF(myPE == PE_IO) THEN
! Construct the error flag.
            lPOS = INDEX(lChemEq,trim(lAlias(L)), BACK=.FALSE.)
            rPOS = INDEX(lChemEq,trim(lAlias(L)), BACK=.TRUE.)
            FLAG = setFlag(20, 1 + int((lPOS + rPOS)/2))
! Write the message to the screen.
            write(*,1000) trim(lAlias(L)), trim(lName)
            write(*,1010) trim(lChemEq), trim(Flag)
            IF(lBA(0)) CALL writeBA()
            write(*,1001)
! Write message to log file.
            write(UNIT_LOG,1000) trim(lAlias(L)), trim(lName)
            write(UNIT_LOG,1010) trim(lChemEq), trim(Flag)
            IF(lBA(0)) CALL writeBA(UNIT_LOG)
            write(UNIT_LOG,1001)
         ENDIF
! Exit MFIX
         CALL MFIX_EXIT(myPE)

      ENDDO ALOOP

      RETURN

 1000 FORMAT(/1X,70('*')/' From: From: setReaction --> mapAliases',/   &
         ' Message: Unable to match species ',A,' in the chemical',    &
         ' equation for ',/' reaction ',A,'.')

 1001 FORMAT(/' Please refer to the Readme file for chemical equation',&
         ' input formats',/' and correct the data file.',/1X,70('*')/)

 1010 FORMAT(/' Chemical Equation: ',A,/1X, A/)

 1011 FORMAT(/' Chemical Equation: ',A,/)


      contains

!......................................................................!
!  Function name: checkMatch                                           !
!                                                                      !
!  Purpose: Takes two species aliases as arguments, converts them to   !
!  uppercase and checks if they match.                                 !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!......................................................................!
      LOGICAL FUNCTION checkMatch(lSA, ceSA)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
      CHARACTER(LEN=32), INTENT(IN) :: lSA, ceSA

! Local Variables:
!---------------------------------------------------------------------//
      CHARACTER(LEN=32) tlSA

! Copy species alias.
      tlSA = lSA
! Remove case sensitivity.
      CALL MAKE_UPPER_CASE (tlSA,32)
! Compare the two strings.
      checkMatch = .FALSE.
      IF(trim(tlSA) == trim(ceSA)) checkMatch = .TRUE.
      RETURN
      END FUNCTION checkMatch

!......................................................................!
!  Function name: updateMap                                            !
!                                                                      !
!  Purpose: Flags that the passed phase is part of the chemical        !
!  reaction. If the phase was not already noted, the number of phases  !
!  in the reaction is increased and the flag set true. Additionally,   !
!  The number of species (either product or reactant) for the phase    !
!  is incremented.                                                     !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!......................................................................!
      SUBROUTINE updateMap(lnP, lpMap, llNoP)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
      INTEGER, INTENT(INOUT) :: lnP
! Map of phases for the reaction.
      LOGICAL, INTENT(INOUT) :: lpMap
      INTEGER, INTENT(INOUT) :: llNoP

! Local Variables:
!---------------------------------------------------------------------//
! None

! Increment the number of reactants/products this phase has involved in
! the current reaction.
      llNoP = llNoP + 1
! If the phase was already identified, return.
      IF(lpMap) RETURN
! If this is the first time the phase is identififed, set the flag to
! true and increment the total number of phases in the reaction.
      lnP = lnP + 1
      lpMap = .TRUE.

      RETURN
      END SUBROUTINE updateMap


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: writeBA()                                            !
!                                                                      !
!  Purpose: Print out which species were not given aliases.            !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE writeBA(FUNIT)

      IMPLICIT NONE

! File unit
      INTEGER, OPTIONAL, INTENT(IN) :: FUNIT

      INTEGER M, N

! Loop counter for continuum and discrete species
      INTEGER Nsp

      IF(.NOT.PRESENT(FUNIT)) THEN
         write(*,1000)
      ELSE
         write(FUNIT,1000)
      ENDIF

! Initialize counters
      Nsp = 0

      DO N = 1, lNg
         Nsp = Nsp + 1
         IF(lBA(Nsp)) THEN
            IF(.NOT.PRESENT(FUNIT)) THEN
               write(*,1001)N
            ELSE
               write(FUNIT,1001) N
            ENDIF
         ENDIF
      ENDDO

! Compare aliaes between solids phases
      DO M = 1, lM
         DO N = 1, lNs(M)
            Nsp = Nsp + 1

            IF(lBA(Nsp)) THEN
               IF(.NOT.PRESENT(FUNIT)) THEN
                  write(*,1002)M, N
               ELSE
                  write(FUNIT,1002)M, N
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      RETURN

 1000 FORMAT(' Species aliases were not provided for the following:')
 1001 FORMAT(3X, ' Gas phase species ',I2)
 1002 FORMAT(3X, ' Solid phase ',I2,' specie ',I2)

      END SUBROUTINE writeBA

      END SUBROUTINE mapAliases


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: setFlag                                              !
!                                                                      !
!  Purpose: Creates a flag pointing to a particular string location.   !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      CHARACTER(len=512) FUNCTION setFlag(fill, flg1, flg2) RESULT(OUT)

      IMPLICIT NONE

! Pass Arguments:
!---------------------------------------------------------------------//
! Number of leading spaces to fill with dashes. This is used to jump
! past any lead-in text.
      INTEGER, INTENT(IN) :: fill
! Location in string to place the pointer.
      INTEGER, INTENT(IN) :: flg1
! Optional - location of second pointer in string
      INTEGER, INTENT(IN), OPTIONAL :: flg2

! Local Variables:
!---------------------------------------------------------------------//
      INTEGER L, FILL1, FILL2

! Create a string with "FILL" dash characters.
      OUT = ''
      DO L = 1, FILL-1
         WRITE(OUT,"(A,A)") trim(OUT), '-'
      ENDDO

! If a second pointer is present, determined the larger of the two.
      IF(PRESENT(flg2)) THEN
         IF(flg1 < flg2) THEN
            FILL1 = flg1 - 1
            FILL2 = (flg2-flg1) - 1
         ELSE
            FILL1 = flg2 - 1
            FILL2 = (flg1-flg2) - 1
         ENDIF
      ELSE
         FILL1 = flg1 - 1
         FILL2 = 0
      ENDIF

! Fill with dashes up to the the first pointer. ----^
      DO L = 1, FILL1
         WRITE(OUT,"(A,A)") trim(OUT), '-'
      ENDDO
      WRITE(OUT,"(A,A)") trim(OUT), '^'
! Fill with dashes up to the second pointer. ----^---^
      IF(FILL2 > 0) THEN
         DO L = 1, FILL2
            WRITE(OUT,"(A,A)") trim(OUT), '-'
         ENDDO
         WRITE(OUT,"(A,A)") trim(OUT), '^'
      ENDIF

      END FUNCTION setFlag

END MODULE parse
