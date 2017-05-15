MODULE RXN_COM

   USE compar
   USE exit, only: mfix_exit
   USE funits
   USE param
   USE param1

! The following data types are used to group chemical reaction data.
!-----------------------------------------------------------------------

! Species belong to PHASE_ associated with a particular reaction.
      TYPE SPECIES_
! A link between the reacting species' arbitrary index and the
! associated phase index in MFiX.
         INTEGER pMap
! A link between the reacting species' arbitrary index and the
! associated species index in MFiX.
         INTEGER sMap
! Stoichiometric coefficient of the species from chemical equation.
         DOUBLE PRECISION Coeff
! Molecular weight
         DOUBLE PRECISION MW
! Fractional mass transfer
         DOUBLE PRECISION xXfr
! Index indicating enthalpy transfer associated with mass transfer.
         INTEGER mXfr
! Molecular weight of speices multiplying the stoichiometric coefficient
         DOUBLE PRECISION MWxStoich
      END TYPE SPECIES_

! Grouping of reaction information.
      TYPE REACTION_BLOCK
! Name of reaction construct from data file.
         CHARACTER(LEN=32) :: Name
! User defined chemical equation from data file.
         CHARACTER(LEN=512) :: ChemEq
! Reaction classification: Homogeneous, Heterogeneous, Catalytic.
         CHARACTER(LEN=16) :: Classification
! Indicates if the automated heat of reaction is to be calculated (T) or
! if the user has supplied a heat of reaction (F).
         LOGICAL Calc_DH
! Number of phases associated with the reaction.
         INTEGER nPhases
! Number of species associated with the reaction.
         INTEGER nSpecies
! User-specified heat of reaction split among phases by fracDH
         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: HoR
! Interphase mass transfer.
         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rPhase
! Reactant/Product information
         TYPE(SPECIES_), DIMENSION(:), ALLOCATABLE :: Species

      END TYPE REACTION_BLOCK

      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: checkDuplicateAlaises                                   !
!                                                                      !
!  Purpose: Loop through species in all phases and ensure that no two  !
!  entries are the same. ***Warning*** Species aliases that were not   !
!  specified are skipped. Non-specified aliases are treated later in   !
!  parse_mod.f/mapAliases.                                             !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE checkDuplicateAliases(lNg, SA_g, lMMx, lNs, SA_s)

      use error_manager

      IMPLICIT NONE

! Number of gas speices
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: SA_g
! Number of solids phases
      INTEGER, INTENT(IN) :: lMMx
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: SA_s

! Loop indices.
      INTEGER M1, N1  ! Phase, Species
      INTEGER M2, N2  ! Phase, Species

      CHARACTER(len=32) SA1, SA2

      CALL INIT_ERR_MSG("RXN_COM --> checkDuplicateAliases")

! Set variables for error messages.
      M1 = 0
      M2 = 0
! Compare gas phase aliases.
      DO N1 = 1, lNg
         SA1 =SA_g(N1)
         IF(len_trim(SA1) == 0) CYCLE
         DO N2=N1+1,lNg
            SA2 = SA_g(N2)
            IF(len_trim(SA2) == 0) CYCLE
            IF(compareAliases(SA1, SA2)) GoTo 100
         ENDDO
! Compare gas and solids phase aliases.
         DO M2 = 1, lMMx
            DO N2 = 1, lNs(M2)
               SA2 = SA_s(M2,N2)
               IF(len_trim(SA2) == 0) CYCLE
               IF(compareAliases(SA1, SA2)) GoTo 100
            ENDDO
         ENDDO
      ENDDO
! Compare aliases between solids phases
      DO M1 = 1, lMMx
         DO N1 = 1, lNs(M1)
            SA1 = SA_s(M1,N1)
            IF(len_trim(SA1) == 0) CYCLE
! Self phase comparison.
            M2 = M1
            DO N2=N1+1, lNs(M2)
               SA2 = SA_s(M2,N2)
               IF(len_trim(SA2) == 0) CYCLE
               IF(compareAliases(SA1, SA2)) GoTo 100
            ENDDO
! Compare with other phases.
            DO M2 = M1+1, lMMx
               DO N2 = 1, lNs(M2)
                  SA2 = SA_s(M2,N2)
                  IF(len_trim(SA2) == 0) CYCLE
                  IF(compareAliases(SA1, SA2)) GoTo 100
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      CALL FINL_ERR_MSG
      RETURN

  100 WRITE(ERR_MSG, 1100) M1, N1, SA1, M2, N2, SA2
      CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1100: Non-unique species aliases detected.',/      &
         3x,'Phase: ',I2,', Species: ',I3,' - Alias: ',A,/             &
         3x,'Phase: ',I2,', Species: ',I3,' - Alias: ',A,//            &
         'Please correct the mfix.dat file.')

      END SUBROUTINE checkDuplicateAliases

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: checkSpeciesInc()                                    !
!                                                                      !
!  Purpose: Loop through the species.inc file and verify that the      !
!  match those provided in the datafile.                               !
!                                                                      !
!  Variables referenced: None                                          !
!                                                                      !
!  Variables modified: None                                            !
!                                                                      !
!  Local variables: None                                               !
!                                                                      !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      SUBROUTINE checkSpeciesInc(lNg, SA_g, lMMx, lNs, SA_s,           &
         lNRxn,  lRNames, lNRxn_DES, lRNames_DES)

      use run, only: REINITIALIZING
      use error_manager
      use toleranc
      USE utilities, ONLY: blank_line, seek_comment

      IMPLICIT NONE

! Number of gas speices
      INTEGER, INTENT(IN) :: lNg
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: SA_g
! Number of solids phases
      INTEGER, INTENT(IN) :: lMMx
! Number of species in each solids phase.
      INTEGER, DIMENSION(DIM_M), INTENT(IN) :: lNs
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: SA_s
! Number of reactions
      INTEGER, INTENT(IN) :: lNRxn
! Reaction Names (aliases)
      CHARACTER(len=32), INTENT(IN) ::  lRNames(DIMENSION_RXN)
! Number of discrete reactions
      INTEGER, INTENT(IN) :: lNRxn_DES
! Reaction Names for discrete solids (aliases)
      CHARACTER(len=32), INTENT(IN) ::  lRNames_DES(DIMENSION_RXN)

! Input/Output status.
      INTEGER :: IOS
! File unit.
      INTEGER, PARAMETER :: FUNIT = 167
! Full path to Burcat and Ruscic database
      CHARACTER(len=255) :: FILENAME
      CHARACTER(len=128) :: INPUT
! Loop counters
      INTEGER :: SRC, M
! Position of interest in string.
      INTEGER :: POS
! Index from species.inc file.
      INTEGER :: lIndex
      CHARACTER(len=64) :: lName
      CHARACTER(len=32) :: tName
! Length of noncomment string
      INTEGER :: LINE_LEN

      CALL INIT_ERR_MSG("RXN_COM --> checkDuplicateAliases")

      SRC = 0

! Loop over possible locations .
      SRC_LP: DO
         SRC = SRC + 1
         SELECT CASE(SRC)

! Check the local run directory.
         CASE(1); FILENAME = 'species.inc'
            OPEN(CONVERT='BIG_ENDIAN',UNIT=FUNIT,FILE=trim(FILENAME),STATUS='OLD',IOSTAT=IOS)
            IF(IOS /= 0) CYCLE SRC_LP
            IF(.NOT.REINITIALIZING)THEN
               WRITE(ERR_MSG, 1000)'species.inc'
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF

 1000 FORMAT(/2X,'Verifying reaction aliases in ',A)

! No species.inc file was located.
         CASE DEFAULT
            IF(.NOT.REINITIALIZING)THEN
               WRITE(ERR_MSG, 1004)
               CALL FLUSH_ERR_MSG
            ENDIF
            EXIT SRC_LP
         END SELECT

 1004 FORMAT('Warning 1004: Unable to locate the species.inc file. No ',&
         'verification',/'of mfix.dat species aliases or reaction ',    &
         'names can be preformed.')

         REWIND(FUNIT)
         READ_LP: DO
            READ(FUNIT,"(A)",IOSTAT=IOS) INPUT

! This is a sanity check because the species.inc file is generated by
! make_mfix and therefore should be the correct format.
            IF(IOS > 0) THEN
               WRITE(ERR_MSG,1200) trim(adjustl(FILENAME))
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
 1200 FORMAT('Error 1200: There was a problem reading file: ',A)

! All entries have been processed.
            ELSEIF(IOS<0)THEN
               EXIT SRC_LP
            ENDIF

! Clean up the input.
            LINE_LEN = SEEK_COMMENT(INPUT,LEN(INPUT)) - 1
            CALL REMOVE_COMMENT(INPUT, LINE_LEN + 1, LEN(INPUT))
            CALL MAKE_UPPER_CASE(INPUT, LINE_LEN)
            CALL REPLACE_TAB(INPUT, LINE_LEN)

! Skip empty entires.
            IF(LINE_LEN <= 0) CYCLE READ_LP
            IF(BLANK_LINE(INPUT)) CYCLE READ_LP

            POS = INDEX(INPUT,"INTEGER, PARAMETER ::")
            IF(POS /= 0) THEN
               INPUT = INPUT((POS + 21):)
            ELSE
               CYCLE READ_LP
            ENDIF

! We only want to process lines that have = as the other are coments.
            POS = INDEX(INPUT,"=")
            IF(POS == 0) CYCLE READ_LP

! Store the species alias.
            WRITE(lName,"(A)") trim(adjustl(INPUT(:(POS-1))))

! Convert the read index from string to integer. Report any errors.
            READ(INPUT((POS+1):),*,IOSTAT=IOS) lIndex
            IF(IOS /= 0) THEN
               WRITE(ERR_MSG,1205) 'index', trim(INPUT)
               CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
            ENDIF

 1205 FORMAT('Error 1205: Unable to obtain alias index from species.', &
         'inc file.',//' INPUT: ',A)

! Match against what was provided in the datafile:
! Gas phase species aliases.
            IF(lIndex <= lNg) THEN
               tName = SA_g(lIndex)
               IF(compareAliases(tName, lName)) CYCLE READ_LP
            ENDIF

! Solids phase species aliases.
            DO M = 1, lMMx
               IF(lIndex <= lNs(M))THEN
                  tName = SA_s(M, lIndex)
                  IF(compareAliases(tName, lName)) CYCLE READ_LP
               ENDIF
            ENDDO

! Reaction Names
            IF(lIndex <= lNRxn)THEN
               tName =  lRNames(lIndex)
               IF(compareAliases(tName, lName)) CYCLE READ_LP
            ENDIF

! Reaction Names for discrete solids
            IF(lIndex <= lNRxn_DES)THEN
               tName =  lRNames_DES(lIndex)
               IF(compareAliases(tName, lName)) CYCLE READ_LP
            ENDIF

            WRITE(ERR_MSG,1300) trim(lName), lIndex
            CALL FLUSH_ERR_MSG

 1300 FORMAT('Error 1300: An entry in the species.inc file does not ', &
         'match any inputs',/'in the mfix.dat file.'/3x,'Name: ',A,4x, &
         'Index: ',I3,/'If the quantity or order of gas species, ',    &
         'solids species, or chemical',/'reactions has changed, then ',&
         'the executable must be re-build. Please',/'see the document',&
         'ation for specifying chemical reactions.')

         ENDDO READ_LP
      ENDDO SRC_LP

      CLOSE(FUNIT)
      CALL FINL_ERR_MSG
      RETURN

      END SUBROUTINE checkSpeciesInc


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Function name: compareAlaises()                                     !
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
      LOGICAL FUNCTION compareAliases(lS1, lS2, N1, N2)

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: lS1, lS2

      INTEGER, OPTIONAL, INTENT(IN) :: N1, N2

      CALL MAKE_UPPER_CASE (lS1, len(lS1))
      CALL MAKE_UPPER_CASE (lS2, len(lS2))

      compareAliases = .FALSE.
      IF(trim(lS1) == trim(lS2)) compareAliases = .TRUE.

      IF(.NOT.compareAliases) RETURN

      IF(PRESENT(N1) .AND. PRESENT(N2)) THEN
         IF(N1 == N2) THEN
            compareAliases = .TRUE.
         ELSE
            compareAliases = .FALSE.
         ENDIF
      ENDIF

      RETURN
      END FUNCTION compareAliases



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: WRITE_RXN_SUMMARY                                       !
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
      SUBROUTINE WRITE_RXN_SUMMARY(RxN, lSAg, lSAs, ABORT, fUNIT)

      USE toleranc

      IMPLICIT NONE

! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(IN) :: RxN

! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs
! Flag to abort the run.
      LOGICAL, INTENT(IN) :: ABORT

! Optional file unit.
      INTEGER, OPTIONAL :: fUNIT

      CHARACTER(LEN=72) :: OUTPUT, full, divided, empty

      CHARACTER(LEN=32) :: lSP

      INTEGER lN, M, N
      INTEGER lS, lE

      INTEGER UNIT_FLAG

      IF(present(fUnit)) THEN
         UNIT_FLAG = fUNIT
      ELSE
         UNIT_FLAG = -1
      ENDIF

      empty = '  '
      CALL WRITE_RS0(empty, UNIT_FLAG)

      full = ''
      WRITE(full,2000)

      divided = ''
      WRITE(divided,2005)

! Lead bar
      CALL WRITE_RS0(full, UNIT_FLAG)
! Reaction Nmae
      OUTPUT = ''
      WRITE(OUTPUT, 2001)trim(RxN%Name)
      OUTPUT(72:72) = '|'
      CALL WRITE_RS0(OUTPUT, UNIT_FLAG)

! Row Divider
      CALL WRITE_RS0(full, UNIT_FLAG)

      OUTPUT = ''
      WRITE(OUTPUT, 2002)trim(RxN%ChemEq(1:54))
      OUTPUT(72:72) = '|'
      CALL WRITE_RS0(OUTPUT, UNIT_FLAG)

      CALL WRITE_RS0(full, UNIT_FLAG)

      IF(RxN%nSpecies > 0) THEN

         OUTPUT = ''
         WRITE(OUTPUT, 2007)trim(RxN%Classification)
         OUTPUT(72:72) = '|'
         CALL WRITE_RS0(OUTPUT, UNIT_FLAG)
! Row Divider
         CALL WRITE_RS0(full, UNIT_FLAG)

         WRITE(OUTPUT,2003); CALL WRITE_RS0(OUTPUT, UNIT_FLAG)
         WRITE(OUTPUT,2004); CALL WRITE_RS0(OUTPUT, UNIT_FLAG)
         CALL WRITE_RS0(divided, UNIT_FLAG)
      ENDIF


      DO lN = 1, RxN%nSpecies

         M = RxN%Species(lN)%pMap
         N = RxN%Species(lN)%sMap

         WRITE(OUTPUT,2006)

         IF(M == 0) THEN
            IF(len_trim(lSAg(N)) > 8) THEN
               lSP = lSAg(N)
               OUTPUT(5:13) = lSP(1:8)
            ELSE
              lS = (9-int(len_trim(lSAg(N))/2))
              lE = lS + len_trim(lSAg(N))
               OUTPUT(lS:lE) = trim(lSAg(N))
            ENDIF
            WRITE(OUTPUT(32:35),"(A)") 'Gas'
         ELSE
            IF(len_trim(lSAs(M,N)) > 8) THEN
               lSP = lSAs(M,N)
               OUTPUT(5:13) = lSP(1:8)
            ELSE
               lS = (9-int(len_trim(lSAs(M,N))/2))
               lE = lS + len_trim(lSAs(M,N))
               OUTPUT(lS:lE) = trim(lSAs(M,N))
            ENDIF
            WRITE(OUTPUT(30:36),"(A,I2)") 'Solid',M
         ENDIF
         WRITE(OUTPUT(43:44),"(I2)") N
         WRITE(OUTPUT(51:60),"(F9.4)") RxN%Species(lN)%MW

         IF(COMPARE(RxN%Species(lN)%Coeff, ZERO)) THEN
            WRITE(OUTPUT(17:26),"(F9.4)") ZERO
            WRITE(OUTPUT(63:71),"(A)") 'Catalyst'
         ELSEIF(RxN%Species(lN)%Coeff < ZERO) THEN
            WRITE(OUTPUT(17:26),"(F9.4)") -RxN%Species(lN)%Coeff
            WRITE(OUTPUT(63:71),"(A)") 'Reactant'
         ELSE
            WRITE(OUTPUT(17:26),"(F9.4)")  RxN%Species(lN)%Coeff
            WRITE(OUTPUT(63:70),"(A)") 'Product'
         ENDIF
         CALL WRITE_RS0(OUTPUT, UNIT_FLAG)
         CALL WRITE_RS0(divided, UNIT_FLAG)

      ENDDO

      CALL WRITE_RS0(empty, UNIT_FLAG)

      IF(ABORT) CALL MFIX_EXIT(myPE)
      RETURN


 2000 FORMAT(2X,'|',68('-'),'|')

 2001 FORMAT(2X,'| Name: ',A)
 2002 FORMAT(2x,'| Chemical Eq: ',A)

 2003 FORMAT('  | Species  |   Stoich    |         | Species |',       &
              ' Molecular  |          |')

 2004 FORMAT('  |  Alias   |   Coeff.    |  Phase  |  Index  |',       &
              '   Weight   |   Type   |')


 2005 FORMAT(2X,'|',10('-'),'|',13('-'),'|',9('-'),'|',9('-'),'|',     &
             12('-'),'|',10('-'),'|')

 2006 FORMAT(2X,'|',10(' '),'|',13(' '),'|',9(' '),'|',9(' '),'|',     &
             12(' '),'|',10(' '),'|')


 2007 FORMAT(2X,'| Classification: ',A)

      contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: WRITE_RS0                                               !
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
      SUBROUTINE WRITE_RS0(LINE, UFLAG)

      use run, only: REINITIALIZING
      use error_manager

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: LINE
      INTEGER, INTENT(IN) :: UFLAG

      CALL INIT_ERR_MSG("WRITE_RXN_SUMMARY --> WRITE_RS0")

      IF(UFLAG == -1)THEN
         IF(.NOT.REINITIALIZING) THEN
            WRITE(ERR_MSG,*) LINE
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
         ENDIF
      ELSE
         WRITE(UFLAG,*) LINE
      ENDIF
      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE WRITE_RS0
      END SUBROUTINE WRITE_RXN_SUMMARY


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: checkThermoReqs                                         !
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
      SUBROUTINE checkThermoReqs(RxN, S_g, S_s, rDB, MWg, MWs, Cpg0, Cps0)

      use error_manager
      use toleranc

      IMPLICIT NONE

! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(INOUT) :: RxN

      CHARACTER(len=18), INTENT(IN) :: S_g(DIM_N_g)
      CHARACTER(len=18), INTENT(in) :: S_s(DIM_M, DIM_N_s)
      LOGICAL, INTENT(inout) :: rDB(0:DIM_M, DIM_N_g)
      DOUBLE PRECISION, INTENT(in) :: Cpg0
      DOUBLE PRECISION, INTENT(in) :: Cps0(DIM_M)
      DOUBLE PRECISION, INTENT(inout) :: MWg(DIM_N_g)
      DOUBLE PRECISION, INTENT(inout) :: MWs(DIM_M, DIM_N_s)

      LOGICAL :: CP_FATAL
      LOGICAL :: CHECK_DATABASE

      INTEGER :: M, N, lN


      CALL INIT_ERR_MSG("RXN_COM --> checkThermoReqs")

      CHECK_DATABASE = .FALSE.
      CP_FATAL = .FALSE.

! Verify that the molecular weights and stoichiometry are consistent and
! determine interphase mass exchanges.
      DO lN = 1, RxN%nSpecies
         M = RxN%Species(lN)%pMap
         N = RxN%Species(lN)%sMap
         IF(M == 0) THEN
            IF(Cpg0 /= UNDEFINED) THEN
               CP_FATAL = .TRUE.
            ELSEIF((RxN%Calc_DH .AND. .NOT.rDB(M,N)) .OR.        &
               (MWg(N) == UNDEFINED)) THEN
               CHECK_DATABASE = .TRUE.
            ENDIF
         ELSE
            IF(Cps0(M) /= UNDEFINED) THEN
               CP_FATAL = .TRUE.
            ELSEIF((RxN%Calc_DH .AND. .NOT.rDB(M,N)) .OR.        &
               (MWs(M,N) == UNDEFINED)) THEN
               CHECK_DATABASE = .TRUE.
            ENDIF
         ENDIF
      ENDDO

! Report errors and messages.
      IF(CP_FATAL) THEN

         WRITE(ERR_MSG, 1100) trim(RxN%Name)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1100: One or more phases associated with ',        &
         'reaction ',A,/'has specified constant specific heat (C_PG0/',&
         'Cps0). This is',/'not permitted for reacting phases. ',     &
         'Please correct the mfix.dat file.')

      ELSEIF(CHECK_DATABASE) THEN

         WRITE(ERR_MSG, 1101) trim(RxN%Name)
         CALL FLUSH_ERR_MSG

 1101 FORMAT('Message 1101: One or more molecular weights and/or ',    &
         'thermochemical data',/'is needed for reaction ',A,'. The ',  &
         'thermochemical database',/'will be used to gather the ',     &
         'necessary data.')

      ENDIF

      IF(CHECK_DATABASE) THEN
         WRITE(ERR_MSG, 1200)
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
      ENDIF

 1200 FORMAT('Message 1200: Searching thermochemical databases for ',&
         'species data.',/'  ')

      DO lN = 1, RxN%nSpecies
         M = RxN%Species(lN)%pMap
         N = RxN%Species(lN)%sMap
         IF(M == 0) THEN
            IF((RxN%Calc_DH .AND. .NOT.rDB(M,N)) .OR.         &
               (MWg(N) == UNDEFINED)) THEN
! Notify the user of the reason the thermochemical database is used.
! Flag that the species name is not provided.
               IF(S_g(N) == UNDEFINED_C) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('SPECIES_g',N))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF

! Update the log files.
               WRITE(ERR_MSG, 3001) N, trim(S_g(N))
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
! Read the database.
               CALL READ_DATABASE(0, N, S_g(N), MWg(N))
! Flag variable to stating that the database was read.
               rDB(0,N) = .TRUE.
            ENDIF
            RxN%Species(lN)%MW = MWg(N)
         ELSE
            IF((RxN%Calc_DH .AND. .NOT.rDB(M,N)) .OR.        &
               (MWs(M,N) == UNDEFINED)) THEN

! Flag that the species name is not provided.
               IF(S_s(M,N) == UNDEFINED_C) THEN
                  WRITE(ERR_MSG,1000) trim(iVar('SPECIES_s',M,N))
                  CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
               ENDIF
! Update the log files.
               WRITE(ERR_MSG, 3001) N, trim(S_s(M,N))
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
               CALL READ_DATABASE(M,N,S_s(M,N),MWs(M,N))
! Flag variable to stating that the database was read.
               rDB(M,N) = .TRUE.
            ENDIF
            RxN%Species(lN)%MW = MWs(M,N)
         ENDIF
      ENDDO
! Finalize the error message.
      IF(CHECK_DATABASE) CALL FLUSH_ERR_MSG(HEADER=.FALSE.)

 3001 FORMAT(/2x,'>',I3,': Species: ',A)

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('Error 1000: Required input not specified: ',A,/'Please ',&
            'correct the mfix.dat file.')

      END SUBROUTINE checkThermoReqs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: checkMassBalance                                        !
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
      SUBROUTINE checkMassBalance(CALLER, RxN, lnMT, IER)

      USE toleranc

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: CALLER

! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(INOUT) :: RxN

      DOUBLE PRECISION, INTENT(OUT) :: lnMT(0:DIM_M)
      INTEGER, INTENT(OUT) :: IER

      INTEGER M, N, lN ! Phase, Species
      DOUBLE PRECISION rSUM, pSUM
      DOUBLE PRECISION MWxStoich

      INTEGER sprCount, sprIndex

      DOUBLE PRECISION, PARAMETER :: massBalanceTol = 1.0d-3

! Initialize variables
      IER = 0
      rSUM = ZERO
      pSUM = ZERO
      lnMT(:) = ZERO
      sprCount = 0

! Verify that the molecular weights and stoichiometry are consistent and
! determine interphase mass exchanges.
      DO lN = 1, RxN%nSpecies
         M = RxN%Species(lN)%pMap
         N = RxN%Species(lN)%sMap

! Multiply the molecular weight and stoichiometric coefficient.
         MWxStoich = RxN%Species(lN)%MW * RxN%Species(lN)%Coeff
         RxN%Species(lN)%MWxStoich = MWxStoich
! Calculate the net mass transfer for phase M.
!  0 : no interphase mass transfder
! >0 : gains mass from anther phase
! <0 : transfers mass to anther phase
         lnMT(M) = lnMT(M) + MWxStoich
! Calculate mass of reactants and products. Used to ensure mass balance.
         IF(MWxStoich < ZERO) THEN
            rSUM = rSUM - MWxStoich
            IF(M /= 0) THEN
               sprCount = sprCount + 1
               IF(sprCount == 1) THEN
                  sprIndex = M
! Verify that there is at most one solids phase fule (reactant).
               ELSEIF( M /= sprIndex) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1002) trim(CALLER), trim(RxN%Name)
                     WRITE(UNIT_LOG,1002) trim(CALLER), trim(RxN%Name)
                     IER = 1
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            pSUM = pSUM + MWxStoich
         ENDIF
      ENDDO
! Verify that the mass of products equals reactants: (Mass Balance)
      IF (.NOT.COMPARE(rSUM,pSUM)) THEN
         IF(DMP_LOG) THEN
            WRITE(*,1001) trim(CALLER), trim(RxN%Name), rSUM, pSUM
            WRITE(UNIT_LOG,1001) trim(CALLER), trim(RxN%Name), rSUM,pSUM
            IER = 1
         ENDIF
      ENDIF

      RETURN

! Error Messages
!---------------------------------------------------------------------//

 1001 FORMAT(/1X,70('*')/' From: ',A,' --> RXN_COM -->',               &
         ' checkMassBalance',/' Error 1001: Stoichiometry is not',     &
         ' consistent with molecular weights',/' for reaction ',A,'.',/&
         ' Mass of reactants: ',F12.4,/' Mass of products:  ',F12.4,/  &
         1X,70('*')/)

 1002 FORMAT(/1X,70('*')/' From: ',A,' --> RXN_COM -->',               &
         ' checkMassBalance',/' Error 1002: More than one solids',     &
         ' phase fules was detected. Unable to',/' determine solids/', &
         'solids heat of reaction unambiguously for',/' reaction ',A,  &
         '.',/1X,70('*')/)

      END SUBROUTINE checkMassBalance

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Subroutine: calcInterphaseTxfr                                      !
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
      SUBROUTINE calcInterphaseTxfr(CALLER, RxN, lnMT, lEEq, lSEq, &
         lSAg, lMMx, lSAs)

         USE exit, only: mfix_exit
         USE toleranc

      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: CALLER

! Data structure for storing reaction data.
      TYPE(REACTION_BLOCK), POINTER, INTENT(INOUT) :: RxN

      DOUBLE PRECISION, INTENT(IN) :: lnMT(0:DIM_M)
! Energy equation flag
      LOGICAL, INTENT(IN) :: lEEq
! Gas/Solids Species Eq Flag
      LOGICAL, INTENT(IN) :: lSEq(0:DIM_M)
! Gas phase species aliases
      CHARACTER(len=32), DIMENSION(DIM_N_g), INTENT(IN) :: lSAg
! Number of solids phases
      INTEGER, INTENT(IN) :: lMMx
! Solids phase speices aliases.
      CHARACTER(len=32), DIMENSION(DIM_M, DIM_N_s), INTENT(IN) :: lSAs

      INTEGER toPhase, toPhaseCount, mCount
      INTEGER fromPhase, fromPhaseCount
      INTEGER catPhase

      INTEGER M, MM, LL
      INTEGER lM, lN

      DOUBLE PRECISION, PARAMETER :: massBalanceTol = 1.0d-3

! Initialize interphase exchange terms.
      IF(Allocated(RxN%rPhase)) RxN%rPhase(:) = ZERO

! If there is only one phase referenced by the reaction, there there
! should be no interphase mass transfer.
      IF(RxN%nPhases == 1) THEN
! Interphase mass transfer is set to zero. Small inconsistencies with
! molecular weights can result in a non-zero value for homogeneous
! reactions. Presumably, the mass balance check caught any major errors.
         RxN%rPhase(:) = ZERO
! Void interphase transfer flags.
         DO lN = 1, RxN%nSpecies
            M = RxN%Species(lN)%pMap
            RxN%Species(lN)%mXfr = M
         ENDDO
         RxN%Classification = "Homogeneous"
! This is a multiphase reaction.
      ELSE
! Initialize.
         toPhaseCount = 0
         fromPhaseCount = 0
         DO M = 0, lMMx
! Determine the number of phases with a net mass gain. Record the index
! of the last phase with a net mass gain.
            IF (lnMT(M) > massBalanceTol) THEN
               toPhaseCount = toPhaseCount + 1
               toPhase = M
! Determine the number of phases with a net mass loss. Record the index
! index of the last phase with a net mass loss.
            ELSEIF(lnMT(M) < -massBalanceTol) THEN
               fromPhaseCount = fromPhaseCount + 1
               fromPhase = M
            ENDIF
         ENDDO

! Only one phase has a net mass gain.
         IF(toPhaseCount == 1) THEN
! Interphase mass transfer flag.
            RxN%Classification = "Heterogeneous"
            DO M = 0, lMMx
               IF(M /= toPhase) THEN
                  IF (toPhase < M) THEN
                     LM = 1 + toPhase + ((M-1)*M)/2
                     RxN%rPhase(LM) = -lnMT(M)
                  ELSE
                     LM = 1 + M + ((toPhase-1)*toPhase)/2
                     RxN%rPhase(LM) = lnMT(M)
                  ENDIF

! Verify that if one phase's species equations are solved, that the
! other phase's species equations are solved.

                  IF(abs(RxN%rPhase(LM)) > SMALL_NUMBER) THEN
                     IF((lSEq(toPhase) .AND. .NOT.lSEq(M)) .OR. &
                        (.NOT.lSEq(toPhase) .AND. lSEq(M))) THEN
                        IF(DMP_LOG) THEN
                           WRITE(*,1001) trim(CALLER)
                           WRITE(UNIT_LOG,1001) trim(CALLER)
                           IF(lSEq(M)) THEN
                              WRITE(*,1101) M, 'Solving'
                              WRITE(*,1101) toPhase, 'Not Solving'
                              WRITE(UNIT_LOG,1101) M, 'Solving'
                              WRITE(UNIT_LOG,1101) toPhase, &
                                 'Not Solving'
                           ELSE
                              WRITE(*,1101) toPhase, 'Solving'
                              WRITE(*,1101) M, 'Not Solving'
                              WRITE(UNIT_LOG,1101) M, 'Solving'
                              WRITE(UNIT_LOG,1101) toPhase, &
                                 'Not Solving'
                           ENDIF
                           WRITE(*,1000)
                           WRITE(UNIT_LOG,1000)
                        ENDIF
                        CALL WRITE_RXN_SUMMARY(RxN, lSAg(:),           &
                           lSAs(:,:), .TRUE.)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

! Set flags for enthalpy transfer associated with mass transfer.
            IF(lEEq .AND. RxN%Calc_DH) THEN
               DO lN = 1, RxN%nSpecies
                  M = RxN%Species(lN)%pMap
! The gas phase is referenced by the reaction.
                  IF(M == 0) THEN
! The gas phase is the destination phase.
                     IF(toPhase == 0) THEN
! Counter for the number of solids phases transfering mass to the
! gas phase.
                        mCount = 0
! Check to see if phase M transfer mass to another solids phase.
                        DO MM = 1, lMMx
                           LM = 1 + (MM-1)*MM/2
! Mass transfer occurs between the gas and solids phase MM.
                           IF(RxN%rPhase(LM) > 0) THEN
! Indicate that phase MM receives mass from phase M.
                              RxN%Species(lN)%mXfr = MM
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                              RxN%Species(lN)%xXfr = ZERO
! Increment the number of phases the gas receives mass from.
                              mCount = mCount + 1
                           ENDIF
                        ENDDO
                        IF(mCount /= 1) THEN
                           IF(DMP_LOG) THEN
                              WRITE(*,1002) trim(CALLER), &
                                 trim(RxN%ChemEq)
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1002) trim(CALLER), &
                                 trim(RxN%ChemEq)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), &
                              lSAs(:,:), .TRUE.)
                        ENDIF

! A solids phase is the destination phase.
                     ELSE
! Since only one phase was detected with a net mass gain and the gas
! phase was detected as a source phase, then all the gas is assigned
! to the destination phase.
                        RxN%Species(lN)%mXfr = toPhase
! This variable is not used for gas/solids reactions.
                        RxN%Species(lN)%xXfr = ZERO
                     ENDIF
! Solids/Solids mass transfer: Enthalpy transfer from mass transfer is
! only calculated from source phase reactants.
                  ELSE
! Check to see if phase M transfer mass to another solids phase.
                     DO LL = 1, lMMx-1
                        DO MM = LL + 1, lMMx
                           IF(M /= LL .AND. M /= MM) CYCLE
                           LM = LL + 1 + (MM-1)*MM/2
                           IF(RxN%rPhase(LM) == ZERO) CYCLE
! Mass transfer occurs between solids phases M and MM, and phase M
! is the source phase.
                           IF( M == LL .AND. &
                              RxN%rPhase(LM) < ZERO) THEN
! Indicate that phase MM receives mass from phase M.
                              RxN%Species(lN)%mXfr = MM
! Calculate the fraction of material consumed from phase M is transfered
! to phase MM.
                              RxN%Species(lN)%xXfr =  &
                                 abs(lnMT(MM) / lnMT(LL))
! Mass transfer occurs between solids phases M and LL, and phase M
! is the source phase.
                           ELSEIF( M == MM .AND. &
                              RxN%rPhase(LM) > ZERO) THEN
! Indicate that phase LL receives mass from phase M.
                              RxN%Species(lN)%mXfr = LL
! Calculate the fraction of material consumed from phase M is transfered
! to phase LL.
                              RxN%Species(lN)%xXfr = &
                                 abs(lnMT(LL) / lnMT(MM))
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF ! Gas or Solids phase.
               ENDDO ! Species Loop
            ENDIF ! Energy Equation
! If there is only one phase with a net mass loss, setup the array for
! interphase mass transfer.
         ELSEIF(fromPhaseCount == 1) THEN
            RxN%Classification = "Heterogeneous"
            DO M = 0, lMMx
               IF (M /= fromPhase) THEN
                  IF(M < fromPhase) THEN
                     LM = 1 + M + ((fromPhase-1)*fromPhase)/2
                     RxN%rPhase(LM) =  lnMT(M)
                  ELSE
                     LM = 1 + fromPhase + ((M-1)*M)/2
                     RxN%rPhase(LM) = -lnMT(M)
                  ENDIF

! Verify that if one phase's species equations are solved, that the
! other phase's species equations are solved.
                  IF(abs(RxN%rPhase(LM)) > SMALL_NUMBER) THEN
                     IF((lSEq(fromPhase) .AND. .NOT.lSEq(M)) .OR.   &
                        (.NOT.lSEq(fromPhase) .AND. lSEq(M))) THEN
                        IF(DMP_LOG) THEN
                           WRITE(*,1001) trim(CALLER)
                           WRITE(UNIT_LOG,1001) trim(CALLER)
                           IF(lSEq(M)) THEN
                              WRITE(*,1101) M, 'Solving'
                              WRITE(*,1101) fromPhase, 'Not Solving'
                              WRITE(UNIT_LOG,1101) M, 'Solving'
                              WRITE(UNIT_LOG,1101) fromPhase, &
                                 'Not Solving'
                           ELSE
                              WRITE(*,1101) toPhase, 'Solving'
                              WRITE(*,1101) M, 'Not Solving'
                              WRITE(UNIT_LOG,1101) fromPhase, 'Solving'
                              WRITE(UNIT_LOG,1101) M, 'Not Solving'
                           ENDIF
                           WRITE(*,1000)
                           WRITE(UNIT_LOG,1000)
                        ENDIF
                        CALL WRITE_RXN_SUMMARY(RxN, lSAg(:),           &
                           lSAs(:,:), .TRUE.)
                     ENDIF
                  ENDIF
               ENDIF
            END DO

! Set flags for enthalpy transfer associated with mass transfer.
            IF(lEEq .AND. RxN%Calc_DH) THEN
               DO lN = 1, RxN%nSpecies
                  M = RxN%Species(lN)%pMap
! Gas/solids reaction: Enthalpy transfer from mass transfer is only
! calculated from gas phase species.
                  IF(M == 0) THEN
! Gas phase is the source phase.
                     IF(fromPhase == 0) THEN
! Counter for the number of solids phases transfering mass to the
! gas phase.
                        mCount = 0
! Check to see if phase M transfer mass to another solids phase.
                        DO MM = 1, lMMx
                           LM = 1 + (MM-1)*MM/2
! Mass transfer occurs between the gas and solids phases MM.
                           IF(RxN%rPhase(LM) < 0) THEN
! Indicate that phase MM receives mass from phase M.
                              RxN%Species(lN)%mXfr = MM
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                              RxN%Species(lN)%xXfr = ZERO
! Increment the number of phases the gas receives mass from.
                              mCount = mCount + 1
                           ENDIF
                        ENDDO
                        IF(mCount /=1 ) THEN
                           IF(DMP_LOG) THEN
                              WRITE(*,1002) trim(CALLER), &
                                 trim(RxN%ChemEq)
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1002) trim(CALLER), &
                                 trim(RxN%ChemEq)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY(RxN, lSAg(:),  &
                              lSAs(:,:), .TRUE.)
                        ENDIF
                     ELSE
! There can be only one solids phase fuel. Store the phase of the
! solids phase reactant.
                        RxN%Species(lN)%mXfr = fromPhase
! Mass fraction of transfered material.
! This variable is not currently used for gas/solids reactions.
                        RxN%Species(lN)%xXfr = ZERO
                     ENDIF
! Solids/Solids mass transfer: Enthalpy transfer from mass transfer is
! only calculated from source phase reactants.
                  ELSE
! Check to see if phase M transfer mass to another solids phase.
                     DO LL = 1, lMMx-1
                        DO MM = LL + 1, lMMx
                           IF(M /= LL .AND. M /= MM) CYCLE
                           LM = LL + 1 + (MM-1)*MM/2
                           IF(RxN%rPhase(LM) == ZERO) CYCLE
! Mass transfer occurs between solids phases M and MM, and phase M
! is the source phase.
                           IF( M == LL .AND. &
                              RxN%rPhase(LM) < ZERO) THEN
! Indicate that phase MM receives mass from phase M.
                              RxN%Species(lN)%mXfr = MM
! Calculate the fraction of material consumed from phase M is transfered
! to phase MM.
                              RxN%Species(lN)%xXfr = &
                                 abs(lnMT(MM) / lnMT(LL))
! Mass transfer occurs between solids phases M and LL, and phase M
! is the source phase.
                           ELSEIF( M == MM .AND. &
                              RxN%rPhase(LM) > ZERO) THEN
! Indicate that phase LL receives mass from phase M.
                              RxN%Species(lN)%mXfr = LL
! Calculate the fraction of material consumed from phase M is transfered
! to phase LL.
                              RxN%Species(lN)%xXfr = &
                                 abs(lnMT(LL) / lnMT(MM))
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF ! Gas or Solids phase.
               ENDDO ! Species Loop
            ENDIF ! Energy Equation

! If there are no phases with a net mass gain/loss, check to see if
! the reaction is turned off.
         ELSEIF(toPhaseCount == 0 .AND. fromPhaseCount == 0) THEN
! If the reaction is active, and there is no interphase mass transfer,
! classify the reaction as catalytic.
            IF(RxN%nPhases > 0) RxN%Classification = "Catalytic"
            RxN%rPhase(:)  = ZERO
! Set flags for enthalpy transfer associated with mass transfer.
            IF(lEEq .AND. RxN%Calc_DH) THEN

! Identify the catalyst phase.
               catPhase = -1
               DO lN= 1, RxN%nSpecies
                  IF(COMPARE(RxN%Species(lN)%Coeff,ZERO)) THEN
                     IF(catPhase /= -1) THEN
                        IF(catPhase /= RxN%Species(lN)%pMap) THEN
                           IF(DMP_LOG) THEN
                              WRITE(*,1002) trim(CALLER), &
                                 trim(RxN%Name)
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1002) trim(CALLER), &
                                 trim(RxN%Name)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY(RxN, lSAg(:),    &
                              lSAs(:,:), .TRUE.)
                        ENDIF
                     ELSE
                        catPhase = RxN%Species(lN)%pMap
                     ENDIF
                  ENDIF
               ENDDO
! Verify that a catalyst phase was found.
               IF(catPhase == -1) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1003) trim(CALLER), 'catalyst', &
                        trim(RxN%Name)
                     WRITE(*,1000)
                     WRITE(UNIT_LOG,1003) trim(CALLER), &
                        'catalyst', trim(RxN%Name)
                     WRITE(UNIT_LOG,1000)
                  ENDIF
                  CALL WRITE_RXN_SUMMARY(RxN, lSAg(:),                 &
                     lSAs(:,:), .TRUE.)
               ENDIF

! Identify the reactant phase.
               toPhase = -1
               DO lN = 1, RxN%nSpecies
                  IF(.NOT.COMPARE(RxN%Species(lN)%Coeff,ZERO)) THEN
                     IF(toPhase /= -1) THEN
                        IF(toPhase /= RxN%Species(lN)%pMap) THEN
                           IF(DMP_LOG) THEN
                              WRITE(*,1002) trim(CALLER), &
                                 trim(RxN%Name)
                              WRITE(*,1000)
                              WRITE(UNIT_LOG,1002) trim(CALLER), &
                                 trim(RxN%Name)
                              WRITE(UNIT_LOG,1000)
                           ENDIF
                           CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), &
                              lSAs(:,:), .TRUE.)
                        ENDIF
                     ELSE
                        toPhase = RxN%Species(lN)%pMap
                     ENDIF
                  ENDIF
               ENDDO
! Verify that a reacting phase was found.
               IF(toPhase == -1) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1003) trim(CALLER), 'reacting', &
                        trim(RxN%Name)
                     WRITE(*,1000)
                     WRITE(UNIT_LOG,1003) trim(CALLER), 'reacting', &
                        trim(RxN%Name)
                     WRITE(UNIT_LOG,1000)
                  ENDIF
                  CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), lSAs(:,:),.TRUE.)
               ENDIF

! Something when wrong.
               IF(catPhase == toPhase) THEN
                  IF(DMP_LOG) THEN
                     WRITE(*,1004) trim(CALLER), trim(RxN%Name)
                     WRITE(*,1000)
                     WRITE(UNIT_LOG,1004) trim(CALLER),trim(RxN%Name)
                     WRITE(UNIT_LOG,1000)
                  ENDIF
                  CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), lSAs(:,:),.TRUE.)
!Gas/solid catalytic reaction:
               ELSEIF(toPhase == 0) THEN
                  DO lN = 1, RxN%nSpecies
                     IF(RxN%Species(lN)%pMap == 0) THEN
! Indicate that phase MM receives mass from phase M.
                        RxN%Species(lN)%mXfr = catPhase
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                        RxN%Species(lN)%xXfr = ZERO
                     ENDIF
                  ENDDO
               ELSEIF(catPhase == 0) THEN
                  DO lN = 1, RxN%nSpecies
                     IF(RxN%Species(lN)%pMap == 0) THEN
! Indicate that phase MM receives mass from phase M.
                        RxN%Species(lN)%mXfr = toPhase
! The fraction of material transfered from phase 0 to phase MM.
! This variable is not currently used for gas/solids reactions.
                        RxN%Species(lN)%xXfr = ZERO
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF ! Energy Equation
         ELSE
! Two or more phases have a net mass loss and two or more phases have
! a net mass gain. Therefore, the interphase mass transfer cannot be
! concluded.
            CALL WRITE_RXN_SUMMARY(RxN, lSAg(:), lSAs(:,:),.FALSE.)
            WRITE(*,1002) trim(CALLER), trim(RxN%ChemEq)
            WRITE(*,1000)
            WRITE(UNIT_LOG,1002) trim(CALLER), trim(RxN%ChemEq)
            WRITE(UNIT_LOG,1000)
            CALL MFiX_EXIT(myPE)
         ENDIF
      ENDIF

      RETURN

! Error Messages
!---------------------------------------------------------------------//

 1000 FORMAT(/' Please refer to the Readme file on the required input',&
         ' format and make',/' the necessary corrections to the data', &
         ' file.',/1X,70('*')//)

 1001 FORMAT(//1X,70('*')/' From: ',A,' --> RXN_COM -->',              &
         ' calcInterphaseTxfr',/' Error 1001: A chemical reaction or', &
         ' phase change was detected between',/' a phases solving',    &
         ' species equations and another phase not solving',/          &
         ' species equations.',/)

 1101 FORMAT(' Phase ',I2,': ',A,' species equations.')

 1002 FORMAT(//1X,70('*')/' From: ',A,' --> RXN_COM -->',              &
         ' calcInterphaseTxfr',/' Error 1002: Reaction complexity',    &
         ' exceeds implementation capabilities.',/' Unable to',        &
         ' determine unambiguously interphase heat or mass transfer.', &
         //' Reaction: ',A,//' Consider splitting the chemical',       &
         ' reaction equation into two or more',/' separate equations.',&
         ' The same reaction rate calculated in usr_rates',/' can be', &
         ' used for the multiple reactions to ensure mass')

 1003 FORMAT(//1X,70('*')/' From: ',A,' --> RXN_COM -->',              &
         ' calcInterphaseTxfr',/' Error 1003: Unable to determine ',A, &
         ' phase for catalytic reaction'/1X,A,'.')

 1004 FORMAT(//1X,70('*')/' From: ',A,' --> RXN_COM -->',              &
         ' calcInterphaseTxfr',/' Error 1004: Unable to distinguish',  &
         ' catalyst phase from reacting phase',/' for catalytic',      &
         ' reaction ',A,'.')

      END SUBROUTINE calcInterphaseTxfr


END MODULE RXN_COM
