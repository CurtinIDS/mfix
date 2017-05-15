!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_CHEMICAL_RXNS                                    !
!  Author: J.Musser                                   Date: 21-MAR-14  !
!                                                                      !
!  Purpose: Check chemical reactions specifications                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CHEMICAL_RXNS

! User defined reaction names from reaction blocks @(RXNS)
      use parse, only: RXN_NAME, DES_RXN_NAME

      use parse, only: RXN_CHEM_EQ, DES_RXN_CHEM_EQ
      use parse, only: usrDH, DES_usrDH
      use parse, only: usrfDH, DES_usrfDH


! Number of continuum reactions and data object
      use rxns, only: NO_OF_RXNS, REACTION
! Number of discrete reactions and data object
      use des_rxns, only: NO_OF_DES_RXNS, DES_REACTION
! User specified species names and aliases:
      use rxns, only: SPECIES_g, SPECIES_ALIAS_g
      use rxns, only: SPECIES_s, SPECIES_ALIAS_s

! Number of continuum solids
      use physprop, only: SMAX
! Number of discrete solids
      use discretelement, only: DES_MMAX
! Number of species comprising each phase
      use physprop, only: NMAX

!
      use param1, only: UNDEFINED_I

      use rxn_com, only: checkSpeciesInc
      use rxn_com, only: checkDuplicateAliases


      use error_manager

      IMPLICIT NONE

! Error flag
      INTEGER :: IER

! Local representation of the number of solids phases.
      INTEGER :: lMMAX

! Undefined indicates that no reaction block was found in the deck file.
      IF(NO_OF_RXNS == UNDEFINED_I) NO_OF_RXNS = 0
      IF(NO_OF_DES_RXNS == UNDEFINED_I) NO_OF_DES_RXNS = 0

! If there are no chemical reactions, then skip this routine.
      IF(NO_OF_RXNS + NO_OF_DES_RXNS == 0) RETURN

      CALL INIT_ERR_MSG('CHECK_CHEMICAL_RXNS')

! Allocate the arrays as a work-around for Intel compiler with debug.
      IF(NO_OF_RXNS == 0) allocate(RXN_NAME(1))
      IF(NO_OF_DES_RXNS == 0) allocate(DES_RXN_NAME(1))

! Initialize the number of solids phases.
      lMMAX = SMAX + DES_MMAX

! Verify that the species aliases are unique.
      CALL checkDuplicateAliases(NMAX(0), SPECIES_ALIAS_g(:), &
         lMMAX, NMAX(1:lMMAX), SPECIES_ALIAS_s(:,:))

! Verify that species aliases in the datafile match those in the
! species.inc file.
      CALL checkSpeciesInc(NMAX(0), SPECIES_ALIAS_g(:), lMMAX,         &
         NMAX(1:lMMAX), SPECIES_ALIAS_s(:,:), NO_OF_RXNS, RXN_NAME(:), &
         NO_OF_DES_RXNS, DES_RXN_NAME)


      IF(NO_OF_RXNS > 0) THEN
         CALL CHECK_CHEMICAL_RXNS_COMMON(NO_OF_RXNS, RXN_NAME,         &
            RXN_CHEM_EQ, usrDH, usrfDH, REACTION)
      ENDIF

      IF(NO_OF_DES_RXNS > 0) THEN
         CALL CHECK_CHEMICAL_RXNS_COMMON(NO_OF_DES_RXNS, DES_RXN_NAME, &
            DES_RXN_CHEM_EQ, DES_usrDH, DES_usrfDH, DES_REACTION)
      ENDIF

      CALL FINL_ERR_MSG

      RETURN

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: CHECK_CHEMICAL_RXNS_COMMON                             !
!  Author: J.Musser                                   Date: 21-MAR-14  !
!                                                                      !
!  Purpose: Check chemical reactions specifications                    !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE CHECK_CHEMICAL_RXNS_COMMON(COUNT, NAME, CHEM_EQ, DH,  &
         fDH, RXN_ARRAY)

! Runtime flags for solving energy and species equations.
      use run, only: ENERGY_EQ, SPECIES_EQ
! Definition of derived data type
      use rxn_com, only: REACTION_BLOCK
! User specified constant specific heats:
      use physprop, only: C_PG0, C_PS0
! Molecular weights:
      use physprop, only: MW_g, MW_s
! Flag marking when the thermochemical database is read.
      use rxns, only: rDatabase

      use param, only: DIM_M, DIMENSION_RXN


      use parse, only: setReaction
      use rxn_com, only: checkThermoReqs
      use rxn_com, only: checkMassBalance
      use rxn_com, only: calcInterphaseTxfr
      use rxn_com, only: WRITE_RXN_SUMMARY



      use error_manager

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: COUNT

! Reaction Names from mfix.dat file:
      CHARACTER(len=32), INTENT(IN) ::  NAME(DIMENSION_RXN)
! Chemical equations:
      CHARACTER(len=512), INTENT(IN) :: CHEM_EQ(DIMENSION_RXN)
! User defined heat of reaction:
      DOUBLE PRECISION, INTENT(IN) :: DH(DIMENSION_RXN)
! User defined heat of reaction partitions.
      DOUBLE PRECISION, INTENT(IN) :: fDH(DIMENSION_RXN,0:DIM_M)

! Array of reaction data objects.
      TYPE(REACTION_BLOCK), TARGET, ALLOCATABLE :: RXN_ARRAY(:)



! loop/variable indices1
      INTEGER :: L

      TYPE(REACTION_BLOCK), POINTER :: This

      DOUBLE PRECISION netMassTransfer(0:DIM_M)




      CALL INIT_ERR_MSG('CHECK_CHEMICAL_RXNS_COMMON')


! Allocate reaction blocks.
      allocate(RXN_ARRAY(COUNT))





! Loop over reaction data pulled from data file.
      DO L=1, COUNT

         This => RXN_ARRAY(L)

! Store the reaction name.
         This%Name = trim(NAME(L))

! This check should not be necessary. Pre-processing by make_mfix and
! reading the data file (PARSE_RXN) should have already caught any
! issues.
         IF(len_trim(This%Name) == 0) THEN
            WRITE(ERR_MSG, 1100) L
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1100 FORMAT('Error 1100: No reaction name identified for reaction ',  &
         I3,'.',/'This should have been caught during the parsing ',   &
         'routines.')

! Store the chemical equation.
         This%ChemEq = trim(CHEM_EQ(L))

! Verify that a chemical equation was given in the data file.
         IF(len_trim(This%ChemEq) == 0) THEN
            WRITE(*,1101) trim(This%Name)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1101 FORMAT('Error 1101: No chemical equation identified for ',       &
         'reaction ',A,'.',/'Please correct the mfix.dat file.')

! Take the data read from the data file and populate the reaction block.
         CALL setReaction(This, NMAX(0), SPECIES_ALIAS_g(:),lMMAX,     &
            NMAX(1:lMMAX), SPECIES_ALIAS_s(:,:), DH(L), fDH(L,:))

! If the energy equations are not being solved and a user provided
! heat of reaction is given, flag error and exit.
         IF(.NOT.ENERGY_EQ .AND. .NOT.This%Calc_DH) THEN
            WRITE(ERR_MSG,1200) trim(This%Name)
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1200 FORMAT('Error 1200: Inconsistent user input. Energy equations ', &
         'are NOT being',/'solved and a user defined heat of reaction',&
         ' was detected',' for chemical',/' reaction ',A,'.',/'Please',&
         ' correct the mfix.dat file.')

! Skip empty reactions.
         IF(This%nSpecies == 0 .AND. This%nPhases == 0) THEN
            CYCLE

! Something went wrong while parsing the reaction. This is a sanity
! check and should never be true.
         ELSEIF((This%nPhases == 0 .AND. This%nSpecies /= 0) .OR.      &
            (This%nPhases /= 0 .AND. This%nSpecies == 0)) THEN

            WRITE(ERR_MSG,1201) trim(This%Name), This%nPhases,         &
               This%nSpecies
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
         ENDIF

 1201 FORMAT('Error 1201: Illogical data returned from setReaction ',  &
         'for chemical',/'reaction',1X,A,'.',//' Number of phases ',   &
         'identified: ',I3,/'Number of species identified: ',I3,//,    &
         'Please check the mfix.dat file.')

! Verify that the necessary information for each species in the reaction
! was defined.
         CALL checkThermoReqs(This, SPECIES_g, SPECIES_s, rDatabase,   &
            MW_G, MW_S, C_PG0, C_PS0)


! Verify Mass Balance (Mass of Reactants = Mass of Products)
!---------------------------------------------------------------------//
         IER = 0
         CALL checkMassBalance('CHECK_CHEMICAL_RXNS', This, &
            netMassTransfer(:), IER)
         IF(IER /= 0) THEN
            CALL WRITE_RXN_SUMMARY(This, SPECIES_ALIAS_g(:), &
               SPECIES_ALIAS_s(:,:), .TRUE.)
         ENDIF

! Determine interphase exchanges
!---------------------------------------------------------------------//
         CALL calcInterphaseTxfr('CHECK_CHEMICAL_RXNS', This,   &
            netMassTransfer(:), ENERGY_EQ, SPECIES_EQ(:), &
            SPECIES_ALIAS_g(:), lMMAX, SPECIES_ALIAS_s(:,:))
      ENDDO

! Write a summary of the chemical reactions
!---------------------------------------------------------------------//
      DO L=1, COUNT
         This => RXN_ARRAY(L)
         CALL WRITE_RXN_SUMMARY(This, SPECIES_ALIAS_g(:), &
            SPECIES_ALIAS_s(:,:), .FALSE.)
      ENDDO

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CHECK_CHEMICAL_RXNS_COMMON




      END SUBROUTINE CHECK_CHEMICAL_RXNS

