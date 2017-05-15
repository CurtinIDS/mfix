!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!     Module name: READ_NAMELIST(POST)                                 !
!     Author: P. Nicoletti                            Date: 25-NOV-91  !
!                                                                      !
!     Purpose: Read in the NAMELIST variables                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE READ_NAMELIST(READ_ACTION, FILENAME)

      USE bc
      USE cdist
      USE compar
      USE constant
      USE cutcell
      USE dashboard
      USE des_bc
      USE des_rxns
      USE des_thermo
      USE discretelement
      USE error_manager
      USE fldvar
      USE funits
      USE geometry
      USE ic
      USE indices
      USE is
      USE iterate, only: max_nit
      USE leqsol
      USE main, only: cmd_line_args, cmd_line_args_count
      USE mfix_pic
      USE output
      USE parallel
      USE param1, only: undefined
      USE particle_filter
      USE physprop
      USE pic_bc
      USE polygon
      USE ps
      USE qmom_kinetic_equation
      USE quadric
      USE residual
      USE run
      USE rxns
      USE scalars
      USE scales
      USE stiff_chem
      USE toleranc
      USE ur_facs
      USE usr
      USE utilities
      USE vtk
      Use stl
      use usr_prop, only: usr_fgs, usr_fss, usr_gama
      use usr_prop, only: usr_rog, usr_cpg, usr_mug, usr_kg, usr_difg
      use usr_prop, only: usr_ros, usr_cps, usr_mus, usr_ks, usr_difs
      use usr_src, only: call_usr_source
      IMPLICIT NONE

! Dummy Arguments:
!------------------------------------------------------------------------//
! Specify how much of the input to process.
      INTEGER, INTENT(IN) :: READ_ACTION

! Filename of the input file
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME

! Local Variables:
!------------------------------------------------------------------------//
! LINE_STRING(1:MAXCOL) has valid input data
      INTEGER, PARAMETER :: MAXCOL = 80
! Holds one line in the input file
      CHARACTER(LEN=512) :: LINE_STRING
! Length of noncomment string
      INTEGER :: LINE_LEN
! Line number
      INTEGER :: LINE_NO
! Coefficient of restitution (old symbol)
      DOUBLE PRECISION   e
! Indicates whether currently reading rxns or rate
      LOGICAL :: RXN_FLAG
! Indicate whether to do a namelist read on the line
      LOGICAL :: READ_FLAG
! Logical to check if file exits.
      LOGICAL :: lEXISTS
! Error flag
      LOGICAL :: ERROR

      CHARACTER(len=256) :: STRING
      INTEGER :: IOS, II

! Flags restricting what data from the mfix.dat to process
      LOGICAL :: READ_LOCKED, READ_FULL

! Local Parameters:
!---------------------------------------------------------------------//
      INTEGER, PARAMETER :: READ_MFIX = 0
      INTEGER, PARAMETER :: READ_POST = 1
      INTEGER, PARAMETER :: READ_INIT = 2

! External Functions
!---------------------------------------------------------------------//
! Returns integer if data past column MAXCOL.
!      INTEGER, EXTERNAL :: LINE_TOO_BIG
! Integer function which returns COMMENT_INDEX
!      INTEGER, EXTERNAL :: SEEK_COMMENT
! Blank line function
!      LOGICAL, EXTERNAL :: BLANK_LINE

      E = UNDEFINED
      RXN_FLAG = .FALSE.
      READ_FLAG = .TRUE.
      NO_OF_RXNS = 0
      LINE_NO = 0

      SELECT CASE(READ_ACTION)
      CASE(READ_MFIX)
         READ_LOCKED = .TRUE.
         READ_FULL = .TRUE.
      CASE(READ_POST)
         READ_LOCKED = .TRUE.
         READ_FULL = .FALSE.
      CASE(READ_INIT)
         READ_LOCKED = .FALSE.
         READ_FULL = .TRUE.
      END SELECT

! Open the mfix.dat file. Report errors if the file is not located or
! there is difficulties opening it.
      inquire(file=filename,exist=lEXISTS)
      IF(.NOT.lEXISTS) THEN
         IF(myPE == PE_IO) WRITE(*,1000)
         CALL MFIX_EXIT(myPE)

 1000 FORMAT(2/,1X,70('*')/' From: READ_NAMELIST',/' Error 1000: ',    &
         'The input data file, mfix.dat, is missing. Aborting.',/1x,   &
         70('*'),2/)

      ELSE
         OPEN(UNIT=UNIT_DAT, FILE=filename, STATUS='OLD', IOSTAT=IOS)
         IF(IOS /= 0) THEN
            IF(myPE == PE_IO) WRITE (*,1001)
            CALL MFIX_EXIT(myPE)
         ENDIF

 1001 FORMAT(2/,1X,70('*')/' From: READ_NAMELIST',/' Error 1001: ',    &
         'Unable to open the mfix.dat file. Aborting.',/1x,70('*'),2/)
      ENDIF


! Loop through the mfix.dat file and process the input data.
      READ_LP: DO
         READ (UNIT_DAT,"(A)",IOSTAT=IOS) LINE_STRING
         IF(IOS < 0) EXIT READ_LP

         LINE_NO = LINE_NO + 1

         LINE_LEN = SEEK_COMMENT(LINE_STRING,LEN(LINE_STRING)) - 1
         CALL REMOVE_COMMENT(LINE_STRING, LINE_LEN+1, LEN(LINE_STRING))

         IF(LINE_LEN <= 0) CYCLE READ_LP           ! comment line
         IF(BLANK_LINE(LINE_STRING)) CYCLE READ_LP ! blank line

         IF(LINE_TOO_BIG(LINE_STRING,LINE_LEN,MAXCOL) > 0) THEN
            WRITE (*, 1100) trim(iVAL(LINE_NO)), trim(ival(MAXCOL)), &
               LINE_STRING(1:MAXCOL)
            CALL MFIX_EXIT(myPE)
         ENDIF

 1100 FORMAT(//1X,70('*')/1x,'From: READ_NAMELIST',/1x,'Error 1100: ', &
         'Line ',A,' in mfix.dat has is too long. Input lines should', &
         /1x,'not pass column ',A,'.',2/3x,A,2/1x,'Please correct ',   &
         'the mfix.dat file.',/1X,70('*'),2/)

! All subsequent lines are thermochemical data
         IF(LINE_STRING(1:11) == 'THERMO DATA') EXIT READ_LP

         CALL SET_KEYWORD(ERROR)
         IF (ERROR) THEN
! At this point, the keyword was not identified therefore it is
! either deprecated or unknown.
            CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))
         ENDIF

      ENDDO READ_LP

      DO II=1, CMD_LINE_ARGS_COUNT
         LINE_STRING = CMD_LINE_ARGS(ii)
         LINE_LEN = len(line_string)
         CALL SET_KEYWORD(ERROR)
         IF (ERROR) THEN
            CALL DEPRECATED_OR_UNKNOWN(LINE_NO, LINE_STRING(1:LINE_LEN))
         ENDIF
      ENDDO

      CLOSE(UNIT=UNIT_DAT)
      IF (E /= UNDEFINED) C_E = E

      RETURN

      CONTAINS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
! Subroutine: SET_KEYWORD(ERROR)                                       !
! Author: P. Nicoletti                                Date: 25-NOV-91  !
!                                                                      !
! Purpose: Process LINE_STRING for MFIX keyword data.                  !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE SET_KEYWORD(ERROR)

      IMPLICIT NONE

      LOGICAL, INTENT(OUT) ::ERROR



! External namelist files:
!---------------------------------------------------------------------//
      INCLUDE 'run_control.inc'
      INCLUDE 'physical_params.inc'
      INCLUDE 'numerical_params.inc'
      INCLUDE 'geometry.inc'
      INCLUDE 'gas_phase.inc'
      INCLUDE 'solids_phase.inc'
      INCLUDE 'tfm_solids.inc'
      INCLUDE 'initial_conditions.inc'
      INCLUDE 'boundary_conditions.inc'
      INCLUDE 'internal_surfaces.inc'
      INCLUDE 'point_sources.inc'
      INCLUDE 'output_control.inc'
      INCLUDE 'usr_hooks.inc'
      INCLUDE 'chem_equations.inc'
      INCLUDE 'dmp_batch_control.inc'
      INCLUDE 'desnamelist.inc'
      INCLUDE 'cartesian_grid_namelist.inc'
      INCLUDE 'qmomknamelist.inc'
      INCLUDE 'legacy.inc'
      INCLUDE 'usrnlst.inc'

      ERROR = .FALSE.

! Make upper case all except species names
      if(index(LINE_STRING,'SPECIES_NAME') == 0 .AND. &
         index(LINE_STRING,'species_name') == 0 .AND. &
         index(LINE_STRING,'Species_Name') == 0 .AND. &
         index(LINE_STRING,'SPECIES_g') == 0 .AND.    &
         index(LINE_STRING,'Species_g') == 0 .AND.    &
         index(LINE_STRING,'species_g') == 0 .AND.    &
         index(LINE_STRING,'SPECIES_s') == 0 .AND.    &
         index(LINE_STRING,'Species_s') == 0 .AND.    &
         index(LINE_STRING,'species_s') == 0)         &
         CALL MAKE_UPPER_CASE (LINE_STRING, LINE_LEN)

      CALL REPLACE_TAB (LINE_STRING, LINE_LEN)
      CALL REMOVE_PAR_BLANKS(LINE_STRING)

! Complete arithmetic operations and expand line
      CALL PARSE_LINE (LINE_STRING, LINE_LEN, RXN_FLAG, READ_FLAG)

! Write the current line to a scratch file
! and read the scratch file in NAMELIST format
      IF(.NOT.READ_FLAG) RETURN


! Run control keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&RUN_CONTROL_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=RUN_CONTROL_LOCKED,  IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&RUN_CONTROL_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=RUN_CONTROL_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Physical parameter keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&PHYSICAL_PARAM_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=PHYSICAL_PARAM_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&PHYSICAL_PARAM_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=PHYSICAL_PARAM_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Numerical parameter keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&NUMERICAL_PARAM_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=NUMERICAL_PARAM_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&NUMERICAL_PARAM_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=NUMERICAL_PARAM_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Geometry and discretization keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&GEOMETRY_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=GEOMETRY_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&GEOMETRY_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=GEOMETRY_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Gas phase keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&GAS_PHASE_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=GAS_PHASE_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&GAS_PHASE_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=GAS_PHASE_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Solidss phase keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&SOLIDS_PHASE_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=SOLIDS_PHASE_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&SOLIDS_PHASE_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=SOLIDS_PHASE_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Two-fluid solids keywords
      STRING=''; STRING = '&TFM_SOLIDS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=TFM_SOLIDS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Initial condtion keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&INITIAL_CONDITIONS_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=INITIAL_CONDITIONS_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&INITIAL_CONDITIONS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=INITIAL_CONDITIONS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Boundary condition keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&BOUNDARY_CONDITIONS_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=BOUNDARY_CONDITIONS_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&BOUNDARY_CONDITIONS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=BOUNDARY_CONDITIONS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Internal surface keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&INTERNAL_SURFACES_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=INTERNAL_SURFACES_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&INTERNAL_SURFACES_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=INTERNAL_SURFACES_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Point source keywords
      STRING=''; STRING = '&POINT_SOURCES_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=POINT_SOURCES_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Output control keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&OUTPUT_CONTROL_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=OUTPUT_CONTROL_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&OUTPUT_CONTROL_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=OUTPUT_CONTROL_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! User hook keywords
      STRING=''; STRING = '&USER_HOOKS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=USER_HOOKS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Chemical equation keywords
      STRING=''; STRING = '&CHEM_EQUATIONS_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=CHEM_EQUATIONS_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! DMP and Batch Queue control keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&DMP_BATCH_CONTROL_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=DMP_BATCH_CONTROL_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      STRING=''; STRING = '&DMP_BATCH_CONTROL_UNLOCKED '//&
         trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
      READ(STRING, NML=DMP_BATCH_CONTROL_UNLOCKED, IOSTAT=IOS)
      IF(IOS == 0)  RETURN


! Legacy keywords
      IF(READ_LOCKED) THEN
         STRING=''; STRING = '&LEGACY_LOCKED '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=LEGACY_LOCKED, IOSTAT=IOS)
         IF(IOS == 0)  RETURN

      ENDIF

! Stop processing keyword inputs if runing POST_MFIX
       IF(.NOT.READ_FULL) RETURN


      IF(READ_LOCKED) THEN

! Discrete Element model input parameters.
         STRING=''; STRING = '&DES_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=DES_INPUT_DATA, IOSTAT=IOS)
         IF(IOS == 0)  RETURN


! User defined input parameters.
         STRING=''; STRING = '&USR_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=USR_INPUT_DATA, IOSTAT=IOS)
         IF(IOS == 0)  RETURN


! Cartesian grid cut-cell input parameters.
         STRING=''; STRING = '&CARTESIAN_GRID_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=CARTESIAN_GRID_INPUT_DATA, IOSTAT=IOS)
         IF(IOS == 0)  RETURN


! QMOMK input parameters.
         STRING=''; STRING = '&QMOMK_INPUT_DATA '//&
            trim(adjustl(LINE_STRING(1:LINE_LEN)))//'/'
         READ(STRING, NML=QMOMK_INPUT_DATA, IOSTAT=IOS)
         IF(IOS == 0)  RETURN
      ENDIF

      IF(READ_LOCKED) ERROR = .TRUE.

      RETURN
      END SUBROUTINE SET_KEYWORD

END SUBROUTINE READ_NAMELIST
