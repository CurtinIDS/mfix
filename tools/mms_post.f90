      MODULE WRITE_UTILITY

      IMPLICIT NONE

! Maximum number of lines a message can have before a flush is needed.
      INTEGER, PARAMETER :: LINE_COUNT  = 32
! Maximum number of characters per line.
      INTEGER, PARAMETER :: LINE_LENGTH = 256
! Output buffer.
      CHARACTER(LEN=LINE_LENGTH) :: OUT(LINE_COUNT)

      INTEGER :: WORKING_CASE
      INTEGER, PARAMETER :: LOWER_CASE = 0
      INTEGER, PARAMETER :: UPPER_CASE = 1 

      INTEGER, PARAMETER :: UNIT_MAKE = 150
      INTEGER, PARAMETER :: UNIT_TEMP = 250

! Flag to include MPIF
      LOGICAL :: INCLUDE_MPIF

! Back slash
      CHARACTER(LEN=1), PARAMETER :: BS = char(92)


      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: SET_WORKING_CASE                                         !
!                                                                      !
! Purpose: Manages flushing the output buffer (OUT) so that if output  !
! needs multiple redirects, it only needs one call to this routine.    !
!......................................................................!
      SUBROUTINE SET_WORKING_CASE

      INTEGER :: ENV_STAT, ENV_LEN
      CHARACTER(LEN=32) :: MODULE_CODE

! Initialize the environment status variable.
      ENV_STAT = 1

! Check to see if the MODULE_CODE environment variable is set.
      CALL GET_ENVIRONMENT_VARIABLE(NAME="MODULE_CODE", &
         VALUE=MODULE_CODE, LENGTH=ENV_LEN, STATUS=ENV_STAT)

      IF(ENV_STAT == 0 .AND. ENV_LEN /=0) THEN
         READ(WORKING_CASE,*) MODULE_CODE
      ELSE
         WORKING_CASE = LOWER_CASE
      ENDIF

! Verify the working case.
      SELECT CASE(WORKING_CASE)
      CASE(LOWER_CASE)
      CASE(UPPER_CASE)
      CASE DEFAULT
         WRITE(*,*) 'Fatal Error: Unknown working case.'
         CALL EXIT(-1)
      END SELECT

      RETURN
      END SUBROUTINE SET_WORKING_CASE

!``````````````````````````````````````````````````````````````````````!
! Subroutine: SET_MPIF_INCLUDE                                         !
!                                                                      !
! Purpose: Manages flushing the output buffer (OUT) so that if output  !
! needs multiple redirects, it only needs one call to this routine.    !
!......................................................................!
      SUBROUTINE SET_MPIF_INCLUDE

      INTEGER :: ENV_STAT, ENV_LEN
      CHARACTER(LEN=256) :: INC_PATH

! Initialize the environment status variable.
      ENV_STAT = 1

! Check to see if the MODULE_CODE environment variable is set.
      CALL GET_ENVIRONMENT_VARIABLE(NAME="MPI_INCLUDE_PATH", &
         VALUE=INC_PATH, LENGTH=ENV_LEN, STATUS=ENV_STAT)

      IF(ENV_STAT == 0 .AND. ENV_LEN /=0) THEN
         INCLUDE_MPIF = .TRUE.
      ELSE
         INCLUDE_MPIF = .FALSE.
      ENDIF

      RETURN
      END SUBROUTINE SET_MPIF_INCLUDE




!``````````````````````````````````````````````````````````````````````!
! Subroutine: FLUSH_OUT_BUFFER                                         !
!                                                                      !
! Purpose: Manages flushing the output buffer (OUT) so that if output  !
! needs multiple redirects, it only needs one call to this routine.    !
!......................................................................!
      SUBROUTINE FLUSH_OUT_BUFFER

! Line Counter
      INTEGER :: LC
! Single line.
      CHARACTER(LEN=LINE_LENGTH) :: LINE
! Line length with trailing space removed.
      INTEGER :: LENGTH
! Index of last line in the message.
      INTEGER :: LAST_LINE, tsize


! Find the end of the message.
      LAST_LINE = 0
      DO LC = 1, LINE_COUNT
         LINE = OUT(LC)
         LENGTH = len_trim(LINE)
         IF(0 < LENGTH .AND. LENGTH < LINE_LENGTH) LAST_LINE = LC
      ENDDO

      DO LC = 1, LAST_LINE
         LINE = OUT(LC)
         LENGTH = len_trim(LINE)
         IF(0 < LENGTH .AND. LENGTH < LINE_LENGTH) THEN
            WRITE(UNIT_MAKE,"(A)") trim(LINE)
         ELSE
            WRITE(UNIT_MAKE,"('  ')")
         ENDIF
      ENDDO

      IF(LAST_LINE == 0) THEN
         WRITE(UNIT_MAKE,"('  ')")
      ENDIF

      OUT = ''

      END SUBROUTINE FLUSH_OUT_BUFFER


!----------------------------------------------------------------------!
!                                                                      !
! Subroutine: FORCE_CASE                                               !
! Author:                                                              !
!                                                                      !
! Purpose: Force the case of a string to upper or lower.               !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE FORCE_CASE(STRING)

! The string modify
      CHARACTER(LEN=*), INTENT(INOUT) :: STRING
! Local variables for the upper/lower bound and differnece
      INTEGER :: LB, UB, DELTA
! Character-to-Integer conversion
      INTEGER :: CHAR_I
! Loop counter.
      INTEGER :: LC

      SELECT CASE(WORKING_CASE)
      CASE(UPPER_CASE)
         LB = ichar('a')
         UB = ichar('z')
         DELTA = ichar('A') - LB

      CASE(LOWER_CASE)
         LB = ichar('A')
         UB = ichar('Z')
         DELTA = ichar('a') - LB
      END SELECT

      DO LC=1, len_trim(STRING)
         CHAR_I = ichar(STRING(LC:LC))
         IF(LB .LE. CHAR_I .AND. CHAR_I .LE. UB) &
            STRING(LC:LC) = char(CHAR_I + DELTA)
      ENDDO

      RETURN
      END SUBROUTINE FORCE_CASE


!----------------------------------------------------------------------!
!                                                                      !
! Subroutine: TRIM_EXT                                                 !
! Author:                                                              !
!                                                                      !
! Purpose: Remove the file extension.                                  !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE TRIM_EXT(STRING, EXT)

      CHARACTER(LEN=*), INTENT(INOUT) :: STRING
      CHARACTER(LEN=*), INTENT(IN) :: EXT

      INTEGER :: LC

      DO LC = index(STRING, EXT, BACK=.TRUE.), len_trim(STRING)
         STRING(LC:LC) = ' '
      ENDDO


      RETURN
      END SUBROUTINE TRIM_EXT

!----------------------------------------------------------------------!
!                                                                      !
! Subroutine: TRIM_EXT                                                 !
! Author:                                                              !
!                                                                      !
! Purpose: Remove the file extension.                                  !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE TRIM_PATH(STRING)

      CHARACTER(LEN=*), INTENT(INOUT) :: STRING

      INTEGER :: LC, FULL, SPLIT

      FULL = len_trim(STRING)
      SPLIT = index(STRING,'/',BACK=.TRUE.)

      DO LC = 1, FULL - SPLIT
         STRING(LC:LC) = STRING(LC+SPLIT:LC+SPLIT)
      ENDDO

      DO LC = FULL-SPLIT+1, FULL
         STRING(LC:LC) = ' '
      ENDDO

      RETURN
      END SUBROUTINE TRIM_PATH


!----------------------------------------------------------------------!
!                                                                      !
! Subroutine: GET_INC_FILES                                            !
! Author:                                                              !
!                                                                      !
! Purpose: Open up FNAME and search for include statements.            !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE GET_INC_FILES(FNAME, LIST, COUNT)

      CHARACTER(LEN=*), INTENT(IN) :: FNAME

      CHARACTER(LEN=*), INTENT(OUT) :: LIST(*)
      INTEGER, INTENT(OUT) :: COUNT

      CHARACTER(LEN=LINE_LENGTH) :: LINE
      CHARACTER(LEN=LINE_LENGTH) :: PATH

      INTEGER :: SPLIT

      SPLIT = index(FNAME,'/',BACK=.TRUE.)
      PATH = FNAME(1:SPLIT)

      COUNT = 0
!      OPEN(UNIT=UNIT_TEMP, FILE=trim(FNAME), STATUS='OLD')

!      DO; READ (UNIT_TEMP,'(A)',END=999) LINE
!         CALL FIND_INC(LINE, LIST, COUNT, PATH)
!      ENDDO

! 999  CLOSE(UNIT_TEMP)
      RETURN
      END SUBROUTINE GET_INC_FILES


!----------------------------------------------------------------------!
!                                                                      !
! Subroutine: FIND_INC                                                 !
! Author:                                                              !
!                                                                      !
! Purpose: Open up FNAME and search for include statements.            !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE FIND_INC(LINE, LIST, COUNT, PATH)

      CHARACTER(LEN=*), INTENT(IN) :: LINE

      INTEGER, INTENT(INOUT) :: COUNT
      CHARACTER(LEN=*), INTENT(INOUT) :: LIST(*)

      CHARACTER(LEN=*), INTENT(IN) :: PATH

      INTEGER :: UC, LC, LEND
      INTEGER :: ldP, lsP
      INTEGER :: lPOS, rPOS, POS

      INTEGER :: COUNT1

      LEND = SEEK_COMMENT(LINE)

      UC = INDEX(LINE(:LEND),'INCLUDE')
      LC = INDEX(LINE(:LEND),'include')

      IF(LC == 0 .AND. UC == 0) RETURN

      IF(LC > 0 .AND. UC == 0) THEN
         POS = LC
      ELSEIF(LC == 0 .AND. UC > 0) THEN
         POS = UC
      ELSE
         write(*,*) 'Error 1 in FIND_INC: ',trim(LINE)
      ENDIF

! Skip lines where 'include' is not the first word.
      IF(.NOT.BLANK_BEFORE(LINE,POS)) RETURN

! Search for quote marks bounding the chemcial equation.
      ldP = POS + INDEX(LINE(POS:LEND),'"')  ! double quote "
      lsP = POS + INDEX(LINE(POS:LEND),"'")  ! single quote '

      IF(ldP .GT. POS .AND. lsP .EQ. POS) THEN
! The include is bounded by double quotes
         lPOS = ldP
! Search for the second quote mark.
         rPOS = lPOS + INDEX(LINE(lPOS+1:LEND),'"')
      ELSEIF(ldP .EQ. POS .AND. lsP .GT. POS) THEN
! The chemical equation is bounded by single quotes
         lPOS = lsP
! Search for the second quote mark.
         rPOS = lPOS + INDEX(LINE(lPOS+1:LEND),"'")

! Flag error and stop.
      ELSE
         write(*,*)'Error 2 in FIND_INC: '
         write(*,*)'LINE: ',trim(LINE)
         CALL EXIT(-1)

      ENDIF

      COUNT1 = COUNT + 1
      LIST(COUNT1) = ''

! Clean up relative references.
      DO LC=lPOS,rPOS-1
         IF(LINE(LC:LC) == '.') CYCLE
         IF(LINE(LC:LC) == '/') CYCLE
         EXIT
      ENDDO

! Clear the entry and cache the informtion.
      IF(LC == lPOS .AND. len_trim(PATH) /= 0) THEN
         LIST(COUNT1) = trim(PATH)//LINE(LC:rPOS-1)
      ELSE
         LIST(COUNT1) = LINE(LC:rPOS-1)
      ENDIF


      IF(.NOT.INCLUDE_MPIF) THEN
         IF(index(LIST(COUNT1),'mpif.h') > 0) RETURN
      ENDIF


! Only keep the entry if it's new.
      DO LC=1, COUNT
         IF(LIST(LC) == LIST(COUNT1)) RETURN
      ENDDO
      COUNT = COUNT1


      RETURN
      END SUBROUTINE FIND_INC 

!----------------------------------------------------------------------!
!                                                                      !
! Subroutine: GET_USE_FILES                                            !
! Author:                                                              !
!                                                                      !
! Purpose: Open up FNAME and search for use statements.                !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE GET_USE_FILES(FNAME, LIST, COUNT)

      CHARACTER(LEN=*), INTENT(IN) :: FNAME

      CHARACTER(LEN=*), INTENT(OUT) :: LIST(*)
      INTEGER, INTENT(OUT) :: COUNT

      CHARACTER(LEN=LINE_LENGTH) :: LINE

      COUNT = 0
      OPEN(UNIT=UNIT_TEMP, FILE=trim(FNAME), STATUS='OLD')

      DO; READ (UNIT_TEMP,'(A)',END=999) LINE
         CALL FIND_USE(LINE, LIST, COUNT)
      ENDDO

 999  CLOSE(UNIT_TEMP)
      RETURN
      END SUBROUTINE GET_USE_FILES

!----------------------------------------------------------------------!
!                                                                      !
! Subroutine: FIND_USE                                                 !
! Author:                                                              !
!                                                                      !
! Purpose: Open up FNAME and search for include statements.            !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE FIND_USE(LINE, LIST, COUNT)

      CHARACTER(LEN=*), INTENT(IN) :: LINE

      INTEGER, INTENT(INOUT) :: COUNT
      CHARACTER(LEN=*), INTENT(INOUT) :: LIST(*)

      INTEGER :: UC, MC, LC, LEND
      INTEGER :: ldP, lsP
      INTEGER :: lPOS, rPOS, POS
      INTEGER :: COUNT1

      CHARACTER(LEN=LINE_LENGTH) :: TEST

      LEND = SEEK_COMMENT(LINE)

      UC = INDEX(LINE(:LEND),'USE ')
      MC = INDEX(LINE(:LEND),'Use ')
      LC = INDEX(LINE(:LEND),'use ')


      IF(LC == 0 .AND. MC == 0 .AND. UC == 0) RETURN

      IF(LC > 0 .AND. MC == 0 .AND. UC == 0) THEN
         POS = LC
      ELSEIF(LC == 0 .AND. MC > 0 .AND. UC == 0) THEN
         POS = MC
      ELSEIF(LC == 0 .AND. MC == 0 .AND. UC > 0) THEN
         POS = UC
      ELSE
         WRITE(*,*)'Error 1 in FIND_USE: ',trim(LINE)
         CALL EXIT(-1)
      ENDIF

! If there are non-blanks this this is not a use statement
      IF(.NOT.BLANK_BEFORE(LINE,POS)) RETURN

      lPOS = POS+4
      DO LC=lPOS, len(LINE)
         IF(LINE(LC:LC) == ' ' .OR. LINE(LC:LC) == ',') THEN
            rPOS = LC
            EXIT
         ENDIF
      ENDDO

      COUNT1 = COUNT+1

! Clear the entry and cache the informtion.
      LIST(COUNT1) = ''
      LIST(COUNT1) = LINE(lPOS:rPOS-1)

! Only keep the entry if it's new.
      DO LC=1, COUNT
         IF(LIST(LC) == LIST(COUNT1)) RETURN
      ENDDO
      COUNT = COUNT1

      RETURN
      END SUBROUTINE FIND_USE

!----------------------------------------------------------------------!
!                                                                      !
! Function: SEEK_COMMENT                                               !
! Author:                                                              !
!                                                                      !
! Purpose: Return the position of a comment.                           !
!                                                                      !
!----------------------------------------------------------------------!
      INTEGER FUNCTION SEEK_COMMENT(STRING)

      CHARACTER(LEN=*), INTENT(IN) :: STRING
      INTEGER :: LC

      SEEK_COMMENT = 0
      DO LC=1, len_trim(STRING)
         IF(STRING(LC:LC) == '!') RETURN
         SEEK_COMMENT = LC
      ENDDO

      END FUNCTION SEEK_COMMENT


!----------------------------------------------------------------------!
!                                                                      !
! Function: BLANK_BEFORE                                               !
! Author:                                                              !
!                                                                      !
! Purpose: Verify that the input before POS in STRING is empty.        !
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION BLANK_BEFORE(STRING, POS)

      CHARACTER(LEN=*), INTENT(IN) :: STRING
      INTEGER, INTENT(IN) :: POS

      INTEGER :: LC

      BLANK_BEFORE = .TRUE.
      IF(POS <= 1) RETURN

      DO LC=1, POS-1
         IF(STRING(LC:LC) == ' ') CYCLE  ! Space
         IF(STRING(LC:LC) == '	') CYCLE ! Tabs
         BLANK_BEFORE = .FALSE.
      ENDDO
      RETURN
      END FUNCTION BLANK_BEFORE

      END MODULE WRITE_UTILITY



!----------------------------------------------------------------------!
!                                                                      !
! Subroutine: MAKEFILE_MFIX                                            !
! Author: P.Nicoletti and M.Syamlal                                    !
!                                                                      !
! Purpose: Create a make file for MFIX.                                !
!                                                                      !
! Assumptions:                                                         !
!   - Module ABC is located in file abc_mod.f                          !
!   - No non-module file contains the string _mod. in the file name.   !
!                                                                      !
!                                                                      !
!                                                                      !
!----------------------------------------------------------------------!
      PROGRAM MAKEFILE_MFIX

      use WRITE_UTILITY

      IMPLICIT NONE

! Max length of file names, modules and includes.
      INTEGER, PARAMETER :: MAX_LEN = 64
! Max number of files.
      INTEGER, PARAMETER :: MAX_LIST = 4000
! List of all files under the mfix/model and all sub dirs
      CHARACTER(LEN=MAX_LEN) :: ALL_FILES(MAX_LIST)
! List of only the module files under mfix/model and all sub dirs.
      CHARACTER(LEN=MAX_LEN) :: MOD_FILES(MAX_LIST)

! List of inlcude and use statements in a file.
      CHARACTER(LEN=MAX_LEN) :: INC_FILES(MAX_LIST)
      CHARACTER(LEN=MAX_LEN) :: USE_FILES(MAX_LIST)

! Counts of various elements.
      INTEGER :: FILE_COUNT  ! number of files under mfix/model
      INTEGER :: MOD_COUNT   ! number of modules under mfix/model
      INTEGER :: INC_COUNT   ! number of include statements (per file)
      INTEGER :: USE_COUNT   ! number of use statements (per file)

! Working verisons of the file names.
      CHARACTER(LEN=MAX_LEN) :: FILENAME ! Full path file name
      CHARACTER(LEN=MAX_LEN) :: BASENAME ! Base name of file

! Loop counters.
      INTEGER :: LC1, LC2
      CHARACTER(LEN=1) :: C1

! Sets the case for module files based on the environment variable. It
! defaults to lower case if the flag is not set.
      CALL SET_WORKING_CASE

! If the MPI include environment variable has been set, then mpif.h
! is included in the Makefile otherwise it is omitted. This removes
! the need for a 'dummy' mpif.h for serial runs as well as issues
! with Cray compilers where includes are not needed.
      CALL SET_MPIF_INCLUDE

! Initialize the file and module count.
      FILE_COUNT = 0
      MOD_COUNT  = 0

! Open the temp make file.
      OPEN(UNIT=UNIT_MAKE,FILE='tmp.make', &
         STATUS='OLD', POSITION='APPEND')

! Open and process the file list information.
      OPEN(UNIT=UNIT_TEMP,FILE='files_post.lis',STATUS='OLD')

      DO; READ(UNIT_TEMP,'(A)',END=100) FILENAME

! Store the file name.
         FILE_COUNT = FILE_COUNT + 1
         ALL_FILES(FILE_COUNT) = FILENAME

! Check for the mod extension.
         IF(index(FILENAME,'_mod.') /= 0) THEN
            MOD_COUNT = MOD_COUNT + 1
            MOD_FILES(MOD_COUNT) = FILENAME
         ENDIF
      ENDDO
 100  CLOSE(UNIT=UNIT_TEMP)


! Write out the start of the Makefile.
      OUT = ''
      WRITE(OUT,5000) BS
      CALL FLUSH_OUT_BUFFER

 5000 FORMAT(/'.$(FORTRAN_EXT).$(OBJ_EXT):',/'	$(FORTRAN_CMD) ',     &
         '$(FORT_FLAGS) $<',2/,'$(EXEC_FILE) : ',A1)

! First write out all the module files. Remove the file extension, path
! information and set the case.
      DO LC1 = 1, MOD_COUNT
         FILENAME = MOD_FILES(LC1)
         CALL TRIM_EXT(FILENAME, '_mod.')
         CALL TRIM_PATH(FILENAME)
         CALL FORCE_CASE(FILENAME)
         WRITE(OUT,"(4X,'$(DPO)',A,'.mod ',A1)") trim(FILENAME), BS
         CALL FLUSH_OUT_BUFFER
      ENDDO

! Now include the object files for the regular source files.
      DO LC1 = 1,FILE_COUNT
         FILENAME = ALL_FILES(LC1)
         IF(index(FILENAME,'_mod.') /= 0) CYCLE
         CALL TRIM_EXT(FILENAME, '.')
         CALL TRIM_PATH(FILENAME)
         WRITE(OUT,"(4X,'$(DPO)',A,'.$(OBJ_EXT) ',A1)") &
            trim(FILENAME), BS
         CALL FLUSH_OUT_BUFFER
      ENDDO


! Set the linking commands.
      WRITE(OUT,"(/'	$(LINK_CMD) $(LINK_FLAGS) ',A1)") BS
      CALL FLUSH_OUT_BUFFER

      DO LC1 = 1,FILE_COUNT
         FILENAME = ALL_FILES(LC1)
         CALL TRIM_EXT(FILENAME, '.')
         CALL TRIM_PATH(FILENAME)
         WRITE(OUT,"(4X,'$(DPO)',A,'.$(OBJ_EXT) ',A1)") &
            trim(FILENAME), BS
         CALL FLUSH_OUT_BUFFER
      ENDDO
      WRITE(OUT,"(' -o $(EXEC_FILE) $(LIB_FLAGS)')")
      CALL FLUSH_OUT_BUFFER

      DO LC1=1, MOD_COUNT

         FILENAME = MOD_FILES(LC1)

         CALL GET_INC_FILES(FILENAME, INC_FILES, INC_COUNT)
         CALL GET_USE_FILES(FILENAME, USE_FILES, USE_COUNT)

         BASENAME = FILENAME
         CALL TRIM_EXT(BASENAME, '_mod')
         CALL TRIM_PATH(BASENAME)
         CALL FORCE_CASE(BASENAME)

         IF((INC_COUNT + USE_COUNT) == 0) THEN
            WRITE(OUT,"(/'$(DPO)',A,'.mod : ',A)")                     &
               trim(BASENAME), trim(FILENAME)
            CALL FLUSH_OUT_BUFFER
         ELSE
            WRITE(OUT,"(/'$(DPO)',A,'.mod : ',A,1X,A1)")               &
               trim(BASENAME), trim(FILENAME), BS
            CALL FLUSH_OUT_BUFFER

            DO LC2=1, USE_COUNT
               C1 = merge(BS,' ',INC_COUNT>0 .OR. LC2/=USE_COUNT)
               CALL FORCE_CASE(USE_FILES(LC2))
               WRITE(OUT,"(11x,'$(DPO)',A,'.mod ',A)") &
                  trim(USE_FILES(LC2)), C1
               CALL FLUSH_OUT_BUFFER
            ENDDO

            DO LC2=1, INC_COUNT
               C1 = merge(BS,' ',LC2 /= INC_COUNT)
               WRITE(OUT,"(11x,A,1x,A)") trim(INC_FILES(LC2)), C1
               CALL FLUSH_OUT_BUFFER
            ENDDO
         ENDIF

         BASENAME = FILENAME
         CALL TRIM_EXT(BASENAME, '_mod')
         CALL TRIM_PATH(BASENAME)

         WRITE(OUT, 5060) trim(FILENAME), trim(BASENAME)
         CALL FLUSH_OUT_BUFFER
     
 5060 FORMAT('	$(FORTRAN_CMD) $(FORT_FLAGS) ',A,' -o $(DPO)',A,&
         '_mod.$(OBJ_EXT) $(MODDIRPREFIX)$(DPO)')

      ENDDO


      DO LC1=1, FILE_COUNT

         FILENAME = ALL_FILES(LC1)
         IF(index(FILENAME,'mod.') /= 0) CYCLE

         CALL GET_INC_FILES(FILENAME, INC_FILES, INC_COUNT)
         CALL GET_USE_FILES(FILENAME, USE_FILES, USE_COUNT)

         BASENAME = FILENAME
         CALL TRIM_EXT(BASENAME, '.')
         CALL TRIM_PATH(BASENAME)

         IF((INC_COUNT + USE_COUNT) == 0) THEN
            WRITE(OUT,"(/'$(DPO)',A,'.$(OBJ_EXT) : ',A)")              &
               trim(BASENAME), trim(FILENAME)
            CALL FLUSH_OUT_BUFFER
         ELSE
            WRITE(OUT,"(/'$(DPO)',A,'.$(OBJ_EXT) : ',A,' 'A1)")        &
               trim(BASENAME), trim(FILENAME), BS
            CALL FLUSH_OUT_BUFFER

            DO LC2=1, USE_COUNT
               C1 = merge(BS,' ',INC_COUNT>0 .OR. LC2/=USE_COUNT)
               CALL FORCE_CASE(USE_FILES(LC2))
               WRITE(OUT,"(11x,'$(DPO)',A,'.mod ',A)") &
                  trim(USE_FILES(LC2)), C1
               CALL FLUSH_OUT_BUFFER
            ENDDO

            DO LC2=1, INC_COUNT
               C1 = merge(BS,' ',LC2 /= INC_COUNT)
               WRITE(OUT,"(11x,A,1x,A)") trim(INC_FILES(LC2)), C1
               CALL FLUSH_OUT_BUFFER
            ENDDO
         ENDIF

         WRITE(OUT, 5061) trim(FILENAME), trim(BASENAME)
         CALL FLUSH_OUT_BUFFER
     
 5061 FORMAT('	$(FORTRAN_CMD) $(FORT_FLAGS) ',A,' -o $(DPO)',A,&
         '.$(OBJ_EXT) $(MODDIRPREFIX)$(DPO)')

      ENDDO

      CLOSE(UNIT_MAKE)

      CALL EXIT(0)
      END PROGRAM MAKEFILE_MFIX

