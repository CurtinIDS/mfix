      MODULE vtp

      use mpi_utility
      use cdist

      use desmpi
      use mpi_comm_des
      use error_manager

      IMPLICIT NONE

      INTEGER, PRIVATE :: GLOBAL_CNT
      INTEGER, PRIVATE :: LOCAL_CNT

      INTEGER :: DES_UNIT = 2000

! file unit for ParaView *.pvd data
      INTEGER, PARAMETER :: PVD_UNIT = 2050

! formatted file name
      CHARACTER(LEN=511) :: FNAME_VTP

      INTERFACE VTP_WRITE_DATA
         MODULE PROCEDURE VTP_WRITE_DP1
         MODULE PROCEDURE VTP_WRITE_DP2
         MODULE PROCEDURE VTP_WRITE_I1
      END INTERFACE

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_DP1                                            !
!                                                                      !
! Purpose: Collect and write 1D double percision arrays to the VTP     !
! file. This routine is designed to collect the data for parallel and  !
! serial runs. This routine also manages the distribted IO case.       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_DP1(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      DOUBLE PRECISION, INTENT(in) :: DATA(:)

      INTEGER :: LC, PC

      IF(bDist_IO) THEN

         WRITE(DES_UNIT,1000) NAME

         PC = 1
         DO LC = 1, MAX_PIP
            IF(PC > PIP) EXIT
            IF(IS_NONEXISTENT(LC)) CYCLE
            PC = PC+1
            IF(IS_GHOST(LC) .OR. IS_ENTERING_GHOST(LC) .OR. IS_EXITING_GHOST(LC)) CYCLE
            WRITE(DES_UNIT, 1001,ADVANCE="NO") real(DATA(LC))
         ENDDO
         WRITE(DES_UNIT,1002)

      ELSE

         allocate (dProcBuf(LOCAL_CNT) )
         allocate (dRootBuf(GLOBAL_CNT))

         CALL DES_GATHER(DATA)

         IF(myPE == PE_IO) THEN
            WRITE(DES_UNIT,1000) NAME
            DO LC=1, GLOBAL_CNT
               WRITE(DES_UNIT,1001,ADVANCE="NO") real(drootbuf(LC))
            ENDDO
            WRITE(DES_UNIT,1002)
         ENDIF

         deallocate(dProcBuf, dRootBuf)

      ENDIF

 1000 FORMAT('<DataArray type="Float32" Name="',A,'" format="ascii">')
 1001 FORMAT(ES14.6,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_DP1

!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_DP2                                            !
!                                                                      !
! Purpose: Collect and write 2D double percision arrays to the VTP     !
! file. This routine is designed to collect the data for parallel and  !
! serial runs. This routine also manages the distribted IO case.       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_DP2(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      DOUBLE PRECISION, INTENT(in) :: DATA(:,:)

      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)

      CHARACTER(len=16) :: NOC
      INTEGER :: LB, UB
      INTEGER :: PC, LC1, LC2

      LB = LBOUND(DATA,2)
      UB = UBOUND(DATA,2)
      NOC=''; WRITE(NOC,*) (UB-LB)+1

      IF(bDist_IO) THEN

         WRITE(DES_UNIT,1000) NAME, trim(adjustl(NOC))

         PC = 1
         DO LC1 = 1, MAX_PIP
            IF(PC > PIP) EXIT
            IF(IS_NONEXISTENT(LC1)) CYCLE
            PC = PC+1
            IF(IS_GHOST(LC1) .OR. IS_ENTERING_GHOST(LC1) .OR. IS_EXITING_GHOST(LC1)) CYCLE
            DO LC2=LB, UB
               WRITE(DES_UNIT,1001,ADVANCE="NO") real(DATA(LC1,LC2))
            ENDDO
         ENDDO
         WRITE(DES_UNIT,1002)

      ELSE

         allocate (dProcBuf(LOCAL_CNT) )
         allocate (dRootBuf(GLOBAL_CNT))
         allocate (ltemp_array((UB-LB)+1,GLOBAL_CNT))

         DO LC1 = LB, UB
            CALL DES_GATHER(DATA(:,LC1))
            ltemp_array(LC1,:) = drootbuf(:)
         ENDDO

         IF(myPE == PE_IO) THEN
            WRITE(DES_UNIT,1000) NAME, trim(adjustl(NOC))
            DO LC1=1, GLOBAL_CNT
               DO LC2=LB, UB
                  WRITE(DES_UNIT,1001,ADVANCE="NO") &
                     real(ltemp_array(LC2,LC1))
               ENDDO
            ENDDO
            WRITE(DES_UNIT,1002)
         ENDIF

         deallocate (dProcBuf, dRootBuf, ltemp_array)

      ENDIF


 1000 FORMAT('<DataArray type="Float32" Name="',A,'" NumberOf',        &
         'Components="',A,'" format="ascii">')
 1001 FORMAT(ES14.6,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_DP2



!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_I1                                             !
!                                                                      !
! Purpose: Collect and write 1D integer arrays to the VTP file. This   !
! routine is designed to collect the data for parallel and serial      !
! runs. This routine also manages the distribted IO case.              !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_I1(NAME, DATA)

      CHARACTER(len=*), INTENT(in) :: NAME
      INTEGER, INTENT(in) :: DATA(:)

      INTEGER :: LC, PC

      IF(bDist_IO) THEN

         WRITE(DES_UNIT,1000) NAME

         PC = 1
         DO LC = 1, MAX_PIP
            IF(PC > PIP) EXIT
            IF(IS_NONEXISTENT(LC)) CYCLE
            PC = PC+1
            IF(IS_GHOST(LC) .OR. IS_ENTERING_GHOST(LC) .OR. IS_EXITING_GHOST(LC)) CYCLE
            WRITE(DES_UNIT, 1001,ADVANCE="NO") DATA(LC)
         ENDDO
         WRITE(DES_UNIT,1002)

      ELSE

         allocate (iProcBuf(LOCAL_CNT) )
         allocate (iRootBuf(GLOBAL_CNT))

         CALL DES_GATHER(DATA)

         IF(myPE == PE_IO) THEN
            WRITE(DES_UNIT,1000) NAME
            DO LC=1, GLOBAL_CNT
               WRITE(DES_UNIT,1001,ADVANCE="NO") irootbuf(LC)
            ENDDO
            WRITE(DES_UNIT,1002)
         ENDIF

         deallocate(iProcBuf, iRootBuf)

      ENDIF

 1000 FORMAT('<DataArray type="Float32" Name="',A,'" format="ascii">')
 1001 FORMAT(I10,1X)
 1002 FORMAT('</DataArray>')

      END SUBROUTINE VTP_WRITE_I1


!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_WRITE_ELEMENT                                        !
!                                                                      !
! Purpose: Write a string to the VTP file. It masks the need to check  !
! the logical before flushing.                                         !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_WRITE_ELEMENT(ELEMENT)

      CHARACTER(len=*), INTENT(in) :: ELEMENT

      IF(bDist_IO .OR. myPE == PE_IO) &
         WRITE(DES_UNIT,"(A)") ELEMENT

      RETURN
      END SUBROUTINE VTP_WRITE_ELEMENT



!``````````````````````````````````````````````````````````````````````!
! Subroutine: VTP_OPEN_FILE                                            !
!                                                                      !
! Purpose: This routine opens the VTP file and calcualtes the offsets  !
! for dmp data collection.                                             !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_OPEN_FILE(NoPc)

! Modules
!-----------------------------------------------
      use discretelement, only: VTP_DIR
      use run, only: RUN_TYPE, RUN_NAME

      IMPLICIT NONE

      CHARACTER(len=*) :: NoPc

      INTEGER :: NumberOfPoints

! Variables related to gather
      integer lgathercnts(0:numpes-1), lproc

! check whether an error occurs in opening a file
      INTEGER :: IOS
! Integer error flag.
      INTEGER :: IER

! logical used for testing is the data file already exists
      LOGICAL :: EXISTS_VTP
! status of the vtp file to be written
      CHARACTER(LEN=8) :: STATUS_VTP

      IF(TRIM(VTP_DIR)/='.') CALL CREATE_DIR(trim(VTP_DIR))

! Initial the global count.
      GLOBAL_CNT = 10
! Calculate the number of 'real' particles on the local process.
      LOCAL_CNT = PIP - iGHOST_CNT

! Distributed IO
      IF(bDIST_IO) THEN
         NumberOfPoints = LOCAL_CNT
         WRITE(NoPc,"(I10.10)") NumberOfPoints

         IF(TRIM(VTP_DIR)/='.') THEN
            WRITE(fname_vtp,'(A,"/",A,"_DES",I4.4,"_",I5.5,".vtp")') &
               trim(VTP_DIR), trim(run_name), vtp_findex, mype
         ELSE
            WRITE(fname_vtp,'(A,"_DES",I4.4,"_",I5.5,".vtp")') &
               trim(run_name), vtp_findex, mype
         ENDIF

! Serial IO
      ELSE

! Calculate the total number of particles system-wide.
         call global_sum(LOCAL_CNT, GLOBAL_CNT)
         NumberOfPoints = GLOBAL_CNT
         WRITE(NoPc,"(I10.10)") NumberOfPoints

! Set the send count from the local process.
         igath_sendcnt = LOCAL_CNT

! Collect the number of particles on each rank.all ranks.
         lgathercnts = 0
         lgathercnts(myPE) = LOCAL_CNT
         call global_sum(lgathercnts,igathercnts)

! Calculate the rank displacements.
         idispls(0) = 0
         DO lPROC = 1,NUMPEs-1
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
         ENDDO

! set the file name and unit number and open file
         IF(TRIM(VTP_DIR)/='.') THEN
            WRITE(fname_vtp,'(A,"/",A,"_DES_",I5.5,".vtp")') &
               trim(VTP_DIR),trim(run_name), vtp_findex
         ELSE
            WRITE(fname_vtp,'(A,"_DES_",I5.5,".vtp")') &
               trim(run_name), vtp_findex
         ENDIF
      ENDIF

      IER = 0
      IF(bDIST_IO .OR. myPE == PE_IO) THEN

! The file should be new but could exist due to restarting.
         STATUS_VTP = 'NEW'
! Check to see if the file already exists.
         INQUIRE(FILE=FNAME_VTP,EXIST=EXISTS_VTP)
! The given file should not exist if the run type is NEW.
         IF(EXISTS_VTP)THEN
! The VTP should never exist for a NEW run.
            IF(RUN_TYPE == 'NEW')THEN
               IER = 1
! The file may exist during a RESTART.
            ELSE
               STATUS_VTP = 'REPLACE'
            ENDIF
         ENDIF

! Open the file and record any erros.
         IF(IER == 0) THEN
            OPEN(CONVERT='BIG_ENDIAN',UNIT=DES_UNIT, FILE=FNAME_VTP,   &
               STATUS=STATUS_VTP, IOSTAT=IOS)
            IF(IOS /= 0) IER = 2
         ENDIF
      ENDIF

      CALL GLOBAL_ALL_MAX(IER)

      IF(IER /= 0) THEN
         CALL INIT_ERR_MSG("VTP_MOD --> OPEN_VTP")
         WRITE(ERR_MSG, 1100) IER
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Unable to open VTP file. This could be ',    &
         'caused by a VTP',/'file with the same file name already ',   &
         'existing. or an error code',/' returned by the OPEN ',       &
         'function.'/'Error code: ',I2,4x,'Aborting.')


      END SUBROUTINE VTP_OPEN_FILE



!......................................................................!
! SUBROUTINE: VTP_CLOSE_FILE                                           !
!                                                                      !
! Purpose: This routine closes the vtp file.                           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE VTP_CLOSE_FILE


      VTP_FINDEX=VTP_FINDEX+1

      IF(bDist_io .OR. (myPE .eq.pe_IO)) CLOSE(des_unit)


      END SUBROUTINE VTP_CLOSE_FILE


!......................................................................!
! SUBROUTINE: ADD_VTP_TO_PVD                                           !
!                                                                      !
! Purpose: This routine opens the pvd file.                            !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE ADD_VTP_TO_PVD

      use discretelement, only: VTP_DIR
      use run, only: RUN_TYPE, RUN_NAME

!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Index position of desired character
      INTEGER IDX_f, IDX_b
! logical used for testing is the data file already exists
      LOGICAL :: EXISTS_PVD
! Generic input limited to 256 characters
      CHARACTER(LEN=256) INPUT

! formatted file name
      CHARACTER(LEN=64) :: FNAME_PVD = ''
! formatted time
      CHARACTER(LEN=64) :: cTIME = ''

      LOGICAL, SAVE :: FIRST_PASS = .TRUE.

! IO Status flag
      INTEGER :: IOS

! Variables related to gather
      integer :: IER

!-----------------------------------------------

      CALL INIT_ERR_MSG('VTP_MOD --> ADD_VTP_TO_PVD')

! Initialize the error flag.
      IER = 0

! Obtain the file name and open the pvd file
      FNAME_PVD = TRIM(RUN_NAME)//'_DES.pvd'

! The PVD file is only written by PE_IO with serial IO.
      IF(myPE == PE_IO .AND. .NOT.bDist_IO) THEN

! Check to see if the file already exists.
         INQUIRE(FILE=FNAME_PVD,EXIST=EXISTS_PVD)

         IF(FIRST_PASS) THEN

! Open the "NEW" file and write the necessary header information.
            IF(RUN_TYPE /= 'RESTART_1')THEN

! The file exists but first_pass is also true so most likely an existing
! file from an earlier/other run is present in the directory. Exit to
! prevent accidently overwriting the existing file.
               IF(EXISTS_PVD) THEN
                  IER = 1
               ELSE
                  OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,STATUS='NEW')
                  WRITE(PVD_UNIT,"(A)")'<?xml version="1.0"?>'
                  WRITE(PVD_UNIT,"(A)")'<VTKFile type="Collection" &
                     &version="0.1" byte_order="LittleEndian">'
                  WRITE(PVD_UNIT,"(3X,'<Collection>')")
               ENDIF

! This is the first pass of a restart run. Extra care is needed to make
! sure that the pvd file is ready to accept new data.
            ELSE ! a restart run
               IF(EXISTS_PVD) THEN
! Open the file at the beginning.
                  OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
                     POSITION="REWIND",STATUS='OLD',IOSTAT=IOS)
                  IF(IOS /= 0) IER = 2
               ELSE ! a pvd file does not exist
                  IER = 3
               ENDIF

               IF(IER == 0) THEN
! Loop over the entries in the PVD file, looking for a match to the
! file that is being written. If no match is found, the data will be
! appended to the end of the pvd file, otherwise, the old data will
! be over-written.
                  DO
! Read in the entires of the PVD file.
                     READ(PVD_UNIT,"(A)",IOSTAT=IOS)INPUT
                     IF(IOS > 0) THEN
                        IER = 4
                        EXIT
                     ELSEIF(IOS<0)THEN
! The end of the pvd file has been reached without finding an entry
! matching the current record. Exit the loop.
                        BACKSPACE(PVD_UNIT)
                        BACKSPACE(PVD_UNIT)
                        BACKSPACE(PVD_UNIT)
                        EXIT
                     ENDIF
! Find the first instances of file=" and "/> in the read data.
                     IDX_f = INDEX(INPUT,'file="')
                     IDX_b = INDEX(INPUT,'"/>')
! Skip rows that do not contain file data
                     IF(IDX_f == 0 .AND. IDX_b == 0) CYCLE
! Truncate the file name from the read data
                     WRITE (INPUT,"(A)") INPUT(IDX_f+6:IDX_b-1)
! If the file name matches the current VTP record, break the loop to
! over-write this record.
                     IF(TRIM(FNAME_VTP) == TRIM(INPUT)) THEN
                        BACKSPACE(PVD_UNIT)
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF ! No errors
            ENDIF ! run_type new or restart

         ELSE ! not FIRST_PASS
            OPEN(UNIT=PVD_UNIT,FILE=FNAME_PVD,&
               POSITION="APPEND",STATUS='OLD',IOSTAT=IOS)
            IF (IOS /= 0) IER = 2
         ENDIF

      ENDIF ! if myPE == PE_IO and not distributed IO


      CAlL GLOBAL_ALL_SUM(IER)
      IF(IER /= 0) THEN
         SELECT CASE(IER)
         CASE(1); WRITE(ERR_MSG,1101) trim(FNAME_PVD)
         CASE(2); WRITE(ERR_MSG,1102) trim(FNAME_PVD)
         CASE(3); WRITE(ERR_MSG,1103) trim(FNAME_PVD)
         CASE(4); WRITE(ERR_MSG,1104) trim(FNAME_PVD)
         CASE DEFAULT; WRITE(ERR_MSG,1105) trim(FNAME_PVD)
         END SELECT
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1101 FORMAT('Error 1101: A PVD file was detected in the run ',        &
         'directory which should',/'not exist for a NEW run.',/        &
         'File: ',A)

 1102 FORMAT('Error 1102: Fatal error status returned while OPENING ', &
         'PVD file.',/'File: ', A)

 1103 FORMAT('Error 1103: PVD file MISSING from run directory.',/      &
         'File: ',A)

 1104 FORMAT('Error 1104: Fatal error status returned while READING ', &
         'PVD file.',/'File: ', A)

 1105 FORMAT('Error 1105:: Fatal unclassified error when processing ', &
         'PVD file.',/'File: ', A)


! If there were no errors, updated the file.
      IF(myPE == PE_IO .AND. .NOT.bDist_IO) THEN

! Remove the last two lines written so that additional data can be added
         IF(.NOT.FIRST_PASS) THEN
            BACKSPACE(PVD_UNIT)
            BACKSPACE(PVD_UNIT)
         ENDIF

         WRITE(cTIME,"(F12.6)") S_TIME
! Write the data to the file
         WRITE(PVD_UNIT,"(6X,A,A,A,A,A,A,A)")&
         '<DataSet timestep="',trim(adjustl(cTIME)),'" ',&
         'group="" part="0" ',& ! necessary file data
         'file="',TRIM(FNAME_VTP),'"/>' ! file name of vtp

! Write the closing tags
         WRITE(PVD_UNIT,"(3X,A)")'</Collection>'
         WRITE(PVD_UNIT,"(A)")'</VTKFile>'

         CLOSE(PVD_UNIT)
      ENDIF
! Identify that the files has been created and opened for next pass
      FIRST_PASS = .FALSE.

      CALL FINL_ERR_MSG

! Return to the calling routine
      RETURN

      END SUBROUTINE ADD_VTP_TO_PVD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VTP_FILE                                         C
!  Purpose: Writes particles data in VTK format (Polydata VTP)         C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_VTP_FILE(LCV,MODE)

      USE vtk, only: DIMENSION_VTK, VTK_DEFINED, FRAME
      USE vtk, only: VTK_REGION,VTK_DEFINED,VTK_DATA
      USE vtk, only: VTK_PART_DIAMETER,VTK_PART_VEL,VTK_PART_USR_VAR,VTK_PART_TEMP
      USE vtk, only: VTK_PART_ANGULAR_VEL,VTK_PART_ORIENTATION
      USE vtk, only: VTK_PART_X_S, VTK_PART_COHESION
      USE vtk, only: TIME_DEPENDENT_FILENAME,VTU_FRAME_UNIT,VTU_FRAME_FILENAME
      USE vtk, only: VTK_DBG_FILE
      USE output, only: FULL_LOG
      use des_thermo, only: DES_T_s
      IMPLICIT NONE
      INTEGER :: L,N,LCV

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2
      INTEGER :: MODE   ! MODE = 0 : Write regular VTK region file
                        ! MODE = 1 : Write debug   VTK region file (VTK_DBG_FILE = .TRUE.)

      VTK_REGION = LCV
! There is nothing to write if we are not in a defined vtk region
      IF(.NOT.VTK_DEFINED(VTK_REGION)) RETURN

      IF(VTK_DATA(LCV)/='P') RETURN
      IF(MODE==0.AND.(VTK_DBG_FILE(LCV))) RETURN
      IF(MODE==1.AND.(.NOT.VTK_DBG_FILE(LCV))) RETURN

      CALL SETUP_VTK_REGION_PARTICLES

      CALL OPEN_VTP_FILE_BIN(MODE)

! Only open pvd file when there are particles in vtk region
      IF(GLOBAL_CNT>0.AND.MODE==0) CALL OPEN_PVD_FILE

! First pass write the data header.
! Second pass writes the data (appended binary format).

      DO PASS=WRITE_HEADER,WRITE_DATA


         CALL WRITE_GEOMETRY_IN_VTP_BIN(PASS)

         IF(VTK_PART_DIAMETER(VTK_REGION)) &
            CALL WRITE_SCALAR_IN_VTP_BIN('Diameter',2.0D0*DES_RADIUS,PASS)

         IF(VTK_PART_VEL(VTK_REGION)) &
            CALL WRITE_VECTOR_IN_VTP_BIN('Velocity',DES_VEL_NEW,PASS)

         IF(VTK_PART_ANGULAR_VEL(VTK_REGION)) &
            CALL WRITE_VECTOR_IN_VTP_BIN('Angular_velocity', OMEGA_NEW,PASS)

         IF(PARTICLE_ORIENTATION) THEN
            IF(VTK_PART_ORIENTATION(VTK_REGION)) &
               CALL WRITE_VECTOR_IN_VTP_BIN('Orientation', ORIENTATION,PASS)
         ENDIF

         DO N=1, DES_USR_VAR_SIZE
            IF(VTK_PART_USR_VAR(VTK_REGION,N)) &
              CALL WRITE_SCALAR_IN_VTP_BIN('User Defined Var '//trim(iVal(N)),DES_USR_VAR(N,:),PASS)
         ENDDO

         IF(ENERGY_EQ.AND.VTK_PART_TEMP(VTK_REGION)) &
           CALL WRITE_SCALAR_IN_VTP_BIN('Temperature', DES_T_s,PASS)

         IF(ANY_SPECIES_EQ) THEN
            DO N=1, DIMENSION_N_S
               IF(VTK_PART_X_s(VTK_REGION,N)) &
                 CALL WRITE_SCALAR_IN_VTP_BIN(trim(iVar('X_s',N)), DES_X_s(:,N),PASS)
            ENDDO
         ENDIF

      IF(USE_COHESION.AND.VTK_PART_COHESION(VTK_REGION)) &
         CALL WRITE_SCALAR_IN_VTP_BIN('CohesiveForce', PostCohesive,PASS)

      ENDDO ! PASS LOOP, EITHER HEADER OR DATA


      CALL CLOSE_VTP_FILE_BIN(MODE)

! Only update pvd file when there are particles in vtk region
      IF(GLOBAL_CNT>0.AND.MODE==0) CALL UPDATE_AND_CLOSE_PVD_FILE

#ifdef MPI
      call MPI_barrier(MPI_COMM_WORLD,mpierr)
#endif

! Update Frames
      IF (myPE == PE_IO.AND.TIME_DEPENDENT_FILENAME) THEN
         OPEN(CONVERT='BIG_ENDIAN',UNIT = VTU_FRAME_UNIT, FILE = TRIM(VTU_FRAME_FILENAME))
         DO L = 1, DIMENSION_VTK
            IF(VTK_DEFINED(L)) WRITE(VTU_FRAME_UNIT,*) L,FRAME(L)
         ENDDO
         CLOSE(VTU_FRAME_UNIT)
      ENDIF

     IF (FULL_LOG.AND.myPE == PE_IO) WRITE(*,20)' DONE.'

20    FORMAT(A,1X/)
      RETURN

      END SUBROUTINE WRITE_VTP_FILE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_VTP_FILE                                          C
!  Purpose: Open a vtp file and writes the header                      C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_VTP_FILE_BIN(MODE)

      USE run, only: TIME
      USE output, only: FULL_LOG
      USE vtk, only: TIME_DEPENDENT_FILENAME, VTU_FRAME_FILENAME, VTU_FRAME_UNIT
      USE vtk, only: RESET_FRAME_AT_TIME_ZERO,PVTU_FILENAME,PVTU_UNIT,BUFFER,END_REC
      USE vtk, only: DIMENSION_VTK, VTK_DEFINED,  FRAME,  VTK_REGION
      USE vtk, only: VTU_FILENAME, VTK_FILEBASE, VTU_DIR, VTU_UNIT
      USE param1, only: ZERO

      IMPLICIT NONE
      LOGICAL :: VTU_FRAME_FILE_EXISTS, NEED_TO_WRITE_VTP
      INTEGER :: ISTAT,BUFF1,BUFF2,L
      INTEGER :: MODE   ! MODE = 0 : Write regular VTK region file
                        ! MODE = 1 : Write debug   VTK region file (VTK_DBG_FILE = .TRUE.)


      IF(BDIST_IO) THEN
         NEED_TO_WRITE_VTP = (LOCAL_CNT>0)
      ELSE
         NEED_TO_WRITE_VTP = (MyPE==0.AND.GLOBAL_CNT>0)
      ENDIF

! Only open the file from head node when not using distributed I/O
      IF (myPE /= PE_IO.AND.(.NOT.BDIST_IO)) RETURN

      IF(TIME_DEPENDENT_FILENAME) THEN
         INQUIRE(FILE=VTU_FRAME_FILENAME,EXIST=VTU_FRAME_FILE_EXISTS)
         IF(VTU_FRAME_FILE_EXISTS) THEN
            OPEN(CONVERT='BIG_ENDIAN',UNIT = VTU_FRAME_UNIT, FILE = TRIM(VTU_FRAME_FILENAME))
            DO L = 1, DIMENSION_VTK
               IF(VTK_DEFINED(L)) THEN
                  READ(VTU_FRAME_UNIT,*)BUFF1,BUFF2
                  FRAME(L)=BUFF2
               ENDIF
            ENDDO
            CLOSE(VTU_FRAME_UNIT)
         ENDIF
         IF(RESET_FRAME_AT_TIME_ZERO.AND.TIME==ZERO) THEN
            DO L = 1, DIMENSION_VTK
               IF(L==VTK_REGION) FRAME(L)=-1
            ENDDO
         ENDIF
         DO L = 1, DIMENSION_VTK
            IF(L==VTK_REGION) FRAME(L) = FRAME(L) + 1
         ENDDO
      ENDIF

! For distributed I/O, define the file name for each processor that owns particles
      IF (BDIST_IO) THEN
         IF (LOCAL_CNT>0) THEN
            IF(TIME_DEPENDENT_FILENAME.AND.MODE==0) THEN
               WRITE(VTU_FILENAME,20) TRIM(VTK_FILEBASE(VTK_REGION)),FRAME(VTK_REGION),MYPE
            ELSE
               WRITE(VTU_FILENAME,25) TRIM(VTK_FILEBASE(VTK_REGION)),MYPE
            ENDIF
         ENDIF
      ELSE
         IF(MYPE.EQ.PE_IO) THEN
            IF(TIME_DEPENDENT_FILENAME.AND.MODE==0) THEN
               WRITE(VTU_FILENAME,30) TRIM(VTK_FILEBASE(VTK_REGION)),FRAME(VTK_REGION)
            ELSE
               WRITE(VTU_FILENAME,35) TRIM(VTK_FILEBASE(VTK_REGION))
            ENDIF
         END IF
      END IF

! Add the VTU directory path if necessary

      IF (NEED_TO_WRITE_VTP) THEN
         IF(TRIM(VTU_DIR)/='.') VTU_FILENAME='./'//TRIM(VTU_DIR)//'/'//VTU_FILENAME
      ENDIF

! Echo
      IF (FULL_LOG) THEN
         IF (.NOT.BDIST_IO) THEN
            WRITE(*,10,ADVANCE='NO')' WRITING VTP FILE : ', TRIM(VTU_FILENAME),' .'
         ELSE
            IF(myPE==PE_IO) WRITE(*,15,ADVANCE='NO')' EACH PROCESOR IS WRITING ITS OWN VTP FILE.'
         ENDIF
      ENDIF

! Open File

      IF (NEED_TO_WRITE_VTP) THEN

         VTU_UNIT = 678
         OPEN(CONVERT='BIG_ENDIAN',UNIT     = VTU_UNIT,           &
              FILE     = TRIM(VTU_FILENAME), &
              FORM     = 'UNFORMATTED',      &  ! works with gfortran 4.3.4 and ifort 10.1 but may not be supported by all compilers
                                                ! use 'BINARY' if 'UNFORMATTED' is not supported
              ACCESS   = 'STREAM',           &  ! works with gfortran 4.3.4 and ifort 10.1 but may not be supported by all compilers
                                                ! use 'SEQUENTIAL' if 'STREAM' is not supported
              ACTION   = 'WRITE',            &
              IOSTAT=ISTAT)


         IF (ISTAT /= 0) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1001) VTU_FILENAME, VTU_UNIT,VTU_DIR
            IF(FULL_LOG.AND.myPE == PE_IO) WRITE(*, 1001) VTU_FILENAME, VTU_UNIT,VTU_DIR
            CALL MFIX_EXIT(myPE)
         ENDIF


1001 FORMAT(/1X,70('*')//, ' From: OPEN_VTP_FILE',/,' Message: ',          &
            'Error opening vtp file. Terminating run.',/10X,          &
            'File name:  ',A,/10X,                                         &
            'VTP_UNIT :  ',i4, /10X,                                       &
            'PLEASE VERIFY THAT VTU_DIR EXISTS: ', A, &
            /1X,70('*')/)


! Write file Header
         BUFFER='<?xml version="1.0"?>'
         WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

         WRITE(BUFFER,*)'<!-- Time =',TIME,' sec. -->'
         WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

         BUFFER='<VTKFile type="PolyData" version="0.1" byte_order="BigEndian">'
         WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

         BUFFER='  <PolyData>'
         WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


      ENDIF
! For distributed I/O, open .p))vtp file that combines all *.vtp files for a given FRAME
! this is a simple ASCII file

      IF (myPE == PE_IO.AND.BDIST_IO.AND.GLOBAL_CNT>0) THEN

         IF(TIME_DEPENDENT_FILENAME.AND.MODE==0) THEN
            WRITE(PVTU_FILENAME,40) TRIM(VTK_FILEBASE(VTK_REGION)),FRAME(VTK_REGION)
         ELSE
            WRITE(PVTU_FILENAME,45) TRIM(VTK_FILEBASE(VTK_REGION))
         ENDIF

         IF(TRIM(VTU_DIR)/='.') PVTU_FILENAME='./'//TRIM(VTU_DIR)//'/'//PVTU_FILENAME

         OPEN(CONVERT='BIG_ENDIAN',UNIT = PVTU_UNIT, FILE = TRIM(PVTU_FILENAME))

         WRITE(PVTU_UNIT,100) '<?xml version="1.0"?>'
         WRITE(PVTU_UNIT,110) '<!-- Time =',TIME,' sec. -->'
         WRITE(PVTU_UNIT,120) '<VTKFile type="PPolyData"',&
                  ' version="0.1" byte_order="BigEndian">'

         WRITE(PVTU_UNIT,100) '  <PPolyData GhostLevel="0">'
         WRITE(PVTU_UNIT,100) '      <PPoints>'
         WRITE(PVTU_UNIT,100) '        <PDataArray type="Float32" Name="coordinates" NumberOfComponents="3" &
              &format="appended" offset=" 0" />'
         WRITE(PVTU_UNIT,100) '      </PPoints>'
         WRITE(PVTU_UNIT,100) ''
         WRITE(PVTU_UNIT,100) '      <PPointData Scalars="Diameter" Vectors="Velocity">'

      ENDIF

100   FORMAT(A)
110   FORMAT(A,E14.7,A)
120   FORMAT(A,A)
10    FORMAT(/1X,3A)
15    FORMAT(/1X,A)
20    FORMAT(A,"_",I4.4,"_",I5.5,".vtp")
25    FORMAT(A,"_",I5.5,".vtp")
30    FORMAT(A,"_",I4.4,".vtp")
35    FORMAT(A,".vtp")
40    FORMAT(A,"_",I4.4,".pvtp")
45    FORMAT(A,".pvtp")

      RETURN

      END SUBROUTINE OPEN_VTP_FILE_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_GEOMETRY_IN_VTP_BIN                              C
!  Purpose: Write Geometry and connectivity in a vtu file              C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_GEOMETRY_IN_VTP_BIN(PASS)

      USE vtk, only: NUMBER_OF_POINTS,BUFFER, VTU_UNIT,END_REC,VTU_OFFSET,BELONGS_TO_VTK_SUBDOMAIN

      IMPLICIT NONE

      REAL(c_float) :: float
      INTEGER(c_int) :: int

      INTEGER ::     nbytes_vector
      INTEGER ::     offset_xyz

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)  ! local
      DOUBLE PRECISION, ALLOCATABLE :: gtemp_array(:,:)  ! global

      INTEGER :: LB, UB
      INTEGER :: PC, LC1, LC2

! Loop through all particles and kee a list of particles belonging to a VTK region

! Since the data is appended (i.e., written after all tags), the
! offset, in number of bytes must be specified.  The offset includes
! the size of the data for each field, plus the size of the integer
! that stores the number of bytes.  this is why the offset of a field
! equals the offset of the previous field plus sizeof(int) plus the
! number of bytes of the field.

! Next, the actual data is written for the geometry (PASS=WRITE_DATA)
! The DATA is converted to single precision to save memory.

      IF (.NOT.BDIST_IO) THEN
! The number of points in the pvd file is the global number of particles
! computed from SETUP_VTK_REGION_PARTICLES

         NUMBER_OF_POINTS = GLOBAL_CNT

! Number of bytes of position field (vector,3 components)
         nbytes_vector       = NUMBER_OF_POINTS * 3 * sizeof(float)

! Offset of each field
         offset_xyz = 0

         IF(PASS==WRITE_HEADER) THEN
            IF(myPE == PE_IO) THEN

               WRITE(BUFFER,*)'    <Piece NumberOfPoints="',NUMBER_OF_POINTS, &
                     '"  NumberOfVerts="0" NumberOfLines ="0" NumberOfStrips="0" NumberOfPolys="0" >'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      <Points>'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'        <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" &
                                       &format="appended" offset="',offset_xyz,'" />'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      </Points>'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'<PointData Scalars="Diameter" Vectors="Velocity"> '!preparing pointData
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

! calculate offset for next field
               VTU_offset = offset_xyz + sizeof(int) + nbytes_vector

            ENDIF

         ELSEIF(PASS==WRITE_DATA) THEN

            IF(myPE == PE_IO) THEN

               WRITE(BUFFER,*)'      </PointData>'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      <Verts> </Verts>'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      <Lines> </Lines>'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      <Strips> </Strips>'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'      <Polys> </Polys>'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'    </Piece>'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'  </PolyData>'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

               WRITE(BUFFER,*)'  <AppendedData encoding="raw">'
               WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC


! Starting raw binary data with an underscore

               WRITE(BUFFER,*)'_'
               WRITE(VTU_UNIT)TRIM(BUFFER)

! Number of bytes for X,Y,Z coordinates
            WRITE(VTU_UNIT) nbytes_vector


         ENDIF

         LB = LBOUND(DES_POS_NEW,2) ! This should always be 1
         UB = UBOUND(DES_POS_NEW,2) ! This should always be 2

         ALLOCATE (dProcBuf(LOCAL_CNT) )
         ALLOCATE (dRootBuf(GLOBAL_CNT))
         ALLOCATE (ltemp_array((UB-LB)+1,LOCAL_CNT))
         ALLOCATE (gtemp_array((UB-LB)+1,GLOBAL_CNT))

! Pack particle coordinates in a temporary local array
         PC = 0
         DO LC1 = 1, MAX_PIP
            IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
               PC =PC + 1
               DO LC2=LB, UB
                  ltemp_array(LC2,PC) = DES_POS_NEW(LC1,LC2)
               ENDDO
            ENDIF
            IF(PC==LOCAL_CNT) EXIT
         ENDDO

! For each coordinate (x,y, and z), gather the local list to global temporary array
         DO LC1 = LB, UB
            dprocbuf(1:LOCAL_CNT)=ltemp_array(LC1,1:LOCAL_CNT)
            CALL desmpi_gatherv(ptype=2)
            gtemp_array(LC1,:) = drootbuf(:)
         ENDDO

! Write the list of coordinates
         IF(myPE == PE_IO) THEN
            DO LC1=1, GLOBAL_CNT
               DO LC2=LB, UB
                  WRITE(VTU_UNIT)  real(gtemp_array(LC2,LC1))
               ENDDO
            ENDDO
         ENDIF

         deallocate (dProcBuf, dRootBuf, ltemp_array,gtemp_array)


         ENDIF


      ELSEIF(BDIST_IO.AND.LOCAL_CNT>0) THEN

         IF(LOCAL_CNT==0) RETURN
! The number of points in the pvd file is the local number of particles
! computed from SETUP_VTK_REGION_PARTICLES

         NUMBER_OF_POINTS = LOCAL_CNT

! Number of bytes of position field (vector,3 components)
         nbytes_vector       = NUMBER_OF_POINTS * 3 * sizeof(float)

! Offset of each field
         offset_xyz = 0

         IF(PASS==WRITE_HEADER) THEN

            WRITE(BUFFER,*)'    <Piece NumberOfPoints="',NUMBER_OF_POINTS, &
                  '"  NumberOfVerts="0" NumberOfLines ="0" NumberOfStrips="0" NumberOfPolys="0" >'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'        <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" &
                                    &format="appended" offset="',offset_xyz,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      </Points>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'<PointData Scalars="Diameter" Vectors="Velocity"> '!preparing pointData
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

! calculate offset for next field
            VTU_offset = offset_xyz + sizeof(int) + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

            WRITE(BUFFER,*)'      </PointData>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Verts> </Verts>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Lines> </Lines>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Strips> </Strips>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'      <Polys> </Polys>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'    </Piece>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'  </PolyData>'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

            WRITE(BUFFER,*)'  <AppendedData encoding="raw">'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

! Starting raw binary data with an underscore

            WRITE(BUFFER,*)'_'
            WRITE(VTU_UNIT)TRIM(BUFFER)

! Number of bytes for X,Y,Z coordinates
            WRITE(VTU_UNIT) nbytes_vector

            LB = LBOUND(DES_POS_NEW,2) ! This should always be 1
            UB = UBOUND(DES_POS_NEW,2) ! This should always be 2

            ALLOCATE (ltemp_array((UB-LB)+1,LOCAL_CNT))

! Pack particle coordinates in a temporary local array
            PC = 0
            DO LC1 = 1, MAX_PIP
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  DO LC2=LB, UB
                     ltemp_array(LC2,PC) = DES_POS_NEW(LC1,LC2)
                  ENDDO
               ENDIF
               IF(PC==LOCAL_CNT) EXIT
            ENDDO


! Write the list of coordinates
            DO LC1=1, LOCAL_CNT
               DO LC2=LB, UB
                  WRITE(VTU_UNIT)  real(ltemp_array(LC2,LC1))
               ENDDO
            ENDDO

            deallocate (ltemp_array)

         ENDIF

      ENDIF

      RETURN

      END SUBROUTINE WRITE_GEOMETRY_IN_VTP_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SCALAR_IN_VTP_BIN                                C
!  Purpose: Write Scalar variable in a vtp file                        C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_SCALAR_IN_VTP_BIN(VAR_NAME,VAR,PASS)

      USE vtk, only: BUFFER,VTU_OFFSET,VTU_UNIT,PVTU_UNIT
      USE vtk, only: END_REC,BELONGS_TO_VTK_SUBDOMAIN
      USE output, only: FULL_LOG

      IMPLICIT NONE
      INTEGER :: I,LC1,PC

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, INTENT(in) :: VAR(:)

      REAL(c_float) :: float

      INTEGER :: nbytes_scalar

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      IF (.NOT.BDIST_IO) THEN

! Number of bytes for each scalar field
         nbytes_scalar = GLOBAL_CNT * sizeof(float)

         IF(PASS==WRITE_HEADER) THEN

! Remove possible white space with underscore
            DO I = 1,LEN_TRIM(VAR_NAME)
               IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
            ENDDO

! For each scalar, write a tag, with corresponding offset
            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

! Prepare the offset for the next field
            VTU_offset = VTU_offset + sizeof(float) + nbytes_scalar


         ELSEIF(PASS==WRITE_DATA) THEN

           allocate (dProcBuf(LOCAL_CNT) )
           allocate (dRootBuf(GLOBAL_CNT))

! Pack scalar list in a local buffer before gathering to root
            PC = 0
            DO LC1 = 1, MAX_PIP
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  dProcBuf(PC) = VAR(LC1)
               ENDIF
               IF(PC==LOCAL_CNT) EXIT
            ENDDO

! Gather local buffer to root
         CALL desmpi_gatherv(ptype=2)

! Write the data, always preceded by its size in number of bytes
! Write root buffer to file
         WRITE(VTU_UNIT) nbytes_scalar

         IF(myPE == PE_IO) THEN
            DO LC1=1, GLOBAL_CNT
               WRITE(VTU_UNIT)  real(drootBuf(LC1))
            ENDDO
         ENDIF

         deallocate (dProcBuf, dRootBuf)


         ENDIF


      ELSEIF(BDIST_IO.AND.LOCAL_CNT>0) THEN

! Number of bytes for each scalar field
         nbytes_scalar = LOCAL_CNT * sizeof(float)

! Remove possible white space with underscore
         DO I = 1,LEN_TRIM(VAR_NAME)
            IF(VAR_NAME(I:I) == ' ') VAR_NAME(I:I) = '_'
         ENDDO

         IF(PASS==WRITE_HEADER) THEN

! For each scalar, write a tag, with corresponding offset
            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

! Prepare the offset for the next field
            VTU_offset = VTU_offset + sizeof(float) + nbytes_scalar


         ELSEIF(PASS==WRITE_DATA) THEN

            allocate (dProcBuf(LOCAL_CNT) )

! Pack scalar list in a local buffer before writing in file
            PC = 0
            DO LC1 = 1, MAX_PIP
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  dProcBuf(PC) = VAR(LC1)
               ENDIF
               IF(PC==LOCAL_CNT) EXIT
            ENDDO

! Write the data, always preceded by its size in number of bytes
! Write root buffer to file
            WRITE(VTU_UNIT) nbytes_scalar

            DO LC1=1, LOCAL_CNT
               WRITE(VTU_UNIT)  real(dProcBuf(LC1))
            ENDDO

            deallocate (dProcBuf)

            IF (myPE == PE_IO) THEN       ! Update pvtu file with variable name
               WRITE(PVTU_UNIT,90) '        <PointArray type="Float32" Name="', &
                    TRIM(VAR_NAME),'" format="appended" offset="',VTU_offset,'" />'
            ENDIF


         ENDIF


      ENDIF


      IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10,ADVANCE='NO')'.'

10    FORMAT(A)
90    FORMAT(A,A,A,I12,A)

      RETURN

      END SUBROUTINE WRITE_SCALAR_IN_VTP_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_VECTOR_IN_VTP                                    C
!  Purpose: Write Vector variable in a vtp file                        C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE WRITE_VECTOR_IN_VTP_BIN(VAR_NAME,VAR,PASS)

      USE vtk, only: BUFFER,VTU_OFFSET,VTU_UNIT,PVTU_UNIT
      USE vtk, only: END_REC,BELONGS_TO_VTK_SUBDOMAIN
      USE output, only: FULL_LOG

      IMPLICIT NONE

      CHARACTER (*) :: VAR_NAME
      DOUBLE PRECISION, INTENT(in) :: VAR(:,:)

      REAL(c_float) :: float

      INTEGER :: nbytes_vector

      INTEGER :: PASS
      INTEGER :: WRITE_HEADER = 1
      INTEGER :: WRITE_DATA   = 2

      DOUBLE PRECISION, ALLOCATABLE :: ltemp_array(:,:)  ! local
      DOUBLE PRECISION, ALLOCATABLE :: gtemp_array(:,:)  ! global

      INTEGER :: LB, UB
      INTEGER :: PC, LC1, LC2

      IF (.NOT.BDIST_IO) THEN

! Number of bytes for each vector field
         nbytes_vector = GLOBAL_CNT * 3 * sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
! For each vector, write a tag, with corresponding offset

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

! Prepare the offset for the next field
            VTU_offset = VTU_offset + sizeof(float) + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

            LB = LBOUND(VAR,2) ! This should always be 1
            UB = UBOUND(VAR,2) ! This should always be 2

            ALLOCATE (dProcBuf(LOCAL_CNT) )
            ALLOCATE (dRootBuf(GLOBAL_CNT))
            ALLOCATE (ltemp_array((UB-LB)+1,LOCAL_CNT))
            ALLOCATE (gtemp_array((UB-LB)+1,GLOBAL_CNT))

! For each vector component, pack component list in a local array
            PC = 0
            DO LC1 = 1, MAX_PIP
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  DO LC2=LB, UB
                     ltemp_array(LC2,PC) = VAR(LC1,LC2)
                  ENDDO
               ENDIF
               IF(PC==LOCAL_CNT) EXIT
            ENDDO


! For each component, gather the local list to global temporary array
         DO LC1 = LB, UB
            dprocbuf(1:LOCAL_CNT)=ltemp_array(LC1,1:LOCAL_CNT)
            CALL desmpi_gatherv(ptype=2)
            gtemp_array(LC1,:) = drootbuf(:)
         ENDDO

! Write the data, always preceded by its size in number of bytes
         IF(myPE == PE_IO) THEN
            WRITE(VTU_UNIT) nbytes_vector
            DO LC1=1, GLOBAL_CNT
               DO LC2=LB, UB
                  WRITE(VTU_UNIT)  real(gtemp_array(LC2,LC1))
               ENDDO
            ENDDO
         ENDIF

         deallocate (dProcBuf, dRootBuf, ltemp_array,gtemp_array)



         ENDIF


      ELSEIF(BDIST_IO.AND.LOCAL_CNT>0) THEN

! Number of bytes for each vector field
         nbytes_vector = LOCAL_CNT * 3 * sizeof(float)

         IF(PASS==WRITE_HEADER) THEN
! For each vector, write a tag, with corresponding offset

            WRITE(BUFFER,90)'        <DataArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
            WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

! Prepare the offset for the next field
            VTU_offset = VTU_offset + sizeof(float) + nbytes_vector


         ELSEIF(PASS==WRITE_DATA) THEN

            LB = LBOUND(VAR,1) ! This should always be 1
            UB = UBOUND(VAR,1) ! This should always be 2

            ALLOCATE (ltemp_array((UB-LB)+1,LOCAL_CNT))

! For each vector component, pack component list in a local array
            PC = 0
            DO LC1 = 1, MAX_PIP
               IF(BELONGS_TO_VTK_SUBDOMAIN(LC1)) THEN
                  PC =PC + 1
                  DO LC2=LB, UB
                     ltemp_array(LC2,PC) = VAR(LC2,LC1)
                  ENDDO
               ENDIF
               IF(PC==LOCAL_CNT) EXIT
            ENDDO


! Write the data, always preceded by its size in number of bytes
            WRITE(VTU_UNIT) nbytes_vector
            DO LC1=1, LOCAL_CNT
               DO LC2=LB, UB
                  WRITE(VTU_UNIT)  real(ltemp_array(LC2,LC1))
               ENDDO
            ENDDO

            deallocate (ltemp_array)


            IF (myPE == PE_IO) THEN       ! Update pvtu file with variable name
               WRITE(PVTU_UNIT,90)'        <PointArray type="Float32" Name="', &
                 TRIM(VAR_NAME),'"  NumberOfComponents="3" format="appended" offset="',VTU_offset,'" />'
            ENDIF


         ENDIF

      ENDIF


      IF (PASS==WRITE_DATA.AND.FULL_LOG.AND.myPE == PE_IO) WRITE(*,10,ADVANCE='NO')'.'

10    FORMAT(A)
90    FORMAT(A,A,A,I12,A)

      RETURN

      END SUBROUTINE WRITE_VECTOR_IN_VTP_BIN

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLOSE_VTP_FILE_BIN                                     C
!  Purpose: Close a vtp file                                           C
!           Binary format                                              C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLOSE_VTP_FILE_BIN(MODE)

      USE vtk, only: BUFFER,VTU_UNIT,END_REC,PVTU_UNIT,TIME_DEPENDENT_FILENAME
      USE vtk, only: VTK_REGION,VTK_FILEBASE,FRAME

      IMPLICIT NONE

      INTEGER:: N
      CHARACTER (LEN=32)  :: VTU_NAME
      INTEGER, DIMENSION(0:numPEs-1) :: ALL_PART_CNT
      INTEGER :: IERR
      INTEGER :: MODE   ! MODE = 0 : Write regular VTK region file
                        ! MODE = 1 : Write debug   VTK region file (VTK_DBG_FILE = .TRUE.)


      IF((myPE == PE_IO.AND.(.NOT.BDIST_IO)).OR.(BDIST_IO.AND.LOCAL_CNT>0)) THEN

! Write last tags and close the vtp file
      WRITE(BUFFER,110)'  </AppendedData>'
      WRITE(VTU_UNIT)END_REC//TRIM(BUFFER)//END_REC

      WRITE(BUFFER,110)'</VTKFile>'
      WRITE(VTU_UNIT)TRIM(BUFFER)//END_REC

      CLOSE(VTU_UNIT)

      ENDIF

! Update pvtu file and close

      IF(BDIST_IO)  THEN
         CALL allgather_1i (LOCAL_CNT,ALL_PART_CNT,IERR)

         IF (myPE == PE_IO.AND.GLOBAL_CNT>0) THEN
            WRITE(PVTU_UNIT,100) '      </PPointData>'

            DO N = 0,NumPEs-1
               IF(ALL_PART_CNT(N)>0) THEN
                  IF(TIME_DEPENDENT_FILENAME.AND.MODE==0) THEN
                     WRITE(VTU_NAME,20) TRIM(VTK_FILEBASE(VTK_REGION)),FRAME(VTK_REGION),N
                  ELSE
                     WRITE(VTU_NAME,25) TRIM(VTK_FILEBASE(VTK_REGION)),N
                  ENDIF

                  WRITE(PVTU_UNIT,110) '      <Piece Source="',TRIM(VTU_NAME),'"/>'
               ENDIF
            ENDDO


            WRITE(PVTU_UNIT,100) '  </PPolyData>'
            WRITE(PVTU_UNIT,100) '</VTKFile>'
            CLOSE(PVTU_UNIT)
         ENDIF
      ENDIF

20    FORMAT(A,"_",I4.4,"_",I5.5,".vtp")
25    FORMAT(A,"_",I5.5,".vtp")

100   FORMAT(A)
110   FORMAT(A,A,A)

      RETURN

      END SUBROUTINE CLOSE_VTP_FILE_BIN


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SETUP_VTK_REGION_PARTICLES                             C
!                                                                      C
!  Purpose: Filter the particles  based on the VTK region bounds and   C
!           set the flag BELONGS_TO_VTK_SUBDOMAIN to .TRUE.            C
!           to keep the particle.                                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 11-Feb-15  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE SETUP_VTK_REGION_PARTICLES

      USE vtk, only: VTK_REGION
      USE vtk, only: VTK_X_E, VTK_X_W, VTK_Y_S, VTK_Y_N, VTK_Z_B, VTK_Z_T
      USE vtk, only: VTK_NXS, VTK_NYS, VTK_NZS
      USE vtk, only: VTK_SLICE_TOL, VTK_SELECT_MODE
      USE vtk, only: BELONGS_TO_VTK_SUBDOMAIN
      USE discretelement, only: MAX_PIP,PIP,DES_POS_NEW

      IMPLICIT NONE

      INTEGER :: PC,LC1
      INTEGER :: NXS,NYS,NZS,NS
      INTEGER :: X_SLICE(DIM_I),Y_SLICE(DIM_J),Z_SLICE(DIM_K)
      DOUBLE PRECISION :: XE,XW,YS,YN,ZB,ZT
      DOUBLE PRECISION :: XP,YP,ZP,XP1,YP1,ZP1,XP2,YP2,ZP2,R

      DOUBLE PRECISION :: SLICE_TOL
      LOGICAL :: KEEP_XDIR,KEEP_YDIR,KEEP_ZDIR

! Variables related to gather
      integer lgathercnts(0:numpes-1), lproc

      CHARACTER(LEN=1) :: SELECT_PARTICLE_BY

! Get VTK region bounds
      XE = VTK_X_E(VTK_REGION)
      XW = VTK_X_W(VTK_REGION)
      YS = VTK_Y_S(VTK_REGION)
      YN = VTK_Y_N(VTK_REGION)
      ZB = VTK_Z_B(VTK_REGION)
      ZT = VTK_Z_T(VTK_REGION)

      NXS = VTK_NXS(VTK_REGION)
      NYS = VTK_NYS(VTK_REGION)
      NZS = VTK_NZS(VTK_REGION)

      SLICE_TOL = VTK_SLICE_TOL(VTK_REGION)

      SELECT_PARTICLE_BY = VTK_SELECT_MODE(VTK_REGION)

! get slice(s) location
      DO NS = 1,NXS
         X_SLICE(NS) = XW + (XE-XW)/(NXS-1)*(NS-1)
      ENDDO

      DO NS = 1,NYS
         Y_SLICE(NS) = YS + (YN-YS)/(NYS-1)*(NS-1)
      ENDDO

      DO NS = 1,NZS
         Z_SLICE(NS) = ZB + (ZT-ZB)/(NZS-1)*(NS-1)
      ENDDO



! Loop through all particles on local rank and keep a list of particles
! belonging to VTK region

      IF(ALLOCATED(BELONGS_TO_VTK_SUBDOMAIN)) DEALLOCATE(BELONGS_TO_VTK_SUBDOMAIN)
      ALLOCATE(BELONGS_TO_VTK_SUBDOMAIN(MAX_PIP))

      BELONGS_TO_VTK_SUBDOMAIN = .FALSE.

      LOCAL_CNT = 0
      PC = 1
      DO LC1 = 1, MAX_PIP
         IF(PC > PIP) EXIT
         IF(IS_NONEXISTENT(LC1)) CYCLE
         PC = PC+1
         IF(IS_GHOST(LC1) .OR. IS_ENTERING_GHOST(LC1) .OR. IS_EXITING_GHOST(LC1)) CYCLE

         SELECT CASE(SELECT_PARTICLE_BY)
            CASE('C')  ! Particle center must be inside vtk region

               XP = DES_POS_NEW(LC1,1)
               YP = DES_POS_NEW(LC1,2)
               ZP = DES_POS_NEW(LC1,3)

! X-direction
               KEEP_XDIR=.FALSE.
               IF(NXS==0) THEN
                  IF(XW<=XP.AND.XP<=XE) KEEP_XDIR=.TRUE.
               ELSE
                  DO NS = 1,NXS
                     IF((X_SLICE(NS)-SLICE_TOL)<=XP.AND.XP<=(X_SLICE(NS)+SLICE_TOL)) KEEP_XDIR=.TRUE.
                  ENDDO
               ENDIF

! Y-direction
               KEEP_YDIR=.FALSE.
               IF(NYS==0) THEN
                  IF(YS<=YP.AND.YP<=YN) KEEP_YDIR=.TRUE.
               ELSE
                  DO NS = 1,NYS
                     IF((Y_SLICE(NS)-SLICE_TOL)<=YP.AND.YP<=(Y_SLICE(NS)+SLICE_TOL)) KEEP_YDIR=.TRUE.
                  ENDDO
               ENDIF

! Z-direction
               KEEP_ZDIR=.FALSE.
               IF(NZS==0) THEN
                  IF(ZB<=ZP.AND.ZP<=ZT) KEEP_ZDIR=.TRUE.
               ELSE
                  DO NS = 1,NZS
                     IF((Z_SLICE(NS)-SLICE_TOL)<=ZP.AND.ZP<=(Z_SLICE(NS)+SLICE_TOL)) KEEP_ZDIR=.TRUE.
                  ENDDO
               ENDIF


            CASE('P')  ! Entire particle must be inside vtk region

               R = DES_RADIUS(LC1)

               XP1 = DES_POS_NEW(LC1,1) - R
               YP1 = DES_POS_NEW(LC1,2) - R
               ZP1 = DES_POS_NEW(LC1,3) - R

               XP2 = DES_POS_NEW(LC1,1) + R
               YP2 = DES_POS_NEW(LC1,2) + R
               ZP2 = DES_POS_NEW(LC1,3) + R

! X-direction
               KEEP_XDIR=.FALSE.
               IF(NXS==0) THEN
                  IF(XW<=XP1.AND.XP2<=XE) KEEP_XDIR=.TRUE.
               ELSE
                  DO NS = 1,NXS
                     IF((X_SLICE(NS)-SLICE_TOL)<=XP1.AND.XP2<=(X_SLICE(NS)+SLICE_TOL)) KEEP_XDIR=.TRUE.
                  ENDDO
               ENDIF

! Y-direction
               KEEP_YDIR=.FALSE.
               IF(NYS==0) THEN
                  IF(YS<=YP1.AND.YP2<=YN) KEEP_YDIR=.TRUE.
               ELSE
                  DO NS = 1,NYS
                     IF((Y_SLICE(NS)-SLICE_TOL)<=YP1.AND.YP2<=(Y_SLICE(NS)+SLICE_TOL)) KEEP_YDIR=.TRUE.
                  ENDDO
               ENDIF

! Z-direction
               KEEP_ZDIR=.FALSE.
               IF(NZS==0) THEN
                  IF(ZB<=ZP1.AND.ZP2<=ZT) KEEP_ZDIR=.TRUE.
               ELSE
                  DO NS = 1,NZS
                     IF((Z_SLICE(NS)-SLICE_TOL)<=ZP1.AND.ZP2<=(Z_SLICE(NS)+SLICE_TOL)) KEEP_ZDIR=.TRUE.
                  ENDDO
               ENDIF


            CASE('I')  ! Particle must be inside or intersect the edge of the vtk region

               R = DES_RADIUS(LC1)

               XP1 = DES_POS_NEW(LC1,1) - R
               YP1 = DES_POS_NEW(LC1,2) - R
               ZP1 = DES_POS_NEW(LC1,3) - R

               XP2 = DES_POS_NEW(LC1,1) + R
               YP2 = DES_POS_NEW(LC1,2) + R
               ZP2 = DES_POS_NEW(LC1,3) + R

! X-direction
               KEEP_XDIR=.FALSE.
               IF(NXS==0) THEN
                  IF(.NOT.(XE<=XP1.OR.XP2<=XW)) KEEP_XDIR=.TRUE.
               ELSE
                  DO NS = 1,NXS
                     IF(.NOT.((X_SLICE(NS)+SLICE_TOL)<=XP1.OR.XP2<=(X_SLICE(NS)-SLICE_TOL))) KEEP_XDIR=.TRUE.
                  ENDDO
               ENDIF

! Y-direction
               KEEP_YDIR=.FALSE.
               IF(NYS==0) THEN
                  IF(.NOT.(YN<=YP1.OR.YP2<=YS)) KEEP_YDIR=.TRUE.
               ELSE
                  DO NS = 1,NYS
                     IF(.NOT.((Y_SLICE(NS)+SLICE_TOL)<=YP1.OR.YP2<=(Y_SLICE(NS)-SLICE_TOL))) KEEP_YDIR=.TRUE.
                  ENDDO
               ENDIF

! Z-direction
               KEEP_ZDIR=.FALSE.
               IF(NZS==0) THEN
                  IF(.NOT.(ZT<=ZP1.OR.ZP2<=ZB)) KEEP_ZDIR=.TRUE.
               ELSE
                  DO NS = 1,NZS
                     IF(.NOT.((Z_SLICE(NS)+SLICE_TOL)<=ZP1.OR.ZP2<=(Z_SLICE(NS)-SLICE_TOL))) KEEP_ZDIR=.TRUE.
                  ENDDO
               ENDIF


            CASE DEFAULT
                  print*,'should not be here'
         END SELECT

! Now combine
         IF(KEEP_XDIR.AND.KEEP_YDIR.AND.KEEP_ZDIR) THEN
            BELONGS_TO_VTK_SUBDOMAIN(LC1) = .TRUE.
            LOCAL_CNT = LOCAL_CNT + 1
         ENDIF
      ENDDO ! particle loop


! Calculate the total number of particles system-wide.
      call global_sum(LOCAL_CNT, GLOBAL_CNT)

! No need to set the send/reccv when using distributed IO
      IF (BDIST_IO) RETURN
! Set the send count from the local process.
      igath_sendcnt = LOCAL_CNT

! Collect the number of particles on each rank.all ranks.
      lgathercnts = 0
      lgathercnts(myPE) = LOCAL_CNT
      call global_sum(lgathercnts,igathercnts)

! Calculate the rank displacements.
      idispls(0) = 0
      DO lPROC = 1,NUMPEs-1
         idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)
      ENDDO

      RETURN

      END SUBROUTINE SETUP_VTK_REGION_PARTICLES

      END MODULE VTP
