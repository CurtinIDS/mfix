MODULE output_man
   CONTAINS
!----------------------------------------------------------------------!
!                                                                      !
!  Subroutine: OUTPUT_MANAGER                                          !
!  Author: J.Musser                                   Date:            !
!                                                                      !
!  Purpose: Relocate calls to write output files (RES, SPx, VTP). This !
!  was done to simplify the time_march code.                           !
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE OUTPUT_MANAGER(EXIT_SIGNAL, FINISHED)

! Global Variables:
!---------------------------------------------------------------------//

      use compar, only: myPE, PE_IO
      use discretelement, only: DISCRETE_ELEMENT
      use machine, only: wall_time
      use output, only: DISK, DISK_TOT
      use output, only: OUT_TIME, OUT_DT
      use output, only: RES_BACKUP_TIME, RES_BACKUP_DT
      use output, only: RES_TIME, RES_DT
      use output, only: SPX_TIME, SPX_DT
      use output, only: USR_TIME, USR_DT
      use param, only: DIMENSION_USR
      use param1, only: N_SPX
      use qmom_kinetic_equation, only: QMOMK
      use run, only: TIME, DT, TSTOP, STEADY_STATE
      use time_cpu, only: CPU_IO
      use vtk, only:    VTK_TIME, VTK_DT
      use vtk, only: DIMENSION_VTK
      use vtk, only: WRITE_VTK_FILES
      use vtp, only: write_vtp_file

      IMPLICIT NONE

! Dummy Arguments:
!---------------------------------------------------------------------//
! Flag that the the user specified batch time (plus buffer) is met.
      LOGICAL, INTENT(IN) :: EXIT_SIGNAL
! Flag that a steady state case is completed.
      LOGICAL, INTENT(IN) :: FINISHED

! Local Variables:
!---------------------------------------------------------------------//
! Loop counter and counter
      INTEGER :: LC, IDX
! Flag to write NetCDF output
      LOGICAL :: bWRITE_NETCDF_FILES
! Flag that the header (time) has not be written.
      LOGICAL :: HDR_MSG
! SPX file extensions.
      CHARACTER(LEN=35) ::  EXT_END
! Wall time at the start of IO operations.
      DOUBLE PRECISION :: WALL_START

!......................................................................!

! Initialize the SPx file extension array.
      EXT_END = '123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
! Initial the header flag.
      HDR_MSG = .TRUE.

! Get the current time before any IO operations begin
      WALL_START = WALL_TIME()

! Create a backup copy of the RES file.
      IF(TIME+0.1d0*DT>=RES_BACKUP_TIME) THEN
         RES_BACKUP_TIME = NEXT_TIME(RES_BACKUP_DT)
         CALL BACKUP_RES
      ENDIF

! Write SPx files, if needed
      IDX = 0
      bWRITE_NETCDF_FILES = .FALSE.

      DO LC=1, N_SPX
         IF(CHECK_TIME(SPX_TIME(LC))) THEN
            SPX_TIME(LC) = NEXT_TIME(SPX_DT(LC))

            CALL WRITE_SPX1(LC, 0)
            CALL NOTIFY_USER('SPx:',EXT_END(LC:LC))

            DISK_TOT = DISK_TOT + DISK(LC)
            IDX = IDX + 1

            bWRITE_NETCDF_FILES = .TRUE.
         ENDIF
      ENDDO
      IF(IDX /=0) CALL FLUSH_LIST


! Write standard output, if needed
      IF(CHECK_TIME(OUT_TIME)) THEN
         OUT_TIME = NEXT_TIME(OUT_DT)
         CALL WRITE_OUT1
         CALL NOTIFY_USER('.OUT;')
      ENDIF

! Write special output, if needed
      IDX = 0
      DO LC = 1, DIMENSION_USR
         IF(CHECK_TIME(USR_TIME(LC))) THEN
            USR_TIME(LC) = NEXT_TIME(USR_DT(LC))
            CALL WRITE_USR1 (LC)
            CALL NOTIFY_USER('.USR:',EXT_END(LC:LC))
            IDX = IDX + 1
         ENDIF
      ENDDO
      IF(IDX /=0) CALL FLUSH_LIST

      CALL FLUSH_NOTIFY_USER

! Write vtk file, if needed
! Only regular (not debug) files are written (second argument is zero)
      IF(WRITE_VTK_FILES) THEN
         DO LC = 1, DIMENSION_VTK
            IF(CHECK_TIME(VTK_TIME(LC))) THEN
               VTK_TIME(LC) = NEXT_TIME(VTK_DT(LC))
               CALL WRITE_VTU_FILE(LC,0)
               IF(DISCRETE_ELEMENT) CALL WRITE_VTP_FILE(LC,0)
            ENDIF
         ENDDO
      ENDIF

! Write NetCDF files.
      IF(bWRITE_NETCDF_FILES) CALL WRITE_NETCDF(0,0,TIME)

! Write restart file, if needed
      IF(CHECK_TIME(RES_TIME) .OR. EXIT_SIGNAL) THEN

         RES_TIME = NEXT_TIME(RES_DT)
         CALL WRITE_RES1
         CALL NOTIFY_USER('.RES;')

         IF(DISCRETE_ELEMENT) THEN
            CALL WRITE_RES0_DES
            CALL NOTIFY_USER('DES.RES;')
         ENDIF

         IF(QMOMK) THEN
            CALL QMOMK_WRITE_RESTART
            CALL NOTIFY_USER('QMOMK.RES;')
         ENDIF

      ENDIF

! Add the amount of time needed for all IO operations to total.
      CPU_IO = CPU_IO + (WALL_TIME() - WALL_START)

      RETURN

      contains

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      LOGICAL FUNCTION CHECK_TIME(lTIME)

      DOUBLE PRECISION, INTENT(IN) :: lTIME

      IF(STEADY_STATE) THEN
         CHECK_TIME = FINISHED
      ELSE
         CHECK_TIME = (TIME+0.1d0*DT>=lTIME).OR.(TIME+0.1d0*DT>=TSTOP)
      ENDIF

      RETURN
      END FUNCTION CHECK_TIME

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      DOUBLE PRECISION FUNCTION NEXT_TIME(lWRITE_DT)

      DOUBLE PRECISION, INTENT(IN) :: lWRITE_DT

      IF (.NOT.STEADY_STATE) THEN
         NEXT_TIME = (INT((TIME + 0.1d0*DT)/lWRITE_DT)+1)*lWRITE_DT
      ELSE
         NEXT_TIME = lWRITE_DT
      ENDIF

      RETURN
      END FUNCTION NEXT_TIME

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE NOTIFY_USER(MSG, EXT)

      use output, only: FULL_LOG
      use funits, only: DMP_LOG
      use funits, only: UNIT_LOG

      CHARACTER(len=*), INTENT(IN) :: MSG
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: EXT


      LOGICAL :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(HDR_MSG) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1000, ADVANCE='NO') TIME
         IF(SCR_LOG) WRITE(*, 1000, ADVANCE='NO') TIME
         HDR_MSG = .FALSE.
      ENDIF

 1000 FORMAT(' ',/' t=',F12.6,' Wrote')

      IF(.NOT.present(EXT)) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG, 1100, ADVANCE='NO') MSG
         IF(SCR_LOG) WRITE(*, 1100, ADVANCE='NO') MSG
      ELSE
         IF(IDX == 0) THEN
            IF(DMP_LOG) WRITE(UNIT_LOG, 1110, ADVANCE='NO') MSG, EXT
            IF(SCR_LOG) WRITE(*, 1110, ADVANCE='NO') MSG, EXT
         ELSE
            IF(DMP_LOG) WRITE(UNIT_LOG, 1120, ADVANCE='NO') EXT
            IF(SCR_LOG) WRITE(*, 1120, ADVANCE='NO') EXT
         ENDIF
      ENDIF

 1100 FORMAT(1X,A)
 1110 FORMAT(1X,A,1x,A)
 1120 FORMAT(',',A)

      RETURN
      END SUBROUTINE NOTIFY_USER

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE FLUSH_LIST

      use output, only: FULL_LOG
      use funits, only: DMP_LOG
      use funits, only: UNIT_LOG

      LOGICAL :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(DMP_LOG) WRITE(UNIT_LOG,1000, ADVANCE='NO')
      IF(SCR_LOG) WRITE(*,1000, ADVANCE='NO')

 1000 FORMAT(';')

      RETURN
      END SUBROUTINE FLUSH_LIST


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      SUBROUTINE FLUSH_NOTIFY_USER

      use discretelement, only: DISCRETE_ELEMENT, DES_CONTINUUM_COUPLED
      use discretelement, only: DTSOLID
      use error_manager
      use funits, only: DMP_LOG
      use funits, only: UNIT_LOG
      use machine, only: wall_time
      use run, only: get_tunit
      use output, only: FULL_LOG
      use output, only: NLOG
      use run, only: TIME, NSTEP, STEADY_STATE
      use time_cpu, only: TIME_START
      use time_cpu, only: WALL_START

      DOUBLE PRECISION :: WALL_ELAP, WALL_LEFT, WALL_NOW
      CHARACTER(LEN=9) :: CHAR_ELAP, CHAR_LEFT
      CHARACTER(LEN=4) :: UNIT_ELAP, UNIT_LEFT

      INTEGER :: TNITS
      LOGICAL :: SCR_LOG

      SCR_LOG = (FULL_LOG .and. myPE.eq.PE_IO)

      IF(.NOT.HDR_MSG) THEN
         IF(DMP_LOG) WRITE(UNIT_LOG,1000)
         IF(SCR_LOG) WRITE(*,1000)
      ENDIF

 1000 FORMAT(' ',/' ')

! Write the elapsed time and estimated remaining time
      IF(MOD(NSTEP,NLOG) == 0) THEN

         IF(DISCRETE_ELEMENT .AND. .NOT.DES_CONTINUUM_COUPLED) THEN
            TNITs = CEILING(real((TSTOP-TIME)/DTSOLID))
            WRITE(ERR_MSG, 1100) TIME, DTSOLID, trim(iVal(TNITs))
            CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE., LOG=.FALSE.)
         ENDIF
 1100 FORMAT(/'Time: ',g12.5,3x,'DT: ',g12.5,3x,'Remaining DEM NITs: ',A)

         WALL_NOW = WALL_TIME()
! Calculate the elapsed wall time.
         WALL_ELAP = WALL_NOW - WALL_START
         CALL GET_TUNIT(WALL_ELAP, UNIT_ELAP)
         CHAR_ELAP=''; WRITE(CHAR_ELAP,"(F9.2)") WALL_ELAP
         CHAR_ELAP = trim(adjustl(CHAR_ELAP))
! Estimate the remaining wall time.
         WALL_LEFT = (WALL_NOW-WALL_START)*(TSTOP-TIME)/               &
            max(TIME-TIME_START,1.0d-6)
         CALL GET_TUNIT(WALL_LEFT, UNIT_LEFT)

         IF (.NOT.STEADY_STATE) THEN
            CHAR_LEFT=''; WRITE(CHAR_LEFT,"(F9.2)") WALL_LEFT
            CHAR_LEFT = trim(adjustl(CHAR_LEFT))
         ELSE
            CHAR_LEFT = '0.0'
            UNIT_LEFT = 's'
         ENDIF

! Notify the user of usage/remaining wall times.
         WRITE(ERR_MSG,2000)                                           &
            'Elapsed:', trim(CHAR_ELAP), trim(UNIT_ELAP),              &
            'Est. Remaining:',trim(CHAR_LEFT), trim(UNIT_LEFT)
         CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
      ENDIF

 2000 FORMAT('Wall Time - ',2(A,1X,A,A,4X))

      RETURN
      END SUBROUTINE FLUSH_NOTIFY_USER

      END SUBROUTINE OUTPUT_MANAGER

!----------------------------------------------------------------------!
! Subroutine: INIT_OUTPUT_VARS                                         !
! Purpose: Initialize variables used for controling ouputs of the      !
! various files.                                                       !
!----------------------------------------------------------------------!
      SUBROUTINE INIT_OUTPUT_VARS

      use geometry, only: IJKMAX2
      use machine, only: wall_time
      use output, only: DISK, DISK_TOT
      use output, only: ONEMEG
      use output, only: OUT_TIME, OUT_DT
      use output, only: RES_TIME, RES_DT
      use output, only: SPX_TIME, SPX_DT
      use output, only: USR_TIME, USR_DT
      use output, only: RES_BACKUP_TIME, RES_BACKUP_DT
      use output, only: RES_BACKUPS
      use param, only: DIMENSION_USR
      use param1, only: N_SPX
      use param1, only: UNDEFINED
      use param1, only: ZERO
      use physprop, only: MMAX, NMAX
      use run, only: K_EPSILON
      use run, only: RUN_TYPE
      use run, only: TIME, DT, STEADY_STATE
      use rxns, only: nRR
      use scalars, only: NScalar
      use time_cpu, only: CPU_IO
      use time_cpu, only: TIME_START
      use time_cpu, only: WALL_START
      use vtk, only:    VTK_TIME, VTK_DT
      use vtk, only: DIMENSION_VTK
      use vtk, only: DIMENSION_VTK
      use vtk, only: VTK_TIME, VTK_DT
      use vtk, only: WRITE_VTK_FILES

      use param1, only:  UNDEFINED_I

      use funits, only: CREATE_DIR

      IMPLICIT NONE

! Disk space needed for one variable and each SPX file
      DOUBLE PRECISION :: DISK_ONE

! Loop counter
      INTEGER :: LC

! Initialize times for writing outputs
      OUT_TIME = merge(TIME, UNDEFINED, OUT_DT /= UNDEFINED)

! Initialize the amount of time spent on IO
      CPU_IO = 0.0d0

! Initialize disk space calculations
      DISK_TOT = ZERO
      DISK_ONE = 4.0*IJKMAX2/ONEMEG

      DISK(1) = 1.0*DISK_ONE                           ! EPg
      DISK(2) = 2.0*DISK_ONE                           ! Pg, Ps
      DISK(3) = 3.0*DISK_ONE                           ! Ug, Vg, Wg
      DISK(4) = 3.0*DISK_ONE*MMAX                      ! Us, Vs, Ws
      DISK(5) = 1.0*DISK_ONE*MMAX                      ! ROPs
      DISK(6) = 1.0*DISK_ONE*(MMAX+1)                  ! Tg, Ts
      DISK(7) = 1.0*DISK_ONE*(sum(NMAX(0:MMAX)))       ! Xg, Xs
      DISK(8) = 1.0*DISK_ONE*MMAX                      ! Theta
      DISK(9) = 1.0*DISK_ONE*NScalar                   ! User Scalars
      DISK(10) = nRR*DISK_ONE                          ! ReactionRates
      DISK(11) = merge(2.0*DISK_ONE, ZERO, K_EPSILON)  ! K-Epsilon


! Initizle RES and SPX_TIME
      IF (RUN_TYPE == 'NEW') THEN
         RES_TIME = TIME
         SPX_TIME(:N_SPX) = TIME
      ELSE
         IF (.NOT. STEADY_STATE) THEN
            RES_TIME = RES_DT *                                        &
               (INT((TIME + 0.1d0*DT)/RES_DT) + 1)
            SPX_TIME(:N_SPX) = SPX_DT(:N_SPX) *                        &
               (INT((TIME + 0.1d0*DT)/SPX_DT(:N_SPX)) + 1)
         ENDIF
      ENDIF

! Initizle RES_BACKUP_TIME
      RES_BACKUP_TIME = UNDEFINED
      IF(RES_BACKUP_DT /= UNDEFINED) RES_BACKUP_TIME =                 &
         RES_BACKUP_DT * (INT((TIME+0.1d0*DT)/RES_BACKUP_DT)+1)

! Initialize USR_TIME
      DO LC = 1, DIMENSION_USR
         USR_TIME(LC) = UNDEFINED
         IF (USR_DT(LC) /= UNDEFINED) THEN
            IF (RUN_TYPE == 'NEW') THEN
               USR_TIME(LC) = TIME
            ELSE
               USR_TIME(LC) = USR_DT(LC) *                             &
                  (INT((TIME+0.1d0*DT)/USR_DT(LC))+1)
            ENDIF
         ENDIF
      ENDDO

! Initialize VTK_TIME

      IF(WRITE_VTK_FILES) THEN
         DO LC = 1, DIMENSION_VTK
            VTK_TIME(LC) = UNDEFINED
            IF (VTK_DT(LC) /= UNDEFINED) THEN
               IF (RUN_TYPE == 'NEW'.OR.RUN_TYPE=='RESTART_2') THEN
                  VTK_TIME(LC) = TIME
               ELSE
                  VTK_TIME(LC) = VTK_DT(LC) *                          &
                     (INT((TIME + 0.1d0*DT)/VTK_DT(LC))+1)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

! Create a subdir for RES backup files.
      IF(RES_BACKUPS /= UNDEFINED_I) CALL CREATE_DIR('BACKUP_RES')

      WALL_START = WALL_TIME()
      TIME_START = TIME

      RETURN
      END SUBROUTINE INIT_OUTPUT_VARS

!----------------------------------------------------------------------!
! Subroutine: BACKUP_RES                                               !
! Purpose: Shift existing RES file backup files by one index, then     !
! create a copy of the current RES file.                               !
!----------------------------------------------------------------------!
      SUBROUTINE BACKUP_RES

      use compar, only: myPE, PE_IO
      use output, only: RES_BACKUPS
      use discretelement, only: DISCRETE_ELEMENT
      use param1, only: UNDEFINED_I

      IMPLICIT NONE

      CHARACTER(len=256) :: FNAME0, FNAME1

      INTEGER :: LC

      IF(myPE /= PE_IO) RETURN
      IF(RES_BACKUPS == UNDEFINED_I) RETURN

! Shift all the existing backups by one.
      DO LC=RES_BACKUPS,2,-1
         CALL SET_FNAME(FNAME0,'.RES', LC-1)
         CALL SET_FNAME(FNAME1,'.RES', LC)
         CALL SHIFT_RES(FNAME0, FNAME1, 'mv')

         IF(DISCRETE_ELEMENT) THEN
            CALL SET_FNAME(FNAME0,'_DES.RES', LC-1)
            CALL SET_FNAME(FNAME1,'_DES.RES', LC)
            CALL SHIFT_RES(FNAME0, FNAME1, 'mv')
         ENDIF
      ENDDO

! Copy RES to RES1
      CALL SET_FNAME(FNAME0, '.RES')
      CALL SET_FNAME(FNAME1, '.RES' ,1)
      CALL SHIFT_RES(FNAME0, FNAME1, 'cp')

      IF(DISCRETE_ELEMENT) THEN
         CALL SET_FNAME(FNAME0, '_DES.RES')
         CALL SET_FNAME(FNAME1, '_DES.RES' ,1)
         CALL SHIFT_RES(FNAME0, FNAME1, 'cp')
      ENDIF

      RETURN

      contains

!----------------------------------------------------------------------!
! Subroutine: SHIFT_RES                                                !
! Purpose: Shift RES(LC-1) to RES(LC)                                  !
!----------------------------------------------------------------------!
      SUBROUTINE SHIFT_RES(pFN0, pFN1, ACT)

      implicit none

      CHARACTER(LEN=*), INTENT(IN) :: pFN0, pFN1, ACT
      CHARACTER(len=1024) :: CMD
      LOGICAL :: EXISTS

      INQUIRE(FILE=trim(pFN0),EXIST=EXISTS)
      IF(EXISTS) THEN
         CMD=''; WRITE(CMD,1000)trim(ACT), trim(pFN0),trim(pFN1)
         CALL SYSTEM(trim(CMD))
      ENDIF

 1000 FORMAT(A,1x,A,1X,A)

      RETURN
      END SUBROUTINE SHIFT_RES

!----------------------------------------------------------------------!
! Subroutine: SET_FNAME                                                !
! Purpose: Set the backup RES file name based on pINDX.                !
!----------------------------------------------------------------------!
      SUBROUTINE SET_FNAME(pFNAME, pEXT, pINDX)

      use run, only: RUN_NAME

      implicit none

      CHARACTER(LEN=*), INTENT(OUT) :: pFNAME
      CHARACTER(LEN=*), INTENT(IN) ::  pEXT
      INTEGER, INTENT(IN), OPTIONAL :: pINDX

! Set the file format for backup copies
      pFNAME=''
      IF(.NOT.PRESENT(pINDX)) THEN
         WRITE(pFNAME,1000) trim(RUN_NAME),pEXT
      ELSE
         IF(RES_BACKUPS < 10) THEN
            WRITE(pFNAME,1001) trim(RUN_NAME), pEXT, pINDX
         ELSEIF(RES_BACKUPS < 100) THEN
            WRITE(pFNAME,1002) trim(RUN_NAME), pEXT, pINDX
         ELSEIF(RES_BACKUPS < 1000) THEN
            WRITE(pFNAME,1003) trim(RUN_NAME), pEXT, pINDX
         ELSEIF(RES_BACKUPS < 10000) THEN
            WRITE(pFNAME,1004) trim(RUN_NAME), pEXT, pINDX
         ELSEIF(RES_BACKUPS < 10000) THEN
            WRITE(pFNAME,1005) trim(RUN_NAME), pEXT, pINDX
         ELSE
            WRITE(pFNAME,1006) trim(RUN_NAME), pEXT, pINDX
         ENDIF
      ENDIF

 1000 FORMAT(2A)
 1001 FORMAT('BACKUP_RES/',2A,I1.1)
 1002 FORMAT('BACKUP_RES/',2A,I2.2)
 1003 FORMAT('BACKUP_RES/',2A,I3.3)
 1004 FORMAT('BACKUP_RES/',2A,I4.4)
 1005 FORMAT('BACKUP_RES/',2A,I5.5)
 1006 FORMAT('BACKUP_RES/',2A,I6.6)

      RETURN
      END SUBROUTINE SET_FNAME

      END SUBROUTINE BACKUP_RES
END MODULE output_man
