      MODULE READ_RES1_DES

      use cdist, only: bDist_IO
      use compar, only: PE_IO
      use compar, only: myPE
      use des_allocate
      use desmpi
      use error_manager
      use mpi_comm_des, only: DESMPI_GATHERV
      use mpi_comm_des, only: DESMPI_SCATTERV

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: INIT_READ_RES_DES
      PUBLIC :: FINL_READ_RES_DES

      PUBLIC :: READ_PAR_POS
      PUBLIC :: READ_PAR_COL

      PUBLIC :: READ_RES_DES
      PUBLIC :: READ_RES_pARRAY
      PUBLIC :: READ_RES_cARRAY

      INTERFACE READ_RES_DES
         MODULE PROCEDURE READ_RES_DES_0I
         MODULE PROCEDURE READ_RES_DES_1I
         MODULE PROCEDURE READ_RES_DES_0D
         MODULE PROCEDURE READ_RES_DES_1D
         MODULE PROCEDURE READ_RES_DES_0L
         MODULE PROCEDURE READ_RES_DES_1L
      END INTERFACE

      INTERFACE READ_RES_pARRAY
         MODULE PROCEDURE READ_RES_pARRAY_1B
         MODULE PROCEDURE READ_RES_pARRAY_1I
         MODULE PROCEDURE READ_RES_pARRAY_1D
         MODULE PROCEDURE READ_RES_pARRAY_1L
      END INTERFACE

      INTERFACE READ_RES_cARRAY
         MODULE PROCEDURE READ_RES_cARRAY_1I
         MODULE PROCEDURE READ_RES_cARRAY_1D
         MODULE PROCEDURE READ_RES_cARRAY_1L
      END INTERFACE


      INTEGER, PARAMETER :: RDES_UNIT = 901

      INTEGER :: pIN_COUNT
      INTEGER :: cIN_COUNT

! Send/Recv parameters for Particle arrays:
      INTEGER :: pROOTCNT, pPROCCNT
      INTEGER :: pRECV
      INTEGER, allocatable :: pSCATTER(:)
      INTEGER, allocatable :: pDISPLS(:)

! Variables used for reading restart file
      INTEGER, ALLOCATABLE :: pRestartMap(:)
      INTEGER, ALLOCATABLE :: cRestartMap(:)

! Send/Recv parameters for Particle arrays:
      INTEGER :: cROOTCNT, cPROCCNT
      INTEGER :: cRECV
      INTEGER, allocatable :: cSCATTER(:)
      INTEGER, allocatable :: cDISPLS(:)

      INTEGER, ALLOCATABLE :: iPAR_COL(:,:)

      CONTAINS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: INIT_READ_RES_DES                                        !
!                                                                      !
! Purpose: Construct the file name and open the DES RES file.          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE INIT_READ_RES_DES(BASE, lVERSION, lNEXT_REC)

      use discretelement, only: MAX_PIP, PIP
      use discretelement, only: iGHOST_CNT
      use discretelement, only: NEIGH_NUM

      use compar, only: numPEs
      use machine, only: OPEN_N1

      use mpi_utility, only: BCAST
      use mpi_utility, only: GLOBAL_ALL_SUM

      implicit none

      CHARACTER(len=*), INTENT(IN)  :: BASE
      DOUBLE PRECISION, INTENT(OUT) :: lVERSION
      INTEGER, INTENT(OUT) :: lNEXT_REC

      CHARACTER(len=32) :: lFNAME

! Integer Error Flag
      INTEGER :: IER


      allocate(pSCATTER(0:numPEs-1))
      allocate(pDISPLS(0:numPEs-1))

      allocate(cSCATTER(0:numPEs-1))
      allocate(cDISPLS(0:numPEs-1))


      IF(bDIST_IO) THEN

         WRITE(lFNAME,'(A,I4.4,A)') BASE//'_DES_',myPE,'.RES'
         OPEN(CONVERT='BIG_ENDIAN',UNIT=RDES_UNIT, FILE=lFNAME,        &
            FORM='UNFORMATTED', STATUS='UNKNOWN', ACCESS='DIRECT',     &
            RECL=OPEN_N1)

         READ(RDES_UNIT, REC=1) lVERSION
         READ(RDES_UNIT, REC=2) pIN_COUNT
         READ(RDES_UNIT, REC=3) iGHOST_CNT
         READ(RDES_UNIT, REC=4) cIN_COUNT

         IF(PIP > MAX_PIP) THEN
            write(*,*) "From des_read_restart:"
            write(*,*) "Error: The pip is grater than current max_pip"
            write(*,*) "pip=" ,pip,"; max_pip =", max_pip

         ENDIF

         PIP = pIN_COUNT
         NEIGH_NUM = cIN_COUNT

         CALL PARTICLE_GROW(NEIGH_NUM)

      ELSE

         IF(myPE == PE_IO) THEN
            WRITE(lFNAME,'(A,A)') BASE//'_DES.RES'
            OPEN(CONVERT='BIG_ENDIAN',UNIT=RDES_UNIT, FILE=lFNAME,     &
               FORM='UNFORMATTED', STATUS='UNKNOWN', ACCESS='DIRECT',  &
               RECL=OPEN_N1)

            READ(RDES_UNIT, REC=1) pIN_COUNT

            READ(RDES_UNIT, REC=1) lVERSION
            READ(RDES_UNIT, REC=2) pIN_COUNT
!           READ(RDES_UNIT, REC=3) -NOTHING-
            READ(RDES_UNIT, REC=4) cIN_COUNT

         ELSE
            pIN_COUNT = 10
         ENDIF

         IER = 0

! Allocate the particle restart map. This is used in determining were
! particle data is sent. Only process zero needs this array.
         allocate( pRestartMap(pIN_COUNT), STAT=IER)
         IF(IER/=0) THEN
            WRITE(ERR_MSG, 1200) 'pRestartMap', trim(iVAL(pIN_COUNT))
            CALL FLUSH_ERR_MSG
         ENDIF

         CALL BCAST(lVERSION, PE_IO)

! Allocate the collision restart map array. All ranks allocatet this
! array so that mapping the collision data can be done in parallel.
         CALL BCAST(cIN_COUNT, PE_IO)
         allocate( cRestartMap(cIN_COUNT), STAT=IER)
         IF(IER/=0) THEN
            WRITE(ERR_MSG, 1200) 'cRestartMap', trim(iVAL(cIN_COUNT))
            CALL FLUSH_ERR_MSG
         ENDIF

 1200 FORMAT('Error 1200: Unable to allocate sufficient memory to ',&
         'read in DES',/'restart file. size(',A,') = ',A)

         CALL GLOBAL_ALL_SUM(IER)
         IF(IER/=0) CALL MFIX_EXIT(myPE)

      ENDIF

      lNEXT_REC = 5

      RETURN
      END SUBROUTINE INIT_READ_RES_DES


!``````````````````````````````````````````````````````````````````````!
! Subroutine: CLOSE_RES_DES                                            !
!                                                                      !
! Purpose: Close the DES RES file.                                     !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE FINL_READ_RES_DES


      IF(bDIST_IO .OR. myPE == PE_IO) close(RDES_UNIT)

      IF(allocated(dPROCBUF)) deallocate(dPROCBUF)
      IF(allocated(dROOTBUF)) deallocate(dROOTBUF)
      IF(allocated(iPROCBUF)) deallocate(iPROCBUF)
      IF(allocated(iROOTBUF)) deallocate(iROOTBUF)

      IF(allocated(pRestartMap)) deallocate(pRestartMap)
      IF(allocated(cRestartMap)) deallocate(cRestartMap)

      IF(allocated(pSCATTER)) deallocate(pSCATTER)
      IF(allocated(pDISPLS)) deallocate(pDISPLS)

      IF(allocated(cSCATTER)) deallocate(cSCATTER)
      IF(allocated(cDISPLS)) deallocate(cDISPLS)




      RETURN
      END SUBROUTINE FINL_READ_RES_DES



!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_PAR_POS                                             !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_PAR_POS(lNEXT_REC)

      use discretelement, only: PIP
      use discretelement, only: DES_POS_NEW
      use geometry, only: NO_K
      use compar, only: numPEs

      use mpi_utility, only: GLOBAL_SUM
      USE in_binary_512

      implicit none

      INTEGER, INTENT(INOUT) :: lNEXT_REC

      INTEGER :: lDIMN
      INTEGER :: LC1, lPROC
      INTEGER :: lScatterCNTS(0:NUMPEs-1)
! The number of particles on each process.
      INTEGER :: PAR_CNT(0:NUMPEs-1)

!-----------------------------------------------

      CALL INIT_ERR_MSG("READ_PAR_POS")

      lDIMN = merge(2,3,NO_K)


! All process read positions for distributed IO restarts.
      IF(bDIST_IO) THEN
         DO LC1 = 1, lDIMN
            CALL READ_RES_DES(lNEXT_REC, DES_POS_NEW(:,LC1))
         ENDDO
         RETURN
      ENDIF

      allocate( dPAR_POS(pIN_COUNT, lDIMN))

! Only the IO proccess reads positions.
      IF(myPE == PE_IO) THEN
         DO LC1=1, merge(2,3,NO_K)
            CALL IN_BIN_512(RDES_UNIT, dPAR_POS(:,LC1),                &
               pIN_COUNT, lNEXT_REC)
         ENDDO
      ENDIF

! Use the particle postions and the domain coverage of each process
! to determine which processor each particle belongs.
      CALL MAP_pARRAY_TO_PROC(PAR_CNT)

! Send the particle position data to the individual ranks.
      CALL SCATTER_PAR_POS(PAR_CNT)

! Set up the read/scatter arrary information.
      pPROCCNT = PIP
      pROOTCNT = pIN_COUNT

! Set the recv count for this process.
      pRECV = PIP

! Construct an array for the Root process that states the number of
! (real) particles on each process.
      lScatterCnts(:) = 0; lScatterCnts(mype) = PIP
      CALL GLOBAL_SUM(lScatterCnts,pSCATTER)

! Calculate the displacements for each process in the global array.
      pDispls(0) = 0
      DO lPROC = 1, NUMPEs-1
         pDispls(lPROC) = pDispls(lPROC-1) + pSCATTER(lPROC-1)
      ENDDO

      IF(allocated(dPAR_POS)) deallocate(dPAR_POS)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE READ_PAR_POS


!``````````````````````````````````````````````````````````````````````!
! Subroutine: MAP_pARRAY_TO_PROC                                       !
!                                                                      !
! Purpose: Use the particle positions to determine which processor     !
! they live on and count the number of particles on each process.      !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE MAP_pARRAY_TO_PROC(lPAR_CNT)

      use discretelement, only: PIP
      use discretelement, only: XE, YN, ZT
      use geometry, only: IMIN1, IMAX1
      use geometry, only: JMIN1, JMAX1
      use geometry, only: KMIN1, KMAX1
      use geometry, only: NO_K, DO_K
      use compar, only: numPEs
      use compar, only: ISTART1_ALL, IEND1_ALL
      use compar, only: JSTART1_ALL, JEND1_ALL
      use compar, only: KSTART1_ALL, KEND1_ALL

      use mpi_utility, only: BCAST
      use mpi_utility, only: GLOBAL_ALL_SUM

      implicit none

      INTEGER, INTENT(OUT) :: lPAR_CNT(0:numPEs-1)

! Data dimensionality flag.
      INTEGER :: lDIMN
! Loop counters.
      INTEGER :: LC1, lPROC
! Error flag.
      INTEGER :: IER(0:numPEs-1)
! The X/Y/Z bounds of the physical space "owned" by each process.
      DOUBLE PRECISION :: lxmin(0:NUMPEs-1), lxmax(0:NUMPEs-1)
      DOUBLE PRECISION :: lymin(0:NUMPEs-1), lymax(0:NUMPEs-1)
      DOUBLE PRECISION :: lzmin(0:NUMPEs-1), lzmax(0:NUMPEs-1)
!-----------------------------------------------

      CALL INIT_ERR_MSG("MAP_pARRAY_TO_PROC")

! Initialize the error flag.
      IER = 0

      lDIMN = merge(2, 3, NO_K)

! set the domain range for each processor
      DO lPROC= 0, NUMPEs-1
         lxmin(lproc) = xe(istart1_all(lproc)-1)
         lxmax(lproc) = xe(iend1_all(lproc))
         lymin(lproc) = yn(jstart1_all(lproc)-1)
         lymax(lproc) = yn(jend1_all(lproc))
         lzmin(lproc) = zt(kstart1_all(lproc)-1)
         lzmax(lproc) = zt(kend1_all(lproc))

! modify the range for mass inlet and outlet, as particles injected
! can lie outside the domain and not ghost particles
         IF(istart1_all(lproc).eq.imin1) &
            lxmin(lproc) = xe(istart1_all(lproc)-2)
         IF(iend1_all(lproc).eq.imax1) &
            lxmax(lproc) = xe(iend1_all(lproc)+1)
         IF(jstart1_all(lproc).eq.jmin1) &
            lymin(lproc) = yn(jstart1_all(lproc)-2)
         IF(jend1_all(lproc).eq.jmax1)  &
            lymax(lproc) = yn(jend1_all(lproc)+1)
         IF(kstart1_all(lproc).eq.kmin1 .AND. DO_K) &
            lzmin(lproc) = zt(kstart1_all(lproc)-2)
         IF(kend1_all(lproc).eq.kmax1 .AND. DO_K) &
            lzmax(lproc) = zt(kend1_all(lproc)+1)
      ENDDO

! build the send buffer in PE_IO proc
! first pass to get the count of particles
      IER = 0
      pRestartMap(:) = -1
      lPAR_CNT(:) = 0
      IF(myPE == PE_IO) THEN
         DO LC1 = 1, pIN_COUNT
            DO lPROC=0, NUMPEs-1
               IF(dPAR_POS(LC1,1) >= lxmin(lproc) .AND. &
                  dPAR_POS(LC1,1) <  lxmax(lproc) .AND. &
                  dPAR_POS(LC1,2) >= lymin(lproc) .AND. &
                  dPAR_POS(LC1,2) <  lymax(lproc)) THEN
                  IF(NO_K)THEN
                     lPAR_CNT(lPROC) = lPAR_CNT(lPROC) + 1
                     pRestartMap(LC1) = lproc
                     EXIT
                  ELSE
                     IF(dPAR_POS(LC1,3) >= lzmin(lproc) .AND. &
                        dPAR_POS(LC1,3) <  lzmax(lproc)) THEN
                        lPAR_CNT(lPROC) = lPAR_CNT(lPROC) + 1
                        pRestartMap(LC1) = lproc
                        EXIT
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO  ! Loop over processes
            IF (pRestartMap(LC1) == -1) then
               IER(myPE) = -1
               WRITE(ERR_MSG,1000) trim(iVal(LC1))
               CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
               IF(NO_K) THEN
                  WRITE(ERR_MSG,1001) dPAR_POS(LC1,1:2)
                  CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
               ELSE
                  WRITE(ERR_MSG,1002) dPAR_POS(LC1,1:3)
                  CALL FLUSH_ERR_MSG(HEADER=.FALSE.)
               ENDIF
            ENDIF
         ENDDO  ! Loop over particles
      ENDIF

 1000 FORMAT('Error 1000: Unable to locate particle inside domain:',/&
         3x,'Particle Number:',A)
 1001 FORMAT(3x,'X POS: ',g12.5,/3x,'Y POS: ',g12.5)
 1002 FORMAT(3x,'X POS: ',g12.5,/3x,'Y POS: ',g12.5,/3x,'Z POS: ',g12.5)

! Send out the error flag and exit if needed.
      CALL BCAST(IER, PE_IO)
      IF(IER(PE_IO) /= 0) CALL MFIX_EXIT(myPE)

! PE_IO sends out the number of particles for each process.
      CALL BCAST(lPAR_CNT(0:NUMPES-1), PE_IO)

! Each process stores the number of particles-on-its-process. The error
! flag is set if that number exceeds the maximum.
      PIP = lPAR_CNT(myPE)
      CALL PARTICLE_GROW(PIP)

! Global collection of error flags to abort it the max was exceeded.
      CALL GLOBAL_ALL_SUM(IER)
      IF(sum(IER) /= 0) THEN
         WRITE(ERR_MSG,1100)
         CALL FLUSH_ERR_MSG(FOOTER=.FALSE.)
         DO LC1=0, numPEs-1
            IF(IER(LC1) /= 0) THEN
               WRITE(ERR_MSG,"(3(2x,I10))")LC1,IER(LC1)-1,lPAR_CNT(LC1)
               CALL FLUSH_ERR_MSG(HEADER=.FALSE., FOOTER=.FALSE.)
            ENDIF
         ENDDO
         WRITE(ERR_MSG,"('Aborting.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE.,ABORT=.TRUE.)
      ENDIF

 1100 FORMAT('Error 1100: Maximum number of particles exceeded.',2/    &
         5x,'Process',5x,'Maximum',7x,'Count')


      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MAP_pARRAY_TO_PROC



!``````````````````````````````````````````````````````````````````````!
! Subroutine: DES_RESTART_MAP                                          !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE SCATTER_PAR_POS(lPAR_CNT)

      use compar, only: numPEs

      use discretelement, only: DES_POS_NEW
      use discretelement, only: PIP
      use functions, only: SET_NORMAL
      use geometry, only: NO_K

      implicit none

! Number of particles on each process.
      INTEGER, INTENT(INOUT) :: lPAR_CNT(0:numPEs-1)
! Dimensionality flag.
      INTEGER :: lDIMN
! Loop counters.
      INTEGER :: LC1, lPROC, lBuf

      lDIMN = merge(2,3,NO_K)

! Set up the recv count and allocate the local process buffer.
      iSCR_RECVCNT = PIP*lDIMN
      allocate (dProcBuf(iscr_recvcnt))

! Allocate the buffer for the root.
      IF (myPE == PE_IO) THEN
         allocate (dRootBuf(pIN_COUNT*lDIMN))
      ELSE
         allocate (dRootBuf(10))
      ENDIF

! The IO processor builds drootbuffer and iDISLS
      IF(myPE == PE_IO) THEN
! Determine the offsets for each process and the amount of data that
! is to be scattered to each.
         iDISPLS(0) = 0
         iScatterCnts(0) = lPAR_CNT(0)*lDIMN
         DO lProc = 1, NUMPES-1
            iDispls(lproc) = iDispls(lproc-1) + iScatterCnts(lproc-1)
            iScatterCnts(lproc) = lPAR_CNT(lProc)*lDIMN
         ENDDO
! Copy the position data into the root buffer, mapped to the owner
! process.
         lPAR_CNT(:) = 0
         DO LC1 = 1,pIN_COUNT
            lPROC = pRestartMap(LC1)
            lbuf = iDispls(lProc) + lPAR_CNT(lProc)*lDIMN+1
            dRootBuf(lBuf:lBuf+lDIMN-1) = dPAR_POS(LC1,1:lDIMN)
            lBuf = lBuf + lDIMN
            lPAR_CNT(lProc) = lPAR_CNT(lProc) + 1
         ENDDO
      ENDIF
      CALL DESMPI_SCATTERV(pTYPE=2)

! Unpack the particle data.
      DO LC1 = 1, PIP
         lBuf = (LC1-1)*lDIMN+1
         DES_POS_NEW(LC1,1:lDIMN) = dProcBuf(lBuf:lBuf+lDIMN-1)
         lBuf = lBuf + lDIMN
         CALL SET_NORMAL(LC1)
      ENDDO

      IF(allocated(dRootBuf)) deallocate(dRootBuf)
      IF(allocated(dProcBuf)) deallocate(dProcBuf)
      IF(allocated(dPAR_POS)) deallocate(dPAR_POS)

      RETURN
      END SUBROUTINE SCATTER_PAR_POS

!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_PAR_COL                                             !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_PAR_COL(lNEXT_REC)

      use discretelement, only: NEIGHBORS, NEIGH_NUM
      use compar, only: numPEs

      use mpi_init_des, only: DES_RESTART_GHOST
      use mpi_utility, only: BCAST
      use mpi_utility, only: GLOBAL_SUM
      use mpi_utility, only: GLOBAL_ALL_SUM
      use in_binary_512i

      implicit none

      INTEGER, INTENT(INOUT) :: lNEXT_REC

      INTEGER :: LC1, lPROC
      INTEGER :: lScatterCNTS(0:NUMPEs-1)
! The number of particles on each process.
      INTEGER :: COL_CNT(0:NUMPEs-1)

!-----------------------------------------------

      CALL INIT_ERR_MSG("READ_PAR_COL")

! All process read positions for distributed IO restarts.
      IF(bDIST_IO) THEN
         CALL READ_RES_DES(lNEXT_REC, NEIGHBORS(:))
      ENDIF

      CALL DES_RESTART_GHOST

      allocate(iPAR_COL(2, cIN_COUNT))
      iPAR_COL = 0

! Only the IO proccess reads positions.
      IF(myPE == PE_IO) THEN
         DO LC1=1, 2
            CALL IN_BIN_512i(RDES_UNIT, iPAR_COL(LC1,:),               &
               cIN_COUNT, lNEXT_REC)
         ENDDO
      ENDIF

! Broadcast collision data to all the other processes.
       CALL GLOBAL_ALL_SUM(iPAR_COL)

! Determine which process owns the neighbor datasets. This is done either
! through matching global ids or a search. The actual method depends
! on the ability to allocate a large enough array.
      CALL MAP_cARRAY_TO_PROC(COL_CNT)

! Send the particle position data to the individual ranks.
      CALL GLOBAL_TO_LOC_COL

! Set up the read/scatter arrary information.
      cPROCCNT = NEIGH_NUM
      cROOTCNT = cIN_COUNT

! Set the recv count for this process.
      cRECV = NEIGH_NUM

! Construct an array for the Root process that states the number of
! (real) particles on each process.
      lScatterCnts(:) = 0; lScatterCnts(mype) = NEIGH_NUM
      CALL GLOBAL_SUM(lScatterCnts,cSCATTER)

! Calculate the displacements for each process in the global array.
      cDispls(0) = 0
      DO lPROC = 1, NUMPEs-1
         cDispls(lPROC) = cDispls(lPROC-1) + cSCATTER(lPROC-1)
      ENDDO

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE READ_PAR_COL


!``````````````````````````````````````````````````````````````````````!
! Subroutine: MAP_cARRAY_TO_PROC                                       !
!                                                                      !
! Purpose: Use the particle positions to determine which processor     !
! they live on and count the number of particles on each process.      !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE MAP_cARRAY_TO_PROC(lCOL_CNT)

      use compar, only: numPEs, myPE
      use discretelement, only: iGLOBAL_ID
      use discretelement, only: PIP
      use discretelement, only: NEIGH_NUM
      use functions, only: IS_GHOST, IS_ENTERING_GHOST, IS_EXITING_GHOST

      use mpi_utility, only: GLOBAL_ALL_SUM
      use mpi_utility, only: GLOBAL_ALL_MAX

      implicit none

      INTEGER, INTENT(OUT) :: lCOL_CNT(0:numPEs-1)

! Loop counters.
      INTEGER :: LC1, LC2
! Error flag.
      INTEGER :: IER
! Max global id.
      INTEGER :: MAX_ID, lSTAT

      INTEGER, ALLOCATABLE :: lGLOBAL_OWNER(:)

!-----------------------------------------------

      CALL INIT_ERR_MSG("MAP_cARRAY_TO_PROC")

! Initialize the error flag.
      IER = 0

      MAX_ID = maxval(IGLOBAL_ID(1:PIP))
      CALL GLOBAL_ALL_MAX(MAX_ID)

      allocate(lGLOBAL_OWNER(MAX_ID), STAT=lSTAT)
      CALL GLOBAL_ALL_SUM(lSTAT)

! All ranks successfully allocated the array. This permits a crude
! but much faster collision owner detection.
      IF(lSTAT == 0) THEN

         WRITE(ERR_MSG,"('Matching DES neighbor data by global owner.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE.,FOOTER=.FALSE.)

         lGLOBAL_OWNER = 0
         DO LC1=1, PIP
            IF(.NOT.IS_GHOST(LC1) .AND. .NOT.IS_ENTERING_GHOST(LC1) &
               .AND. .NOT.IS_EXITING_GHOST(LC1)) &
               lGLOBAL_OWNER(iGLOBAL_ID(LC1)) = myPE + 1
         ENDDO

! Loop over the neighbor list and match the read global ID to
! one of the global IDs.
         lCOL_CNT = 0
         cRestartMap = 0
         DO LC1=1, cIN_COUNT
            IF(lGLOBAL_OWNER(iPAR_COL(1,LC1)) == myPE + 1) THEN
               cRestartMap(LC1) = myPE + 1
               lCOL_CNT(myPE) = lCOL_CNT(myPE) + 1
            ENDIF
         ENDDO
! One or more ranks could not allocate the memory needed to do the
! quick and dirty match so do a search instead.
      ELSE

         WRITE(ERR_MSG,"('Matching DES neighbor data by search.')")
         CALL FLUSH_ERR_MSG(HEADER=.FALSE.,FOOTER=.FALSE.)

! Loop over the neighbor list and match the read global ID to
! one of the global IDs.
         lCOL_CNT = 0
         cRestartMap = 0
         LC1_LP: DO LC1=1, cIN_COUNT
            DO LC2=1, PIP!-iGHOST_CNT
               IF(iPAR_COL(1,LC1) == iGLOBAL_ID(LC2)) THEN
                  cRestartMap(LC1) = myPE + 1
                  lCOL_CNT(myPE) = lCOL_CNT(myPE) + 1
                  CYCLE LC1_LP
               ENDIF
            ENDDO
         ENDDO LC1_LP

      ENDIF

! Clean up the large array as it is no longer needed.
      IF(allocated(lGLOBAL_OWNER)) deallocate(lGLOBAL_OWNER)

! Calculate the number of matched collisions over all processes. Throw
! and error if it doesn't match the number of read collisions.
      CALL GLOBAL_ALL_SUM(lCOL_CNT)
      IF(sum(lCOL_CNT) /= cIN_COUNT) THEN
         WRITE(ERR_MSG,1000) cIN_COUNT, sum(lCOL_CNT)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

1000 FORMAT('Error 1000: Unable to establish the own of all read ',    &
         'collision data.',/3x,'Number of Collisions: ',I10,/3x,       &
         'Matched Collisions:   ',I10)

! Sync the collision restart map arcross all ranks.
      CALL GLOBAL_ALL_SUM(cRestartMap)

! Error checking and cleanup.
      DO LC1 = 1, cIN_COUNT
! Verify that each collision is owned by a rank.
         IF (cRestartMap(LC1) == 0) THEN
            IER = -1
            WRITE(ERR_MSG,1100) trim(iVal(LC1)), trim(iVal(            &
               iPAR_COL(1,LC1))), trim(iVal(iPAR_COL(2,LC1)))
            CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1100 FORMAT('Error 1100: Unable to locate process neighbor owner:',/  &
         3x,'Neighbor Number:',A,/3x,'Particles: ',A,' and ',A)

         ELSEIF(cRestartMap(LC1) > numPEs) THEN

            IER = -1
            WRITE(ERR_MSG,1101) trim(iVal(LC1)), trim(iVal(            &
              iPAR_COL(1,LC1))), trim(iVal(iPAR_COL(2,LC1)))
             CALL FLUSH_ERR_MSG(ABORT=.TRUE.)

 1101 FORMAT('Error 1101: More than one process neighbor owner:',/     &
         3x,'Neighbor Number:',A,/3x,'Particles: ',A,' and ',A)

! Shift the rank ID to the correct value.
         ELSE
            cRestartMap(LC1) = cRestartMap(LC1) - 1
         ENDIF
      ENDDO

! Send out the error flag and exit if needed.
      CALL GLOBAL_ALL_SUM(IER, PE_IO)
      IF(IER /= 0) CALL MFIX_EXIT(myPE)

! Each process stores the number of particles-on-its-process. The error
! flag is set if that number exceeds the maximum.
      NEIGH_NUM = lCOL_CNT(myPE)

      CALL NEIGHBOR_GROW(NEIGH_NUM)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE MAP_cARRAY_TO_PROC


!``````````````````````````````````````````````````````````````````````!
! Subroutine: GLOBAL_TO_LOC_COL                                        !
!                                                                      !
! Purpose: Generates the mapping used by the scatter routines to send  !
! read data to the correct rank.                                       !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE GLOBAL_TO_LOC_COL

      use discretelement, only: iGLOBAL_ID
      use discretelement, only: PIP

      use mpi_utility, only: GLOBAL_ALL_SUM
      use mpi_utility, only: GLOBAL_ALL_MAX

      use funits, only: DMP_LOG

      use error_manager

      implicit none

! Loop counters.
      INTEGER :: LC1, LC2, LC3, IER
      INTEGER :: UNMATCHED
      INTEGER, ALLOCATABLE :: iLOCAL_ID(:)

! Max global id.
      INTEGER :: MAX_ID, lSTAT
! Debug flags.
      LOGICAL :: dFlag
      LOGICAL, parameter :: setDBG = .FALSE.

      CALL INIT_ERR_MSG("GLOBAL_TO_LOC_COL")

! Initialize the error flag.
      IER = 0

! Set the local debug flag.
      dFlag = (DMP_LOG .AND. setDBG)

      MAX_ID = maxval(IGLOBAL_ID(1:PIP))
      CALL GLOBAL_ALL_MAX(MAX_ID)

      allocate(iLOCAL_ID(MAX_ID), STAT=lSTAT)
      CALL GLOBAL_ALL_SUM(lSTAT)

! All ranks successfully allocated the array. This permits a crude
! but much faster collision owner detection.
      IF(lSTAT /= 0) THEN
         WRITE(ERR_MSG,1000)
         CALL FLUSH_ERR_MSG(ABORT=.TRUE.)
      ENDIF

 1000 FORMAT('Error 1000: Unable to allocate sufficient memory to ',&
         'generate the',/'map from global to local particle IDs.')

      iLOCAL_ID = 0
      DO LC1=1, PIP
         iLOCAL_ID(iGLOBAL_ID(LC1)) = LC1
      ENDDO

! Store the particle data.
      LC3 = 1
      LC2 = 0
      UNMATCHED = 0

! FIXME Fix Restart
!      LP1: DO LC1 = 1, cIN_COUNT
!          IF(cRestartMap(LC1) == myPE) THEN
!             LC2 = LC2 + 1
!             NEIGHBORS(LC2) = iLOCAL_ID(iPAR_COL(1,LC1))
!             NEIGHBORS(2,LC2) = iLOCAL_ID(iPAR_COL(2,LC1))
! ! Verify that the local indices are valid. If they do not match it is
! ! likely because one of the neighbor was removed via an outlet at the time
! ! the RES file was written but the ghost data wasn't updated.
!             IF(NEIGHBORS(1,LC2) == 0 .OR. NEIGHBORS(2,LC2) == 0) THEN
!                UNMATCHED = UNMATCHED + 1
!                IF(dFLAG) THEN
!                   WRITE(ERR_MSG,1100) iPAR_COL(1,LC1), NEIGHBORS(1,LC2),   &
!                      iPAR_COL(2,LC1), NEIGHBORS(2,LC2)
!                   CALL FLUSH_ERR_MSG(ABORT=.FALSE.)
!                ENDIF
!                DO WHILE(PEA(LC3,1))
!                   LC3 = LC3 + 1
!                ENDDO
!                NEIGHBORS(2,LC2) = LC3
!             ENDIF
!          ENDIF
!       ENDDO LP1

! 1100 FORMAT('Error 1100: Particle neighbor local indices are invalid.',/  &
!         5x,'Global-ID    Local-ID',/' 1:  ',2(3x,I9),/' 2:  ',2(3x,I9))

      CALL GLOBAL_ALL_SUM(UNMATCHED)
      IF(UNMATCHED /= 0) THEN
         WRITE(ERR_MSG,1101) trim(iVal(UNMATCHED))
         CALL FLUSH_ERR_MSG
      ENDIF

 1101 FORMAT(' Warning: 1101: ',A,' particle neighbor datasets were ',&
         'not matched',/' during restart.')

      IF(allocated(iLOCAL_ID)) deallocate(iLOCAL_ID)

      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE GLOBAL_TO_LOC_COL



!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_0I(lNEXT_REC, INPUT_I)

      use mpi_utility, only: BCAST

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: INPUT_I

      IF(bDIST_IO) THEN
         READ(RDES_UNIT, REC=lNEXT_REC) INPUT_I
      ELSE
         IF(myPE == PE_IO) READ(RDES_UNIT, REC=lNEXT_REC) INPUT_I
         CALL BCAST(INPUT_I, PE_IO)
      ENDIF

      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE READ_RES_DES_0I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_1I                                              !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1I(lNEXT_REC, INPUT_I)

      use mpi_utility, only: BCAST
      USE in_binary_512i

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: INPUT_I(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_I)

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) &
            CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
         CALL BCAST(INPUT_I, PE_IO)
      ENDIF


      RETURN
      END SUBROUTINE READ_RES_DES_1I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0D                                          !
!                                                                      !
! Purpose: Write scalar double percision values to RES file.           !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_0D(lNEXT_REC, INPUT_D)

      use mpi_utility, only: BCAST

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: INPUT_D

      IF(bDIST_IO) THEN
         READ(RDES_UNIT, REC=lNEXT_REC) INPUT_D
      ELSE
         IF(myPE == PE_IO) READ(RDES_UNIT, REC=lNEXT_REC) INPUT_D
         CALL BCAST(INPUT_D, PE_IO)
      ENDIF
      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE READ_RES_DES_0D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1D                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1D(lNEXT_REC, INPUT_D)

      use mpi_utility, only: BCAST
      USE in_binary_512

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: INPUT_D(:)

      INTEGER :: lSIZE

      lSIZE = size(INPUT_D)

      IF(bDIST_IO) THEN
         CALL IN_BIN_512(RDES_UNIT, INPUT_D, lSIZE, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) &
            CALL IN_BIN_512(RDES_UNIT, INPUT_D, lSIZE, lNEXT_REC)
         CALL BCAST(INPUT_D, PE_IO)
      ENDIF


      RETURN
      END SUBROUTINE READ_RES_DES_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_0L                                          !
!                                                                      !
! Purpose: Write scalar logical values to RES file.                    !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_0L(lNEXT_REC, OUTPUT_L)

      use mpi_utility, only: BCAST

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: OUTPUT_L

      INTEGER :: OUTPUT_I

      OUTPUT_L = .TRUE.

      IF(bDIST_IO)THEN
         READ(RDES_UNIT, REC=lNEXT_REC) OUTPUT_I
      ELSE
         IF(myPE == PE_IO) READ(RDES_UNIT, REC=lNEXT_REC) OUTPUT_I
         CALL BCAST(OUTPUT_I, PE_IO)
      ENDIF

      IF(OUTPUT_I == 1) OUTPUT_L = .TRUE.
      lNEXT_REC = lNEXT_REC + 1

      RETURN
      END SUBROUTINE READ_RES_DES_0L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1L                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_DES_1L(lNEXT_REC, INPUT_L)

      use mpi_utility, only: BCAST
      USE in_binary_512i

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: INPUT_L(:)

      INTEGER, ALLOCATABLE :: INPUT_I(:)

      INTEGER :: lSIZE, LC1

      lSIZE = size(INPUT_I)
      ALLOCATE( INPUT_I(lSIZE))

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) &
            CALL IN_BIN_512i(RDES_UNIT, INPUT_I, lSIZE, lNEXT_REC)
         CALL BCAST(INPUT_I, PE_IO)
      ENDIF

      DO LC1=1, LSIZE
         IF(INPUT_I(LC1) == 1) THEN
            INPUT_L(LC1) = .TRUE.
         ELSE
            INPUT_L(LC1) = .FALSE.
         ENDIF
      ENDDO

      IF(allocated(INPUT_I)) deallocate(INPUT_I)

      RETURN
      END SUBROUTINE READ_RES_DES_1L

!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1B                                          !
!                                                                      !
! Purpose: Write scalar bytes to RES file.                             !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1B(lNEXT_REC, OUTPUT_B)

      use discretelement, only: PIP

      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf

      use compar, only: numPEs
      USE in_binary_512i

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER(KIND=1), INTENT(OUT) :: OUTPUT_B(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: OUTPUT_I(:)
      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)

      allocate(iPROCBUF(pPROCCNT))
      allocate(iROOTBUF(pROOTCNT))


      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      allocate(output_i(size(output_b)))
      OUTPUT_I(:) = OUTPUT_B(:)

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, OUTPUT_I, pIN_COUNT, lNEXT_REC)
         OUTPUT_B(:) = OUTPUT_I(:)
      ELSE

         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(pIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, pIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, pIN_COUNT
               lPROC = pRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               iRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_I(LC1)
            ENDDO

            deallocate(lBUF_I)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=1)
         DO LC1=1, PIP
            OUTPUT_B(LC1) = iProcBuf(LC1)
         ENDDO

      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)
      deallocate(output_i)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1B

!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1I(lNEXT_REC, OUTPUT_I)

      use discretelement, only: PIP

      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf

      use compar, only: numPEs
      USE in_binary_512i

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: OUTPUT_I(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(iPROCBUF(pPROCCNT))
      allocate(iROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, OUTPUT_I, pIN_COUNT, lNEXT_REC)
      ELSE

         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(pIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, pIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, pIN_COUNT
               lPROC = pRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               iRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_I(LC1)
            ENDDO

            deallocate(lBUF_I)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=1)
         DO LC1=1, PIP
            OUTPUT_I(LC1) = iProcBuf(LC1)
         ENDDO

      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1I



!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_pARRAY_1D                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1D(lNEXT_REC, OUTPUT_D)

      use discretelement, only: PIP
      use desmpi, only: dRootBuf
      use desmpi, only: dProcBuf
      use compar, only: numPEs
      USE in_binary_512

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT_D(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      DOUBLE PRECISION, ALLOCATABLE :: lBUF_D(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(dPROCBUF(pPROCCNT))
      allocate(dROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      IF(bDIST_IO) THEN
         CALL IN_BIN_512(RDES_UNIT, OUTPUT_D, pIN_COUNT, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_D(pIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512(RDES_UNIT, lBUF_D, pIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, pIN_COUNT
               lPROC = pRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               dRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_D(LC1)
            ENDDO

            deallocate(lBUF_D)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=2)
         DO LC1=1, PIP
            OUTPUT_D(LC1) = dProcBuf(LC1)
         ENDDO
      ENDIF

      deallocate(dPROCBUF)
      deallocate(dROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_pARRAY_1L                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_pARRAY_1L(lNEXT_REC, OUTPUT_L)

      use discretelement, only: PIP
      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf
      use compar, only: numPEs
      USE in_binary_512i

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: OUTPUT_L(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)

      allocate(iPROCBUF(pPROCCNT))
      allocate(iROOTBUF(pROOTCNT))

      iDISPLS = pDISPLS
      iScr_RecvCNT = pRECV
      iScatterCNTS = pSCATTER

      IF(bDIST_IO) THEN
         allocate(lBUF_I(pIN_COUNT))
         CALL IN_BIN_512i(RDES_UNIT, lBUF_I, pIN_COUNT, lNEXT_REC)
         DO LC1=1,pIN_COUNT
            IF(lBUF_I(LC1) == 1) THEN
               OUTPUT_L(LC1) = .TRUE.
            ELSE
               OUTPUT_L(LC1) = .FALSE.
            ENDIF
         ENDDO
         deallocate(lBUF_I)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(pIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, pIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, pIN_COUNT
               lPROC = pRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               iRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_I(LC1)
            ENDDO

            deallocate(lBUF_I)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=1)
         DO LC1=1, PIP
            IF(iProcBuf(LC1) == 1) THEN
               OUTPUT_L(LC1) = .TRUE.
            ELSE
               OUTPUT_L(LC1) = .FALSE.
            ENDIF
         ENDDO
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_pARRAY_1L


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_DES_1I                                          !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_cARRAY_1I(lNEXT_REC, OUTPUT_I)

      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf
      use compar, only: numPEs
      use discretelement, only: NEIGH_NUM
      USE in_binary_512i

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      INTEGER, INTENT(OUT) :: OUTPUT_I(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(iPROCBUF(cPROCCNT))
      allocate(iROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iScr_RecvCNT = cRECV
      iScatterCNTS = cSCATTER

      IF(bDIST_IO) THEN
         CALL IN_BIN_512i(RDES_UNIT, OUTPUT_I, cIN_COUNT, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(cIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, cIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, cIN_COUNT
               lPROC = cRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               iRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_I(LC1)
            ENDDO

            deallocate(lBUF_I)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=1)
         DO LC1=1, NEIGH_NUM
            OUTPUT_I(LC1) = iProcBuf(LC1)
         ENDDO
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_cARRAY_1I


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_cARRAY_1D                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_cARRAY_1D(lNEXT_REC, OUTPUT_D)

      use compar, only: numPEs
      use discretelement, only: NEIGH_NUM
      use desmpi, only: dRootBuf
      use desmpi, only: dProcBuf
      USE in_binary_512

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT_D(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      DOUBLE PRECISION, ALLOCATABLE :: lBUF_D(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)


      allocate(dPROCBUF(cPROCCNT))
      allocate(dROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iScr_RecvCNT = cRECV
      iScatterCNTS = cSCATTER


      IF(bDIST_IO) THEN
         CALL IN_BIN_512(RDES_UNIT, OUTPUT_D, cIN_COUNT, lNEXT_REC)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_D(cIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512(RDES_UNIT, lBUF_D, cIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, cIN_COUNT
               lPROC = cRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               dRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_D(LC1)
            ENDDO

            deallocate(lBUF_D)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=2)
         DO LC1=1, NEIGH_NUM
            OUTPUT_D(LC1) = dProcBuf(LC1)
         ENDDO
      ENDIF

      deallocate(dPROCBUF)
      deallocate(dROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_cARRAY_1D


!``````````````````````````````````````````````````````````````````````!
! Subroutine: READ_RES_pARRAY_1L                                       !
!                                                                      !
! Purpose: Write scalar integers to RES file.                          !
!``````````````````````````````````````````````````````````````````````!
      SUBROUTINE READ_RES_cARRAY_1L(lNEXT_REC, OUTPUT_L)

      use compar, only: numPEs
      use discretelement, only: NEIGH_NUM
      use desmpi, only: iRootBuf
      use desmpi, only: iProcBuf
      USE in_binary_512i

      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: lNEXT_REC
      LOGICAL, INTENT(OUT) :: OUTPUT_L(:)

! Loop counters
      INTEGER :: LC1

      INTEGER :: lPROC

      INTEGER, ALLOCATABLE :: lBUF_I(:)
      INTEGER, ALLOCATABLE :: lCOUNT(:)

      allocate(iPROCBUF(cPROCCNT))
      allocate(iROOTBUF(cROOTCNT))

      iDISPLS = cDISPLS
      iScr_RecvCNT = cRECV
      iScatterCNTS = cSCATTER

      IF(bDIST_IO) THEN
         allocate(lBUF_I(cIN_COUNT))
         CALL IN_BIN_512i(RDES_UNIT, lBUF_I, cIN_COUNT, lNEXT_REC)
         DO LC1=1,cIN_COUNT
            IF(lBUF_I(LC1) == 1) THEN
               OUTPUT_L(LC1) = .TRUE.
            ELSE
               OUTPUT_L(LC1) = .FALSE.
            ENDIF
         ENDDO
         deallocate(lBUF_I)
      ELSE
         IF(myPE == PE_IO) THEN
            allocate(lBUF_I(cIN_COUNT))
            allocate(lCOUNT(0:NUMPEs-1))

            CALL IN_BIN_512i(RDES_UNIT, lBUF_I, cIN_COUNT, lNEXT_REC)

            lCOUNT = 0
            DO LC1=1, cIN_COUNT
               lPROC = cRestartMap(LC1)
               lCOUNT(lPROC) = lCOUNT(lPROC) + 1
               iRootBuf(iDispls(lPROC) + lCOUNT(lPROC)) = lBUF_I(LC1)
            ENDDO

            deallocate(lBUF_I)
            deallocate(lCOUNT)
         ENDIF
         CALL DESMPI_SCATTERV(ptype=1)
         DO LC1=1, NEIGH_NUM
            IF(iProcBuf(LC1) == 1) THEN
               OUTPUT_L(LC1) = .TRUE.
            ELSE
               OUTPUT_L(LC1) = .FALSE.
            ENDIF
         ENDDO
      ENDIF

      deallocate(iPROCBUF)
      deallocate(iROOTBUF)

      RETURN
      END SUBROUTINE READ_RES_cARRAY_1L

      END MODULE READ_RES1_DES
