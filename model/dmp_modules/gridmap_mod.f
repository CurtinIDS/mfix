!--------------------------------------------------------------------
! Purpose:
! Contains following subroutines:
!    partition, gridmap_init
!--------------------------------------------------------------------
       module gridmap

!-----------------------------------------------
! Modules
!-----------------------------------------------
        use mpi_utility
        use parallel_mpi
        use geometry
        use sendrecv
        use compar
        use run
        use indices

        use error_manager

        implicit none

        contains

!----------------------------------------------------------------------!
! Purpose: Routine to partition the grid. It works for 1-d, 2-d        !
! decomposition in the current implementation                          !
!----------------------------------------------------------------------!
      SUBROUTINE PARTITION(CYC_XLL, CYC_YLL, CYC_ZLL)

      implicit none

! DUMMY Arguments
!---------------------------------------------------------------------//
! Local flags for cyclic boundarys
       LOGICAL :: CYC_XLL, CYC_YLL, CYC_ZLL

! Local variables
!---------------------------------------------------------------------//
      INTEGER, dimension(0:nodesi-1) :: isize1_all
      INTEGER, dimension(0:nodesj-1) :: jsize1_all
      INTEGER, dimension(0:nodesk-1) :: ksize1_all

      INTEGER :: ip, iproc, isize, iremain
      INTEGER :: jp, jproc, jsize, jremain
      INTEGER :: kp, kproc, ksize, kremain

      INTEGER :: ijkproc

      LOGICAL :: PRESENT

!......................................................................!


! Initialize the error manager.
      CALL INIT_ERR_MSG("PARTITION")

! Set the number of layers for BICGSTAB
      IF(NODESI .NE. 1 .AND. CYC_XLL) nlayers_bicgs = 2
      IF(NODESJ .NE. 1 .AND. CYC_YLL) nlayers_bicgs = 2
      IF(NODESK .NE. 1 .AND. CYC_ZLL) nlayers_bicgs = 2

! Flag that the current setup may not be efficient.
      IF(NODESJ .NE. 1) THEN
         WRITE(ERR_MSG,1000)
         CALL FLUSH_ERR_MSG
      ENDIF

! Get Domain size from ADJUST_IJK_SIZE Subroutine
      IF(DOMAIN_SIZE_ADJUSTED) THEN
         isize1_all = isize_all
         jsize1_all = jsize_all
         ksize1_all = ksize_all
      ELSE
! Determine the size in i direction and add the remainder sequentially
         isize = (imax1-imin1+1)/nodesi
         isize1_all(0:nodesi-1) = isize
         iremain = (imax1-imin1+1) - nodesi*isize
         IF (iremain.ge.1) isize1_all( 0:(iremain-1) ) = isize + 1

! Determine the size in j direction and add the remainder sequentially
         jsize = (jmax1-jmin1+1)/nodesj
         jsize1_all(0:nodesj-1) = jsize
         jremain = (jmax1-jmin1+1) - nodesj*jsize
         IF (jremain.ge.1) jsize1_all( 0:(jremain-1) ) = jsize + 1

! Determine the size in k direction and add the remainder sequentially
         ksize = (kmax1-kmin1+1)/nodesk
         ksize1_all(0:nodesk-1) = ksize
         kremain = (kmax1-kmin1+1) - nodesk*ksize
         IF (kremain.ge.1) ksize1_all( 0:(kremain-1) ) = ksize + 1
      ENDIF


! Get Domain size from gridmap.dat
! This works only in the j-direction   <-------------------------
      IF(.NOT.DOMAIN_SIZE_ADJUSTED) THEN
         INQUIRE(FILE='gridmap.dat',EXIST=PRESENT)
         IF(PRESENT) THEN
            IF(MyPE == PE_IO) THEN
               WRITE(*,*)'Reading gridmap from grimap.dat...'
               OPEN(CONVERT='BIG_ENDIAN',UNIT=777, FILE='gridmap.dat', STATUS='OLD')

               READ (777, *) NODESI,NODESJ,NODESK
               DO IPROC = 0,NODESI-1
                  READ(777,*) jPROC,Isize1_all(IPROC)
               ENDDO
               DO IPROC = 0,NODESJ-1
                  READ(777,*) jPROC,Jsize1_all(IPROC)
               ENDDO
               DO IPROC = 0,NODESK-1
                  READ(777,*) jPROC,Ksize1_all(IPROC)
               ENDDO

               CLOSE(777)
            ENDIF
            CALL BCAST(ISIZE1_ALL)
            CALL BCAST(JSIZE1_ALL)
            CALL BCAST(KSIZE1_ALL)
            allocate( ISIZE_ALL(0:NODESI-1))
            allocate( JSIZE_ALL(0:NODESJ-1))
            allocate( KSIZE_ALL(0:NODESK-1))
            isize_all = isize1_all
            jsize_all = jsize1_all
            ksize_all = ksize1_all
         ENDIF
      ENDIF

! The following is general for 1-d or 2-d or 3-d decompostion
! Determining  istart, jstart and kstart for all the processors
      ijkproc = 0
      kp = kmin1
      do kproc=0,nodesk-1
         jp = jmin1
         do jproc=0,nodesj-1
            ip = imin1
            do iproc=0,nodesi-1

               istart1_all(ijkproc) = ip + sum(isize1_all(0:iproc-1))
               iend1_all(ijkproc) = istart1_all(ijkproc) + isize1_all(iproc)-1

               jstart1_all(ijkproc) = jp + sum(jsize1_all(0:jproc-1))
               jend1_all(ijkproc) = jstart1_all(ijkproc) + jsize1_all(jproc)-1

               kstart1_all(ijkproc) = kp + sum(ksize1_all(0:kproc-1))
               kend1_all(ijkproc) = kstart1_all(ijkproc) + ksize1_all(kproc)-1

               ijkproc = ijkproc+1

            ENDDO
         ENDDO
      ENDDO

      CALL FINL_ERR_MSG

      RETURN

 1000 FORMAT('WARNING 1000: The preconditioner for the linear solver',/&
         'MIGHT NOT be very efficient with DMP partitions in the y-',  &
         'axis.')

      END SUBROUTINE PARTITION


!----------------------------------------------------------------------!
! Purpose: Initializing all the variables from the information         !
! obtained in the above routine.                                       !
!----------------------------------------------------------------------!
        SUBROUTINE GRIDMAP_INIT

        use functions
        use toleranc

        implicit none

! Local variables
!---------------------------------------------------------------------//
! Loop indicies
      integer :: iproc, ii, jj, kk
! Local flags for cyclic boundarys
      LOGICAL :: CYC_XL, CYC_YL, CYC_ZL
! Communicator (MPI_COMM_WORLD)
      INTEGER :: COMM
! Amount of load imbalance
      INTEGER :: IMBALANCE
! Theoritical speedup (based on load imbalance)
      CHARACTER(len=32) :: AMDAHL_SPEEDUP

!......................................................................!

! Initialize the error manager.
      CALL INIT_ERR_MSG("GRIDMAP_INIT")

! Set local flags for cyclic boundaries.
      CYC_XL = (CYCLIC_X .OR. CYCLIC_X_PD)
      CYC_YL = (CYCLIC_Y .OR. CYCLIC_Y_PD)
      CYC_ZL = (CYCLIC_Z .OR. CYCLIC_Z_PD)
      IF(DO_K .AND. COMPARE(ZLENGTH,8.D0*ATAN(ONE))  .AND. &
         (COORDINATES == 'CYLINDRICAL')) CYC_ZL = .TRUE.


      IF(.NOT.ALLOCATED(ijksize3_all))   allocate( ijksize3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(ijkstart3_all))  allocate( ijkstart3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(ijkend3_all))    allocate( ijkend3_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(ijksize4_all))   allocate( ijksize4_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(ijkstart4_all))  allocate( ijkstart4_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(ijkend4_all))    allocate( ijkend4_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(istart_all))     allocate( istart_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jstart_all))     allocate( jstart_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kstart_all))     allocate( kstart_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(istart1_all))    allocate( istart1_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jstart1_all))    allocate( jstart1_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kstart1_all))    allocate( kstart1_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(istart2_all))    allocate( istart2_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jstart2_all))    allocate( jstart2_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kstart2_all))    allocate( kstart2_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(istart3_all))    allocate( istart3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jstart3_all))    allocate( jstart3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kstart3_all))    allocate( kstart3_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(istart4_all))    allocate( istart4_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jstart4_all))    allocate( jstart4_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kstart4_all))    allocate( kstart4_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(iend_all))       allocate( iend_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jend_all))       allocate( jend_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kend_all))       allocate( kend_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(iend1_all))      allocate( iend1_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jend1_all))      allocate( jend1_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kend1_all))      allocate( kend1_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(iend2_all))      allocate( iend2_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jend2_all))      allocate( jend2_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kend2_all))      allocate( kend2_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(iend3_all))      allocate( iend3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jend3_all))      allocate( jend3_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kend3_all))      allocate( kend3_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(iend4_all))      allocate( iend4_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(jend4_all))      allocate( jend4_all(0:numPEs-1) )
      IF(.NOT.ALLOCATED(kend4_all))      allocate( kend4_all(0:numPEs-1) )

      IF(.NOT.ALLOCATED(displs))         allocate( displs(0:numPEs-1) )


      CALL PARTITION(CYC_XL, CYC_YL, CYC_ZL)

! The upper and lower bounds are prescribed such that two ghost
! layers are allowed at the physical boundaries - this is consistent
! with our present approach - need to be generalized if only one ghost
! layer is needed
      do iproc=0,numPEs-1
         istart2_all(iproc) = max(imin1-1,min(imax1+1,istart1_all(iproc)-1))
         if(nodesi.ne.1) then
            istart3_all(iproc) = max(imin1-2,min(imax1+2,istart2_all(iproc)-1))
            istart4_all(iproc) = max(imin1-3,min(imax1+3,istart3_all(iproc)-1))
         else
            istart3_all(iproc) = istart2_all(iproc)
            istart4_all(iproc) = istart3_all(iproc)
         endif

         jstart2_all(iproc) = max(jmin1-1,min(jmax1+1,jstart1_all(iproc)-1))
         if(nodesj.ne.1) then
            jstart3_all(iproc) = max(jmin1-2,min(jmax1+2,jstart2_all(iproc)-1))
            jstart4_all(iproc) = max(jmin1-3,min(jmax1+3,jstart3_all(iproc)-1))
         else
            jstart3_all(iproc) = jstart2_all(iproc)
            jstart4_all(iproc) = jstart3_all(iproc)
         endif

         if(no_k) then
            kstart2_all(iproc) = kstart1_all(iproc)
            kstart3_all(iproc) = kstart1_all(iproc)
            kstart4_all(iproc) = kstart1_all(iproc)
         else
            kstart2_all(iproc) = max(kmin1-1,min(kmax1+1,kstart1_all(iproc)-1))
            if(nodesk.ne.1) then
               kstart3_all(iproc) = max(kmin1-2,min(kmax1+2,kstart2_all(iproc)-1))
               kstart4_all(iproc) = max(kmin1-3,min(kmax1+3,kstart3_all(iproc)-1))
            else
               kstart3_all(iproc) =  kstart2_all(iproc)
               kstart4_all(iproc) =  kstart3_all(iproc)
            endif
         endif

         iend2_all(iproc) = max(imin1-1,min(imax1+1,iend1_all(iproc)+1))
         if(nodesi.ne.1) then
            iend3_all(iproc) = max(imin1-2,min(imax1+2,iend2_all(iproc)+1))
            iend4_all(iproc) = max(imin1-3,min(imax1+3,iend3_all(iproc)+1))
         else
            iend3_all(iproc) = iend2_all(iproc)
            iend4_all(iproc) = iend3_all(iproc)
         endif

         jend2_all(iproc) = max(jmin1-1,min(jmax1+1,jend1_all(iproc)+1))
         if(nodesj.ne.1) then
            jend3_all(iproc) = max(jmin1-2,min(jmax1+2,jend2_all(iproc)+1))
            jend4_all(iproc) = max(jmin1-3,min(jmax1+3,jend3_all(iproc)+1))
         else
            jend3_all(iproc) = jend2_all(iproc)
            jend4_all(iproc) = jend3_all(iproc)
         endif

         if(no_k) then
            kend2_all(iproc) = kend1_all(iproc)
            kend3_all(iproc) = kend1_all(iproc)
            kend4_all(iproc) = kend1_all(iproc)
         else
            kend2_all(iproc) = max(kmin1-1,min(kmax1+1,kend1_all(iproc)+1))
            if(nodesk.ne.1) then
               kend3_all(iproc) = max(kmin1-2,min(kmax1+2,kend2_all(iproc)+1))
               kend4_all(iproc) = max(kmin1-3,min(kmax1+3,kend3_all(iproc)+1))
            else
               kend3_all(iproc) = kend2_all(iproc)
               kend4_all(iproc) = kend3_all(iproc)
            endif
         endif
      enddo

! for higher order methods
      if(.not.fpfoi) then
         do iproc=0,numPEs-1
            istart4_all(iproc) = istart3_all(iproc)
            jstart4_all(iproc) = jstart3_all(iproc)
            kstart4_all(iproc) = kstart3_all(iproc)
            iend4_all(iproc)   = iend3_all(iproc)
            jend4_all(iproc)   = jend3_all(iproc)
            kend4_all(iproc)   = kend3_all(iproc)
         enddo
      endif

      do iproc=0,numPEs-1
         ijkstart3_all(iproc) = 1
         ijkend3_all(iproc) =  1 + (iend3_all(iproc) - istart3_all(iproc)) &
           + (jend3_all(iproc)-jstart3_all(iproc))*(iend3_all(iproc)-istart3_all(iproc)+1) &
           + (kend3_all(iproc)-kstart3_all(iproc))*(jend3_all(iproc)-jstart3_all(iproc)+1)* &
             (iend3_all(iproc)-istart3_all(iproc)+1)

         ijkstart4_all(iproc) = 1
         ijkend4_all(iproc) =  1 + (iend4_all(iproc) - istart4_all(iproc)) &
            + (jend4_all(iproc)-jstart4_all(iproc))*(iend4_all(iproc)-istart4_all(iproc)+1) &
            + (kend4_all(iproc)-kstart4_all(iproc))*(jend4_all(iproc)-jstart4_all(iproc)+1)* &
              (iend4_all(iproc)-istart4_all(iproc)+1)
      enddo

      do iproc=0,numPEs-1
         ijksize3_all(iproc) = ijkend3_all(iproc) - ijkstart3_all(iproc) + 1
         ijksize4_all(iproc) = ijkend4_all(iproc) - ijkstart4_all(iproc) + 1
      enddo

      displs(0) = 0
      do iproc=1,numPEs-1
         displs(iproc) = displs(iproc-1)+ijksize3_all(iproc-1)
!       write(*,*) 'displ',displs(iproc),iproc, ijksize3_all(iproc)
      enddo


      ijkstart3 = ijkstart3_all(myPE)
      ijkend3   = ijkend3_all(myPE)
      ijksize3  = ijksize3_all(myPE)

      ijkstart4 = ijkstart4_all(myPE)
      ijkend4   = ijkend4_all(myPE)
      ijksize4  = ijksize4_all(myPE)

      istart1   = istart1_all(myPE)
      iend1     = iend1_all(myPE)
      jstart1   = jstart1_all(myPE)
      jend1     = jend1_all(myPE)
      kstart1   = kstart1_all(myPE)
      kend1     = kend1_all(myPE)

      istart2   = istart2_all(myPE)
      iend2     = iend2_all(myPE)
      jstart2   = jstart2_all(myPE)
      jend2     = jend2_all(myPE)
      kstart2   = kstart2_all(myPE)
      kend2     = kend2_all(myPE)

      istart3   = istart3_all(myPE)
      iend3     = iend3_all(myPE)
      jstart3   = jstart3_all(myPE)
      jend3     = jend3_all(myPE)
      kstart3   = kstart3_all(myPE)
      kend3     = kend3_all(myPE)

      istart4   = istart4_all(myPE)
      iend4     = iend4_all(myPE)
      jstart4   = jstart4_all(myPE)
      jend4     = jend4_all(myPE)
      kstart4   = kstart4_all(myPE)
      kend4     = kend4_all(myPE)

      IF(.not.allocated(NCPP_UNIFORM)) allocate( NCPP_UNIFORM(0:NumPEs-1))

      IF(.NOT.NCPP_UNIFORM_BACKED_UP) THEN
         NCPP_UNIFORM(MyPE) = ijksize3_all(MyPE)
      ENDIF
      NCPP_UNIFORM_BACKED_UP = .TRUE.

      IF(SHORT_GRIDMAP_INIT) THEN
!        do iproc=0,numPEs-1
!           NCPP_UNIFORM(iproc) = ijksize3_all(iproc)
!        enddo
         RETURN
      ENDIF

! Setup mapping to take care of cyclic boundary conditions
! ---------------------------------------------------------------->>>
! consider cyclic boundary condition using the imap(:),jmap(:),kmap(:)
! indirection arrays

      allocate( imap( imin4:imax4 ) )
      allocate( jmap( jmin4:jmax4 ) )
      allocate( kmap( kmin4:kmax4 ) )

      allocate( imap_c( imin4:imax4 ) )
      allocate( jmap_c( jmin4:jmax4 ) )
      allocate( kmap_c( kmin4:kmax4 ) )

      do kk=kmin4,kmax4
        kmap(kk) = kk
      enddo

      do jj=jmin4,jmax4
        jmap(jj) = jj
      enddo

      do ii=imin4,imax4
        imap(ii) = ii
      enddo

      if (CYC_ZL) then
         kmap( kmax2 ) = kmin1
         kmap( kmin2 ) = kmax1
         if (kmax3.gt.kmax2) kmap(kmax3) = kmap(kmax2)+1
         if (kmin3.lt.kmin2) kmap(kmin3) = kmap(kmin2)-1
         if (kmax4.gt.kmax3) kmap(kmax4) = kmap(kmax3)+1
         if (kmin4.lt.kmin3) kmap(kmin4) = kmap(kmin3)-1
      endif

      if (CYC_YL) then
         jmap( jmax2 ) = jmin1
         jmap( jmin2 ) = jmax1
         if (jmax3.gt.jmax2) jmap(jmax3) = jmap(jmax2)+1
         if (jmin3.lt.jmin2) jmap(jmin3) = jmap(jmin2)-1
         if (jmax4.gt.jmax3) jmap(jmax4) = jmap(jmax3)+1
         if (jmin4.lt.jmin3) jmap(jmin4) = jmap(jmin3)-1
      endif

      if (CYC_XL) then
         imap( imax2 ) = imin1
         imap( imin2 ) = imax1
         if (imax3.gt.imax2) imap(imax3) = imap(imax2)+1
         if (imin3.lt.imin2) imap(imin3) = imap(imin2)-1
         if (imax4.gt.imax3) imap(imax4) = imap(imax3)+1
         if (imin4.lt.imin3) imap(imin4) = imap(imin3)-1
      endif

      do kk=kmin4,kmax4
        kmap_c(kk) = kk
      enddo

      do jj=jmin4,jmax4
        jmap_c(jj) = jj
      enddo

      do ii=imin4,imax4
        imap_c(ii) = ii
      enddo

      if (CYC_ZL.and.nodesk.eq.1) then
         kmap_c( kmax2 ) = kmin1
         kmap_c( kmin2 ) = kmax1
         if (kmax3.gt.kmax2) kmap_c(kmax3) = kmap_c(kmax2)+1
         if (kmin3.lt.kmin2) kmap_c(kmin3) = kmap_c(kmin2)-1
         if (kmax4.gt.kmax3) kmap_c(kmax4) = kmap_c(kmax3)+1
         if (kmin4.lt.kmin3) kmap_c(kmin4) = kmap_c(kmin3)-1
      endif

      if (CYC_YL.and.nodesj.eq.1) then
         jmap_c( jmax2 ) = jmin1
         jmap_c( jmin2 ) = jmax1
         if (jmax3.gt.jmax2) jmap_c(jmax3) = jmap_c(jmax2)+1
         if (jmin3.lt.jmin2) jmap_c(jmin3) = jmap_c(jmin2)-1
         if (jmax4.gt.jmax3) jmap_c(jmax4) = jmap_c(jmax3)+1
         if (jmin4.lt.jmin3) jmap_c(jmin4) = jmap_c(jmin3)-1
      endif

      if (CYC_XL.and.nodesi.eq.1) then
         imap_c( imax2 ) = imin1
         imap_c( imin2 ) = imax1
         if (imax3.gt.imax2) imap_c(imax3) = imap_c(imax2)+1
         if (imin3.lt.imin2) imap_c(imin3) = imap_c(imin2)-1
         if (imax4.gt.imax3) imap_c(imax4) = imap_c(imax3)+1
         if (imin4.lt.imin3) imap_c(imin4) = imap_c(imin3)-1
      endif
! End setup mapping to take care of cyclic boundary conditions
! ----------------------------------------------------------------<<<


! Defining new set of varaibles to define upper and lower bound of the
! indices to include actual physical boundaries of the problem
      do iproc = 0, numPEs-1
         istart = istart1
         iend = iend1
         jstart = jstart1
         jend = jend1
         kstart = kstart1
         kend = kend1

         if(istart3.eq.imin3) istart = istart2
         if(iend3.eq.imax3) iend = iend2
         if(jstart3.eq.jmin3) jstart = jstart2
         if(jend3.eq.jmax3) jend = jend2
         if(kstart3.eq.kmin3) kstart = kstart2
         if(kend3.eq.kmax3) kend = kend2

         istart_all(iproc) = istart
         iend_all(iproc)   = iend
         jstart_all(iproc) = jstart
         jend_all(iproc)   = jend
         kstart_all(iproc) = kstart
         kend_all(iproc)   = kend
      enddo

      IF(numPEs .GT. 1) THEN
! Calculate any load imbalance.
         IMBALANCE = DBLE(maxval(ijksize3_all)- minval(ijksize3_all)) /    &
            minval(ijksize3_all)*100.0
! Calculate potential speedup based on Amdahl's Law
         IF(IMBALANCE == 0)THEN
            AMDAHL_SPEEDUP='+Inf'
         ELSE
            AMDAHL_SPEEDUP=''
            WRITE(AMDAHL_SPEEDUP,*)1.0/dble(IMBALANCE)
         ENDIF
! Construct a message for the user telling them the grid partition info.
         WRITE(ERR_MSG,1000)maxval(ijksize3_all), maxloc(ijksize3_all),&
            minval(ijksize3_all), minloc(ijksize3_all),                &
            sum(ijksize3_all)/numPEs, trim(AMDAHL_SPEEDUP)
         CALL FLUSH_ERR_MSG
      ENDIF

! Setup coefficients of FUINIJK
        c0 = 1 - jstart3_all(myPE)
        c1 = (jend3_all(myPE)-jstart3_all(myPE)+1)
        c2 = (jend3_all(myPE)-jstart3_all(myPE)+1)* (iend3_all(myPE)-istart3_all(myPE)+1)
        c0 =  c0  - c1*istart3_all(myPE) - c2*kstart3_all(myPE)

! Setup coefficients of FUINIJK3
        c0_3 = 1 - jstart4_all(myPE)
        c1_3 = (jend4_all(myPE)-jstart4_all(myPE)+1)
        c2_3 = (jend4_all(myPE)-jstart4_all(myPE)+1)* (iend4_all(myPE)-istart4_all(myPE)+1)
        c0_3 =  c0_3  - c1_3*istart4_all(myPE) - c2_3*kstart4_all(myPE)

!   Initialize Array mapping (I,J,K) to IJK
        INCREMENT_ARRAYS_ALLOCATED = .FALSE.

! These arrays could already be allocated in post_mfix
! when interpolating old data to new grid

        if(allocated(IJK_ARRAY_OF)) deallocate(IJK_ARRAY_OF)
        if(allocated(FUNIJK_MAP_C)) deallocate(FUNIJK_MAP_C)
        if(allocated(DEAD_CELL_AT)) deallocate(DEAD_CELL_AT)

        ! Must extend range such that neighbors (IM,JP etc...) stay in bound
        allocate(IJK_ARRAY_OF(istart3-1:iend3+1,jstart3-1:jend3+1,kstart3-1:kend3+1))
        allocate(FUNIJK_MAP_C(istart3-1:iend3+1,jstart3-1:jend3+1,kstart3-1:kend3+1))
        allocate(DEAD_CELL_AT(imin3-1:imax3+1,jmin3-1:jmax3+1,kmin3-1:kmax3+1))

        DEAD_CELL_AT = .FALSE.

! Save IJK value of (I,J,K) cell in an array
! IJK_ARRAY_OF(I,J,K) will replace the use of FUNIJK(I,J,K)
        DO ii = istart3,iend3
           DO jj = jstart3,jend3
              DO kk = kstart3,kend3
                 IJK_ARRAY_OF(ii,jj,kk)=FUNIJK_0(ii,jj,kk)
              ENDDO
           ENDDO
        ENDDO

        DO ii = istart3,iend3
           DO jj = jstart3,jend3
              DO kk = kstart3,kend3
                 FUNIJK_MAP_C(ii,jj,kk)=IJK_ARRAY_OF(IMAP_C(ii),JMAP_C(jj),KMAP_C(kk))
              ENDDO
           ENDDO
        ENDDO


! Call to sendrecv_init to set all the communication pattern
#ifdef MPI
      COMM = MPI_COMM_WORLD
      CALL SENDRECV_INIT(COMM, CYC_XL, CYC_YL, CYC_ZL, idebug=0)
#endif

      CALL FINL_ERR_MSG
      RETURN

 1000 FORMAT('Parallel load balancing statistics:',2/,13x,'Comp. cells',&
      4X,'Processor',/3X,'maximum   ',I11,4X,I9,/3X,'minimum   ',I11,&
      4X,I9,/3X,'average   ',I11,6X,'-N/A-',2/,3X,'Maximum speedup ',&
      '(Amdahls Law) = ',A)

      END SUBROUTINE GRIDMAP_INIT

      END MODULE GRIDMAP
