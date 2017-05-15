!----------------------------------------------------------------------!
!  Module: MPI_PACK_DES                                                !
!  Author: Pradeep Gopalakrishnan, J.Musser                            !
!                                                                      !
!  Purpose: Contains routines for passing data interpolted to ghost    !
!     cells to the owner cell via summation.                           !
!----------------------------------------------------------------------!
      module sendrecvnode

      PRIVATE

      PUBLIC :: INIT_DES_COLLECT_gDATA, DES_COLLECT_gDATA
      PUBLIC :: DES_SETNODEINDICES, DES_EXCHANGENODE

      interface DES_COLLECT_gDATA
        module procedure DES_COLLECT_gDATA_db1
        module procedure DES_COLLECT_gDATA_db2
      end interface DES_COLLECT_gDATA

      integer :: itotalneigh, itotalindx

      integer, allocatable :: itoproc(:)
      integer, allocatable :: istartsend(:)
      integer, allocatable :: istartrecv(:)

! Following variables are used to exchange grid index values when
! des_interp_on is true
      integer, allocatable :: isendnodes(:)
      integer, allocatable :: irecvnodes(:)

      double precision, allocatable :: dsendnodebuf(:)
      double precision, allocatable :: drecvnodebuf(:)

      integer, allocatable :: irecvreqnode(:)
      integer, allocatable :: isendreqnode(:)

      contains

!----------------------------------------------------------------------!
!  Subroutine: INIT_DES_COLLECT_gDATA                                  !
!  Author: J.Musser                                                    !
!                                                                      !
! Purpose: Setup the send/recv schedules for summing ghost cell data   !
!    back into the owner cell for DES interpolation.                   !
!----------------------------------------------------------------------!
      subroutine INIT_DES_COLLECT_gDATA

      use compar, only: istart1, istart2, iend1, iend2
      use compar, only: jstart1, jstart2, jend1, jend2
      use compar, only: kstart1, kstart2, kend1, kend2
      use compar, only: mype, nodesi, nodesj, nodesk, dead_cell_at
      use desgrid, only: IofPROC, JofPROC, KofPROC
      use desgrid, only: procIJK
      use discretelement, only: des_periodic_walls_x, des_periodic_walls_y, des_periodic_walls_z
      use functions, only: funijk, wall_at

      implicit none

! Local variables
!-----------------------------------------------
      integer :: lijkproc,liproc,ljproc,lkproc
      integer :: li,lj,lk
      integer :: li2,lj2,lk2

      integer :: liproc_start, liproc_end
      integer :: ljproc_start, ljproc_end
      integer :: lkproc_start, lkproc_end

      integer :: lci,lcj,lck,lproc,lcount
      integer :: linode_start,linode_end, linode
      integer :: ljnode_start,ljnode_end, ljnode
      integer :: lknode_start,lknode_end, lknode
      logical :: lpresent

      integer, allocatable :: iprocsumindx(:)

!-----------------------------------------------

! set flags for interprocessor boundaries
      liproc = iofproc(mype)
      ljproc = jofproc(mype)
      lkproc = kofproc(mype)

! if not periodic then limit the processor
      if(des_periodic_walls_x .and. nodesi > 1) then
         liproc_start=liproc-1
         liproc_end=liproc+1
      else
         liproc_start =max(liproc-1,0)
         liproc_end=min(liproc+1,nodesi-1)
      end if

      if(des_periodic_walls_y .and. nodesj > 1) then
         ljproc_start=ljproc-1
         ljproc_end=ljproc+1
      else
         ljproc_start =max(ljproc-1,0)
         ljproc_end=min(ljproc+1,nodesj-1)
      end if

      if(des_periodic_walls_z .and. nodesk > 1) then
         lkproc_start=lkproc-1
         lkproc_end=lkproc+1
      else
         lkproc_start =max(lkproc-1,0)
         lkproc_end=min(lkproc+1,nodesk-1)
      end if

      itotalneigh = (liproc_end-liproc_start+1)*&
         (ljproc_end-ljproc_start+1)*(lkproc_end-lkproc_start+1)-1

! allocate the variables
      allocate (itoproc(itotalneigh))
      allocate (iprocsumindx(itotalneigh))
      allocate (istartsend(itotalneigh+1))
      allocate (istartrecv(itotalneigh+1))
      allocate (isendreqnode(itotalneigh))
      allocate (irecvreqnode(itotalneigh))

! First loop to count the total index for each processor and count the
! neighbour processor
      itotalneigh = 0
      itoproc(:)=-1
      iprocsumindx(:) =0

      do lk = lkproc_start, lkproc_end
      do lj = ljproc_start, ljproc_end
      do li = liproc_start, liproc_end

         li2 = mod(li,nodesi); if(li2 < 0) li2 = nodesi-1
         lj2 = mod(lj,nodesj); if(lj2 < 0) lj2 = nodesj-1
         lk2 = mod(lk,nodesk); if(lk2 < 0) lk2 = nodesk-1

         lijkproc = procijk(li2,lj2,lk2)

         if (lijkproc.eq.mype) cycle

! check if the processor exits in the previous list
         lpresent = .false.
         do lproc = 1,itotalneigh
            if (lijkproc .eq.itoproc(lproc)) then
               lpresent = .true.
               exit
            end if
         end do
         if(.not.lpresent) then
            itotalneigh = itotalneigh + 1
            lproc = itotalneigh
         end if

         itoproc(lproc) = lijkproc

         lci=(liproc-li)
         if(lci == 1) then
            linode_start = iStart2
            linode_end = iStart2
         elseif(lci == -1) then
            linode_start = iEnd2
            linode_end = iEnd2
         else
            linode_start = iStart1
            linode_end = iEnd1
         endif

         lcj=(ljproc-lj)
         if(lcj == 1) then
            ljnode_start = jStart2
            ljnode_end = jStart2
         elseif(lcj == -1) then
            ljnode_start = jEnd2
            ljnode_end = jEnd2
         else
            ljnode_start = jStart1
            ljnode_end=jEnd1
         endif

         lck=(lkproc-lk)
         if(lck == 1) then
            lknode_start = kStart2
            lknode_end = kStart2
         elseif(lck == -1) then
            lknode_start = kEnd2
            lknode_end = kEnd2
         else
            lknode_start = kStart1
            lknode_end=kEnd1
         endif

         do lknode = lknode_start,lknode_end
         do linode = linode_start,linode_end
         do ljnode = ljnode_start,ljnode_end
            IF(DEAD_CELL_AT(linode,ljnode,lknode)) CYCLE
            IF(WALL_AT(FUNIJK(linode,ljnode,lknode))) CYCLE
            iprocsumindx(lproc) = iprocsumindx(lproc) + 1
         end do
         end do
         end do
      end do
      end do
      end do


!assign the start index
      do lproc =1,itotalneigh+1
         istartsend(lproc)=sum(iprocsumindx(1:lproc-1))+1
      end do
      itotalindx=istartsend(itotalneigh+1)-1

! allocate the variables
      allocate(isendnodes(itotalindx))
      allocate(dsendnodebuf(itotalindx))

! second loop to assign actual index for send map
      iprocsumindx(:)=0
      do lk = lkproc_start,lkproc_end
      do lj = ljproc_start,ljproc_end
      do li = liproc_start,liproc_end
         li2 = mod(li,nodesi);if(li2.lt.0)li2=nodesi-1
         lj2 = mod(lj,nodesj);if(lj2.lt.0)lj2=nodesj-1
         lk2 = mod(lk,nodesk);if(lk2.lt.0)lk2=nodesk-1
         lijkproc = procijk(li2,lj2,lk2)
         if (lijkproc.eq.mype) cycle
! find the index of the neighbour
         do lproc =1,itotalneigh
            if(lijkproc.eq.itoproc(lproc)) then
               exit
            end if
         end do

         lci=(liproc-li)
         if(lci == 1) then
            linode_start = iStart2
            linode_end = iStart2
         elseif(lci == -1) then
            linode_start = iEnd2
            linode_end = iEnd2
         else
            linode_start = iStart1
            linode_end = iEnd1
         endif

         lcj=(ljproc-lj)
         if(lcj == 1) then
            ljnode_start = jStart2
            ljnode_end = jStart2
         elseif(lcj == -1) then
            ljnode_start = jEnd2
            ljnode_end = jEnd2
         else
            ljnode_start = jStart1
            ljnode_end=jEnd1
         endif

         lck=(lkproc-lk)
         if(lck == 1) then
            lknode_start = kStart2
            lknode_end = kStart2
         elseif(lck == -1) then
            lknode_start = kEnd2
            lknode_end = kEnd2
         else
            lknode_start = kStart1
            lknode_end=kEnd1
         endif

         lcount = istartsend(lproc)+iprocsumindx(lproc)
         do lknode = lknode_start,lknode_end
         do linode = linode_start,linode_end
         do ljnode = ljnode_start,ljnode_end
            IF(DEAD_CELL_AT(linode,ljnode,lknode)) CYCLE
            IF(WALL_AT(FUNIJK(linode,ljnode,lknode))) CYCLE
            isendnodes(lcount)=funijk(linode,ljnode,lknode)
            iprocsumindx(lproc)=iprocsumindx(lproc)+1
            lcount = lcount+1
         end do
         end do
         end do

      end do
      end do
      end do

! Build the recv schedule

      iprocsumindx(:) =0
      do lk = lkproc_start, lkproc_end
      do lj = ljproc_start, ljproc_end
      do li = liproc_start, liproc_end

         li2 = mod(li,nodesi); if(li2 < 0) li2 = nodesi-1
         lj2 = mod(lj,nodesj); if(lj2 < 0) lj2 = nodesj-1
         lk2 = mod(lk,nodesk); if(lk2 < 0) lk2 = nodesk-1

         lijkproc = procijk(li2,lj2,lk2)

         if (lijkproc.eq.mype) cycle

! check if the processor exits in the previous list
         do lproc = 1,itotalneigh
            if(lijkproc .eq. itoproc(lproc)) exit
         end do

         lci=(liproc-li);lcj=(ljproc-lj);lck=(lkproc-lk)

         linode_start = istart1; linode_end=iend1
         ljnode_start = jstart1; ljnode_end=jend1
         lknode_start = kstart1; lknode_end=kend1
         if(lci.eq. 1) linode_end = iStart1
         if(lci.eq.-1) linode_start = iEnd1
         if(lcj.eq. 1) ljnode_end = jStart1
         if(lcj.eq.-1) ljnode_start = jEnd1
         if(lck.eq. 1) lknode_end = kStart1
         if(lck.eq.-1) lknode_start = kEnd1

         do lknode = lknode_start,lknode_end
         do linode = linode_start,linode_end
         do ljnode = ljnode_start,ljnode_end
            IF(DEAD_CELL_AT(linode,ljnode,lknode)) CYCLE
            IF(WALL_AT(FUNIJK(linode,ljnode,lknode))) CYCLE
            iprocsumindx(lproc) = iprocsumindx(lproc) + 1
         end do
         end do
         end do
      end do
      end do
      end do

!assign the start index
      do lproc =1,itotalneigh+1
         istartrecv(lproc)=sum(iprocsumindx(1:lproc-1))+1
      end do
      itotalindx=istartrecv(itotalneigh+1)-1

      allocate(irecvnodes(itotalindx))
      allocate(drecvnodebuf(itotalindx))

! second loop to assign actual index
      iprocsumindx(:)=0
      do lk = lkproc_start,lkproc_end
      do lj = ljproc_start,ljproc_end
      do li = liproc_start,liproc_end

         li2 = mod(li,nodesi);if(li2.lt.0)li2=nodesi-1
         lj2 = mod(lj,nodesj);if(lj2.lt.0)lj2=nodesj-1
         lk2 = mod(lk,nodesk);if(lk2.lt.0)lk2=nodesk-1
         lijkproc = procijk(li2,lj2,lk2)

         if (lijkproc.eq.mype) cycle

! find the index of the neighbour
         do lproc =1,itotalneigh
            if(lijkproc.eq.itoproc(lproc)) then
               exit
            end if
         end do

         lci=(liproc-li);lcj=(ljproc-lj);lck=(lkproc-lk)

! Set up the receive map
         linode_start = istart1; linode_end=iend1
         ljnode_start = jstart1; ljnode_end=jend1
         lknode_start = kstart1; lknode_end=kend1

         if(lci.eq. 1) linode_end = iStart1
         if(lci.eq.-1) linode_start = iEnd1
         if(lcj.eq. 1) ljnode_end = jStart1
         if(lcj.eq.-1) ljnode_start = jEnd1
         if(lck.eq. 1) lknode_end = kStart1
         if(lck.eq.-1) lknode_start = kEnd1

         lcount = istartrecv(lproc)+iprocsumindx(lproc)
         do lknode = lknode_start,lknode_end
         do linode = linode_start,linode_end
         do ljnode = ljnode_start,ljnode_end
            IF(DEAD_CELL_AT(linode,ljnode,lknode)) CYCLE
            IF(WALL_AT(FUNIJK(linode,ljnode,lknode))) CYCLE
            irecvnodes(lcount)=funijk(linode,ljnode,lknode)
            iprocsumindx(lproc)=iprocsumindx(lproc)+1
            lcount = lcount+1
         end do
         end do
         end do

      end do
      end do
      end do

!      call des_dbgnodesr()

      RETURN
      end subroutine INIT_DES_COLLECT_gDATA


!----------------------------------------------------------------------!
!  Subroutine: DES_COLLECT_gDATA_db1                                   !
!  Author: J.Musser                                                    !
!                                                                      !
! Purpose: Conduct a reverse halo type exchange were data that was     !
!    stored in ghost cells is summed into the 'real' cell. This is     !
!    needed for DES interpolation where data is interpolted into ghost !
!    cells.                                                            !
!----------------------------------------------------------------------!
      subroutine DES_COLLECT_gDATA_db1(pvar)

      use desmpi_wrapper, only: DES_MPI_WAIT
      use desmpi_wrapper, only: DES_MPI_iSEND
      use desmpi_wrapper, only: DES_MPI_iRECV
      use parallel_mpi, only: MPI_Check
      use compar, only: myPE, numPEs
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      double precision, intent(inout) :: pvar(:)
!-----------------------------------------------
! local variables
!-----------------------------------------------
      character(len=80), parameter :: name = 'des_exchangenode'
      integer :: lindx,lcount,lcount2,lneigh,ltag,lerr
      integer :: lstart,lend,ltotal
!-----------------------------------------------

! steps pack the buffer call isend and irecv
      do lcount = 1,itotalneigh
         lneigh = itoproc(lcount)
         lstart = istartsend(lcount);lend=istartsend(lcount+1)-1
         do lcount2 = lstart,lend
            dsendnodebuf(lcount2) = pvar(isendnodes(lcount2))
         end do

         ltag = message_tag(lneigh,mype)
         lstart = istartrecv(lcount);lend=istartrecv(lcount+1)-1
         ltotal = lend-lstart+1
         call des_mpi_irecv(drecvnodebuf(lstart:lend),ltotal, &
                            lneigh,ltag,irecvreqnode(lcount),lerr)
         call mpi_check( name //':mpi_irecv ', lerr )

         ltag = message_tag(mype,lneigh)
         lstart = istartsend(lcount);lend=istartsend(lcount+1)-1
         ltotal = lend-lstart+1
         call des_mpi_isend(dsendnodebuf(lstart:lend),ltotal, &
                            lneigh,ltag,isendreqnode(lcount),lerr)
         call mpi_check( name //':mpi_irecv ', lerr )
      end do
! call mpi wait to complete the exchange
      do lcount = 1,itotalneigh
         call des_mpi_wait(isendreqnode(lcount),lerr)
         call mpi_check( name //':mpi_wait-send', lerr )
         call des_mpi_wait(irecvreqnode(lcount),lerr)
         call mpi_check( name //':mpi_wait-recv', lerr )
      end do

! after receiving the buffer the values are either added or
! replaced based on the flag
      do lcount = 1,itotalindx
         lindx = irecvnodes(lcount)
         pvar(lindx) = pvar(lindx) + drecvnodebuf(lcount)
      end do

      return

      contains

      integer function message_tag(lsource,ldest)
         implicit none
         integer, intent(in) :: lsource,ldest
         message_tag = lsource+numpes*ldest+200
      end function message_tag

      end subroutine DES_COLLECT_gDATA_db1


!----------------------------------------------------------------------!
!  Subroutine: DES_COLLECT_gDATA_db2                                   !
!  Author: J.Musser                                                    !
!                                                                      !
! Purpose: Wrapper for 2D arrays. See DES_COLLECT_gDATA_db1.           !
!----------------------------------------------------------------------!
      subroutine DES_COLLECT_gDATA_db2(pvar)

      implicit none

! dummy arguments 
!-----------------------------------------------
      double precision, intent(inout) :: pvar(:,:)
      integer :: lc

      do lc=lbound(pVAR,2), ubound(pVAR,2)
         call des_collect_gDATA_db1(pVAR(:,lc))
      enddo
      return
      end subroutine des_collect_gdata_db2



!!############################################################################################
!!############################################################################################
!!############################################################################################
!!############################################################################################

!------------------------------------------------------------------------
! Subroutine       : des_setnodesindices
! Purpose          : allocates and initializes the variables related
!                    to send and recv for grid node values
!------------------------------------------------------------------------
      subroutine des_setnodeindices
      use compar, only: mype, nodesi, nodesj, nodesk, dead_cell_at
      use desgrid, only: IofPROC, JofPROC, KofPROC
      use desgrid, only: procIJK
      use discretelement, only: des_periodic_walls_x, des_periodic_walls_y, des_periodic_walls_z
      use compar, only: istart2, iend1
      use compar, only: jstart2, jend1
      use compar, only: kstart2, kend1
      use functions, only: funijk
!-----------------------------------------------
      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer :: lijkproc,liproc,ljproc,lkproc
      integer :: li,lj,lk
      integer :: li2,lj2,lk2
      integer :: liproc_start,liproc_end
      integer :: ljproc_start,ljproc_end
      integer :: lkproc_start,lkproc_end
      integer :: lci,lcj,lck,lproc,lcount
      integer :: linode_start,linode_end, linode
      integer :: ljnode_start,ljnode_end, ljnode
      integer :: lknode_start,lknode_end, lknode
      logical :: lpresent
      integer, allocatable :: iprocsumindx(:)
!-----------------------------------------------

! set flags for interprocessor boundaries and set the corresponding to proc
      liproc = iofproc(mype)
      ljproc = jofproc(mype)
      lkproc = kofproc(mype)

! if not periodic then limit the processor
      if(des_periodic_walls_x.and.nodesi.gt.1) then
         liproc_start=liproc-1
         liproc_end=liproc+1
      else
         liproc_start =max(liproc-1,0)
         liproc_end=min(liproc+1,nodesi-1)
      end if
      if(des_periodic_walls_y.and.nodesj.gt.1) then
         ljproc_start=ljproc-1
         ljproc_end=ljproc+1
      else
         ljproc_start =max(ljproc-1,0)
         ljproc_end=min(ljproc+1,nodesj-1)
      end if
      if(des_periodic_walls_z.and.nodesk.gt.1) then
         lkproc_start=lkproc-1
         lkproc_end=lkproc+1
      else
         lkproc_start =max(lkproc-1,0)
         lkproc_end=min(lkproc+1,nodesk-1)
      end if
      itotalneigh = (liproc_end-liproc_start+1)*&
         (ljproc_end-ljproc_start+1)*(lkproc_end-lkproc_start+1)-1

! allocate the variables
      allocate (itoproc(itotalneigh))
      allocate (iprocsumindx(itotalneigh))
      allocate (istartsend(itotalneigh+1))
      allocate (irecvreqnode(itotalneigh))
      allocate (isendreqnode(itotalneigh))

! First loop to count the total index for each processor and count the
! neighbour processor
      itotalneigh = 0
      itoproc(:)=-1
      iprocsumindx(:) =0
      do lk = lkproc_start,lkproc_end
      do lj = ljproc_start,ljproc_end
      do li = liproc_start,liproc_end
         li2 = mod(li,nodesi);if(li2.lt.0)li2=nodesi-1
         lj2 = mod(lj,nodesj);if(lj2.lt.0)lj2=nodesj-1
         lk2 = mod(lk,nodesk);if(lk2.lt.0)lk2=nodesk-1
         lijkproc = procijk(li2,lj2,lk2)
         if (lijkproc.eq.mype) cycle
! check if the processor exits in the previous list
         lpresent = .false.
         do lproc = 1,itotalneigh
            if (lijkproc .eq.itoproc(lproc)) then
               lpresent = .true.
               exit
            end if
         end do
         if(.not.lpresent) then
            itotalneigh = itotalneigh + 1
            lproc = itotalneigh
         end if
         itoproc(lproc) = lijkproc
         lci=(liproc-li);lcj=(ljproc-lj);lck=(lkproc-lk)
         linode_start = istart2; linode_end=iend1
         ljnode_start = jstart2; ljnode_end=jend1
         lknode_start = kstart2; lknode_end=kend1
         if(lci.eq.1) linode_end = istart2
         if(lci.eq.-1)  linode_start = iend1
         if(lcj.eq.1) ljnode_end = jstart2
         if(lcj.eq.-1)  ljnode_start = jend1
         if(lck.eq.1) lknode_end = kstart2
         if(lck.eq.-1)  lknode_start = kend1
         do lknode = lknode_start,lknode_end
         do linode = linode_start,linode_end
         do ljnode = ljnode_start,ljnode_end
            IF(DEAD_CELL_AT(linode,ljnode,lknode)) CYCLE
            iprocsumindx(lproc) = iprocsumindx(lproc) + 1
         end do
         end do
         end do
      end do
      end do
      end do
!assign the start index
      do lproc =1,itotalneigh+1
         istartsend(lproc)=sum(iprocsumindx(1:lproc-1))+1
      end do
      itotalindx=istartsend(itotalneigh+1)-1

! allocate the variables
      allocate (isendnodes(itotalindx))
      allocate (dsendnodebuf(itotalindx))
      allocate (drecvnodebuf(itotalindx))

! second loop to assign actual index
      iprocsumindx(:)=0
      do lk = lkproc_start,lkproc_end
      do lj = ljproc_start,ljproc_end
      do li = liproc_start,liproc_end
         li2 = mod(li,nodesi);if(li2.lt.0)li2=nodesi-1
         lj2 = mod(lj,nodesj);if(lj2.lt.0)lj2=nodesj-1
         lk2 = mod(lk,nodesk);if(lk2.lt.0)lk2=nodesk-1
         lijkproc = procijk(li2,lj2,lk2)
         if (lijkproc.eq.mype) cycle
! find the index of the neighbour
         do lproc =1,itotalneigh
            if(lijkproc.eq.itoproc(lproc)) then
               exit
            end if
         end do
         lci=(liproc-li);lcj=(ljproc-lj);lck=(lkproc-lk)
         linode_start = istart2; linode_end=iend1
         ljnode_start = jstart2; ljnode_end=jend1
         lknode_start = kstart2; lknode_end=kend1
         if(lci.eq.1) linode_end = istart2
         if(lci.eq.-1)  linode_start = iend1
         if(lcj.eq.1) ljnode_end = jstart2
         if(lcj.eq.-1)  ljnode_start = jend1
         if(lck.eq.1) lknode_end = kstart2
         if(lck.eq.-1)  lknode_start = kend1
         lcount = istartsend(lproc)+iprocsumindx(lproc)
         do lknode = lknode_start,lknode_end
         do linode = linode_start,linode_end
         do ljnode = ljnode_start,ljnode_end
            IF(DEAD_CELL_AT(linode,ljnode,lknode)) CYCLE
            isendnodes(lcount)=funijk(linode,ljnode,lknode)
            iprocsumindx(lproc)=iprocsumindx(lproc)+1
            lcount = lcount+1
         end do
         end do
         end do
      end do
      end do
      end do

!     call  des_dbgnodesr()
      end subroutine des_setnodeindices

!------------------------------------------------------------------------
! Subroutine       : des_exchangenode
! Purpose          : calls send and recv to exchange the node values and
!                    adds based on the flag
!                    to send and recv for grid node values
! Parameters       : pvar - variable that has to be exchanged
!                    padd - if true node values will be added
!                           else node values will be replaced
!------------------------------------------------------------------------
      subroutine des_exchangenode(pvar,padd)

      use desmpi_wrapper, only: DES_MPI_WAIT
      use desmpi_wrapper, only: DES_MPI_iSEND
      use desmpi_wrapper, only: DES_MPI_iRECV
      use parallel_mpi, only: MPI_Check
      use compar, only: myPE, numPEs

!-----------------------------------------------
      implicit none
!-----------------------------------------------
! dummy variables
!-----------------------------------------------
      double precision,dimension(:),intent(inout) ::pvar
      logical:: padd
!-----------------------------------------------
! local variables
!-----------------------------------------------
      character(len=80), parameter :: name = 'des_exchangenode'
      integer :: lindx,lcount,lcount2,lneigh,ltag,lerr
      integer :: lstart,lend,ltotal
!-----------------------------------------------

! steps pack the buffer call isend and irecv
      do lcount = 1,itotalneigh
         lneigh = itoproc(lcount)
         lstart = istartsend(lcount);lend=istartsend(lcount+1)-1
         do lcount2 = lstart,lend
            dsendnodebuf(lcount2) = pvar(isendnodes(lcount2))
         end do
         ltag = message_tag(lneigh,mype)
         ltotal = lend-lstart+1
         call des_mpi_irecv(drecvnodebuf(lstart:lend),ltotal, &
                            lneigh,ltag,irecvreqnode(lcount),lerr)
         call mpi_check( name //':mpi_irecv ', lerr )
         ltag = message_tag(mype,lneigh)
         call des_mpi_isend(dsendnodebuf(lstart:lend),ltotal, &
                            lneigh,ltag,isendreqnode(lcount),lerr)
         call mpi_check( name //':mpi_irecv ', lerr )
      end do
! call mpi wait to complete the exchange
      do lcount = 1,itotalneigh
         call des_mpi_wait(isendreqnode(lcount),lerr)
         call mpi_check( name //':mpi_wait-send', lerr )
         call des_mpi_wait(irecvreqnode(lcount),lerr)
         call mpi_check( name //':mpi_wait-recv', lerr )
      end do
! after receiving the buffer the values are either added or
! replaced based on the flag
      if (padd) then
         do lcount = 1,itotalindx
            lindx = isendnodes(lcount)
            pvar(lindx) = pvar(lindx) + drecvnodebuf(lcount)
         end do
      else
         do lcount = 1,itotalindx
            lindx = isendnodes(lcount)
            pvar(lindx) = drecvnodebuf(lcount)
         end do
      end if
      return

      contains

        integer function message_tag(lsource,ldest)
          implicit none
          integer, intent(in) :: lsource,ldest
          message_tag = lsource+numpes*ldest+200
        end function message_tag

      end subroutine des_exchangenode

!------------------------------------------------------------------------
! Subroutine       : des_dbgnodesr
! Purpose          : For debugging prints the indices
!------------------------------------------------------------------------
      subroutine des_dbgnodesr()
!-----------------------------------------------

      use indices, only: i_of, j_of, k_of
      use compar, only: mype

      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      character (255) filename
      integer ijk
      integer lcount,lcount2,lstart,lend
!-----------------------------------------------

! pradeep remove print the flags
      write(filename,'("dbg_nodesr",I4.4,".dat")') mype
      open(44,file=filename,convert='big_endian')
      do lcount = 1,itotalneigh
         lstart = istartsend(lcount);lend=istartsend(lcount+1)-1
         write(44,"(2/,72('*'))")
         write(44,1100) myPE, itoproc(lcount)
         write(44,"(/2x,'Start:',I6)") lstart
         write(44,"( 2x,'End:  ',I6)") lend
         write(44,"(72('-'))")
         do lcount2 = lstart,lend
            ijk = isendnodes(lcount2)
            write(44,1000)'SEND', i_of(ijk),j_of(ijk),k_of(ijk),ijk
         end do
         write(44,"(72('-'))")
      end do

      if(allocated(irecvnodes)) then
         do lcount = 1,itotalneigh
            lstart = istartrecv(lcount);lend=istartrecv(lcount+1)-1
            write(44,"(2/,72('*'))")
            write(44,1100) itoproc(lcount), myPE
            write(44,"(/2x,'Start:',I6)") lstart
            write(44,"( 2x,'End:  ',I6)") lend
            write(44,"(72('-'))")
            do lcount2 = lstart,lend
               ijk = irecvnodes(lcount2)
               write(44,1000) 'RECV', i_of(ijk),j_of(ijk),k_of(ijk),ijk
            end do
            write(44,"(72('-'))")
         end do
      end if

      close (44)

 1100 FORMAT(2x,'Send Proc ',I2,'   -->   Recv Proc' I2)
 1000 FORMAT(3x,A,': (',I3,',',I3,',',I3,') :: ',I7)

      end subroutine des_dbgnodesr

      end module


