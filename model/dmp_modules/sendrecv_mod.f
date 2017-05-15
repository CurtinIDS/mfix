!--------------------------------------------------------------------
! Purpose:
! Contains following subroutines:
!    ijk_of, ijk_of_gl, sendrecv_init
!    sendrecv_begin_1d, sendrecv_begin_1i, sendrecv_begin_1c
!    sendrecv_end_1dm, sendrecv_end_1c, sendrecv_end_1i
!    send_recv_1c, send_recv_1d, send_recv_2d, send_recv_3d
!    send_recv_1i
!--------------------------------------------------------------------
module sendrecv

  !-----------------------------------------------
  ! Modules
  !-----------------------------------------------
  use compar
  use debug
  use exit, only: mfix_exit
  use functions
  use geometry
  use indices
  use parallel_mpi
  implicit none
  !-----------------------------------------------

  logical,parameter :: localfunc=.false.
  logical,parameter :: use_persistent_message=.true.


  integer, pointer, dimension(:) :: &
       recvproc1, recvtag1, xrecv1, recvijk1, &
       sendproc1, sendtag1, xsend1, sendijk1, &
       recvproc2, recvtag2, xrecv2, recvijk2, &
       sendproc2, sendtag2, xsend2, sendijk2

  integer,pointer, dimension(:) :: &
       send_persistent_request, recv_persistent_request,     &
       send_persistent_request1, send_persistent_request2,   &
       recv_persistent_request1, recv_persistent_request2

  integer :: nrecv1,nsend1, nrecv2,nsend2


  double precision, dimension(:), pointer :: &
       dsendbuffer, drecvbuffer
  integer, dimension(:), pointer :: &
       isendbuffer, irecvbuffer
  character, dimension(:), pointer :: &
       csendbuffer, crecvbuffer

  integer :: nrecv,nsend
  integer, pointer, dimension(:) :: &
       recvrequest, sendrequest, &
       xrecv,recvproc, recvijk, recvtag, &
       xsend,sendproc, sendijk, sendtag

  integer :: communicator

  ! -----------------
  ! generic interface
  ! -----------------
  interface sendrecv_begin
     module procedure &
          sendrecv_begin_1d, &
          sendrecv_begin_1i, &
          sendrecv_begin_1c
  end interface sendrecv_begin

  interface sendrecv_end
     module procedure &
          sendrecv_end_1d, &
          sendrecv_end_1i, &
          sendrecv_end_1c
  end interface sendrecv_end

  interface send_recv
     module procedure &
          send_recv_1d, send_recv_2d, send_recv_3d, &
          send_recv_1i, &
          send_recv_1c
  end interface send_recv

contains

  !--------------------------------------------------------------------
  ! Purpose:
  !--------------------------------------------------------------------
  subroutine ijk_of( ijkp, i,j,k )

    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    integer, intent(in) :: ijkp
    integer, intent(out) :: i,j,k

    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    integer :: k1,k2, j1,j2, i1,i2, &
         ijk, isize,jsize,ksize, gijk

    character(len=32), parameter :: name = "ijk_of"
    logical :: isok_k, isok_j, isok_i, is_same, isok
    !-----------------------------------------------

    ijk = ijkp

    i1 = istart3_all(myPE)
    i2 = iend3_all(myPE)
    j1 = jstart3_all(myPE)
    j2 = jend3_all(myPE)
    k1 = kstart3_all(myPE)
    k2 = kend3_all(myPE)

    ksize = (k2-k1+1)
    jsize = (j2-j1+1)
    isize = (i2-i1+1)


    if (mod(ijk,isize*jsize).ne.0) then
       k = int( ijk/(isize*jsize) ) + k1
    else
       k = int( ijk/(isize*jsize) ) + k1 -1
    endif
    ijk = ijk - (k-k1)*(isize*jsize)

    if (mod(ijk,isize).ne.0) then
       j = int( ijk/isize ) + j1
    else
       j = int( ijk/isize ) + j1 - 1
    endif
    ijk = ijk - (j-j1)*isize

    i = (ijk-1) + i1

    ! double check
    isok_i = (i1 <= i) .and. (i <= i2)
    isok_j = (j1 <= j) .and. (j <= j2)
    isok_k = (k1 <= k) .and. (k <= k2)
    gijk = 1 + (i-i1) + (j-j1)*(i2-i1+1) + &
         (k-k1)*(j2-j1+1)*(i2-i1+1)
    is_same = (gijk .eq. ijkp)
    isok = isok_i .and. isok_j .and. isok_k .and. is_same
    if (.not.isok) then
       call write_debug( name, 'i,j,k ', i,j,k )
       call write_debug( name, 'ijkp, gijk ', ijkp, gijk )
    endif

    return
  end subroutine ijk_of

  !--------------------------------------------------------------------
  ! Purpose:
  !--------------------------------------------------------------------
  subroutine ijk_of_gl( ijkp, i,j,k )

    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    integer, intent(in) :: ijkp
    integer, intent(out) :: i,j,k

    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    integer :: k1,k2, j1,j2, i1,i2, &
         ijk, isize,jsize,ksize, gijk

    character(len=32), parameter :: name = "ijk_of_gl"
    logical :: isok_k, isok_j, isok_i, is_same, isok
    !-----------------------------------------------

    ijk = ijkp

    k1 = minval( kstart3_all(:) )
    k2 = maxval( kend3_all(:) )
    j1 = minval( jstart3_all(:) )
    j2 = maxval( jend3_all(:) )
    i1 = minval( istart3_all(:) )
    i2 = maxval( iend3_all(:) )

    ksize = (k2-k1+1)
    jsize = (j2-j1+1)
    isize = (i2-i1+1)

    if (mod(ijk,isize*jsize).ne.0) then
       k = int( ijk/(isize*jsize) ) + k1
    else
       k = int( ijk/(isize*jsize) ) + k1 -1
    endif
    ijk = ijk - (k-k1)*(isize*jsize)

    if (mod(ijk,isize).ne.0) then
       j = int( ijk/isize ) + j1
    else
       j = int( ijk/isize ) + j1 - 1
    endif
    ijk = ijk - (j-j1)*isize

    i = (ijk-1) + i1

    ! double check
    isok_i = (i1 <= i) .and. (i <= i2)
    isok_j = (j1 <= j) .and. (j <= j2)
    isok_k = (k1 <= k) .and. (k <= k2)
    gijk = 1 + (i-i1) + (j-j1)*(i2-i1+1) + &
         (k-k1)*(j2-j1+1)*(i2-i1+1)
    is_same = (gijk .eq. ijkp)
    isok = isok_i .and. isok_j .and. isok_k .and. is_same
    if (.not.isok) then
       call write_debug( name, 'i,j,k ', i,j,k )
       call write_debug( name, 'ijkp, gijk ', ijkp, gijk )
    endif

    return
  end subroutine ijk_of_gl

  !--------------------------------------------------------------------
  ! Purpose:
  ! set up tables and data structures for exchanging ghost regions
  !--------------------------------------------------------------------

  subroutine sendrecv_init(comm, &
       cyclic_i,cyclic_j,cyclic_k, idebug )

    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    integer, intent(in) :: comm
    logical,intent(in) :: cyclic_i,cyclic_j,cyclic_k
    integer, intent(in), optional :: idebug

#ifdef MPI
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    logical, parameter :: jfastest = .true.
    integer, parameter :: message_tag_offset = 11

    character(len=80), parameter :: name = 'sendrecv_init'

    character(len=80), pointer, dimension(:) :: line
    integer :: ip, lmax

    integer :: layer,request, source, tag, datatype

    integer :: lidebug
    integer :: isize,jsize,ksize, ijksize
    integer :: recvsize1, recvsize2, &
         sendsize1, sendsize2

    integer :: iter, i,j,k, ii, jj,kk, &
         ntotal, icount,ipos, &
         ilayer, i1,i2, j1,j2, k1,k2,  &
         ijk, ijk2, iproc, jproc, src,dest, &
         ierror

    logical :: isok, isvalid, ismine, is_halobc

    integer, dimension(:,:,:), pointer :: ijk2proc
    integer, pointer, dimension(:) :: &
         istartx,iendx, jstartx,jendx, kstartx,kendx, &
         ncount, &
         recvproc, recvtag, xrecv, recvijk,  &
         sendproc, sendtag, xsend, sendijk
    !-----------------------------------------------
    ! Inline functions
    !-----------------------------------------------
    integer :: message_tag
    !-----------------------------------------------

    message_tag(src,dest) = message_tag_offset + (1+src + dest*numPEs)

    nullify( &
         recvproc1, recvtag1, xrecv1, recvijk1, &
         sendproc1, sendtag1, xsend1, sendijk1, &
         recvproc2, recvtag2, xrecv2, recvijk2, &
         sendproc2, sendtag2, xsend2, sendijk2)

    nullify( &
         send_persistent_request, recv_persistent_request,   &
         send_persistent_request1, send_persistent_request2, &
         recv_persistent_request1, recv_persistent_request2 )

    nullify( dsendbuffer, drecvbuffer )
    nullify( isendbuffer, irecvbuffer )
    nullify( csendbuffer, crecvbuffer )

    nullify( &
         recvrequest, sendrequest, &
         xrecv,recvproc, recvijk, recvtag, &
         xsend,sendproc, sendijk, sendtag )

    ! initialize variables
    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    communicator = comm
    call MPI_COMM_SIZE( comm, numPEs, ierror )
    call MPI_Check( 'sendrecv_init:MPI_COMM_SIZE ', ierror )

    call MPI_COMM_RANK( comm, myPE, ierror )
    call MPI_Check( 'sendrecv_init:MPI_COMM_RANK ', ierror )

    ! check obtain bounds of domain
    ! check bounds for k-axis
    call assert( kmin1 .eq. minval( kstart1_all(:) ), &
         '** sendrecv_init: invalid kmin1, ' // &
         ' kmin1, minval(kstart1_all(:)) ', &
         kmin1, minval(kstart1_all(:)) )

    call assert( kmin2 .eq. minval( kstart2_all(:) ), &
         '** sendrecv_init: invalid kmin2, ' // &
         ' kmin2, minval(kstart2_all(:)) ', &
         kmin2, minval(kstart2_all(:)) )

    call assert( kmin3 .eq. minval( kstart3_all(:) ), &
         '** sendrecv_init: invalid kmin3, ' // &
         ' kmin3, minval(kstart3_all(:)) ', &
         kmin3, minval(kstart3_all(:)) )

    call assert( kmax1 .eq. maxval( kend1_all(:) ), &
         '** sendrecv_init: invalid kmax1, ' // &
         ' kmax1, maxval(kend1_all(:)) ', &
         kmax1, maxval(kend1_all(:)) )

    call assert( kmax2 .eq. maxval( kend2_all(:) ), &
         '** sendrecv_init: invalid kmax2, ' // &
         ' kmax2, maxval(kend2_all(:)) ', &
         kmax2, maxval(kend2_all(:)) )

    call assert( kmax3 .eq. maxval( kend3_all(:) ), &
         '** sendrecv_init: invalid kmax3, ' // &
         ' kmax3, maxval(kend3_all(:)) ', &
         kmax3, maxval(kend3_all(:)) )

    ! check bounds for j-axis
    call assert( jmin1 .eq. minval( jstart1_all(:) ), &
         '** sendrecv_init: invalid jmin1, ' // &
         ' jmin1, minval(jstart1_all(:)) ', &
         jmin1, minval(jstart1_all(:)) )

    call assert( jmin2 .eq. minval( jstart2_all(:) ), &
         '** sendrecv_init: invalid jmin2, ' // &
         ' jmin2, minval(jstart2_all(:)) ', &
         jmin2, minval(jstart2_all(:)) )

    call assert( jmin3 .eq. minval( jstart3_all(:) ), &
         '** sendrecv_init: invalid jmin3, ' // &
         ' jmin3, minval(jstart3_all(:)) ', &
         jmin3, minval(jstart3_all(:)) )

    call assert( jmax1 .eq. maxval( jend1_all(:) ), &
         '** sendrecv_init: invalid jmax1, ' // &
         ' jmax1, maxval(jend1_all(:)) ', &
         jmax1, maxval(jend1_all(:)) )

    call assert( jmax2 .eq. maxval( jend2_all(:) ), &
         '** sendrecv_init: invalid jmax2, ' // &
         ' jmax2, maxval(jend2_all(:)) ', &
         jmax2, maxval(jend2_all(:)) )

    call assert( jmax3 .eq. maxval( jend3_all(:) ), &
         '** sendrecv_init: invalid jmax3, ' // &
         ' jmax3, maxval(jend3_all(:)) ', &
         jmax3, maxval(jend3_all(:)) )

    ! check bounds for i-axis
    call assert( imin1 .eq. minval( istart1_all(:) ), &
         '** sendrecv_init: invalid imin1, ' // &
         ' imin1, minval(istart1_all(:)) ', &
         imin1, minval(istart1_all(:)) )

    call assert( imin2 .eq. minval( istart2_all(:) ), &
         '** sendrecv_init: invalid imin2, ' // &
         ' imin2, minval(istart2_all(:)) ', &
         imin2, minval(istart2_all(:)) )

    call assert( imin3 .eq. minval( istart3_all(:) ), &
         '** sendrecv_init: invalid imin3, ' // &
         ' imin3, minval(istart3_all(:)) ', &
         imin3, minval(istart3_all(:)) )

    call assert( imax1 .eq. maxval( iend1_all(:) ), &
         '** sendrecv_init: invalid imax1, ' // &
         ' imax1, maxval(iend1_all(:)) ', &
         imax1, maxval(iend1_all(:)) )

    call assert( imax2 .eq. maxval( iend2_all(:) ), &
         '** sendrecv_init: invalid imax2, ' // &
         ' imax2, maxval(iend2_all(:)) ', &
         imax2, maxval(iend2_all(:)) )

    call assert( imax3 .eq. maxval( iend3_all(:) ), &
         '** sendrecv_init: invalid imax3, ' // &
         ' imax3, maxval(iend3_all(:)) ', &
         imax3, maxval(iend3_all(:)) )



    call assert( jmin1 .le. jmax1, &
         '** sendrecv_init: jmin1,jmax1 ', jmin1,jmax1 )
    call assert( jmin2 .le. jmax2, &
         '** sendrecv_init: jmin2,jmax2 ', jmin2,jmax2 )
    call assert( jmin3 .le. jmax3, &
         '** sendrecv_init: jmin3,jmax3 ', jmin3,jmax3 )

    call assert( kmin1 .le. kmax1, &
         '** sendrecv_init: kmin1,kmax1 ', kmin1,kmax1 )
    call assert( kmin2 .le. kmax2, &
         '** sendrecv_init: kmin2,kmax2 ', kmin2,kmax2 )
    call assert( kmin3 .le. kmax3, &
         '** sendrecv_init: kmin3,kmax3 ', kmin3,kmax3 )

    call assert( imin1 .le. imax1, &
         '** sendrecv_init: imin1,imax1 ', imin1,imax1 )
    call assert( imin2 .le. imax2, &
         '** sendrecv_init: imin2,imax2 ', imin2,imax2 )
    call assert( imin3 .le. imax3, &
         '** sendrecv_init: imin3,imax3 ', imin3,imax3 )



    k1 = min( kmin1, min(kmin2, kmin3) )
    k2 = max( kmax1, max(kmax2, kmax3) )
    j1 = min( jmin1, min(jmin2, jmin3) )
    j2 = max( jmax1, max(jmax2, jmax3) )
    i1 = min( imin1, min(imin2, imin3) )
    i2 = max( imax1, max(imax2, imax3) )

    allocate( ijk2proc( i1:i2, j1:j2, k1:k2 ) )

    if(localfunc) then
       ! double check ijk_of()
       do k=kstart3_all(myPE),kend3_all(myPE)
          do j=jstart3_all(myPE),jend3_all(myPE)
             do i=istart3_all(myPE),iend3_all(myPE)
                ijk = funijk(i,j,k)
                call ijk_of(ijk, ii,jj,kk)
                ijk2 = funijk( ii,jj,kk)

                isvalid = (ii.eq.i).and.(jj.eq.j).and.(kk.eq.k).and.(ijk.eq.ijk2)
                if (.not.isvalid) then
                   call write_debug( name, 'error with ijk_of ')

                   call write_debug( name, &
                        'istart3_all(myPE),iend3_all(myPE) ', &
                        istart3_all(myPE),iend3_all(myPE) )
                   call write_debug( name, &
                        'jstart3_all(myPE),jend3_all(myPE) ', &
                        jstart3_all(myPE),jend3_all(myPE) )
                   call write_debug( name, &
                        'kstart3_all(myPE),kend3_all(myPE) ', &
                        kstart3_all(myPE),kend3_all(myPE) )

                   call write_debug( name, 'i,j,k, ijk ', i,j,k, ijk )
                   call write_debug( name, 'ii,jj,kk,  ijk2 ',&
                        ii,jj,kk,ijk2 )

                endif
             enddo
          enddo
       enddo
    endif ! Local Function

    if (lidebug.ge.1) then
       call write_debug( name, 'imap ', imap )
       call write_debug( name, 'jmap ', jmap )
       call write_debug( name, 'kmap ', kmap )
    endif


    ! ----------------------------
    ! set up table ijk2proc(:,:,:)
    !
    ! ijk2proc(i,j,k) maps (i,j,k) index to
    ! unique processor that 'owns' that node.
    ! ----------------------------

    ijk2proc( :,:,: ) = 0

    ! --------------------------------------------------
    ! double check domain decomposition that
    ! each interior node is assigned to UNIQUE processor
    ! --------------------------------------------------
    do iproc=0,numPEs-1
       i1 = istart1_all(iproc)
       i2 = iend1_all(iproc)
       j1 = jstart1_all(iproc)
       j2 = jend1_all(iproc)
       k1 = kstart1_all(iproc)
       k2 = kend1_all(iproc)
       if(istart3_all(iproc).eq.imin3) i1 = istart3_all(iproc)
       if(iend3_all(iproc).eq.imax3) i2 = iend3_all(iproc)
       if(jstart3_all(iproc).eq.jmin3) j1 = jstart3_all(iproc)
       if(jend3_all(iproc).eq.jmax3) j2 = jend3_all(iproc)
       if(kstart3_all(iproc).eq.kmin3) k1 = kstart3_all(iproc)
       if(kend3_all(iproc).eq.kmax3) k2 = kend3_all(iproc)
       do k=k1,k2
          do j=j1,j2
             do i=i1,i2
                ijk2proc(i,j,k) = ijk2proc(i,j,k) + 1
             enddo
          enddo
       enddo
    enddo

    do k=lbound(ijk2proc,3),ubound(ijk2proc,3)
       do j=lbound(ijk2proc,2),ubound(ijk2proc,2)
          do i=lbound(ijk2proc,1),ubound(ijk2proc,1)
             isvalid = (ijk2proc(i,j,k) .eq. 1)
             if (.not.isvalid) then
                call write_debug(name, ' invalid decomposition ')
                call write_debug(name, 'i,j,k ',i,j,k )
                call write_debug(name, 'ijk2proc(i,j,k) ', &
                     ijk2proc(i,j,k))
                call mfix_exit( myPE )
             endif
          enddo
       enddo
    enddo

    ijk2proc(:,:,:) = -1
    do iproc=0,numPEs-1
       i1 = istart1_all(iproc)
       i2 = iend1_all(iproc)
       j1 = jstart1_all(iproc)
       j2 = jend1_all(iproc)
       k1 = kstart1_all(iproc)
       k2 = kend1_all(iproc)
       if(istart3_all(iproc).eq.imin3) i1 = istart3_all(iproc)
       if(iend3_all(iproc).eq.imax3) i2 = iend3_all(iproc)
       if(jstart3_all(iproc).eq.jmin3) j1 = jstart3_all(iproc)
       if(jend3_all(iproc).eq.jmax3) j2 = jend3_all(iproc)
       if(kstart3_all(iproc).eq.kmin3) k1 = kstart3_all(iproc)
       if(kend3_all(iproc).eq.kmax3) k2 = kend3_all(iproc)
       do k=k1,k2
          do j=j1,j2
             do i=i1,i2
                ijk2proc(i,j,k) = iproc
             enddo
          enddo
       enddo
    enddo


    allocate( ncount(0:numPEs-1) )
    allocate( istartx(0:numPEs-1) )
    allocate( jstartx(0:numPEs-1) )
    allocate( kstartx(0:numPEs-1) )
    allocate( iendx(0:numPEs-1) )
    allocate( jendx(0:numPEs-1) )
    allocate( kendx(0:numPEs-1) )


    do ilayer=1,2
       if (ilayer.eq.1) then
          kstartx(:) = kstart2_all(:)
          kendx(:) = kend2_all(:)
          jstartx(:) = jstart2_all(:)
          jendx(:) = jend2_all(:)
          istartx(:) = istart2_all(:)
          iendx(:) = iend2_all(:)
       else
          kstartx(:) = kstart3_all(:)
          kendx(:) = kend3_all(:)
          jstartx(:) = jstart3_all(:)
          jendx(:) = jend3_all(:)
          istartx(:) = istart3_all(:)
          iendx(:) = iend3_all(:)
       endif

       if (lidebug.ge.1) then
          call write_debug(name, 'determine send schedule ', myPE )
       endif

       ! -----------------------
       ! determine send schedule
       ! examine all neighboring processors to see if they need my data
       !  -----------------------

       ! first pass to determine array sizes
       ! ---------------------------------------------------------------->>>
       ncount(:) = 0
       do iproc=0,numPEs-1
          if (iproc.ne.myPE) then
             k1 = lbound(ijk2proc,3)
             k2 = ubound(ijk2proc,3)
             j1 = lbound(ijk2proc,2)
             j2 = ubound(ijk2proc,2)
             i1 = lbound(ijk2proc,1)
             i2 = ubound(ijk2proc,1)

             do k=kstartx(iproc),kendx(iproc)
                do j=jstartx(iproc),jendx(iproc)
                   do i=istartx(iproc),iendx(iproc)
                      ii = imap(i)
                      jj = jmap(j)
                      kk = kmap(k)

                      isvalid  = (k1.le.kk).and.(kk.le.k2)
                      call assert( isvalid, '** sendrecv_init: invalid kk ', kk )
                      isvalid  = (j1.le.jj).and.(jj.le.j2)
                      call assert( isvalid, '** sendrecv_init: invalid jj ', jj )
                      isvalid  = (i1.le.ii).and.(ii.le.i2)
                      call assert( isvalid, '** sendrecv_init: invalid ii ', ii )
                      jproc = ijk2proc( ii,jj,kk )

                      ismine = (jproc .eq. myPE)
                      if (ismine) then
                         ncount(iproc) = ncount(iproc) + 1
                      endif
                   enddo
                enddo
             enddo
          endif
       enddo   ! end do (iproc=0,numPEs-1)
       ! ----------------------------------------------------------------<<<

       ! prepare arrays
       ! ---------------------------------------------------------------->>>
       ntotal = 0
       nsend = 0
       do iproc=0,numPEs-1
          ntotal = ntotal + ncount(iproc)
          if (ncount(iproc).ge.1) then
             nsend = nsend + 1
          endif
       enddo

       if (lidebug.ge.1) then
          call write_debug( name, 'ncount = ', ncount )
          call write_debug( name, 'nsend, ntotal ', nsend, ntotal )
       endif

       allocate( xsend(nsend+1) )
       allocate( sendijk( max(1,ntotal) ) )
       allocate( sendproc(max(1,nsend)) )

       nsend = 0
       do iproc=0,numPEs-1
          if (ncount(iproc).ne.0) then
             nsend = nsend + 1
             sendproc(nsend) = iproc
          endif
       enddo

       xsend(1) = 1
       do i=1,nsend
          iproc = sendproc(i)
          xsend(i+1) = xsend(i) + ncount(iproc)
       enddo

       allocate( sendtag( max(1,nsend) ) )
       do ii=1,nsend
          iproc = sendproc(ii)
          src = myPE
          dest = iproc
          sendtag(ii) = message_tag( src, dest )
       enddo
       ! ----------------------------------------------------------------<<<


       ! second pass to fill in arrays
       ! ---------------------------------------------------------------->>>
       ipos = 1
       do iter=1,nsend
          iproc = sendproc(iter)
          icount = 0
          do k=kstartx(iproc),kendx(iproc)

             if (jfastest) then
                do i=istartx(iproc),iendx(iproc)
                   do j=jstartx(iproc),jendx(iproc)
                      ii = imap(i)
                      jj = jmap(j)
                      kk = kmap(k)
                      jproc = ijk2proc(ii,jj,kk)
                      ismine = (jproc.eq.myPE)
                      if (ismine) then
                         icount = icount + 1
                         ijk = funijk(ii,jj,kk)
                         ipos = xsend(iter)-1 + icount
                         sendijk( ipos ) = ijk
                      endif
                   enddo
                enddo
             else
                do j=jstartx(iproc),jendx(iproc)
                   do i=istartx(iproc),iendx(iproc)
                      ii = imap(i)
                      jj = jmap(j)
                      kk = kmap(k)
                      jproc = ijk2proc(ii,jj,kk)
                      ismine = (jproc.eq.myPE)
                      if (ismine) then
                         icount = icount + 1
                         ijk = funijk(ii,jj,kk)
                         ipos = xsend(iter)-1 + icount
                         sendijk( ipos ) = ijk
                      endif
                   enddo
                enddo
             endif
          enddo   ! end do (k=kstartx,kendx)
          isvalid = (icount .eq. ncount(iproc))
          call assert( isvalid, &
               '** sendrecv_init: icount != ncount(iproc) ', iproc)
       enddo   ! end do (iter=1,nsend)

       if (lidebug.ge.1) then
          call write_debug(name, 'determine recv schedule ', myPE )
       endif
       ! ----------------------------------------------------------------<<<


       ! ---------------------------
       ! determine recv schedule
       ! examine nodes in my ghost region and see what data is needed from
       ! my neighbors
       ! ---------------------------


       ! first pass to determine array sizes
       ! ---------------------------------------------------------------->>>

       ncount(:) = 0
       k1 = lbound(ijk2proc,3)
       k2 = ubound(ijk2proc,3)
       j1 = lbound(ijk2proc,2)
       j2 = ubound(ijk2proc,2)
       i1 = lbound(ijk2proc,1)
       i2 = ubound(ijk2proc,1)

       do k=kstartx(myPE),kendx(myPE)
          do j=jstartx(myPE),jendx(myPE)
             do i=istartx(myPE),iendx(myPE)
                ii = imap(i)
                jj = jmap(j)
                kk = kmap(k)

                isvalid  = (k1.le.kk).and.(kk.le.k2)
                call assert( isvalid, '** sendrecv_init: invalid kk ', kk )

                isvalid  = (j1.le.jj).and.(jj.le.j2)
                call assert( isvalid, '** sendrecv_init: invalid jj ', jj )

                isvalid  = (i1.le.ii).and.(ii.le.i2)
                call assert( isvalid, '** sendrecv_init: invalid ii ', ii )


                iproc = ijk2proc(ii,jj,kk)
                is_halobc = (iproc.eq.-1)
                ismine = (iproc.eq.myPE)
                if (.not.ismine) then
                   isvalid = (0 .le. iproc) .and. &
                        (iproc.le.numPEs-1) .and. &
                        (iproc.ne.myPE)
                   call assert( isvalid, &
                        '** sendrecv_init: invalid iproc ',iproc)

                   ncount(iproc) = ncount(iproc) + 1
                endif
             enddo
          enddo
       enddo

       ncount(myPE) = 0

       ntotal = 0
       do iproc=0,numPEs-1
          ntotal = ntotal + ncount(iproc)
       enddo

       nrecv = count( ncount(:) .ne. 0)

       allocate( recvproc( max(1,nrecv) ) )

       nrecv = 0
       do iproc=0,numPEs-1
          if (ncount(iproc).ne.0) then
             nrecv = nrecv + 1
             recvproc(nrecv) = iproc
          endif
       enddo

       allocate( xrecv(nrecv+1) )
       allocate( recvijk(max(1,ntotal)) )

       xrecv(1) = 1
       do iter=1,nrecv
          iproc = recvproc(iter)
          xrecv(iter+1) = xrecv(iter) + ncount(iproc)
       enddo

       allocate( recvtag( max(1,nrecv) ) )

       do iter=1,nrecv
          iproc = recvproc(iter)
          src = iproc
          dest = myPE
          recvtag(iter) = message_tag( src, dest )
       enddo
       ! ----------------------------------------------------------------<<<


       ! second pass to fill in array
       ! ---------------------------------------------------------------->>>
       if (lidebug.ge.1) then
          call write_debug( name, 'recv second pass ', myPE )
       endif

       ipos = 1

       do iter=1,nrecv
          jproc = recvproc(iter)
          do k=kstartx(myPE),kendx(myPE)

             if (jfastest) then
                do i=istartx(myPE),iendx(myPE)
                   do j=jstartx(myPE),jendx(myPE)
                      ii = imap(i)
                      jj = jmap(j)
                      kk = kmap(k)

                      iproc = ijk2proc(ii,jj,kk)
                      is_halobc = (iproc.eq.-1)
                      ismine = (iproc.eq.myPE)

                      if ((.not.ismine) .and. (iproc.eq.jproc)) then
                         ijk = funijk( i,j,k)
                         recvijk( ipos ) = ijk
                         ipos = ipos + 1
                      endif
                   enddo
                enddo

             else
                do j=jstartx(myPE),jendx(myPE)
                   do i=istartx(myPE),iendx(myPE)
                      ii = imap(i)
                      jj = jmap(j)
                      kk = kmap(k)

                      iproc = ijk2proc(ii,jj,kk)
                      is_halobc = (iproc.eq.-1)
                      ismine = (iproc.eq.myPE)

                      if ((.not.ismine) .and. (iproc.eq.jproc)) then
                         ijk = funijk( i,j,k)
                         recvijk( ipos ) = ijk
                         ipos = ipos + 1
                      endif
                   enddo
                enddo
             endif   ! end if/else (if(jfastest))
          enddo   ! end do (k=kstartx,kendx)
       enddo   ! end do (iter=1,nrecv)
       ! ----------------------------------------------------------------<<<

       if (ilayer.eq.1) then
          nsend1 = nsend
          xsend1 => xsend
          sendijk1 => sendijk
          sendproc1 => sendproc
          sendtag1 => sendtag

          nrecv1 = nrecv
          xrecv1 => xrecv
          recvijk1 => recvijk
          recvproc1 => recvproc
          recvtag1 => recvtag
       else
          nsend2 = nsend
          xsend2 => xsend
          sendijk2 => sendijk
          sendproc2 => sendproc
          sendtag2 => sendtag

          nrecv2 = nrecv
          xrecv2 => xrecv
          recvijk2 => recvijk
          recvproc2 => recvproc
          recvtag2 => recvtag
       endif


       nullify( xsend )
       nullify( sendijk )
       nullify( sendproc )
       nullify( sendtag )
       nullify( xrecv )
       nullify( recvijk )
       nullify( recvproc )
       nullify( recvtag )

    enddo ! end do (ilayer=1,2)


    deallocate( ncount )
    deallocate( ijk2proc )

    deallocate( istartx )
    deallocate( jstartx )
    deallocate( kstartx )
    deallocate( iendx )
    deallocate( jendx )
    deallocate( kendx )

    nullify( ncount )
    nullify( ijk2proc )

    nullify( istartx )
    nullify( jstartx )
    nullify( kstartx )
    nullify( iendx )
    nullify( jendx )
    nullify( kendx )


    ! ---------------------------------------------------------------->>>
    if (lidebug.ge.1) then
       call write_debug( name, ' allocate message buffers ' )
       call write_debug( name, 'nrecv1 ', nrecv1 )
       call write_debug( name, 'recvproc1 ', recvproc1 )
       call write_debug( name, 'recvtag1 ', recvtag1 )
       call write_debug( name, 'xrecv1 ', xrecv1 )

       lmax = size(recvijk1)
       allocate( line(lmax) )
       line(:) = " "
       ip = 1
       do ii=lbound(recvijk1,1),ubound(recvijk1,1)
          ijk = recvijk1(ii)
          if(localfunc) then
             call ijk_of(ijk,i,j,k)
          else
             i = i_of(ijk)
             j = j_of(ijk)
             k = k_of(ijk)
          endif

          write(line(ip),9001) ii,ijk, i,j,k
9001      format('recvijk1( ', i6,') = ', &
               i6, '( ', i6,',',i6,',',i6,') ')

          ip = ip + 1
       enddo
       call write_error( name, line, lmax )
       deallocate( line )
       nullify( line )


       lmax = size(recvijk2)
       allocate( line(lmax) )
       line(:) = " "
       ip = 1
       do ii=lbound(recvijk2,1),ubound(recvijk2,1)
          ijk = recvijk2(ii)
          if(localfunc) then
             call ijk_of(ijk,i,j,k)
          else
             i = i_of(ijk)
             j = j_of(ijk)
             k = k_of(ijk)
          endif

          write(line(ip),9101) ii,ijk, i,j,k
9101      format('recvijk2( ', i6,') = ', &
               i6, '( ', i6,',',i6,',',i6,') ')

          ip = ip + 1
       enddo
       call write_error( name, line, lmax )
       deallocate( line )
       nullify( line )

       call write_debug( name, ' allocate message buffers ' )
       call write_debug( name, 'nsend1 ', nsend1 )
       call write_debug( name, 'sendproc1 ', sendproc1 )
       call write_debug( name, 'sendtag1 ', sendtag1 )
       call write_debug( name, 'xsend1 ', xsend1 )

       lmax = size(sendijk1)
       allocate(line(lmax))
       line(:) = " "
       ip = 1
       do ii=lbound(sendijk1,1),ubound(sendijk1,1)
          ijk = sendijk1(ii)
          if(localfunc) then
             call ijk_of(ijk,i,j,k)
          else
             i = i_of(ijk)
             j = j_of(ijk)
             k = k_of(ijk)
          endif

          write(line(ip),9002) ii,ijk,   i,j,k
9002      format('sendijk1( ', i6,') = ', &
               i6, '( ', i6,',',i6,',',i6,') ')

          ip = ip + 1
       enddo
       call write_error( name, line, lmax )
       deallocate( line )
       nullify( line )

       lmax = size(sendijk2)
       allocate(line(lmax))
       line(:) = " "
       ip = 1
       do ii=lbound(sendijk2,1),ubound(sendijk2,1)
          ijk = sendijk2(ii)
          if(localfunc) then
             call ijk_of(ijk,i,j,k)
          else
             i = i_of(ijk)
             j = j_of(ijk)
             k = k_of(ijk)
          endif

          write(line(ip),9102) ii,ijk,   i,j,k
9102      format('sendijk2( ', i6,') = ', &
               i6, '( ', i6,',',i6,',',i6,') ')

          ip = ip + 1
       enddo
       call write_error( name, line, lmax )
       deallocate( line )
       nullify( line )
    endif   ! end if (lidebug.ge.1)
    ! ----------------------------------------------------------------<<<



    ! allocate message buffers
    isize = max(1,max(nsend1,nsend2))
    allocate( sendrequest( isize ) )
    allocate( send_persistent_request1( isize ) )
    allocate( send_persistent_request2( isize ) )

    isize = max(1,max(nrecv1,nrecv2))
    allocate( recvrequest( isize ) )
    allocate( recv_persistent_request1( isize ) )
    allocate( recv_persistent_request2( isize ) )


    ! preallocate buffers for common case
    recvsize1 = xrecv1( nrecv1+1)-1
    recvsize2 = xrecv2( nrecv2+1)-1

    isize = max(1,max(recvsize1,recvsize2))
    allocate( drecvbuffer( isize ) )

    sendsize1 = xsend1( nsend1+1)-1
    sendsize2 = xsend2( nsend2+1)-1

    isize = max(1,max(sendsize1,sendsize2))
    allocate( dsendbuffer( isize ) )


    if (use_persistent_message) then
       datatype = MPI_DOUBLE_PRECISION

       do layer=1,2
          if (layer.eq.1) then
             nrecv = nrecv1
             recvtag =>recvtag1
             recvproc => recvproc1
             recvijk => recvijk1
             xrecv => xrecv1

             nsend = nsend1
             sendtag => sendtag1
             sendproc => sendproc1
             sendijk => sendijk1
             xsend => xsend1

             send_persistent_request => send_persistent_request1
             recv_persistent_request => recv_persistent_request1

          else
             nrecv = nrecv2
             recvtag =>recvtag2
             recvproc => recvproc2
             recvijk => recvijk2
             xrecv => xrecv2

             nsend = nsend2
             sendtag => sendtag2
             sendproc => sendproc2
             sendijk => sendijk2
             xsend => xsend2

             send_persistent_request => send_persistent_request2
             recv_persistent_request => recv_persistent_request2
          endif   ! end if/else (layer.eq.1)

          do ii=1,nrecv
             j1 = xrecv(ii)
             j2 = xrecv(ii+1)-1
             icount = j2-j1+1
             source = recvproc( ii )
             tag = recvtag( ii )

             if (lidebug.ge.2) then
                call write_debug(name, 'mpi_recv_init: ii,j1,j2 ', &
                     ii,j1,j2 )
                call write_debug(name, 'icount, source, tag ', &
                     icount,source,tag )
             endif

             call MPI_RECV_INIT( drecvbuffer(j1), icount, datatype, &
                  source, tag, comm, request, ierror )
             call MPI_Check( 'sendrecv_begin_1d:MPI_IRECV ', ierror )

             recv_persistent_request(ii) = request
          enddo   ! end do (ii=1,nrecv)

          do ii=1,nsend
             j1 = xsend(ii)
             j2 = xsend(ii+1)-1
             dest = sendproc( ii )
             tag = sendtag( ii )
             icount = j2-j1+1

             if (lidebug.ge.2) then
                call write_debug(name, 'mpi_send_init: ii,j1,j2 ', &
                     ii,j1,j2)
                call write_debug(name, 'icount, dest, tag ', &
                     icount,dest,tag )
             endif

             call MPI_SEND_INIT( dsendbuffer(j1), icount, datatype, &
                  dest, tag, comm, request, ierror )
             call MPI_Check( 'sendrecv_begin_1d:MPI_SEND_INIT ', &
                  ierror)

             send_persistent_request( ii ) = request
          enddo   ! end do (ii=1,nsend)
       enddo   ! end do (layer=1,2)

    endif   ! end if (use_persistent_message)

    if (lidebug.ge.1) then
       call write_debug(name, ' end of sendrecv_init ', myPE )
    endif
#endif

  end subroutine sendrecv_init


  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine sendrecv_begin_1d( XX, ilayer, idebug )

    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    integer, intent(in),optional :: ilayer
    double precision, intent(inout), dimension(:) :: XX
    integer, intent(in), optional :: idebug

#ifdef MPI
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    character(len=80), parameter :: name = 'sendrecv_begin_1d'
    integer :: lidebug
    integer :: layer, datatype, comm, recvsize, sendsize, &
         request, count, source,dest, tag, ierror
    integer :: ijk, jj, j1, j2, ii

    !-----------------------------------------------
    lidebug = 0

    if (present(idebug)) then
       lidebug = idebug
    endif

    layer = 1
    if (present(ilayer)) then
       layer = ilayer
    endif

    if (layer.eq.1) then
       nrecv = nrecv1
       recvtag =>recvtag1
       recvproc => recvproc1
       recvijk => recvijk1
       xrecv => xrecv1

       nsend = nsend1
       sendtag => sendtag1
       sendproc => sendproc1
       sendijk => sendijk1
       xsend => xsend1

       send_persistent_request => send_persistent_request1
       recv_persistent_request => recv_persistent_request1

    else
       nrecv = nrecv2
       recvtag =>recvtag2
       recvproc => recvproc2
       recvijk => recvijk2
       xrecv => xrecv2

       nsend = nsend2
       sendtag => sendtag2
       sendproc => sendproc2
       sendijk => sendijk2
       xsend => xsend2

       send_persistent_request => send_persistent_request2
       recv_persistent_request => recv_persistent_request2

    endif   ! end if/else (layer.eq.1)


    ! post asynchronous receives
    ! ---------------------------------------------------------------->>>
    if (lidebug.ge.1) then
       call write_debug(name, 'post asynchronous receives, nrecv = ', &
            nrecv)
    endif

    if (nrecv.ge.1) then
       recvsize = xrecv( nrecv+1)-1

       if (lidebug.ge.1) then
          call write_debug( name, 'recvsize, ubound(drecvbuffer,1) ',&
               recvsize, ubound(drecvbuffer,1) )

          call write_debug( name, 'ubound(xrecv,1) ', &
               ubound(xrecv,1) )
          call write_debug( name, 'ubound(recvproc,1) ', &
               ubound(recvproc,1) )
          call write_debug( name, 'ubound(recvtag,1) ', &
               ubound(recvtag,1) )
       endif

       ! post receives
       datatype = MPI_DOUBLE_PRECISION
       comm = communicator

       if (use_persistent_message) then
          ! persistent request already established
          if (lidebug.ge.2) then
             call write_debug( name,'before startall for recv ',&
                  recv_persistent_request)
          endif

          call MPI_STARTALL( nrecv, recv_persistent_request, ierror )

          if (lidebug.ge.2) then
             call write_debug( name,'after startall for recv, ierror',&
                  ierror)
          endif

          call MPI_Check( 'sendrecv_begin: MPI_STARTALL ', ierror )

       else
          ! use irecv
          do ii=1,nrecv
             j1 = xrecv(ii)
             j2 = xrecv(ii+1)-1
             count = j2-j1+1
             source = recvproc( ii )
             tag = recvtag( ii )

             if (lidebug.ge.2) then
                call write_debug(name, 'mpi_irecv: ii,j1,j2 ', &
                     ii, j1, j2)
                call write_debug(name, 'count, source, tag ', &
                     count,source,tag )
             endif

             call MPI_IRECV( drecvbuffer(j1), count, datatype, &
                  source, tag, comm, request, ierror )

             call MPI_Check( 'sendrecv_begin_1d:MPI_IRECV ', ierror )

             recvrequest( ii ) = request
          enddo
       endif   ! end if/else (use_persistent_message)
    endif   ! end if (nrecv.ge.1)
    ! ----------------------------------------------------------------<<<

    ! post asynchronous sends
    ! ---------------------------------------------------------------->>>
    if (lidebug.ge.1) then
       call write_debug(name, 'post asynchronous sends ')
    endif

    if (nsend.ge.1) then
       sendsize = xsend( nsend+1)-1

       if (lidebug.ge.1) then
          call write_debug( name, &
               'sendsize, ubound(dsendbuffer,1) ', &
               sendsize, ubound(dsendbuffer,1) )

          call write_debug( name, 'ubound(xsend,1) ', &
               ubound(xsend,1) )
          call write_debug( name, 'ubound(sendproc,1) ', &
               ubound(sendproc,1) )
          call write_debug( name, 'ubound(sendtag,1) ', &
               ubound(sendtag,1) )
       endif

       ! perform sends
       datatype = MPI_DOUBLE_PRECISION
       comm = communicator

       if (use_persistent_message) then
          ! persistent request already established

          ! perform copy into dsendbuffer
          j1 = xsend(1)
          j2 = xsend(nsend+1)-1

          do jj=j1,j2
             ijk = sendijk( jj )
             dsendbuffer( jj )  = XX(ijk)
          enddo

          if (lidebug.ge.2) then
             call write_debug(name,'before mpi_startall send ',&
                  send_persistent_request )
          endif

          call MPI_STARTALL( nsend, send_persistent_request, ierror )

          if (lidebug .ge.2) then
             call write_debug(name,'after mpi_startall send ',&
                  send_persistent_request )
          endif

          call MPI_Check( 'sendrecv_begin_1d:MPI_STARTALL ', ierror )

       else

          do ii=1,nsend
             ! perform copy into dsendbuffer
             j1 = xsend(ii)
             j2 = xsend(ii+1)-1
             count = j2-j1+1

             do jj=j1,j2
                ijk = sendijk( jj )
                dsendbuffer(jj) = XX(ijk)
             enddo

             dest = sendproc( ii )
             tag = sendtag( ii )

             if (lidebug.ge.2) then
                call write_debug(name, 'mpi_isend: ii,j1,j2 ', &
                     ii,j1,j2)
                call write_debug(name, 'count, dest, tag ', &
                     count,dest,tag)
             endif

             call MPI_ISEND( dsendbuffer(j1), count, datatype, dest, &
                  tag, comm, request, ierror )
             call MPI_Check( 'sendrecv_begin_1d:MPI_ISEND ', ierror )

             sendrequest( ii ) = request
          enddo   ! end do (ii=1,nsend)
       endif   ! end if/else (use_persistent_message)
    endif   ! end if  (nsend.ge.1)
    ! ----------------------------------------------------------------<<<
#endif

    return
  end subroutine sendrecv_begin_1d

  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine sendrecv_begin_1i( XX, ilayer, idebug )

    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    integer, intent(in),optional :: ilayer
    integer, intent(inout), dimension(:) :: XX
    integer, intent(in), optional :: idebug

#ifdef MPI
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    character(len=80), parameter :: name = 'sendrecv_begin_1i'
    integer :: lidebug
    integer :: layer, datatype, comm, recvsize, sendsize, &
         request, count, source, dest, tag, ierror
    integer :: ijk, jj, j1, j2, ii

    !-----------------------------------------------

    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    layer = 1
    if (present(ilayer)) then
       layer = ilayer
    endif

    if (layer.eq.1) then
       nrecv = nrecv1
       recvtag =>recvtag1
       recvproc => recvproc1
       recvijk => recvijk1
       xrecv => xrecv1

       nsend = nsend1
       sendtag => sendtag1
       sendproc => sendproc1
       sendijk => sendijk1
       xsend => xsend1
    else
       nrecv = nrecv2
       recvtag =>recvtag2
       recvproc => recvproc2
       recvijk => recvijk2
       xrecv => xrecv2

       nsend = nsend2
       sendtag => sendtag2
       sendproc => sendproc2
       sendijk => sendijk2
       xsend => xsend2
    endif   ! end if/else (layer.eq.1)


    ! post asynchronous receives
    ! ---------------------------------------------------------------->>>

    if (lidebug.ge.1) then
       call write_debug(name, &
            'post asynchronous receives, nrecv = ', nrecv )
    endif

    if (nrecv.ge.1) then
       recvsize = xrecv( nrecv+1)-1
       allocate( irecvbuffer( recvsize ) )

       if (lidebug.ge.1) then
          call write_debug( name, &
               'recvsize, ubound(irecvbuffer,1) ', &
               recvsize, ubound(irecvbuffer,1) )
          call write_debug( name, 'ubound(xrecv,1) ', &
               ubound(xrecv,1) )
          call write_debug( name, 'ubound(recvproc,1) ', &
               ubound(recvproc,1) )
          call write_debug( name, 'ubound(recvtag,1) ', &
               ubound(recvtag,1) )
       endif

       ! post receives
       datatype = MPI_INTEGER
       comm = communicator

       do ii=1,nrecv
          j1 = xrecv(ii)
          j2 = xrecv(ii+1)-1
          count = j2-j1+1
          source = recvproc( ii )
          tag = recvtag( ii )

          if (lidebug.ge.2) then
             call write_debug(name, 'mpi_irecv: ii,j1,j2 ', ii,j1,j2 )
             call write_debug(name, 'count, source, tag ', &
                  count,source,tag )
          endif

          call MPI_IRECV( irecvbuffer(j1), count, datatype, &
               source, tag, comm, request, ierror )
          call MPI_Check( 'sendrecv_begin_1i:MPI_IRECV ', ierror )

          recvrequest( ii ) = request
       enddo   ! end do (ii=1,nrecv)
    endif   ! end if (nrecv.ge.1)
    ! ----------------------------------------------------------------<<<



    !  post asynchronous sends
    ! ---------------------------------------------------------------->>>
    if (lidebug.ge.1) then
       call write_debug(name, 'post asynchronous sends ')
    endif

    if (nsend.ge.1) then
       sendsize = xsend( nsend+1)-1
       allocate( isendbuffer( sendsize ) )

       if (lidebug.ge.1) then
          call write_debug( name, 'sendsize, ubound(isendbuffer,1) ',&
               sendsize, ubound(isendbuffer,1) )
          call write_debug( name, 'ubound(xsend,1) ', &
               ubound(xsend,1) )
          call write_debug( name, 'ubound(sendproc,1) ', &
               ubound(sendproc,1) )
          call write_debug( name, 'ubound(sendtag,1) ', &
               ubound(sendtag,1) )
       endif

       ! perform sends
       datatype = MPI_INTEGER
       comm = communicator

       do ii=1,nsend
          ! perform copy into sendbuffer
          j1 = xsend(ii)
          j2 = xsend(ii+1)-1
          count = j2-j1+1

          do jj=j1,j2
             ijk = sendijk( jj )
             isendbuffer(jj) = XX(ijk)
          enddo

          dest = sendproc( ii )
          tag = sendtag( ii )

          if (lidebug.ge.2) then
             call write_debug(name, 'mpi_isend: ii,j1,j2 ', ii,j1,j2)
             call write_debug(name, 'count, dest, tag ', count, &
                  dest, tag)
          endif

          call MPI_ISEND( isendbuffer(j1), count, datatype, dest, &
               tag, comm, request, ierror )
          call MPI_Check( 'sendrecv_begin_1i:MPI_ISEND ', ierror )

          sendrequest( ii ) = request
       enddo   ! end do (ii=1,nsend)
    endif   ! end if (nsend.ge.1)
    ! ----------------------------------------------------------------<<<
#endif

    return
  end subroutine sendrecv_begin_1i


  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine sendrecv_begin_1c( XX, ilayer, idebug )

    use functions

    implicit none

    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    integer, intent(in),optional :: ilayer
    character(len=*), intent(inout), dimension(:) :: XX
    integer, intent(in), optional :: idebug

#ifdef MPI
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    character(len=80), parameter :: name = 'sendrecv_begin_1c'
    integer :: lidebug
    integer :: layer, datatype, comm, recvsize, sendsize, &
         request, count, source, dest, tag, ierror
    integer :: ijk, jj, j1, j2, ii
    integer :: ic, clen, jpos

    !-----------------------------------------------

    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    layer = 1
    if (present(ilayer)) then
       layer = ilayer
    endif

    jpos = lbound(XX,1)
    clen = len( XX( jpos ) )

    if (layer.eq.1) then
       nrecv = nrecv1
       recvtag =>recvtag1
       recvproc => recvproc1
       recvijk => recvijk1
       xrecv => xrecv1

       nsend = nsend1
       sendtag => sendtag1
       sendproc => sendproc1
       sendijk => sendijk1
       xsend => xsend1
    else
       nrecv = nrecv2
       recvtag =>recvtag2
       recvproc => recvproc2
       recvijk => recvijk2
       xrecv => xrecv2

       nsend = nsend2
       sendtag => sendtag2
       sendproc => sendproc2
       sendijk => sendijk2
       xsend => xsend2
    endif   ! end if/else (layer.eq.1)


    ! post asynchronous receives
    ! ---------------------------------------------------------------->>>
    if (lidebug.ge.1) then
       call write_debug(name, 'post asynchronous receives, nrecv = ',&
            nrecv )
    endif

    if (nrecv.ge.1) then
       recvsize = xrecv( nrecv+1)-1

       allocate( crecvbuffer( recvsize*clen ) )

       if (lidebug.ge.1) then
          call write_debug( name, 'recvsize, ubound(crecvbuffer,1) ', &
               recvsize, ubound(crecvbuffer,1) )
          call write_debug( name, 'ubound(xrecv,1) ', &
               ubound(xrecv,1) )
          call write_debug( name, 'ubound(recvproc,1) ', &
               ubound(recvproc,1) )
          call write_debug( name, 'ubound(recvtag,1) ', &
               ubound(recvtag,1) )
       endif

       ! post receives
       datatype = MPI_CHARACTER
       comm = communicator

       do ii=1,nrecv
          j1 = xrecv(ii)
          j2 = xrecv(ii+1)-1

          count = j2-j1+1
          count = count*clen

          source = recvproc( ii )
          tag = recvtag( ii )

          if (lidebug.ge.2) then
             call write_debug(name, 'mpi_irecv: ii,j1,j2 ', ii,j1,j2 )
             call write_debug(name, 'count, source, tag ', &
                  count,source,tag )
          endif

          jpos = 1 + (j1-1)*clen
          call MPI_IRECV( crecvbuffer(jpos), count, datatype, source, &
               tag, comm, request, ierror )
          call MPI_Check( 'sendrecv_begin_1c:MPI_IRECV ', ierror )

          recvrequest( ii ) = request
       enddo   ! end do (ii=1,nrecv)
    endif   ! end if (nrecv.ge.1)
    ! ----------------------------------------------------------------<<<


    ! post asynchronous sends
    ! ---------------------------------------------------------------->>>
    if (lidebug.ge.1) then
       call write_debug(name, 'post asynchronous sends ')
    endif

    if (nsend.ge.1) then
       sendsize = xsend( nsend+1)-1

       allocate( csendbuffer( sendsize*clen ) )

       if (lidebug.ge.1) then
          call write_debug( name, 'sendsize, ubound(csendbuffer,1) ', &
               sendsize, ubound(csendbuffer,1) )
          call write_debug( name, 'ubound(xsend,1) ', &
               ubound(xsend,1) )
          call write_debug( name, 'ubound(sendproc,1) ', &
               ubound(sendproc,1) )
          call write_debug( name, 'ubound(sendtag,1) ', &
               ubound(sendtag,1) )
       endif

       ! perform sends
       datatype = MPI_CHARACTER
       comm = communicator

       do ii=1,nsend
          ! perform copy into sendbuffer
          j1 = xsend(ii)
          j2 = xsend(ii+1)-1

          count = j2-j1+1
          count = count*clen

          do jj=j1,j2
             ijk = sendijk( jj )
             do ic=1,clen
                jpos = (jj-1)*clen + ic
                csendbuffer(jpos) = XX(ijk)(ic:ic)
             enddo
          enddo

          dest = sendproc( ii )
          tag = sendtag( ii )

          if (lidebug.ge.2) then
             call write_debug(name, 'mpi_isend: ii,j1,j2 ', ii,j1,j2)
             call write_debug(name, 'count, dest, tag ', count, &
                  dest, tag )
          endif

          jpos = (j1-1)*clen + 1
          call MPI_ISEND( csendbuffer(jpos), count, datatype, dest, &
               tag, comm, request, ierror )
          call MPI_Check( 'sendrecv_begin_1c:MPI_ISEND ', ierror )
          sendrequest( ii ) = request
       enddo   ! end do (ii=1,nsend)
    endif   ! end if (nsend.ge.1)
    ! ----------------------------------------------------------------<<<
#endif

    return
  end subroutine sendrecv_begin_1c

  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine sendrecv_end_1d( XX, idebug )

    use functions

    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    double precision, intent(inout), dimension(:) :: XX
    integer, intent(in), optional :: idebug

#ifdef MPI
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    character(len=80), parameter :: name = 'sendrecv_end_1d'
    logical, parameter :: use_waitany = .false.
    integer :: lidebug
    integer :: jj, ijk, jindex, ii, j1, j2, ierror
    integer, dimension(MPI_STATUS_SIZE) :: recv_status_any
    integer, dimension(:,:), pointer :: recv_status
    integer, dimension(:,:), pointer :: send_status
    !-----------------------------------------------


    ! wait for sends to complete
    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    if (nsend.ge.1) then
       if (lidebug.ge.1) then
          call write_debug(name, &
               'waiting for sends to complete, nsend  = ', nsend )
       endif

       allocate( send_status(MPI_STATUS_SIZE,nsend))

       if (use_persistent_message) then
          call MPI_WAITALL( nsend, send_persistent_request, &
               send_status, ierror )
       else
          call MPI_WAITALL( nsend, sendrequest, send_status, ierror )
       endif

       call MPI_Check( 'sendrecv_end_1d:MPI_WAITALL ', ierror )

       deallocate( send_status )
       nullify( send_status )
    endif   ! end if (nsend.ge.1)


    ! wait for recvs to complete
    if (nrecv.ge.1) then
       if (lidebug.ge.1) then
          call write_debug( name, &
               'waiting for receives to complete, nrecv =  ', nrecv )
       endif

       if (use_waitany) then
          do ii=1,nrecv
             if (use_persistent_message) then
                call MPI_WAITANY( nrecv, recv_persistent_request, &
                     jindex, recv_status_any, ierror )
             else
                call MPI_WAITANY( nrecv, recvrequest,   &
                     jindex, recv_status_any, ierror )
             endif

             call MPI_Check( 'sendrecv_end_1d:MPI_WAITANY ', ierror )

             j1 = xrecv( jindex )
             j2 = xrecv( jindex + 1)-1

             if (lidebug.ge.2) then
                call write_debug(name, 'jindex, j1,j2 ', jindex,j1,j2)
             endif

             do jj=j1,j2
                ijk = recvijk( jj )
                XX(ijk) = drecvbuffer(jj)
             enddo
          enddo   ! end do (ii=nrecv)

       else

          allocate( recv_status(MPI_STATUS_SIZE,nrecv))
          if (use_persistent_message) then
             call MPI_WAITALL( nrecv, recv_persistent_request, &
                  recv_status, ierror )
          else
             call MPI_WAITALL( nrecv, recvrequest, recv_status,&
                  ierror )
          endif

          call MPI_Check( 'sendrecv_end_1d:MPI_WAITALL recv ', &
               ierror )
          deallocate( recv_status )
          nullify( recv_status )

          j1 = xrecv(1)
          j2 = xrecv( nrecv +1)-1
          do jj=j1,j2
             ijk = recvijk( jj )
             XX(ijk) = drecvbuffer(jj)
          enddo
       endif   ! end if/else (use_waitany)
    endif   ! end if (nrecv.ge.1)
#endif

    return
  end subroutine sendrecv_end_1d

  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine sendrecv_end_1c( XX, idebug )

    use functions

    implicit none

    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    character(len=*), intent(inout), dimension(:) :: XX
    integer, intent(in), optional :: idebug

#ifdef MPI
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    character(len=80), parameter :: name = 'sendrecv_end_1c'
    integer :: ic, clen, jpos
    logical, parameter :: use_waitany = .false.
    integer :: lidebug
    integer :: jj, ijk, jindex, ii, j1, j2, ierror
    integer, dimension(MPI_STATUS_SIZE) :: recv_status_any
    integer, dimension(:,:), pointer :: recv_status
    integer, dimension(:,:), pointer :: send_status
    !-----------------------------------------------


    ! wait for sends to complete
    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    jpos = lbound(XX,1)
    clen = len(XX(jpos))

    if (nsend.ge.1) then
       if (lidebug.ge.1) then
          call write_debug(name, &
               'waiting for sends to complete, nsend  = ', nsend )
       endif

       allocate( send_status(MPI_STATUS_SIZE,nsend))

       call MPI_WAITALL( nsend, sendrequest, send_status, ierror )
       call MPI_Check( 'sendrecv_end_1c:MPI_WAITALL ', ierror )

       deallocate( send_status )
       nullify( send_status )

       deallocate( csendbuffer )
       nullify( csendbuffer )

    endif   ! end if (nsend.ge.1)


    ! wait for recvs to complete
    if (nrecv.ge.1) then
       if (lidebug.ge.1) then
          call write_debug( name, &
               'waiting for receives to complete, nrecv =  ', nrecv )
       endif

       if (use_waitany) then
          do ii=1,nrecv
             call MPI_WAITANY( nrecv, recvrequest, jindex, &
                  recv_status_any, ierror )
             call MPI_Check( 'sendrecv_end_1c:MPI_WAITANY ', ierror )

             j1 = xrecv( jindex )
             j2 = xrecv( jindex + 1)-1

             if (lidebug.ge.2) then
                call write_debug(name, 'jindex, j1,j2 ', jindex,j1,j2 )
             endif

             do jj=j1,j2
                ijk = recvijk( jj )

                do ic=1,clen
                   jpos = (jj-1)*clen + ic
                   XX(ijk)(ic:ic) = crecvbuffer(jpos)
                enddo
             enddo
          enddo   ! end do (ii=1,nrecv)

       else

          allocate( recv_status(MPI_STATUS_SIZE,nrecv))
          call MPI_WAITALL( nrecv, recvrequest, recv_status, ierror )
          call MPI_Check( 'sendrecv_end_1c:MPI_WAITALL recv ', ierror )

          deallocate( recv_status )
          nullify( recv_status )

          j1 = xrecv(1)
          j2 = xrecv( nrecv +1)-1
          do jj=j1,j2
             ijk = recvijk( jj )
             do ic=1,clen
                jpos = (jj-1)*clen + ic
                XX(ijk)(ic:ic) = crecvbuffer(jpos)
             enddo
          enddo
       endif   ! end if/else (use_waitany)

       deallocate( crecvbuffer )
       nullify( crecvbuffer )

    endif   ! end if (nrecv.ge.1)
#endif

    return
  end subroutine sendrecv_end_1c

  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine sendrecv_end_1i( XX, idebug )

    use functions

    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    integer, intent(inout), dimension(:) :: XX
    integer, intent(in), optional :: idebug
    !-----------------------------------------------

#ifdef MPI
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    character(len=80), parameter :: name = 'sendrecv_end_1i'
    logical, parameter :: use_waitany = .false.
    integer :: lidebug
    integer :: jj, ijk, jindex, ii, j1, j2, ierror
    integer, dimension(MPI_STATUS_SIZE) :: recv_status_any
    integer, dimension(:,:), pointer :: recv_status
    integer, dimension(:,:), pointer :: send_status
    !-----------------------------------------------

    ! wait for sends to complete
    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    if (nsend.ge.1) then
       if (lidebug.ge.1) then
          call write_debug(name, &
               'waiting for sends to complete, nsend  = ', nsend )
       endif

       allocate( send_status(MPI_STATUS_SIZE,nsend))

       call MPI_WAITALL( nsend, sendrequest, send_status, ierror )
       call MPI_Check( 'sendrecv_end_1i:MPI_WAITALL ', ierror )

       deallocate( send_status )
       nullify( send_status )

       deallocate( isendbuffer )
       nullify( isendbuffer )
    endif   ! end if (nsend.ge.1)

    ! wait for recvs to complete
    if (nrecv.ge.1) then
       if (lidebug.ge.1) then
          call write_debug( name, &
               'waiting for receives to complete, nrecv =  ', nrecv )
       endif

       if (use_waitany) then
          do ii=1,nrecv
             call MPI_WAITANY( nrecv, recvrequest, jindex, &
                  recv_status_any, ierror )
             call MPI_Check( 'sendrecv_end_1i:MPI_WAITANY ', ierror )

             j1 = xrecv( jindex )
             j2 = xrecv( jindex + 1)-1

             if (lidebug.ge.2) then
                call write_debug(name, 'jindex, j1,j2 ', jindex,j1,j2 )
             endif

             do jj=j1,j2
                ijk = recvijk( jj )
                XX(ijk) = irecvbuffer(jj)
             enddo
          enddo    ! end do (ii=1,nrecv)
       else
          allocate( recv_status(MPI_STATUS_SIZE,nrecv))
          call MPI_WAITALL( nrecv, recvrequest, recv_status, ierror )
          call MPI_Check( 'sendrecv_end_1i:MPI_WAITALL recv ', ierror )
          deallocate( recv_status )
          nullify( recv_status )

          j1 = xrecv(1)
          j2 = xrecv( nrecv +1)-1
          do jj=j1,j2
             ijk = recvijk( jj )
             XX(ijk) = irecvbuffer(jj)
          enddo
       endif   ! end if/else (use_waitany)

       deallocate( irecvbuffer )
       nullify( irecvbuffer )

    endif   ! end if (nrecv.ge.1)
#endif

    return
  end subroutine sendrecv_end_1i



  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine send_recv_1c( XX, ilayer, idebug )
    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    character(len=*),  dimension(:), intent(inout) :: XX
    integer, intent(in), optional :: ilayer,idebug
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    integer :: lidebug, layer
    !-----------------------------------------------

#ifdef MPI
    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    layer = 1
    if (present(ilayer)) then
       layer = ilayer
    endif

    call sendrecv_begin(XX,layer,lidebug)
    call sendrecv_end( XX, lidebug )
#endif

    return
  end subroutine send_recv_1c


  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine send_recv_1d( XX, ilayer, idebug )
    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    double precision,  dimension(:), intent(inout) :: XX
    integer, intent(in), optional :: ilayer,idebug
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    integer :: lidebug, layer
    !-----------------------------------------------

#ifdef MPI
    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    layer = 1
    if (present(ilayer)) then
       layer = ilayer
    endif

    call sendrecv_begin(XX,layer,lidebug)
    call sendrecv_end( XX, lidebug )
#endif

    return
  end subroutine send_recv_1d

  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine send_recv_2d( XX, ilayer, idebug )
    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    double precision,  dimension(:,:), intent(inout) :: XX
    integer, intent(in), optional :: ilayer,idebug
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    integer :: lidebug, layer
    integer :: j
    !-----------------------------------------------

#ifdef MPI
    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    layer = 1
    if (present(ilayer)) then
       layer = ilayer
    endif

    do j=lbound(XX,2),ubound(XX,2)
       call sendrecv_begin(XX(:,j),layer,lidebug)
       call sendrecv_end( XX(:,j), lidebug )
    enddo
#endif

    return
  end subroutine send_recv_2d

  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine send_recv_3d( XX, ilayer, idebug )
    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    double precision,  dimension(:,:,:), intent(inout) :: XX
    integer, intent(in), optional :: ilayer,idebug
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    integer :: lidebug, layer
    integer :: j,k
    !-----------------------------------------------
#ifdef MPI
    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    layer = 1
    if (present(ilayer)) then
       layer = ilayer
    endif

    do k=lbound(XX,3),ubound(XX,3)
       do j=lbound(XX,2),ubound(XX,2)
          call sendrecv_begin(XX(:,j,k),layer,lidebug)
          call sendrecv_end( XX(:,j,k), lidebug )
       enddo
    enddo
#endif

    return
  end subroutine send_recv_3d


  !--------------------------------------------------------------------
  ! Purpose:
  !
  !--------------------------------------------------------------------
  subroutine send_recv_1i( XX, ilayer, idebug )
    implicit none
    !-----------------------------------------------
    ! Dummy arguments
    !-----------------------------------------------
    integer,  dimension(:), intent(inout) :: XX
    integer, intent(in), optional :: ilayer,idebug
    !-----------------------------------------------
    ! Local variables
    !-----------------------------------------------
    integer :: lidebug, layer
    !-----------------------------------------------
#ifdef MPI
    lidebug = 0
    if (present(idebug)) then
       lidebug = idebug
    endif

    layer = 1
    if (present(ilayer)) then
       layer = ilayer
    endif

    call sendrecv_begin(XX,layer,lidebug)
    call sendrecv_end( XX, lidebug )
#endif

    return
  end subroutine send_recv_1i



  ! Re-initialize send/receive after re-indexing
  subroutine sendrecv_re_init_after_re_indexing(comm, idebug )

    use functions

    implicit none

    integer, intent(in) :: comm

    integer, intent(in), optional :: idebug

    !       -------------------------------------
    !       set up tables and data structures for
    !       exchanging ghost regions
    !       -------------------------------------

    !       ---------------
    !       local variables
    !       ---------------
    character(len=80), parameter :: name = 'sendrecv_init'

    character(len=80), pointer, dimension(:) :: line
    integer :: ip, lmax

    integer :: layer,request, source, tag, datatype

    integer :: lidebug
    integer :: isize,jsize,ksize, ijksize
    integer :: recvsize1, recvsize2, &
         sendsize1, sendsize2

    integer :: iter, i,j,k, ii, jj,kk, &
         ntotal, icount,ipos, &
         ilayer,        i1,i2,  j1,j2, k1,k2,  &
         ijk, ijk2, iproc, jproc, src,dest, &
         ierror

    logical :: isok, isvalid, ismine, is_halobc

    integer, dimension(:,:,:), pointer :: ijk2proc
    integer, pointer, dimension(:) :: &
         istartx,iendx, jstartx,jendx, kstartx,kendx, &
         ncount, &
         recvproc, recvtag, xrecv, recvijk,  &
         sendproc, sendtag, xsend, sendijk

    logical, parameter :: jfastest = .true.


    integer, parameter :: message_tag_offset = 11


    !       ----------------
    !       inline functions
    !       ----------------
    integer :: message_tag

#ifdef MPI
    !  NEW SEND_RECV INIT HERE
    if (use_persistent_message) then

       datatype = MPI_DOUBLE_PRECISION

       do layer=1,2


          if (layer.eq.1) then
             nrecv = nrecv1
             recvtag =>recvtag1
             recvproc => recvproc1
             recvijk => recvijk1
             xrecv => xrecv1

             nsend = nsend1
             sendtag => sendtag1
             sendproc => sendproc1
             sendijk => sendijk1
             xsend => xsend1

             send_persistent_request => send_persistent_request1
             recv_persistent_request => recv_persistent_request1

          else
             nrecv = nrecv2
             recvtag =>recvtag2
             recvproc => recvproc2
             recvijk => recvijk2
             xrecv => xrecv2

             nsend = nsend2
             sendtag => sendtag2
             sendproc => sendproc2
             sendijk => sendijk2
             xsend => xsend2

             send_persistent_request => send_persistent_request2
             recv_persistent_request => recv_persistent_request2

          endif



          do ii=1,nrecv
             j1 = xrecv(ii)
             j2 = xrecv(ii+1)-1
             icount = j2-j1+1
             source = recvproc( ii )
             tag = recvtag( ii )



             if (lidebug.ge.2) then

                !                  call write_debug(name, 'mpi_recv_init: ii,j1,j2 ', &
                !                                             ii,j1,j2 )
                !                  call write_debug(name, 'icount, source, tag ', &
                !                                             icount,source,tag )
             endif


             call MPI_RECV_INIT( drecvbuffer(j1), icount, datatype, &
                  source, tag, comm, request, ierror )
             call MPI_Check( 'sendrecv_begin_1d:MPI_IRECV ', ierror )

             recv_persistent_request(ii) = request
          enddo


          do ii=1,nsend
             j1 = xsend(ii)
             j2 = xsend(ii+1)-1
             dest = sendproc( ii )
             tag = sendtag( ii )
             icount = j2-j1+1

             if (lidebug.ge.2) then

                !                  call write_debug(name, 'mpi_send_init: ii,j1,j2 ', &
                !                                         ii,j1,j2)
                !                  call write_debug(name, 'icount, dest, tag ', &
                !                                         icount,dest,tag )
             endif


             call MPI_SEND_INIT( dsendbuffer(j1), icount, datatype, &
                  dest, tag, &
                  comm, request, ierror )
             call MPI_Check( 'sendrecv_begin_1d:MPI_SEND_INIT ', ierror )

             send_persistent_request( ii ) = request
          enddo

       enddo  ! layers

    endif  ! use_persistent_message
#endif

    return
  end subroutine sendrecv_re_init_after_re_indexing


end module sendrecv
