!       MPI Modules written at ORNL by Ed and Sreekanth for MFIX
!       under joint effort with FETC - 06/08/99.

#include "version.inc"

module mpi_utility

  !       module to perform most of the mpi functionalities like scatter,
  !       gather, bcast, globalsum and so on.

  use geometry
  use compar
  use parallel_mpi
  use debug
  use indices
  implicit none

  !       Object-oriented approach to direct to the correct procedure
  !       depending on the argument type. i stands for integer, r for real
  !       and d for double precision. 0 for scalar, 1 for vector, 2 for
  !       2-D array and similarly 3.

  !==============================================================================
  !  JFD: Interfaces used for vtk file writting (Cartesian grid):
  !==============================================================================

  interface allgather
     module  procedure allgather_1i
  end interface allgather

  interface gatherv
     module  procedure gatherv_1d
  end interface gatherv

  !==============================================================================
  !  JFD: End of Interfaces used for vtk file writting (Cartesian grid):
  !==============================================================================

  interface scatter
     module  procedure scatter_1i, scatter_2i, scatter_3i, &
          scatter_1r, scatter_2r, scatter_3r, &
          scatter_1d, scatter_2d, scatter_3d, &
          scatter_1c,scatter_1l
  end interface scatter

  interface gather
     module  procedure gather_1i, gather_2i, gather_3i, &
          gather_1r, gather_2r, gather_3r, &
          gather_1d, gather_2d, gather_3d, &
          gather_1c, gather_1l
  end interface gather

  interface bcast
     module  procedure bcast_0i, bcast_1i, bcast_2i, bcast_3i, &
          bcast_0r, bcast_1r, bcast_2r, bcast_3r, &
          bcast_0d, bcast_1d, bcast_2d, bcast_3d, &
          bcast_0l, bcast_1l, bcast_0c, bcast_1c
  end interface bcast

  interface global_sum
     module  procedure global_sum_0i, global_sum_1i, global_sum_2i, global_sum_3i, &
          global_sum_0r, global_sum_1r, global_sum_2r, global_sum_3r, &
          global_sum_0d, global_sum_1d, global_sum_2d, global_sum_3d
  end interface global_sum

  interface global_all_sum
     module  procedure &
          global_all_sum_0i, global_all_sum_1i, &
          global_all_sum_2i, global_all_sum_3i, &
          global_all_sum_0r, global_all_sum_1r, &
          global_all_sum_2r, global_all_sum_3r, &
          global_all_sum_0d, global_all_sum_1d, &
          global_all_sum_2d, global_all_sum_3d, &
          global_all_sum_onevar_0i, global_all_sum_onevar_1i, &
          global_all_sum_onevar_2i, global_all_sum_onevar_3i, &
          global_all_sum_onevar_0r, global_all_sum_onevar_1r, &
          global_all_sum_onevar_2r, global_all_sum_onevar_3r, &
          global_all_sum_onevar_0d, global_all_sum_onevar_1d, &
          global_all_sum_onevar_2d, global_all_sum_onevar_3d
  end interface global_all_sum

  interface global_min
     module  procedure global_min_0i, global_min_1i, global_min_2i, global_min_3i, &
          global_min_0r, global_min_1r, global_min_2r, global_min_3r, &
          global_min_0d, global_min_1d, global_min_2d, global_min_3d
  end interface global_min

  interface global_all_min
     module  procedure &
          global_all_min_0i, global_all_min_1i, &
          global_all_min_2i, global_all_min_3i, &
          global_all_min_0r, global_all_min_1r, &
          global_all_min_2r, global_all_min_3r, &
          global_all_min_0d, global_all_min_1d, &
          global_all_min_2d, global_all_min_3d, &
          global_all_min_onevar_0i, global_all_min_onevar_1i, &
          global_all_min_onevar_2i, global_all_min_onevar_3i, &
          global_all_min_onevar_0r, global_all_min_onevar_1r, &
          global_all_min_onevar_2r, global_all_min_onevar_3r, &
          global_all_min_onevar_0d, global_all_min_onevar_1d, &
          global_all_min_onevar_2d, global_all_min_onevar_3d
  end interface global_all_min

  interface global_max
     module  procedure global_max_0i, global_max_1i, global_max_2i, global_max_3i, &
          global_max_0r, global_max_1r, global_max_2r, global_max_3r, &
          global_max_0d, global_max_1d, global_max_2d, global_max_3d
  end interface global_max

  interface global_all_max
     module  procedure &
          global_all_max_0i, global_all_max_1i, &
          global_all_max_2i, global_all_max_3i, &
          global_all_max_0r, global_all_max_1r, &
          global_all_max_2r, global_all_max_3r, &
          global_all_max_0d, global_all_max_1d, &
          global_all_max_2d, global_all_max_3d, &
          global_all_max_onevar_0i, global_all_max_onevar_1i, &
          global_all_max_onevar_2i, global_all_max_onevar_3i, &
          global_all_max_onevar_0r, global_all_max_onevar_1r, &
          global_all_max_onevar_2r, global_all_max_onevar_3r, &
          global_all_max_onevar_0d, global_all_max_onevar_1d, &
          global_all_max_onevar_2d, global_all_max_onevar_3d
  end interface global_all_max

  interface global_all_and
     module procedure &
          global_all_and_0d, global_all_and_1d, &
          global_all_and_onevar_0d, global_all_and_onevar_1d
  end interface global_all_and

  interface global_all_or
     module procedure &
          global_all_or_0d, global_all_or_1d, &
          global_all_or_onevar_0d, global_all_or_onevar_1d
  end interface global_all_or

contains


  !==============================================================================
  !  JFD: Subroutines used for vtk file writting (Cartesian grid):
  !==============================================================================

  subroutine allgather_1i( lbuf, gbuf, idebug )
    integer, intent(in) :: lbuf
    integer, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) ::  idebug

#ifdef MPI
    integer :: sendtype,recvtype,sendcnt,recvcnt,ierr,lidebug,mpierr

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = 1
    recvcnt = sendcnt

    CALL MPI_ALLGATHER(lbuf,sendcnt,sendtype,  &
         gbuf,recvcnt,recvtype,MPI_COMM_WORLD, IERR)
#else
    gbuf = 0
#endif

    return
  end subroutine allgather_1i

  subroutine allgather_1d( lbuf, gbuf, idebug )
    double precision, intent(in) :: lbuf
    double precision, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) ::  idebug

#ifdef MPI
    integer :: sendtype,recvtype,sendcnt,recvcnt,ierr,lidebug,mpierr

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = MPI_DOUBLE_PRECISION

    sendcnt = 1
    recvcnt = sendcnt

    CALL MPI_ALLGATHER(lbuf,sendcnt,sendtype,  &
         gbuf,recvcnt,recvtype,MPI_COMM_WORLD, IERR)
#else
    gbuf = 0
#endif

    return
  end subroutine allgather_1d

  subroutine gatherv_1i( lbuf, sendcnt, gbuf, rcount, disp, mroot, idebug )
    integer, intent(in), dimension(:) :: lbuf
    integer, intent(in), dimension(:) :: rcount
    integer, intent(in), dimension(:) :: disp
    integer, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug
    integer :: sendtype,recvtype,sendcnt,recvcnt,lroot,ierr,lidebug

#ifdef MPI
    !       check to see whether there is root

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = MPI_INTEGER

    CALL MPI_GATHERV(lbuf,sendcnt,sendtype,  &
         gbuf,rcount,disp,recvtype, &
         lroot,MPI_COMM_WORLD, IERR)
#else
    gbuf = lbuf
#endif

    return
  end subroutine gatherv_1i

  subroutine gatherv_1d( lbuf, sendcnt, gbuf, rcount, disp, mroot, idebug )
    double precision, intent(in), dimension(:) :: lbuf
    integer, intent(in), dimension(:) :: rcount
    integer, intent(in), dimension(:) :: disp
    double precision, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug
    integer :: sendtype,recvtype,sendcnt,recvcnt,lroot,ierr,lidebug

#ifdef MPI
    !       check to see whether there is root

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = MPI_DOUBLE_PRECISION

    CALL MPI_GATHERV(lbuf,sendcnt,sendtype,  &
         gbuf,rcount,disp,recvtype, &
         lroot,MPI_COMM_WORLD, IERR)
#else
    gbuf = lbuf
#endif

    return
  end subroutine gatherv_1d

  !==============================================================================
  !  JFD: End of Subroutines used for vtk file writting (Cartesian grid):
  !==============================================================================


  !       Routine to scatter gbuf available on root to all the processors

  subroutine scatter_1i( lbuf, gbuf, mroot, idebug )

    use functions

    implicit none

    integer, intent(in), dimension(:) :: gbuf
    integer, intent(out), dimension(:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

    integer, allocatable, dimension(:) :: gbuf_pack

    integer :: sendtype, recvtype, ijk1, ijk2, recvcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk

#ifdef MPI
    !       check to see whether there is root

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(sum(ijksize3_all(:))))
    else
       allocate(gbuf_pack(10))
    endif

    if( myPE.eq.lroot) then
       ioffset = 0
       do iproc = 0,numPEs-1
          ibuffer = 0
          do k = kstart3_all(iproc), kend3_all(iproc)
             do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                   ibuffer = funijk_proc(i,j,k,iproc) + ioffset
                   gbuf_pack(ibuffer) = gbuf(funijk_gl(i,j,k))

                enddo
             enddo
          enddo
          ioffset = ibuffer
       enddo
    endif

    sendtype = MPI_INTEGER
    recvtype = sendtype

    ijk1 = ijkstart3
    ijk2 = ijkend3

    recvcnt = ijk2-ijk1+1

    !       Call MPI routines

    call MPI_Scatterv( gbuf_pack, ijksize3_all, displs, sendtype, &
         lbuf, recvcnt, recvtype,  &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'scatter_1i:MPI_Scatterv', ierr )

    deallocate(gbuf_pack)
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_1i

  subroutine scatter_2i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:) :: gbuf
    integer, intent(out), dimension(:,:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** scatter_2i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call scatter_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_2i

  subroutine scatter_3i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:,:) :: gbuf
    integer, intent(out), dimension(:,:,:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** scatter_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** scatter_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call scatter_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_3i

  subroutine scatter_1r( lbuf, gbuf, mroot, idebug )

    use functions

    implicit none

    real, intent(in), dimension(:) :: gbuf
    real, intent(out), dimension(:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    real, allocatable, dimension(:) :: gbuf_pack

    integer :: sendtype, recvtype, ijk1, ijk2, recvcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(sum(ijksize3_all(:))))
    else
       allocate(gbuf_pack(10))
    endif

    if( myPE.eq.lroot) then
       ioffset = 0
       do iproc = 0,numPEs-1
          ibuffer = 0
          do k = kstart3_all(iproc), kend3_all(iproc)
             do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                   ibuffer = funijk_proc(i,j,k,iproc) + ioffset
                   gbuf_pack(ibuffer) = gbuf(funijk_gl(i,j,k))

                enddo
             enddo
          enddo
          ioffset = ibuffer
       enddo
    endif

    sendtype = MPI_REAL
    recvtype = sendtype

    ijk1 = ijkstart3
    ijk2 = ijkend3

    recvcnt = ijk2-ijk1+1

    call MPI_Scatterv( gbuf_pack, ijksize3_all, displs, sendtype, &
         lbuf, recvcnt, recvtype,  &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'scatter_1r:MPI_Scatterv', ierr )

    deallocate(gbuf_pack)
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_1r


  subroutine scatter_2r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:) :: gbuf
    real, intent(out), dimension(:,:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** scatter_2r: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call scatter_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_2r

  subroutine scatter_3r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:,:) :: gbuf
    real, intent(out), dimension(:,:,:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** scatter_3r: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** scatter_3r: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call scatter_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_3r


  subroutine scatter_1d( lbuf, gbuf, mroot, idebug )

    use functions
    implicit none

    double precision, intent(in), dimension(:) :: gbuf
    double precision, intent(out), dimension(:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    double precision, allocatable, dimension(:) :: gbuf_pack

    integer :: sendtype, recvtype, ijk1,ijk2,recvcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(sum(ijksize3_all(:))))
    else
       allocate(gbuf_pack(10))
    endif

    if( myPE.eq.lroot) then
       ioffset = 0
       do iproc = 0,numPEs-1
          ibuffer = 0
          do k = kstart3_all(iproc), kend3_all(iproc)
             do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                   ibuffer = funijk_proc(i,j,k,iproc) + ioffset
                   gbuf_pack(ibuffer) = gbuf(funijk_gl(i,j,k))

                enddo
             enddo
          enddo
          ioffset = ibuffer
       enddo
    endif

    sendtype = MPI_DOUBLE_PRECISION
    recvtype = sendtype

    ijk1 = ijkstart3
    ijk2 = ijkend3

    recvcnt = ijk2-ijk1+1

    call MPI_Scatterv( gbuf_pack, ijksize3_all, displs, sendtype, &
         lbuf, recvcnt, recvtype,  &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'scatter_1d:MPI_Scatterv', ierr )

    deallocate(gbuf_pack)
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_1d


  subroutine scatter_2d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:) :: gbuf
    double precision, intent(out), dimension(:,:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** scatter_2d: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call scatter_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_2d

  subroutine scatter_3d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:,:) :: gbuf
    double precision, intent(out), dimension(:,:,:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** scatter_3d: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** scatter_3d: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call scatter_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_3d

  subroutine scatter_1c( lbuf, gbuf, mroot, idebug )

    use functions
    implicit none

    character(len=*), intent(in), dimension(:) :: gbuf
    character(len=*), intent(out), dimension(:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer, allocatable, dimension(:,:) :: gbuf_pack,lbuf1
    character(len=80) :: string

    integer :: sendtype, recvtype, ijk1, ijk2, recvcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk
    integer :: lenchar, icount

    !       check to see whether there is root

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    ijk1 = ijkstart3
    ijk2 = ijkend3

    lenchar = len(gbuf(1))

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(ijkmax3,lenchar))
    else
       allocate(gbuf_pack(10,lenchar))
    endif

    allocate(lbuf1(ijk1:ijk2,lenchar))

    if(myPE.eq.lroot) then
       do i = 1,ijkmax3
          do j = 1,lenchar

             string = gbuf(i)(1:lenchar)
             gbuf_pack(i,j) = ichar(string(j:j))

          enddo
       enddo
    endif

    call scatter_2i(lbuf1,gbuf_pack)

    do i = ijk1, ijk2
       do j = 1,lenchar

          lbuf(i)(j:j) = char(lbuf1(i,j))

       enddo
    enddo

    deallocate(gbuf_pack)
    deallocate(lbuf1)
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_1c


  subroutine scatter_1l( lbuf, gbuf, mroot, idebug )

    use functions
    implicit none

    logical, intent(in), dimension(:) :: gbuf
    logical, intent(out), dimension(:) :: lbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    logical, allocatable, dimension(:) :: gbuf_pack

    integer :: sendtype, recvtype, ijk1, ijk2, recvcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk

    !       check to see whether there is root

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(sum(ijksize3_all(:))))
    else
       allocate(gbuf_pack(10))
    endif

    if( myPE.eq.lroot) then
       ioffset = 0
       do iproc = 0,numPEs-1
          ibuffer = 0
          do k = kstart3_all(iproc), kend3_all(iproc)
             do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                   ibuffer = funijk_proc(i,j,k,iproc) + ioffset
                   gbuf_pack(ibuffer) = gbuf(funijk_gl(i,j,k))

                enddo
             enddo
          enddo
          ioffset = ibuffer
       enddo
    endif

    sendtype = MPI_LOGICAL
    recvtype = sendtype

    ijk1 = ijkstart3
    ijk2 = ijkend3

    recvcnt = ijk2-ijk1+1

    !       Call MPI routines

    call MPI_Scatterv( gbuf_pack, ijksize3_all, displs, sendtype, &
         lbuf, recvcnt, recvtype,  &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'scatter_1l:MPI_Scatterv', ierr )

    deallocate(gbuf_pack)
#else
    lbuf = gbuf
#endif

    return
  end subroutine scatter_1l


  !       Routines to gather lbuf from individual processors and put it on
  !       processor root in gbuf
  !       Logic is similar to the scatter routines above.

  subroutine gather_1i( lbuf, gbuf, mroot, idebug )

    use functions
    implicit none

    integer, intent(in), dimension(:) :: lbuf
    integer, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer, allocatable, dimension(:) :: gbuf_pack

    integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk, ijk_gl
    integer :: istartl, iendl, jstartl, jendl, kstartl, kendl
    logical :: isok_k,isok_j,isok_i, isinterior
    logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy

    !       check to see whether there is root

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(sum(ijksize3_all(:))))
    else
       allocate(gbuf_pack(10))
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    ijk1 = ijkstart3
    !        ijk2 = ijkend3
    ijk2 = max(ijkend3,BACKGROUND_IJKEND3)   !  For cell re-indexing

    sendcnt = ijk2-ijk1+1

    call MPI_Gatherv( lbuf, sendcnt, sendtype,  &
         gbuf_pack, ijksize3_all, displs, recvtype, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'gather_1i:MPI_Gatherv', ierr )

    if( myPE.eq.lroot) then
       ioffset = 0
       do iproc = 0,numPEs-1
          ibuffer = 0
          istartl = istart1_all(iproc)
          iendl = iend1_all(iproc)
          jstartl = jstart1_all(iproc)
          jendl = jend1_all(iproc)
          kstartl = kstart1_all(iproc)
          kendl = kend1_all(iproc)

          if(istart3_all(iproc).eq.imin3) istartl = istart3_all(iproc)
          if(iend3_all(iproc).eq.imax3) iendl = iend3_all(iproc)
          if(jstart3_all(iproc).eq.jmin3) jstartl = jstart3_all(iproc)
          if(jend3_all(iproc).eq.jmax3) jendl = jend3_all(iproc)
          if(kstart3_all(iproc).eq.kmin3) kstartl = kstart3_all(iproc)
          if(kend3_all(iproc).eq.kmax3) kendl = kend3_all(iproc)

          do k = kstart3_all(iproc), kend3_all(iproc)
             do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                   ibuffer = funijk_proc(i,j,k,iproc) + ioffset
                   isok_k = (kstartl <= k) .and. (k <=kendl)
                   isok_j = (jstartl <= j) .and. (j <=jendl)
                   isok_i = (istartl <= i) .and. (i <=iendl)

                   need_copy = isok_k .and. isok_j .and. isok_i

                   if (need_copy) then
                      ijk_gl = funijk_gl(i,j,k)
                      gbuf( ijk_gl ) = gbuf_pack(ibuffer)
                   endif

                enddo
             enddo
          enddo
          ioffset = ibuffer
       enddo
    endif

    deallocate(gbuf_pack)
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_1i


  subroutine gather_2i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:) :: lbuf
    integer, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** gather_2i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call gather_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_2i

  subroutine gather_3i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:,:) :: lbuf
    integer, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** gather_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** gather_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call gather_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_3i

  subroutine gather_1r( lbuf, gbuf, mroot, idebug )

    use functions
    implicit none

    real, intent(in), dimension(:) :: lbuf
    real, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    real, allocatable, dimension(:) :: gbuf_pack

    integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk, ijk_gl
    integer :: istartl, iendl, jstartl, jendl, kstartl, kendl
    logical :: isok_k,isok_j,isok_i, isinterior
    logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(sum(ijksize3_all(:))))
    else
       allocate(gbuf_pack(10))
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    ijk1 = ijkstart3
    !        ijk2 = ijkend3
    ijk2 = max(ijkend3,BACKGROUND_IJKEND3)   !  For cell re-indexing

    sendcnt = ijk2-ijk1+1

    call MPI_Gatherv( lbuf, sendcnt, sendtype,  &
         gbuf_pack, ijksize3_all, displs, recvtype, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'gather_1r:MPI_Gatherv', ierr )

    if( myPE.eq.lroot) then
       ioffset = 0
       do iproc = 0,numPEs-1
          ibuffer = 0
          istartl = istart1_all(iproc)
          iendl = iend1_all(iproc)
          jstartl = jstart1_all(iproc)
          jendl = jend1_all(iproc)
          kstartl = kstart1_all(iproc)
          kendl = kend1_all(iproc)

          if(istart3_all(iproc).eq.imin3) istartl = istart3_all(iproc)
          if(iend3_all(iproc).eq.imax3) iendl = iend3_all(iproc)
          if(jstart3_all(iproc).eq.jmin3) jstartl = jstart3_all(iproc)
          if(jend3_all(iproc).eq.jmax3) jendl = jend3_all(iproc)
          if(kstart3_all(iproc).eq.kmin3) kstartl = kstart3_all(iproc)
          if(kend3_all(iproc).eq.kmax3) kendl = kend3_all(iproc)

          do k = kstart3_all(iproc), kend3_all(iproc)
             do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                   ibuffer = funijk_proc(i,j,k,iproc) + ioffset
                   isok_k = (kstartl <= k) .and. (k <=kendl)
                   isok_j = (jstartl <= j) .and. (j <=jendl)
                   isok_i = (istartl <= i) .and. (i <=iendl)

                   need_copy = isok_k .and. isok_j .and. isok_i

                   if (need_copy) then
                      ijk_gl = funijk_gl(i,j,k)
                      gbuf( ijk_gl ) = gbuf_pack(ibuffer)
                   endif

                enddo
             enddo
          enddo
          ioffset = ibuffer
       enddo
    endif

    deallocate(gbuf_pack)
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_1r

  subroutine gather_2r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:) :: lbuf
    real, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** gather_2r: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call gather_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_2r

  subroutine gather_3r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:,:) :: lbuf
    real, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** gather_3r: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** gather_3r: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call gather_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_3r


  subroutine gather_1d( lbuf, gbuf, mroot, idebug )

    use functions
    implicit none

    double precision, intent(in), dimension(:) :: lbuf
    double precision, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    double precision, allocatable, dimension(:) :: gbuf_pack

    integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk, ijk_gl
    logical :: isok_k,isok_j,isok_i, isinterior
    logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy
    integer :: istartl, iendl, jstartl, jendl, kstartl, kendl

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(sum(ijksize3_all(:))))
    else
       allocate(gbuf_pack(10))
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    ijk1 = ijkstart3
    !        ijk2 = ijkend3
    ijk2 = max(ijkend3,BACKGROUND_IJKEND3)   !  For cell re-indexing

    sendcnt = ijk2-ijk1+1

    call MPI_Gatherv( lbuf, sendcnt, sendtype,  &
         gbuf_pack, ijksize3_all, displs, recvtype, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'gather_1d:MPI_Gatherv', ierr )

    if( myPE.eq.lroot) then
       ioffset = 0
       do iproc = 0,numPEs-1
          ibuffer = 0
          istartl = istart1_all(iproc)
          iendl = iend1_all(iproc)
          jstartl = jstart1_all(iproc)
          jendl = jend1_all(iproc)
          kstartl = kstart1_all(iproc)
          kendl = kend1_all(iproc)

          if(istart3_all(iproc).eq.imin3) istartl = istart3_all(iproc)
          if(iend3_all(iproc).eq.imax3) iendl = iend3_all(iproc)
          if(jstart3_all(iproc).eq.jmin3) jstartl = jstart3_all(iproc)
          if(jend3_all(iproc).eq.jmax3) jendl = jend3_all(iproc)
          if(kstart3_all(iproc).eq.kmin3) kstartl = kstart3_all(iproc)
          if(kend3_all(iproc).eq.kmax3) kendl = kend3_all(iproc)

          do k = kstart3_all(iproc), kend3_all(iproc)
             do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                   ibuffer = funijk_proc(i,j,k,iproc) + ioffset
                   isok_k = (kstartl <= k) .and. (k <=kendl)
                   isok_j = (jstartl <= j) .and. (j <=jendl)
                   isok_i = (istartl <= i) .and. (i <=iendl)

                   need_copy = isok_k .and. isok_j .and. isok_i

                   if (need_copy) then
                      ijk_gl = funijk_gl(i,j,k)
                      gbuf( ijk_gl ) = gbuf_pack(ibuffer)
                   endif

                enddo
             enddo
          enddo
          ioffset = ibuffer
       enddo
    endif

    deallocate(gbuf_pack)
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_1d

  subroutine gather_2d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:) :: lbuf
    double precision, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** gather_2d: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call gather_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_2d

  subroutine gather_3d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:,:) :: lbuf
    double precision, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** gather_3d: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** gather_3d: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call gather_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_3d


  subroutine gather_1c( lbuf, gbuf, mroot, idebug )

    use functions
    implicit none

    character(len=*), intent(in), dimension(:) :: lbuf
    character(len=*), intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer, allocatable, dimension(:,:) :: gbuf_pack,lbuf1
    character(len=80) :: string

    integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk, ijk_gl
    integer :: istartl, iendl, jstartl, jendl, kstartl, kendl
    integer :: lenchar, icount
    logical :: isok_k,isok_j,isok_i, isinterior
    logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy

    !       check to see whether there is root

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif


    ijk1 = ijkstart3
    !        ijk2 = ijkend3
    ijk2 = max(ijkend3,BACKGROUND_IJKEND3)   !  For cell re-indexing

    lenchar = len(lbuf(1))

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(ijkmax3,lenchar))
    else
       allocate(gbuf_pack(10,lenchar))
    endif

    allocate(lbuf1(ijk1:ijk2,lenchar))

    do i = ijk1,ijk2
       string = lbuf(i)(1:lenchar)
       do j = 1,lenchar
          lbuf1(i,j) = ichar(string(j:j))
       enddo
    enddo

    call gather_2i(lbuf1, gbuf_pack)

    if(myPE.eq.lroot) then
       do i = 1,ijkmax3
          do j = 1,lenchar

             string(j:j) = char(gbuf_pack(i,j))

          enddo
          gbuf(i)(1:lenchar) = string(1:lenchar)

       enddo
    endif

    deallocate(gbuf_pack)
    deallocate(lbuf1)
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_1c


  subroutine gather_1l( lbuf, gbuf, mroot, idebug )

    use functions
    implicit none

    logical, intent(in), dimension(:) :: lbuf
    logical, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    logical, allocatable, dimension(:) :: gbuf_pack

    integer :: recvtype, sendtype, ijk1,ijk2,sendcnt, ierr,lroot, lidebug
    integer :: i,j,k,ibuffer,iproc, ioffset
    integer :: ijk, ijk_gl
    integer :: istartl, iendl, jstartl, jendl, kstartl, kendl
    logical :: isok_k,isok_j,isok_i, isinterior
    logical :: isbc_k,isbc_j,isbc_i, isboundary, need_copy

    !       check to see whether there is root

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       allocate(gbuf_pack(sum(ijksize3_all(:))))
    else
       allocate(gbuf_pack(10))
    endif

    recvtype = MPI_LOGICAL
    sendtype = recvtype

    ijk1 = ijkstart3
    !        ijk2 = ijkend3
    ijk2 = max(ijkend3,BACKGROUND_IJKEND3)   !  For cell re-indexing

    sendcnt = ijk2-ijk1+1

    call MPI_Gatherv( lbuf, sendcnt, sendtype,  &
         gbuf_pack, ijksize3_all, displs, recvtype, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'gather_1l:MPI_Gatherv', ierr )

    if( myPE.eq.lroot) then
       ioffset = 0
       do iproc = 0,numPEs-1
          ibuffer = 0
          istartl = istart1_all(iproc)
          iendl = iend1_all(iproc)
          jstartl = jstart1_all(iproc)
          jendl = jend1_all(iproc)
          kstartl = kstart1_all(iproc)
          kendl = kend1_all(iproc)

          if(istart3_all(iproc).eq.imin3) istartl = istart3_all(iproc)
          if(iend3_all(iproc).eq.imax3) iendl = iend3_all(iproc)
          if(jstart3_all(iproc).eq.jmin3) jstartl = jstart3_all(iproc)
          if(jend3_all(iproc).eq.jmax3) jendl = jend3_all(iproc)
          if(kstart3_all(iproc).eq.kmin3) kstartl = kstart3_all(iproc)
          if(kend3_all(iproc).eq.kmax3) kendl = kend3_all(iproc)

          do k = kstart3_all(iproc), kend3_all(iproc)
             do j = jstart3_all(iproc), jend3_all(iproc)
                do i = istart3_all(iproc), iend3_all(iproc)

                   ibuffer = funijk_proc(i,j,k,iproc) + ioffset
                   isok_k = (kstartl <= k) .and. (k <=kendl)
                   isok_j = (jstartl <= j) .and. (j <=jendl)
                   isok_i = (istartl <= i) .and. (i <=iendl)

                   need_copy = isok_k .and. isok_j .and. isok_i

                   if (need_copy) then
                      ijk_gl = funijk_gl(i,j,k)
                      gbuf( ijk_gl ) = gbuf_pack(ibuffer)
                   endif

                enddo
             enddo
          enddo
          ioffset = ibuffer
       enddo
    endif

    deallocate(gbuf_pack)
#else
    gbuf = lbuf
#endif

    return
  end subroutine gather_1l


  !       Routines to broadcast information from processor 0 in buffer to all
  !       the processors

  subroutine bcast_0i( buffer, mroot, idebug )
    integer, intent(inout) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    datatype = MPI_INTEGER

    count = 1

    call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_0i:MPI_Bcast', ierr )
#endif

    return
  end subroutine bcast_0i


  subroutine bcast_1i( buffer, mroot, idebug )
    integer, intent(inout), dimension(:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    datatype = MPI_INTEGER

    count = size(buffer,1)

    call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_1i:MPI_Bcast', ierr )
#endif

    return
  end subroutine bcast_1i


  subroutine bcast_2i( buffer, mroot, idebug )
    integer, intent(inout), dimension(:,:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    do j=lbound(buffer,2),ubound(buffer,2)
       call bcast_1i( buffer(:,j), lroot, lidebug )
    enddo
#endif

    return
  end subroutine bcast_2i

  subroutine bcast_3i( buffer, mroot, idebug )
    integer, intent(inout), dimension(:,:,:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    do k=lbound(buffer,3),ubound(buffer,3)
       do j=lbound(buffer,2),ubound(buffer,2)
          call bcast_1i( buffer(:,j,k), lroot, lidebug )
       enddo
    enddo
#endif

    return
  end subroutine bcast_3i

  subroutine bcast_0r( buffer, mroot, idebug )
    real, intent(inout) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    datatype = MPI_REAL

    count = 1

    call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_0r:MPI_Bcast', ierr )
#endif

    return
  end subroutine bcast_0r

  subroutine bcast_1r( buffer, mroot, idebug )
    real, intent(inout), dimension(:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    datatype = MPI_REAL

    count = size(buffer,1)

    call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_1r:MPI_Bcast', ierr )
#endif

    return
  end subroutine bcast_1r

  subroutine bcast_2r( buffer, mroot, idebug )
    real, intent(inout), dimension(:,:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    do j=lbound(buffer,2),ubound(buffer,2)
       call bcast_1r( buffer(:,j), lroot, lidebug )
    enddo
#endif

    return
  end subroutine bcast_2r

  subroutine bcast_3r( buffer, mroot, idebug )
    real, intent(inout), dimension(:,:,:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    do k=lbound(buffer,3),ubound(buffer,3)
       do j=lbound(buffer,2),ubound(buffer,2)
          call bcast_1r( buffer(:,j,k), lroot, lidebug )
       enddo
    enddo
#endif

    return
  end subroutine bcast_3r

  subroutine bcast_0d( buffer, mroot, idebug )
    double precision, intent(inout) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    datatype = MPI_DOUBLE_PRECISION

    count = 1

    call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_0d:MPI_Bcast', ierr )
#endif

    return
  end subroutine bcast_0d


  subroutine bcast_1d( buffer, mroot, idebug )
    double precision, intent(inout), dimension(:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    datatype = MPI_DOUBLE_PRECISION

    count = size(buffer,1)

    call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_1d:MPI_Bcast', ierr )
#endif

    return
  end subroutine bcast_1d

  subroutine bcast_2d( buffer, mroot, idebug )
    double precision, intent(inout), dimension(:,:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    do j=lbound(buffer,2),ubound(buffer,2)
       call bcast_1d( buffer(:,j), lroot, lidebug )
    enddo
#endif

    return
  end subroutine bcast_2d

  subroutine bcast_3d( buffer, mroot, idebug )
    double precision, intent(inout), dimension(:,:,:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    do k=lbound(buffer,3),ubound(buffer,3)
       do j=lbound(buffer,2),ubound(buffer,2)
          call bcast_1d( buffer(:,j,k), lroot, lidebug )
       enddo
    enddo
#endif

    return
  end subroutine bcast_3d

  subroutine bcast_0c( buffer, mroot, idebug )
    character(len=*), intent(inout) :: buffer
    integer, optional, intent(in) :: mroot, idebug
    character, allocatable, dimension(:) :: buffer1

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug
    integer :: lenchar,icount, i, j

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    lenchar = len(buffer)

    allocate(buffer1(lenchar))

    icount = 0
    do j = 1,lenchar

       icount = icount+1
       buffer1(icount) = buffer(j:j)

    enddo

    datatype = MPI_CHARACTER

    count = 1

    call MPI_Bcast( buffer1, count*lenchar, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_0c:MPI_Bcast', ierr )

    icount = 0
    do j = 1,lenchar

       icount = icount+1
       buffer(j:j) = buffer1(icount)

    enddo

    deallocate(buffer1)
#endif

    return
  end subroutine bcast_0c


  subroutine bcast_1c( buffer, mroot, idebug )
    character(len=*), intent(inout), dimension(:) :: buffer
    integer, optional, intent(in) :: mroot, idebug
    character, allocatable, dimension(:) :: buffer1

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug
    integer :: lenchar,icount, i, j
    character(len=len(buffer(1))) :: string

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    lenchar = len(buffer(1))

    allocate(buffer1(size(buffer)*lenchar))

    icount = 0
    do i = 1,size(buffer)
       string = buffer(i)(1:lenchar)
       do j = 1,lenchar

          icount = icount+1
          buffer1(icount) = string(j:j)

       enddo
    enddo

    datatype = MPI_CHARACTER

    count = size(buffer,1)

    call MPI_Bcast( buffer1, count*lenchar, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_1c:MPI_Bcast', ierr )

    icount = 0
    do i = 1,size(buffer)
       do j = 1,lenchar

          icount = icount+1
          string(j:j) = buffer1(icount)

       enddo
       buffer(i) = string
    enddo

    deallocate(buffer1)
#endif

    return
  end subroutine bcast_1c

  subroutine bcast_0l( buffer, mroot, idebug )
    logical, intent(inout) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    datatype = MPI_LOGICAL

    count = 1

    call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_0l:MPI_Bcast', ierr )
#endif

    return
  end subroutine bcast_0l


  subroutine bcast_1l( buffer, mroot, idebug )
    logical, intent(inout), dimension(:) :: buffer
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: datatype, count, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    datatype = MPI_LOGICAL

    count = size(buffer,1)

    call MPI_Bcast( buffer, count, datatype, lroot, MPI_COMM_WORLD, ierr)
    call MPI_Check( 'bcast_1l:MPI_Bcast', ierr )
#endif

    return
  end subroutine bcast_1l

  !       Procedures to do global operations (Sum, Min, Max). _all_ routines
  !       send the information to all the processors otherwise they are
  !       kept on processor 0.

  subroutine global_sum_0i( lbuf, gbuf, mroot, idebug )
    integer, intent(in) :: lbuf
    integer, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = 1

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_sum_0i:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_0i


  subroutine global_sum_1i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:) :: lbuf
    integer, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_sum_1i:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_1i

  subroutine global_sum_2i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:) :: lbuf
    integer, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_sum_2i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_sum_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_2i

  subroutine global_sum_3i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:,:) :: lbuf
    integer, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** global_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_sum_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_3i

  subroutine global_sum_0r( lbuf, gbuf, mroot, idebug )
    real, intent(in) :: lbuf
    real, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = 1

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_sum_0r:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_0r


  subroutine global_sum_1r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:) :: lbuf
    real, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_sum_1r:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_1r

  subroutine global_sum_2r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:) :: lbuf
    real, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_sum_2r: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_sum_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_2r

  subroutine global_sum_3r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:,:) :: lbuf
    real, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** global_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_sum_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_3r

  subroutine global_sum_0d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in) :: lbuf
    double precision, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = 1

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_sum_0d:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_0d


  subroutine global_sum_1d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:) :: lbuf
    double precision, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_sum_1d:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_1d

  subroutine global_sum_2d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:) :: lbuf
    double precision, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_sum_2d: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_sum_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_2d

  subroutine global_sum_3d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:,:) :: lbuf
    double precision, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** global_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_sum_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_sum_3d

  subroutine global_all_sum_onevar_0d( gbuf )
    doubleprecision, intent(inout) :: gbuf
    doubleprecision :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_0d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_0d


  subroutine global_all_sum_onevar_1d( gbuf )
    doubleprecision, dimension(:), intent(inout) :: gbuf
    doubleprecision, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_1d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_1d

  subroutine global_all_sum_onevar_2d( gbuf )
    doubleprecision, dimension(:,:), intent(inout) :: gbuf
    doubleprecision, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_2d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_2d


  subroutine global_all_sum_onevar_3d( gbuf )
    doubleprecision, dimension(:,:,:), intent(inout) :: gbuf
    doubleprecision, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_3d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_3d



  subroutine global_all_sum_onevar_0i( gbuf )
    integer, intent(inout) :: gbuf
    integer :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_0i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_0i

  subroutine global_all_sum_onevar_1i( gbuf )
    integer, dimension(:), intent(inout) :: gbuf
    integer, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_1i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_1i

  subroutine global_all_sum_onevar_2i( gbuf )
    integer, dimension(:,:), intent(inout) :: gbuf
    integer, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_2i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_2i


  subroutine global_all_sum_onevar_3i( gbuf )
    integer, dimension(:,:,:), intent(inout) :: gbuf
    integer, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_3i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_3i

  subroutine global_all_sum_onevar_0r( gbuf )
    real, intent(inout) :: gbuf
    real :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_0r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_0r


  subroutine global_all_sum_onevar_1r( gbuf )
    real, dimension(:), intent(inout) :: gbuf
    real, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_1r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_1r

  subroutine global_all_sum_onevar_2r( gbuf )
    real, dimension(:,:), intent(inout) :: gbuf
    real, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_2r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_2r


  subroutine global_all_sum_onevar_3r( gbuf )
    real, dimension(:,:,:), intent(inout) :: gbuf
    real, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_sum_3r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_sum_onevar_3r


  subroutine global_all_sum_0i( lbuf, gbuf, mroot, idebug )
    integer, intent(in) :: lbuf
    integer, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = 1

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_sum_0i:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_0i


  subroutine global_all_sum_1i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:) :: lbuf
    integer, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_sum_1i:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_1i

  subroutine global_all_sum_2i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:) :: lbuf
    integer, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_sum_2i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_all_sum_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_2i

  subroutine global_all_sum_3i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:,:) :: lbuf
    integer, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    call assert( size(lbuf,3).eq.size(gbuf,3),  &
         '** global_all_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
         size(lbuf,3), size(gbuf,3) )

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_sum_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_3i

  subroutine global_all_sum_0r( lbuf, gbuf, mroot, idebug )
    real, intent(in) :: lbuf
    real, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = 1

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_sum_0r:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_0r


  subroutine global_all_sum_1r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:) :: lbuf
    real, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_sum_1r:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_1r

  subroutine global_all_sum_2r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:) :: lbuf
    real, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_sum_2r: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_all_sum_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_2r

  subroutine global_all_sum_3r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:,:) :: lbuf
    real, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    call assert( size(lbuf,3).eq.size(gbuf,3),  &
         '** global_all_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
         size(lbuf,3), size(gbuf,3) )

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_sum_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_3r

  subroutine global_all_sum_0d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in) :: lbuf
    double precision, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = 1

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_SUM, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_sum_0d:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_0d


  subroutine global_all_sum_1d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:) :: lbuf
    double precision, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype,  MPI_SUM, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_sum_1d:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_1d

  subroutine global_all_sum_2d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:) :: lbuf
    double precision, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_sum_2d: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_all_sum_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_2d

  subroutine global_all_sum_3d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:,:) :: lbuf
    double precision, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_sum_3i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    call assert( size(lbuf,3).eq.size(gbuf,3),  &
         '** global_all_sum_3i: size(lbuf,3).ne.size(gbuf,3) ', &
         size(lbuf,3), size(gbuf,3) )

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_sum_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_sum_3d

  subroutine global_min_0i( lbuf, gbuf, mroot, idebug )
    integer, intent(in) :: lbuf
    integer, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = 1

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_min_0i:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_0i


  subroutine global_min_1i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:) :: lbuf
    integer, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_min_1i:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_1i

  subroutine global_min_2i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:) :: lbuf
    integer, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_min_2i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_min_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_2i

  subroutine global_min_3i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:,:) :: lbuf
    integer, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** global_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_min_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_3i

  subroutine global_min_0r( lbuf, gbuf, mroot, idebug )
    real, intent(in) :: lbuf
    real, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = 1

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_min_0r:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_0r


  subroutine global_min_1r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:) :: lbuf
    real, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_min_1r:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_1r

  subroutine global_min_2r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:) :: lbuf
    real, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_min_2r: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_min_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_2r

  subroutine global_min_3r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:,:) :: lbuf
    real, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** global_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_min_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_3r

  subroutine global_min_0d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in) :: lbuf
    double precision, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = 1

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_min_0d:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_0d


  subroutine global_min_1d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:) :: lbuf
    double precision, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_min_1d:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_1d

  subroutine global_min_2d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:) :: lbuf
    double precision, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_min_2d: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_min_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_2d

  subroutine global_min_3d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:,:) :: lbuf
    double precision, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** global_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_min_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_min_3d

  subroutine global_all_min_onevar_0d( gbuf )
    doubleprecision, intent(inout) :: gbuf
    doubleprecision :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_0d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_0d


  subroutine global_all_min_onevar_1d( gbuf )
    doubleprecision, dimension(:), intent(inout) :: gbuf
    doubleprecision, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_1d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_1d

  subroutine global_all_min_onevar_2d( gbuf )
    doubleprecision, dimension(:,:), intent(inout) :: gbuf
    doubleprecision, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_2d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_2d


  subroutine global_all_min_onevar_3d( gbuf )
    doubleprecision, dimension(:,:,:), intent(inout) :: gbuf
    doubleprecision, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_3d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_3d




  subroutine global_all_min_onevar_0i( gbuf )
    integer, intent(inout) :: gbuf
    integer :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_0i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_0i


  subroutine global_all_min_onevar_1i( gbuf )
    integer, dimension(:), intent(inout) :: gbuf
    integer, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_1i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_1i

  subroutine global_all_min_onevar_2i( gbuf )
    integer, dimension(:,:), intent(inout) :: gbuf
    integer, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_2i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_2i


  subroutine global_all_min_onevar_3i( gbuf )
    integer, dimension(:,:,:), intent(inout) :: gbuf
    integer, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_3i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_3i

  subroutine global_all_min_onevar_0r( gbuf )
    real, intent(inout) :: gbuf
    real :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_0r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_0r


  subroutine global_all_min_onevar_1r( gbuf )
    real, dimension(:), intent(inout) :: gbuf
    real, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_1r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_1r

  subroutine global_all_min_onevar_2r( gbuf )
    real, dimension(:,:), intent(inout) :: gbuf
    real, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_2r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_2r



  subroutine global_all_min_onevar_3r( gbuf )
    real, dimension(:,:,:), intent(inout) :: gbuf
    real, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_min_3r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_min_onevar_3r


  subroutine global_all_min_0i( lbuf, gbuf, mroot, idebug )
    integer, intent(in) :: lbuf
    integer, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = 1

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_min_0i:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_0i


  subroutine global_all_min_1i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:) :: lbuf
    integer, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_min_1i:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_1i

  subroutine global_all_min_2i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:) :: lbuf
    integer, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_min_2i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_all_min_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_2i

  subroutine global_all_min_3i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:,:) :: lbuf
    integer, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    call assert( size(lbuf,3).eq.size(gbuf,3),  &
         '** global_all_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
         size(lbuf,3), size(gbuf,3) )

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_min_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_3i

  subroutine global_all_min_0r( lbuf, gbuf, mroot, idebug )
    real, intent(in) :: lbuf
    real, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = 1

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_min_0r:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_0r


  subroutine global_all_min_1r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:) :: lbuf
    real, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_min_1r:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_1r

  subroutine global_all_min_2r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:) :: lbuf
    real, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_min_2r: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_all_min_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_2r

  subroutine global_all_min_3r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:,:) :: lbuf
    real, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    call assert( size(lbuf,3).eq.size(gbuf,3),  &
         '** global_all_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
         size(lbuf,3), size(gbuf,3) )

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_min_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_3r

  subroutine global_all_min_0d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in) :: lbuf
    double precision, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = 1

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MIN, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_min_0d:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_0d


  subroutine global_all_min_1d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:) :: lbuf
    double precision, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MIN, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_min_1d:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_1d

  subroutine global_all_min_2d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:) :: lbuf
    double precision, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_min_2d: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_all_min_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_2d

  subroutine global_all_min_3d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:,:) :: lbuf
    double precision, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_min_3i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    call assert( size(lbuf,3).eq.size(gbuf,3),  &
         '** global_all_min_3i: size(lbuf,3).ne.size(gbuf,3) ', &
         size(lbuf,3), size(gbuf,3) )

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_min_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_min_3d

  subroutine global_max_0i( lbuf, gbuf, mroot, idebug )
    integer, intent(in) :: lbuf
    integer, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = 1

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_max_0i:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_0i


  subroutine global_max_1i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:) :: lbuf
    integer, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_max_1i:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_1i

  subroutine global_max_2i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:) :: lbuf
    integer, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_max_2i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_max_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_2i

  subroutine global_max_3i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:,:) :: lbuf
    integer, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** global_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_max_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_3i

  subroutine global_max_0r( lbuf, gbuf, mroot, idebug )
    real, intent(in) :: lbuf
    real, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = 1

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_max_0r:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_0r


  subroutine global_max_1r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:) :: lbuf
    real, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_max_1r:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_1r

  subroutine global_max_2r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:) :: lbuf
    real, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_max_2r: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_max_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_2r

  subroutine global_max_3r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:,:) :: lbuf
    real, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** global_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_max_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_3r

  subroutine global_max_0d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in) :: lbuf
    double precision, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = 1

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_max_0d:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_0d


  subroutine global_max_1d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:) :: lbuf
    double precision, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Reduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
         lroot, MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_max_1d:MPI_Reduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_1d

  subroutine global_max_2d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:) :: lbuf
    double precision, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_max_2d: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )
    endif

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_max_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_2d

  subroutine global_max_3d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:,:) :: lbuf
    double precision, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    if(myPE.eq.lroot) then
       call assert( size(lbuf,2).eq.size(gbuf,2),  &
            '** global_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
            size(lbuf,2), size(gbuf,2) )

       call assert( size(lbuf,3).eq.size(gbuf,3),  &
            '** global_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
            size(lbuf,3), size(gbuf,3) )
    endif

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_max_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_max_3d

  subroutine global_all_max_onevar_0d( gbuf )
    doubleprecision, intent(inout) :: gbuf
    doubleprecision :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_0d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_0d


  subroutine global_all_max_onevar_1d( gbuf )
    doubleprecision, dimension(:), intent(inout) :: gbuf
    doubleprecision, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_1d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_1d

  subroutine global_all_max_onevar_2d( gbuf )
    doubleprecision, dimension(:,:), intent(inout) :: gbuf
    doubleprecision, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_2d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_2d


  subroutine global_all_max_onevar_3d( gbuf )
    doubleprecision, dimension(:,:,:), intent(inout) :: gbuf
    doubleprecision, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_3d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_3d




  subroutine global_all_max_onevar_0i( gbuf )
    integer, intent(inout) :: gbuf
    integer :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_0i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_0i


  subroutine global_all_max_onevar_1i( gbuf )
    integer, dimension(:), intent(inout) :: gbuf
    integer, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_1i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_1i

  subroutine global_all_max_onevar_2i( gbuf )
    integer, dimension(:,:), intent(inout) :: gbuf
    integer, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_2i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_2i


  subroutine global_all_max_onevar_3i( gbuf )
    integer, dimension(:,:,:), intent(inout) :: gbuf
    integer, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_3i( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_3i

  subroutine global_all_max_onevar_0r( gbuf )
    real, intent(inout) :: gbuf
    real :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_0r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_0r

  subroutine global_all_max_onevar_1r( gbuf )
    real, dimension(:), intent(inout) :: gbuf
    real, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_1r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_1r

  subroutine global_all_max_onevar_2r( gbuf )
    real, dimension(:,:), intent(inout) :: gbuf
    real, dimension(size(gbuf,1),size(gbuf,2)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_2r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_2r


  subroutine global_all_max_onevar_3r( gbuf )
    real, dimension(:,:,:), intent(inout) :: gbuf
    real, dimension(size(gbuf,1),size(gbuf,2),size(gbuf,3)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_max_3r( lbuf, gbuf )
#endif
    return
  end subroutine global_all_max_onevar_3r



  subroutine global_all_max_0i( lbuf, gbuf, mroot, idebug )
    integer, intent(in) :: lbuf
    integer, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = 1

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_max_0i:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_0i


  subroutine global_all_max_1i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:) :: lbuf
    integer, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_INTEGER
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_max_1i:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_1i

  subroutine global_all_max_2i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:) :: lbuf
    integer, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_max_2i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_all_max_1i( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_2i

  subroutine global_all_max_3i( lbuf, gbuf, mroot, idebug )
    integer, intent(in), dimension(:,:,:) :: lbuf
    integer, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    call assert( size(lbuf,3).eq.size(gbuf,3),  &
         '** global_all_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
         size(lbuf,3), size(gbuf,3) )

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_max_1i( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_3i

  subroutine global_all_max_0r( lbuf, gbuf, mroot, idebug )
    real, intent(in) :: lbuf
    real, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = 1

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_max_0r:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_0r


  subroutine global_all_max_1r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:) :: lbuf
    real, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_REAL
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_max_1r:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_1r

  subroutine global_all_max_2r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:) :: lbuf
    real, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_max_2r: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_all_max_1r( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_2r

  subroutine global_all_max_3r( lbuf, gbuf, mroot, idebug )
    real, intent(in), dimension(:,:,:) :: lbuf
    real, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    call assert( size(lbuf,3).eq.size(gbuf,3),  &
         '** global_all_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
         size(lbuf,3), size(gbuf,3) )

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_max_1r( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_3r

  subroutine global_all_max_0d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in) :: lbuf
    double precision, intent(out) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = 1

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype, MPI_MAX, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_max_0d:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_0d


  subroutine global_all_max_1d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:) :: lbuf
    double precision, intent(out), dimension(:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: recvtype, sendtype, sendcnt, ierr,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    recvtype = MPI_DOUBLE_PRECISION
    sendtype = recvtype

    sendcnt = size(lbuf)

    call MPI_Allreduce( lbuf, gbuf, sendcnt, sendtype,  MPI_MAX, &
         MPI_COMM_WORLD, ierr )
    call MPI_Check( 'global_all_max_1d:MPI_Allreduce', ierr )
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_1d

  subroutine global_all_max_2d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:) :: lbuf
    double precision, intent(out), dimension(:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: i,j,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_max_2d: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    do j=lbound(lbuf,2),ubound(lbuf,2)
       call global_all_max_1d( lbuf(:,j), gbuf(:,j), lroot, lidebug )
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_2d

  subroutine global_all_max_3d( lbuf, gbuf, mroot, idebug )
    double precision, intent(in), dimension(:,:,:) :: lbuf
    double precision, intent(out), dimension(:,:,:) :: gbuf
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    integer :: j,k,lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    call assert( size(lbuf,2).eq.size(gbuf,2),  &
         '** global_all_max_3i: size(lbuf,2).ne.size(gbuf,2) ', &
         size(lbuf,2), size(gbuf,2) )

    call assert( size(lbuf,3).eq.size(gbuf,3),  &
         '** global_all_max_3i: size(lbuf,3).ne.size(gbuf,3) ', &
         size(lbuf,3), size(gbuf,3) )

    do k=lbound(lbuf,3),ubound(lbuf,3)
       do j=lbound(lbuf,2),ubound(lbuf,2)
          call global_all_max_1d( lbuf(:,j,k), gbuf(:,j,k), lroot, lidebug )
       enddo
    enddo
#else
    gbuf = lbuf
#endif

    return
  end subroutine global_all_max_3d



  subroutine global_all_or_onevar_0d( gbuf )
    logical, intent(inout) :: gbuf
    logical :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_or_0d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_or_onevar_0d

  subroutine global_all_or_onevar_1d( gbuf )
    logical, dimension(:), intent(inout) :: gbuf
    logical, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_or_1d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_or_onevar_1d

  subroutine global_all_and_onevar_0d( gbuf )
    logical, intent(inout) :: gbuf
    logical :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_and_0d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_and_onevar_0d

  subroutine global_all_and_onevar_1d( gbuf )
    logical, dimension(:), intent(inout) :: gbuf
    logical, dimension(size(gbuf)) :: lbuf

#ifdef MPI
    lbuf = gbuf
    call global_all_and_1d( lbuf, gbuf )
#endif
    return
  end subroutine global_all_and_onevar_1d


  subroutine global_all_and_0d( lvalue, gvalue, mroot, idebug )
    logical, intent(in) :: lvalue
    logical, intent(out) :: gvalue
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    !       ---------------
    !       local variables
    !       ---------------
    integer :: ierror, icount
    integer :: lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif

    icount = 1

    call MPI_Allreduce( lvalue, gvalue, icount, MPI_LOGICAL, &
         MPI_LAND, MPI_COMM_WORLD, ierror )

    call MPI_Check( 'global_all_and_0d ', ierror )
#else
    gvalue = lvalue
#endif
    return
  end subroutine  global_all_and_0d


  subroutine global_all_and_1d( lvalue, gvalue, mroot, idebug )
    logical, intent(in), dimension(:) :: lvalue
    logical, intent(out), dimension(:) :: gvalue
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    !       ---------------
    !       local variables
    !       ---------------
    integer :: ierror, icount
    integer :: lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif


    icount = size( lvalue )

    call MPI_Allreduce( lvalue, gvalue, icount, MPI_LOGICAL, &
         MPI_LAND, MPI_COMM_WORLD, ierror )

    call MPI_Check( 'global_all_and_1d ', ierror )
#else
    gvalue = lvalue
#endif
    return
  end subroutine global_all_and_1d



  subroutine global_all_or_0d( lvalue, gvalue, mroot, idebug )
    logical, intent(in) :: lvalue
    logical, intent(out) :: gvalue
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    !       ---------------
    !       local variables
    !       ---------------
    integer :: ierror, icount
    integer :: lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif


    icount = 1

    call MPI_Allreduce( lvalue, gvalue, icount, MPI_LOGICAL, &
         MPI_LOR, MPI_COMM_WORLD, ierror )

    call MPI_Check( 'global_all_or_0d ', ierror )
#else
    gvalue = lvalue
#endif
    return
  end subroutine global_all_or_0d


  subroutine global_all_or_1d( lvalue, gvalue, mroot, idebug )
    logical, intent(in), dimension(:) :: lvalue
    logical, intent(out), dimension(:) :: gvalue
    integer, optional, intent(in) :: mroot, idebug

#ifdef MPI
    !       ---------------
    !       local variables
    !       ---------------
    integer :: ierror, icount
    integer :: lroot, lidebug

    if (.not. present(mroot)) then
       lroot = 0
    else
       lroot = mroot
    endif

    if (.not. present(idebug)) then
       lidebug = 0
    else
       lidebug = idebug
    endif


    icount = size( lvalue )

    call MPI_Allreduce( lvalue, gvalue, icount, MPI_LOGICAL, &
         MPI_LOR, MPI_COMM_WORLD, ierror )

    call MPI_Check( 'global_all_or_1d ', ierror )
#else
    gvalue = lvalue
#endif
    return
  end subroutine global_all_or_1d


  !``````````````````````````````````````````````````````````````````````!
  ! Subroutine: ExitMPI                                                  !
  !                                                                      !
  ! Purpose: Clean abort from a parallel run. This routine is invoked by !
  ! by calling MFIX_EXIT                                                 !
  !......................................................................!
  SUBROUTINE ExitMPI(myid)

    USE funits, only: UNIT_LOG
    USE funits, only: DMP_LOG

    implicit none

    INTEGER, optional, intent(in) :: MyID

#ifdef MPI
    INTEGER :: MyID_l
    INTEGER :: ERRORCODE

    ! Flag to call MPI_ABORT and bypass the call to MPI_Finalize.
    ! This is only needed if debugging a 'deadlocked' run.
    LOGICAL, PARAMETER :: FORCED_ABORT = .FALSE.

    ! Process ID (myPE) converted to a character string.
    CHARACTER(len=64) :: myID_c


    ! Set the ID of the caller.
    myID_l= merge(MyID, myPE, PRESENT(MyID))
    myID_c=''; WRITE(myID_c,*) myID_l

    ! Hard abort. If you need this functionality, then you need to figure
    ! out why the code has deadlocked. Most likely, a call to MFIX_EXIT
    ! was put inside of a logical branch that only a few ranks execute.
    ! DON'T JUST USE A FORCED ABORT --> FIX THE CODE CAUSE DEADLOCK <--
    IF(FORCED_ABORT) THEN
       ERRORCODE = 100 + myPE
       CALL MPI_ABORT(MPI_COMM_WORLD, ERRORCODE, MPIERR)
       WRITE(*,2000) myID_c, MPIERR

       ! Calls to ExitMPI (via MFIX_EXIT) should be made by all processes
       ! and therefore calling MPI_Finalize should be sufficient to exit
       ! a failed run. However, a FORCED_ABORT can be issued if deadlock
       ! is an issue.
    ELSE
       CALL MPI_BARRIER(MPI_COMM_WORLD, MPIERR)
       CALL MPI_Finalize(MPIERR)
    ENDIF

    ! Notify that MPI was cleanly terminated. This point will not be
    ! reached if MPI is aborted.
    IF(myPE == PE_IO) WRITE(*,1000)
#endif

    RETURN

1000 FORMAT(2/,1X,'MPI Terminated.')
2000 FORMAT(2/,1X,'Rank ',A,' :: MPI_ABORT CODE = ',I4)

  END SUBROUTINE exitMPI

end module mpi_utility
