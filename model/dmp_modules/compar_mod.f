! -*- f90 -*-
!-------------------------------------------------------------------!
!                                                                   !
! Purpose:                                                          !
!    Variables to be declared for parallel information. Removed     !
!    from geometry_mod and put in COMPAR module with some           !
!    additional variables used by AEOLUS                            !
!                                                                   !
! Added by Ed and Sreekanth on 06/22/99.                            !
!-------------------------------------------------------------------!

      MODULE compar

!-----------------------------------------------
! Modules
!-----------------------------------------------
#ifdef MPI
        USE mpi ! ignore-depcomp
#endif
!-----------------------------------------------

! myPE - my processor id (it varies from 0 to nproc-1)
! numPEs - total number of nodes
      integer :: myPE, numPEs

! mpierr - used by AEOLUS for error checking
      INTEGER :: mpierr

! specify the rank of the PE to be used for I/O
      INTEGER :: PE_IO = 0

! nodesi, nodesj and nodesk represent the number of nodes
! in i, j, k directions respectively.
! nodesj = 1 (No decomposition along j-direction)
! For 1-D decomposition, nodesk = nproc for a 3d problem and
! nodesi = nproc for a 2D problem.
      integer :: nodesi, nodesj, nodesk

! root represents the 'root' processor. For now it is defaulted to
! zero
      integer :: root
      data root /0/

! nlayers_bicgs - Number of layers for send_recv in bicgs
      integer :: nlayers_bicgs = 1


! -istart1_all contains the starting i value for all the processors
!  excluding the ghost regions. istart2_all is for one extra ghost
!  layer and istart3_all is for two ghost layers. Similarly
!  iend1_all, iend2_all and iend3_all contain the ending values.
!  Similarly for j and k, jstart..., kstart.... are  prescribed.
! -All the variables without the '_all' represent that processor
!  values. So ijkstart3 denotes the starting value of ijk, which
!  belongs to the processor = funijk(istart3_all(myid),
!  jstart3_all(myid), kstart3_all(myid) for a 1-d decompostion of a
!  3D problem. For more details see gridmap_mod.f90. Similarly the
!  end values are denoted by ijkend3_all
! -displs has the necessary shift information to position the buffer
!  in the scatterv and gatherv routines.
! -ijksize3 is the size of the element owned by each processor plus
!  the ghost regions.
! -'_all' has information about all the processor mapping based on
!  above convention
      integer, allocatable,dimension(:) ::  &
                ijkstart3_all,ijkend3_all,    &
                istart_all,istart1_all,istart2_all,istart3_all, &
                jstart_all,jstart1_all,jstart2_all,jstart3_all, &
                kstart_all,kstart1_all,kstart2_all,kstart3_all, &
                iend_all,iend1_all,iend2_all,iend3_all, &
                jend_all,jend1_all,jend2_all,jend3_all, &
                kend_all,kend1_all,kend2_all,kend3_all, &
                ijksize3_all, displs

! Variables used for mapping i, j, k to ii, jj, kk to take care of
! of cyclic conditions...
      integer, allocatable,dimension(:) :: imap, jmap, kmap
      integer, allocatable,dimension(:) :: imap_c, jmap_c, kmap_c

      integer :: &
                ijksize3, ijkstart3,ijkend3, &
                istart3, iend3, jstart3, jend3, &
                kstart3, kend3, istart2, iend2, jstart2, jend2, &
                kstart2, kend2, istart1, iend1, jstart1, jend1, &
                kstart1, kend1

      integer :: istart, iend, jstart, jend, kstart, kend

! Variables used for fourth order methods
      integer, allocatable,dimension(:) ::  &
                ijkstart4_all,ijkend4_all, ijksize4_all,&
                istart4_all, jstart4_all, kstart4_all, &
                iend4_all, jend4_all, kend4_all
      integer :: &
                istart4, jstart4, kstart4, &
                iend4, jend4, kend4, &
                ijkstart4,ijkend4,ijksize4

! declaration for storing filebasename, e.g. mfix00000.dat
      CHARACTER(len=5) :: fbname
      INTEGER :: idbg = 1

! Funijk coefficients
      integer :: c0, c1, c2

! Funijk3 coefficients
      integer :: c0_3, c1_3, c2_3


!       Integer Array of IJK values at each (I,J,K) cell

        integer, allocatable, dimension(:,:,:) :: IJK_ARRAY_OF,FUNIJK_MAP_C

!        integer, allocatable, dimension(:,:,:) :: funijk

!       Integer Array of neighbor cells

        integer, allocatable, dimension(:)     :: WEST_ARRAY_OF,EAST_ARRAY_OF
        integer, allocatable, dimension(:)     :: SOUTH_ARRAY_OF,NORTH_ARRAY_OF
        integer, allocatable, dimension(:)     :: BOTTOM_ARRAY_OF,TOP_ARRAY_OF
        integer, allocatable, dimension(:)     :: IM_ARRAY_OF,IP_ARRAY_OF
        integer, allocatable, dimension(:)     :: JM_ARRAY_OF,JP_ARRAY_OF
        integer, allocatable, dimension(:)     :: KM_ARRAY_OF,KP_ARRAY_OF

!       Flag to identify dead (unused cells)

        LOGICAL, allocatable, dimension(:,:,:) :: DEAD_CELL_AT

!       Flag to know if above neighbor arrays have been allocated

        LOGICAL :: INCREMENT_ARRAYS_ALLOCATED

!       Number of Ghost Cells

        INTEGER :: NGC_EAST
        INTEGER :: NGC_WEST
        INTEGER :: NGC_NORTH
        INTEGER :: NGC_SOUTH
        INTEGER :: NGC_TOP
        INTEGER :: NGC_BOTTOM

!       List of Ghost Cells

        INTEGER, ALLOCATABLE, DIMENSION(:) ::  LGC_EAST
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  LGC_WEST
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  LGC_NORTH
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  LGC_SOUTH
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  LGC_TOP
        INTEGER, ALLOCATABLE, DIMENSION(:) ::  LGC_BOTTOM

!       Domain size of each processor

        INTEGER, ALLOCATABLE, DIMENSION(:) :: ISIZE_ALL,JSIZE_ALL,KSIZE_ALL

        LOGICAL :: DOMAIN_SIZE_ADJUSTED = .FALSE.

        INTEGER, ALLOCATABLE, DIMENSION(:) :: NCPP_UNIFORM

        LOGICAL :: NCPP_UNIFORM_BACKED_UP = .FALSE.

        integer, allocatable,dimension(:) ::  new_ijksize3_all

!       Flag to exit gridmap_init after domain size is assigned
        LOGICAL :: SHORT_GRIDMAP_INIT = .FALSE.


      END MODULE compar

