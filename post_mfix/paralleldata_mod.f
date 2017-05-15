! *************************  module ParallelData *********************

      module ParallelData

      type CellMap
         integer :: proc
         integer :: ijk
      end type CellMap


      ! given IJK_IO (from 1 to ijkmax2) , then the value for that
      ! cell is found
      !
      ! on processor : cell_map(ijk_io)%proc
      ! and ijk      : cell_map(ijk_io)%ijk
      TYPE (CellMap) , dimension (:) , allocatable :: cell_map

      integer , allocatable :: cell_map_v2(:,:)



      ! the following are defined for each processor
      integer , dimension(:), allocatable  :: is3       ! istart3
      integer , dimension(:), allocatable  :: ie3       ! iend3
      integer , dimension(:), allocatable  :: js3       ! jstart3
      integer , dimension(:), allocatable  :: je3       ! jend3
      integer , dimension(:), allocatable  :: ks3       ! kstart3
      integer , dimension(:), allocatable  :: ke3       ! kend3
      integer , dimension(:), allocatable  :: n_cells   ! number of cells

      ! current record for each processor file
      integer , dimension(:), allocatable  :: cr


      ! the array used to hold the processor IO data
      real , dimension(:), allocatable    :: r_tmp


      integer :: np  ! mumber of processors for this run


      character(LEN=80) :: fname_scav, fname_dist

      character(LEN=35) :: ext

      character(LEN=512) :: pbuffer

      integer , allocatable :: cellcount(:)


      end module ParallelData

