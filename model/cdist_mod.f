! -*- f90 -*-
      module cdist

      logical :: bDist_IO
      logical :: bStart_with_one_RES
      logical :: bDoing_postmfix

      integer :: netCDF_file_index

      logical :: bWrite_netcdf(20)

      logical :: bFirst_netcdf_write = .true.

      logical :: bGlobalNetcdf = .true.

      end module cdist
