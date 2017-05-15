module MFIX_netcdf

#ifdef NETCDF
  include 'netcdf.inc'
  USE netcdf ! ignore-depcomp
#else
  integer NF90_64BIT_OFFSET

  integer NF90_DOUBLE
  integer NF90_INT
  integer NF90_NOWRITE
  integer NF90_CLOBBER

  parameter (NF90_64BIT_OFFSET = 0)
  parameter (NF90_DOUBLE = 0)
  parameter (NF90_INT = 0)
  parameter (NF90_NOWRITE = 0)
  parameter (NF90_CLOBBER = 0)
#endif

  ! Overloaded variable functions
  interface MFIX_nf90_def_var
     module procedure MFIX_nf90_def_var_Scalar, MFIX_nf90_def_var_oneDim, MFIX_nf90_def_var_ManyDims
  end interface MFIX_nf90_def_var ! MFIX_nf90_def_var

  interface MFIX_nf90_put_var
     module procedure MFIX_nf90_put_var_text,                                   &
          MFIX_nf90_put_var_FourByteInt,  &
          MFIX_nf90_put_var_FourByteReal, MFIX_nf90_put_var_EightByteReal
     module procedure MFIX_nf90_put_var_1D_text,                                      &
          MFIX_nf90_put_var_1D_FourByteInt, &
          MFIX_nf90_put_var_1D_FourByteReal, MFIX_nf90_put_var_1D_EightByteReal
     module procedure MFIX_nf90_put_var_2D_text,                                       &
          MFIX_nf90_put_var_2D_FourByteInt,   &
          MFIX_nf90_put_var_2D_FourByteReal, MFIX_nf90_put_var_2D_EightByteReal
     module procedure MFIX_nf90_put_var_3D_text,                                       &
          MFIX_nf90_put_var_3D_FourByteInt, &
          MFIX_nf90_put_var_3D_FourByteReal, MFIX_nf90_put_var_3D_EightByteReal
     module procedure MFIX_nf90_put_var_4D_FourByteInt , MFIX_nf90_put_var_4D_EightByteReal
  end interface MFIX_nf90_put_var ! MFIX_nf90_put_var

  interface MFIX_nf90_get_var
     module procedure MFIX_nf90_get_var_text,                                   &
          MFIX_nf90_get_var_FourByteInt, & ! MFIX_nf90_get_var_EightByteInt, &
          MFIX_nf90_get_var_FourByteReal, MFIX_nf90_get_var_EightByteReal
     module procedure MFIX_nf90_get_var_1D_text,                                      &
          MFIX_nf90_get_var_1D_FourByteInt , &
          MFIX_nf90_get_var_1D_FourByteReal, MFIX_nf90_get_var_1D_EightByteReal
     module procedure MFIX_nf90_get_var_2D_text,                                      &
          MFIX_nf90_get_var_2D_FourByteReal, MFIX_nf90_get_var_2D_EightByteReal
     module procedure MFIX_nf90_get_var_3D_text,                                      &
          MFIX_nf90_get_var_3D_FourByteInt, &
          MFIX_nf90_get_var_3D_FourByteReal, MFIX_nf90_get_var_3D_EightByteReal
     module procedure MFIX_nf90_get_var_4D_text
     module procedure MFIX_nf90_get_var_5D_text
     module procedure MFIX_nf90_get_var_6D_text
     module procedure MFIX_nf90_get_var_7D_text
  end interface MFIX_nf90_get_var ! MFIX_nf90_get_var

contains
  function MFIX_nf90_def_var_Scalar(ncid, name, xtype, varid)
    integer,               intent( in) :: ncid
    character (len = *),   intent( in) :: name
    integer,               intent( in) :: xtype
    integer,               intent(inout) :: varid
    integer                            :: MFIX_nf90_def_var_Scalar

    ! Dummy - shouldn't get used
    integer, dimension(1) :: dimids

#ifdef NETCDF
    MFIX_nf90_def_var_Scalar = nf90_def_var(ncid, name, xtype, varid)
#else
    MFIX_nf90_def_var_Scalar = 1
#endif
  end function MFIX_nf90_def_var_Scalar
  ! -----
  function MFIX_nf90_def_var_oneDim(ncid, name, xtype, dimids, varid)
    integer,               intent( in) :: ncid
    character (len = *),   intent( in) :: name
    integer,               intent( in) :: xtype
    integer,               intent( in) :: dimids
    integer,               intent(inout) :: varid
    integer                            :: MFIX_nf90_def_var_oneDim

    integer, dimension(1) :: dimidsA
    dimidsA(1) = dimids
#ifdef NETCDF
    MFIX_nf90_def_var_oneDim = nf90_def_var(ncid, name, xtype, dimids, varid)
#else
    MFIX_nf90_def_var_oneDim = 1
#endif
  end function MFIX_nf90_def_var_oneDim
  ! -----
  function MFIX_nf90_def_var_ManyDims(ncid, name, xtype, dimids, varid)
    integer,               intent( in) :: ncid
    character (len = *),   intent( in) :: name
    integer,               intent( in) :: xtype
    integer, dimension(:), intent( in) :: dimids
    integer,               intent(inout) :: varid
    integer                            :: MFIX_nf90_def_var_ManyDims

#ifdef NETCDF
    MFIX_nf90_def_var_ManyDims = nf90_def_var(ncid, name, xtype, dimids, varid)
#else
    MFIX_nf90_def_var_ManyDims = 1
#endif
  end function MFIX_nf90_def_var_ManyDims


  function MFIX_nf90_put_var_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *),             intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_text


#ifdef NETCDF
    MFIX_nf90_put_var_text = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_text = 1
#endif
  end function MFIX_nf90_put_var_text

  function MFIX_nf90_get_var_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *),             intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_text


#ifdef NETCDF
    MFIX_nf90_get_var_text = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_text = 1
#endif
  end function MFIX_nf90_get_var_text


  function MFIX_nf90_put_var_1D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_1D_text

    integer, parameter                    :: numDims = 1

#ifdef NETCDF
    MFIX_nf90_put_var_1D_text = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_1D_text = 1
#endif

  end function MFIX_nf90_put_var_1D_text


  function MFIX_nf90_put_var_2D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:, :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_2D_text

    integer, parameter                    :: numDims = 2

#ifdef NETCDF
    MFIX_nf90_put_var_2D_text = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_2D_text = 1
#endif

  end function MFIX_nf90_put_var_2D_text


  function MFIX_nf90_put_var_3D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_3D_text

    integer, parameter                    :: numDims = 3
#ifdef NETCDF
    MFIX_nf90_put_var_3D_text = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_3D_text = 1
#endif

  end function MFIX_nf90_put_var_3D_text



  function MFIX_nf90_get_var_1D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_1D_text

    integer, parameter                  :: numDims = 1
#ifdef NETCDF
    MFIX_nf90_get_var_1D_text = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_1D_text = 1
#endif
  end function MFIX_nf90_get_var_1D_text


  function MFIX_nf90_get_var_2D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_2D_text

    integer, parameter                  :: numDims = 2
#ifdef NETCDF
    MFIX_nf90_get_var_2D_text = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_2D_text = 1
#endif
  end function MFIX_nf90_get_var_2D_text


  function MFIX_nf90_get_var_3D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_3D_text

    integer, parameter                  :: numDims = 3
#ifdef NETCDF
    MFIX_nf90_get_var_3D_text = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_3D_text = 1
#endif
  end function MFIX_nf90_get_var_3D_text


  function MFIX_nf90_get_var_4D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_4D_text

    integer, parameter                  :: numDims = 4
#ifdef NETCDF
    MFIX_nf90_get_var_4D_text = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_4D_text = 1
#endif
  end function MFIX_nf90_get_var_4D_text


  function MFIX_nf90_get_var_5D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_5D_text

    integer, parameter                  :: numDims = 5

#ifdef NETCDF
    MFIX_nf90_get_var_5D_text = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_5D_text = 1
#endif
  end function MFIX_nf90_get_var_5D_text


  function MFIX_nf90_get_var_6D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :, :, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_6D_text

    integer, parameter                  :: numDims = 6
#ifdef NETCDF
    MFIX_nf90_get_var_6D_text = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_6D_text = 1
#endif
  end function MFIX_nf90_get_var_6D_text


  function MFIX_nf90_get_var_7D_text(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :, :, :, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_7D_text

    integer, parameter                  :: numDims = 7
#ifdef NETCDF
    MFIX_nf90_get_var_7D_text = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_7D_text = 1
#endif
  end function MFIX_nf90_get_var_7D_text

  function MFIX_nf90_put_var_FourByteInt(ncid, varid, values, start)
    integer,                         intent( in) :: ncid, varid
    integer , intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer                                      :: MFIX_nf90_put_var_FourByteInt

#ifdef NETCDF
    MFIX_nf90_put_var_FourByteInt = nf90_put_var(ncid, varid, values, start)
#else
    MFIX_nf90_put_var_FourByteInt = 1
#endif
  end function MFIX_nf90_put_var_FourByteInt

  function MFIX_nf90_put_var_FourByteReal(ncid, varid, values, start)
    integer,                         intent( in) :: ncid, varid
    real , intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer                                      :: MFIX_nf90_put_var_FourByteReal


#ifdef NETCDF
    MFIX_nf90_put_var_FourByteReal = nf90_put_var(ncid, varid, values, start)
#else
    MFIX_nf90_put_var_FourByteReal = 1
#endif
  end function MFIX_nf90_put_var_FourByteReal


  function MFIX_nf90_put_var_EightByteReal(ncid, varid, values, start)
    integer,                         intent( in) :: ncid, varid
    double precision, intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer                                      :: MFIX_nf90_put_var_EightByteReal


#ifdef NETCDF
    MFIX_nf90_put_var_EightByteReal = nf90_put_var(ncid, varid, values, start)
#else
    MFIX_nf90_put_var_EightByteReal = 1
#endif
  end function MFIX_nf90_put_var_EightByteReal



  function MFIX_nf90_get_var_TwoByteInt(ncid, varid, values, start)
    integer,                         intent( in) :: ncid, varid
    integer , intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer                                      :: MFIX_nf90_get_var_TwoByteInt

    integer                               :: counter


#ifdef NETCDF
    MFIX_nf90_get_var_TwoByteInt = nf90_get_var(ncid, varid, values, start)
#else
    MFIX_nf90_get_var_TwoByteInt = 1
#endif
  end function MFIX_nf90_get_var_TwoByteInt


  function MFIX_nf90_get_var_FourByteInt(ncid, varid, values, start)
    integer,                         intent( in) :: ncid, varid
    integer , intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer                                      :: MFIX_nf90_get_var_FourByteInt

    integer                               :: counter
    integer                               :: defaultInteger

#ifdef NETCDF
    MFIX_nf90_get_var_FourByteInt = nf90_get_var(ncid, varid, values, start)
#else
    MFIX_nf90_get_var_FourByteInt = 1
#endif
  end function MFIX_nf90_get_var_FourByteInt


  !   function MFIX_nf90_get_var_EightByteInt(ncid, varid, values, start)
  !     integer,                         intent( in) :: ncid, varid
  !     integer , intent(inout) :: values
  !     integer, dimension(:), optional, intent( in) :: start
  !     integer                                      :: MFIX_nf90_get_var_EightByteInt
  !
  !
#ifdef NETCDF
  !     MFIX_nf90_get_var_EightByteInt = nf90_get_var(ncid, varid, values, start)
#else
  !     MFIX_nf90_get_var_EightByteInt = 1
#endif
  !   end function MFIX_nf90_get_var_EightByteInt


  function MFIX_nf90_get_var_FourByteReal(ncid, varid, values, start)
    integer,                         intent( in) :: ncid, varid
    real , intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer                                      :: MFIX_nf90_get_var_FourByteReal


#ifdef NETCDF
    MFIX_nf90_get_var_FourByteReal = nf90_get_var(ncid, varid, values, start)
#else
    MFIX_nf90_get_var_FourByteReal = 1
#endif
  end function MFIX_nf90_get_var_FourByteReal


  function MFIX_nf90_get_var_EightByteReal(ncid, varid, values, start)
    integer,                         intent( in) :: ncid, varid
    double precision, intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer                                      :: MFIX_nf90_get_var_EightByteReal

#ifdef NETCDF
    MFIX_nf90_get_var_EightByteReal = nf90_get_var(ncid, varid, values, start)
#else
    MFIX_nf90_get_var_EightByteReal = 1
#endif
  end function MFIX_nf90_get_var_EightByteReal


  function MFIX_nf90_get_var_1D_TwoByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_1D_TwoByteInt

#ifdef NETCDF
    MFIX_nf90_get_var_1D_TwoByteInt = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_1D_TwoByteInt = 1
#endif

  end function MFIX_nf90_get_var_1D_TwoByteInt


  function MFIX_nf90_get_var_2D_TwoByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_2D_TwoByteInt

#ifdef NETCDF
    MFIX_nf90_get_var_2D_TwoByteInt = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_2D_TwoByteInt = 1
#endif
  end function MFIX_nf90_get_var_2D_TwoByteInt


  function MFIX_nf90_get_var_3D_TwoByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:, :, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_3D_TwoByteInt

#ifdef NETCDF
    MFIX_nf90_get_var_3D_TwoByteInt = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_3D_TwoByteInt = 1
#endif
  end function MFIX_nf90_get_var_3D_TwoByteInt

  function MFIX_nf90_get_var_1D_FourByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_1D_FourByteInt
#ifdef NETCDF
    MFIX_nf90_get_var_1D_FourByteInt = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_1D_FourByteInt = 1
#endif
  end function MFIX_nf90_get_var_1D_FourByteInt


  function MFIX_nf90_get_var_2D_FourByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_2D_FourByteInt

#ifdef NETCDF
    MFIX_nf90_get_var_2D_FourByteInt = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_2D_FourByteInt = 1
#endif
  end function MFIX_nf90_get_var_2D_FourByteInt


  function MFIX_nf90_get_var_3D_FourByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:, :, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_3D_FourByteInt
#ifdef NETCDF
    MFIX_nf90_get_var_3D_FourByteInt = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_3D_FourByteInt = 1
#endif
  end function MFIX_nf90_get_var_3D_FourByteInt



  function MFIX_nf90_get_var_1D_FourByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    real , dimension(:), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_1D_FourByteReal
#ifdef NETCDF
    MFIX_nf90_get_var_1D_FourByteReal = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_1D_FourByteReal = 1
#endif
  end function MFIX_nf90_get_var_1D_FourByteReal


  function MFIX_nf90_get_var_2D_FourByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    real , dimension(:, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_2D_FourByteReal

#ifdef NETCDF
    MFIX_nf90_get_var_2D_FourByteReal = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_2D_FourByteReal = 1
#endif
  end function MFIX_nf90_get_var_2D_FourByteReal


  function MFIX_nf90_get_var_3D_FourByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    real , dimension(:, :, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_3D_FourByteReal
#ifdef NETCDF
    MFIX_nf90_get_var_3D_FourByteReal = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_3D_FourByteReal = 1
#endif
  end function MFIX_nf90_get_var_3D_FourByteReal

  function MFIX_nf90_get_var_1D_EightByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    double precision, dimension(:), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_1D_EightByteReal

#ifdef NETCDF
    MFIX_nf90_get_var_1D_EightByteReal = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_1D_EightByteReal = 1
#endif
  end function MFIX_nf90_get_var_1D_EightByteReal


  function MFIX_nf90_get_var_2D_EightByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    double precision, dimension(:, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_2D_EightByteReal

#ifdef NETCDF
    MFIX_nf90_get_var_2D_EightByteReal = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_2D_EightByteReal = 1
#endif
  end function MFIX_nf90_get_var_2D_EightByteReal


  function MFIX_nf90_get_var_3D_EightByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    double precision, dimension(:, :, :), &
         intent(inout) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_get_var_3D_EightByteReal

#ifdef NETCDF
    MFIX_nf90_get_var_3D_EightByteReal = nf90_get_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_get_var_3D_EightByteReal = 1
#endif
  end function MFIX_nf90_get_var_3D_EightByteReal

  function MFIX_nf90_inquire_dimension(ncid, dimid, name, len)
    integer,                       intent( in) :: ncid, dimid
    character (len = *), optional, intent(inout) :: name
    integer,             optional, intent(inout) :: len
    integer                                    :: MFIX_nf90_inquire_dimension

#ifdef NETCDF
    MFIX_nf90_inquire_dimension = nf90_inquire_dimension(ncid, dimid, name, len)
#else
    MFIX_nf90_inquire_dimension = 1
#endif
  end function MFIX_nf90_inquire_dimension

  function MFIX_nf90_put_var_1D_FourByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_1D_FourByteInt

#ifdef NETCDF
    MFIX_nf90_put_var_1D_FourByteInt = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_1D_FourByteInt = 1
#endif
  end function MFIX_nf90_put_var_1D_FourByteInt


  function MFIX_nf90_put_var_2D_FourByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:, :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_2D_FourByteInt

#ifdef NETCDF
    MFIX_nf90_put_var_2D_FourByteInt = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_2D_FourByteInt = 1
#endif
  end function MFIX_nf90_put_var_2D_FourByteInt


  function MFIX_nf90_put_var_3D_FourByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:, :, :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_3D_FourByteInt
#ifdef NETCDF
    MFIX_nf90_put_var_3D_FourByteInt = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_3D_FourByteInt = 1
#endif
  end function MFIX_nf90_put_var_3D_FourByteInt

  function MFIX_nf90_put_var_4D_FourByteInt(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    integer , dimension(:, :, : , :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_4D_FourByteInt
#ifdef NETCDF
    MFIX_nf90_put_var_4D_FourByteInt = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_4D_FourByteInt = 1
#endif
  end function MFIX_nf90_put_var_4D_FourByteInt

  function MFIX_nf90_put_var_1D_FourByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    real , dimension(:), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_1D_FourByteReal

#ifdef NETCDF
    MFIX_nf90_put_var_1D_FourByteReal = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_1D_FourByteReal = 1
#endif
  end function MFIX_nf90_put_var_1D_FourByteReal


  function MFIX_nf90_put_var_2D_FourByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    real , dimension(:, :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_2D_FourByteReal

#ifdef NETCDF
    MFIX_nf90_put_var_2D_FourByteReal = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_2D_FourByteReal = 1
#endif
  end function MFIX_nf90_put_var_2D_FourByteReal


  function MFIX_nf90_put_var_3D_FourByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    real , dimension(:, :, :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_3D_FourByteReal
#ifdef NETCDF
    MFIX_nf90_put_var_3D_FourByteReal = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_3D_FourByteReal = 1
#endif
  end function MFIX_nf90_put_var_3D_FourByteReal


  function MFIX_nf90_put_var_1D_EightByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    double precision, dimension(:), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_1D_EightByteReal

#ifdef NETCDF
    MFIX_nf90_put_var_1D_EightByteReal = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_1D_EightByteReal = 1
#endif
  end function MFIX_nf90_put_var_1D_EightByteReal


  function MFIX_nf90_put_var_2D_EightByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    double precision, dimension(:, :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_2D_EightByteReal
#ifdef NETCDF
    MFIX_nf90_put_var_2D_EightByteReal = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_2D_EightByteReal = 1
#endif
  end function MFIX_nf90_put_var_2D_EightByteReal


  function MFIX_nf90_put_var_3D_EightByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    double precision, dimension(:, :, :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_3D_EightByteReal
#ifdef NETCDF
    MFIX_nf90_put_var_3D_EightByteReal = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_3D_EightByteReal = 1
#endif
  end function MFIX_nf90_put_var_3D_EightByteReal


  function MFIX_nf90_put_var_4D_EightByteReal(ncid, varid, values, start, count, stride, map)
    integer,                         intent( in) :: ncid, varid
    double precision, dimension(:, :, : , :), &
         intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer                                      :: MFIX_nf90_put_var_4D_EightByteReal
#ifdef NETCDF
    MFIX_nf90_put_var_4D_EightByteReal = nf90_put_var(ncid, varid, values, start, count, stride, map)
#else
    MFIX_nf90_put_var_4D_EightByteReal = 1
#endif
  end function MFIX_nf90_put_var_4D_EightByteReal

  function MFIX_usingNETCDF()
    logical :: MFIX_usingNETCDF
#ifdef NETCDF
    MFIX_usingNETCDF = .true.
#else
    MFIX_usingNETCDF = .false.
#endif
  end function MFIX_usingNETCDF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                          check_netcdf                             !
  !
  subroutine MFIX_check_netcdf( status )
    implicit none
    integer, intent ( in) :: status

#ifdef NETCDF
    if (status /= nf90_noerr) then
       write (*,*) ' ******************************************'
       write (*,*) trim(nf90_strerror(status))
       write (*,*) ' ******************************************'
    end if
#endif
  end subroutine MFIX_check_netcdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function MFIX_nf90_create(path, cmode, ncid, initialsize, chunksize)
    character (len = *), intent(in   ) :: path
    integer,             intent(in   ) :: cmode
    integer,             intent(inout) :: ncid
    integer, optional,   intent(in   ) :: initialsize
    integer, optional,   intent(inout) :: chunksize
    integer                            :: MFIX_nf90_create

    integer :: fileSize

#ifdef NETCDF
    MFIX_nf90_create = nf90_create(path, cmode, ncid, initialsize, chunksize)
#else
    MFIX_nf90_create = 0
#endif
    return
  end function MFIX_nf90_create


  function MFIX_nf90_def_dim(ncid, name, len, dimid)
    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent( in) :: len
    integer,             intent(inout) :: dimid
    integer                          :: MFIX_nf90_def_dim

#ifdef NETCDF
    MFIX_nf90_def_dim = nf90_def_dim(ncid, name, len, dimid)
#else
    MFIX_nf90_def_dim = 1
#endif
  end function MFIX_nf90_def_dim


  function MFIX_nf90_enddef(ncid, h_minfree, v_align, v_minfree, r_align)
    integer,           intent( in) :: ncid
    integer, optional, intent( in) :: h_minfree, v_align, v_minfree, r_align
    integer                        :: MFIX_nf90_enddef

#ifdef NETCDF
    MFIX_nf90_enddef = nf90_enddef(ncid, h_minfree, v_align, v_minfree, r_align)
#else
    MFIX_nf90_enddef = 1
#endif
  end function MFIX_nf90_enddef


  function MFIX_nf90_open(path, mode, ncid, chunksize)
    character (len = *), intent(in   ) :: path
    integer,             intent(in   ) :: mode
    integer,             intent(inout) :: ncid
    integer, optional,   intent(inout) :: chunksize
    integer                            :: MFIX_nf90_open
#ifdef NETCDF
    MFIX_nf90_open = nf90_open(path, mode, ncid, chunksize)
#else
    MFIX_nf90_open = 1
#endif
  end function MFIX_nf90_open


  function MFIX_nf90_inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimId, formatNum)
    integer,           intent( in) :: ncid
    integer, optional, intent(inout) :: nDimensions, nVariables, nAttributes, unlimitedDimId, formatNum
    integer                        :: MFIX_nf90_inquire

    integer :: nDims, nVars, nGAtts, unlimDimId, frmt

#ifdef NETCDF
    MFIX_nf90_inquire = nf90_inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimId, formatNum)
#else
    MFIX_nf90_inquire = 1
#endif
  end function MFIX_nf90_inquire


  function MFIX_nf90_inq_dimid(ncid, name, dimid)
    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent(inout) :: dimid
    integer                          :: MFIX_nf90_inq_dimid

#ifdef NETCDF
    MFIX_nf90_inq_dimid = nf90_inq_dimid(ncid, name, dimid)
#else
    MFIX_nf90_inq_dimid = 1
#endif
  end function MFIX_nf90_inq_dimid

  function MFIX_nf90_inq_varid(ncid, name, varid)
    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent(inout) :: varid
    integer                          :: MFIX_nf90_inq_varid

#ifdef NETCDF
    MFIX_nf90_inq_varid = nf90_inq_varid(ncid, name, varid)
#else
    MFIX_nf90_inq_varid = 1
#endif
  end function MFIX_nf90_inq_varid

  function MFIX_nf90_close(ncid)
    integer,             intent( in) :: ncid
    integer                          :: MFIX_nf90_close

#ifdef NETCDF
    MFIX_nf90_close = nf90_close(ncid)
#else
    MFIX_nf90_close = 1
#endif
  end function MFIX_nf90_close

  subroutine MFIX_ncvinq(ncid,varid,varnam,vartyp,nvdims,vdims,nvatts,rcode)
    implicit none

    integer       :: ncid , varid , vartyp , nvdims , nvatts , rcode
    integer       :: vdims(*)
    character(len=*)    varnam

#ifdef NETCDF
    call ncvinq(ncid,varid,varnam,vartyp,nvdims,vdims,nvatts,rcode)
#endif
    return
  end subroutine MFIX_ncvinq


end module MFIX_netcdf
